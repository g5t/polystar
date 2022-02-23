/* This file is part of brille

Copyright © 2019-2022 Greg Tucker <gregory.tucker@ess.eu>

brille is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

brille is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with brille. If not, see <https://www.gnu.org/licenses/>.            */

#ifndef BRILLE_GEOMETRY_H_
#define BRILLE_GEOMETRY_H_
/*! \file
    \author Greg Tucker
    \brief A class representing convex polyhedra in three dimensions
*/
#include <random>
#include <utility>
#include "tetgen.h"
#include "array_.hpp" // defines bArray
#include "utilities.hpp"
namespace brille {

  template<class T>
  std::tuple<bArray<T>, bArray<T>> plane_points_from_normal(const bArray<T> &n, const bArray<T> &p) {
    bArray<T> a(n.shape(), T(1)), b;
    auto n_norm = norm(n); // this *should* be all 1, but maybe it's not
    auto nn = n / n_norm;
    auto a_dot_n = dot(a, nn) * nn;
    auto is_n = a_dot_n.row_is(brille::cmp::eq, nn * norm(a));
    if (is_n.any()) {
      for (ind_t i = 0; i < a.size(0); ++i)
        if (is_n.val(i)) {
          // set to (2, 0, 1) so that subtracting (1, 1, 1) gives (1, -1, 0)
          a[{i, 0}] = T(2);
          a[{i, 1}] = T(0);
          for (ind_t j = 2; j < a.size(1); ++j) a[{i, j}] = T(1);
        }
    }
    a -= a_dot_n;
    a /= norm(a);
    b = cross(nn, a);  // so that the points (p, a, b) define a plane with normal n
    // the returned a and b should have the property (x - p) . n == 0
    return std::make_tuple(a + p, b + p);
  }

  template<class T>
  bArray<T> three_point_normal(const bArray<T> &a, const bArray<T> &b, const bArray<T> &c) {
    auto n = cross(b - a, c - b);
    return n / norm(n);
  }

  template<class T, class I>
  bArray<T> three_point_normal(const bArray<T> &p, const I a, const I b, const I c) {
    return three_point_normal(p.view(a), p.view(b), p.view(c));
  }

  template<class T, class I>
  bArray<T> three_point_normal(const bArray<T> &p, const std::vector<I> &f) {
    return three_point_normal(p.view(f[0]), p.view(f[1]), p.view(f[2]));
  }

  template<class T>
  bArray<T> three_point_normal(const bArray<T> &p, const bArray<ind_t> &t) {
    bArray<T> out(t.size(0), 3u);
    for (ind_t i = 0; i < t.size(0); ++i)
      out.set(i, three_point_normal(p, t[{i, 0}], t[{i, 1}], t[{i, 2}]));
    return out;
  }


//! Determine if a facet of a Polyhedron is part of a convex volume
  template<typename T>
  static bool is_not_dangling(const std::vector<size_t> &counts, const std::vector<T> &face) {
    return std::all_of(face.begin(), face.end(), [counts](const T &x) { return counts[x] > 2u; });
  }

  template<class T> static void on_plane_vector(const bArray<T>& x, const bArray<T>& y, const bArray<T>& v3, T* v2){
    v2[0] = dot(v3, x).sum();
    v2[1] = dot(v3, y).sum();
  }

//  //! Determine the winding angles for the vertices of a Polyhedron facet
//  template<class T>
//  std::vector<T> bare_winding_angles(const bArray<T> &vecs, const ind_t i, const bArray<T> &n) {
//    if (vecs.ndim() != 2 || vecs.size(1) != 3)
//      throw std::runtime_error("Finding a winding angle requires the cross product, which is only defined in 3D");
//    // vecs should be normalized already
//    std::vector<T> angles(vecs.size(0), 0.);
//    T dotij, y_len, angij;
//    bArray<T> x(1u, 3u), y(1u, 3u); // ensure all have memory allocated
//    T crsij[3]; // to hold the cross product
//    for (ind_t j = 0; j < vecs.size(0); ++j)
//      if (i != j) {
//        dotij = vecs.dot(i, j);
//        brille::utils::vector_cross(crsij, vecs.ptr(i), vecs.ptr(j));
//        x = dotij * vecs.view(i);
//        y = vecs.view(j) - x;
//        y_len = y.norm(0) * (std::signbit(brille::utils::vector_dot(crsij, n.ptr(0))) ? -1 : 1);
//        angij = std::atan2(y_len, dotij);
//        angles[j] = angij < 0 ? angij + 2 * brille::pi : angij;
//      }
//    return angles;
//  }

  template<class T>
  [[nodiscard]] std::vector<T> bare_winding_angles(const bArray<T> & vecs, const ind_t i, const bArray<T> & n){
    T centre[2]{0, 0}, a[2]{1, 0}, b[2];
    // double check that we're dealing with vectors in a plane:
    auto nn = norm(n);
    auto pv = vecs - dot(vecs, n / nn) * (n / nn);
    auto x = pv.view(i) / norm(pv.view(i));
    auto y = cross(n / nn, x);
    std::vector<T> angles(vecs.size(0), 0.);
    for (ind_t j=0; j<vecs.size(0); ++j) if(i != j) {
      on_plane_vector(x, y, pv.view(j), b);
      auto angij = std::atan2(b[1], b[0]);
      angles[j] = angij < 0 ? angij + 2 * brille::pi : angij;
    }
    return angles;
  }


  //! Check whether the points defining a Polyhedron facet encompass a finite area
  template<class T>
  bool face_has_area(const bArray<T> &points, const bool strict = false) {
    // first verify that all points are coplanar
    // pick the first three points to define a plane, then ensure all points are in it
    if (points.size(0) < 3) return -2; // can't be a face
    // move to the 2-D face coordinate system so that we can use orient2d
    auto centre = points.sum(0) / static_cast<T>(points.size(0));
    auto facet = points - centre;
    // pick the first on-face vector to be our x-axis
    auto x = facet.view(0).decouple(); // decouple since we're going to normalise this
    auto z = three_point_normal(facet, 0, 1, 2);
    auto y = cross(z, x);
    x /= norm(x);
    y /= norm(y);
    T a[2] = {1, 0}, b[2], c[2];
    T s{0};
    for (ind_t i = 1; i < facet.size(0) - 1; ++i) {
      on_plane_vector(x, y, facet.view(i), b);
      on_plane_vector(x, y, facet.view(i+1), c);
      auto piece = orient2d(a, b, c);
      if (!strict) piece = std::abs(piece);
      s += piece;
    }
    return (s > T(0));
  }

  /*! \brief Check if a point is 'behind' a plane

  \param plane_a the first point defining the plane
  \param plane_b the second point defining the plane
  \param plane_c the third point defining the plane
  \param x the point or points to check
  \returns whether `x` is on or behind the plane

  \note Uses the geometry predicates orient3d to perform the check within machine precision
  */
  template<class T>
  std::vector<bool>
  point_inside_plane(const bArray<T> &plane_a, const bArray<T> &plane_b, const bArray<T> &plane_c, const bArray<T> &x) {
    assert(plane_a.numel() == 3 && plane_b.numel() == 3 && plane_c.numel() == 3);
    assert(plane_a.is_contiguous() && plane_b.is_contiguous() && plane_c.is_contiguous());
    assert(x.is_contiguous() && x.is_row_ordered() && x.size(1u) == 3);
    std::vector<bool> pip;
    for (ind_t i = 0; i < x.size(0u); ++i) {
      auto o3d = orient3d(plane_a.ptr(0), plane_b.ptr(0), plane_c.ptr(0), x.ptr(i));
      pip.push_back(o3d >= 0.);
    }
    return pip;
  }

  template<class I>
  std::vector<I> two_in_one_indexes(const std::vector<std::vector<I>> &one, const std::vector<I> &two) {
    std::vector<I> indexes;
    indexes.reserve(one.size());
    auto has = [one](const size_t &i, const I &x) {
      return std::find(one[i].begin(), one[i].end(), x) != one[i].end();
    };
    for (size_t index = 0; index < one.size(); ++index) {
      if (std::any_of(two.begin(), two.end(), [index, has](const I &x) { return has(index, x); }))
        indexes.push_back(static_cast<I>(index));
    }
    return indexes;
  }

  template<class I>
  std::vector<I> keep_to_cut_list(const std::vector<bool> &keep, const std::vector<std::vector<I>> &faces) {
    std::vector<int> to_remove;
    for (size_t j = 0; j < keep.size(); ++j) if (!keep[j]) to_remove.push_back(static_cast<I>(j));
    return two_in_one_indexes(faces, to_remove);
  }


  /*! \brief Find the common-edge vertex indices between two facets

  \param  one the first facet vertex indexes
  \param  two the second facet vertex indexes
  \param  strict whether the two vertex indexes must be neighbours and in opposite order on the two facets
  \returns A tuple containing a success value and the two indexes; if there is no common edge the success value is false.
  */
  template<class I>
  std::tuple<bool, I, I> common_edge(const std::vector<I> &one, const std::vector<I> &two, const bool strict = false) {
    bool ok{false};
    I a{0}, b{0};
    if (strict) {
      for (size_t i = 0; !ok && i < one.size(); ++i) {
        a = one[i];
        b = one[(i + 1) % one.size()];
        for (size_t j = 0; !ok && j < two.size(); ++j) {
          ok = b == two[j] && a == two[(j + 1) % two.size()];
        }
      }
    } else {
      for (size_t ia = 0; !ok && ia < one.size(); ++ia) {
        a = one[ia];
        for (size_t ib = ia + 1; !ok && ib < one.size() + 1; ++ib) {
          b = one[ib % one.size()];
          for (size_t j = 0; !ok && j < two.size(); ++j) {
            const auto &ta{two[j]};
            for (size_t k = j + 1; !ok && k < two.size() + 1; ++k) {
              const auto &tb{two[k % two.size()]};
              ok = (ta == a && tb == b) || (tb == a && ta == b);
            }
          }
        }
      }
    }
    return std::make_tuple(ok, a, b);
  }

  /*! \brief Find the intersection of the edge between two facets and a plane defined by three points
  \param v      the vertices of the facets
  \param vpf    the vertex indices for each known facet
  \param j      the first facet index (accesses vpf)
  \param k      the second facet index (accesses vpf)
  \param a      the first point on the plane
  \param b      the second point on the plane
  \param c      the third point on the plane
  \param strict whether the facet vertex indices are sorted correctly
  \returns A bArray containing the intersection point, or an empty bArray if there is no intersection
  */
  template<class T, class I>
  bArray<T> edge_plane_intersection(
    const bArray<T> &v, const std::vector<I> &one, const std::vector<I> &two,
    const bArray<T> &a, const bArray<T> &b, const bArray<T> &c, const bool strict = false) {
    // find the correct pair of vertices which form the edge:
    auto[ok, j, k] = common_edge(one, two, strict);
    if (ok) {
      // check whether they're on opposite sides of the plane (a, b, c)
      auto o3dj = orient3d(a.ptr(0), b.ptr(0), c.ptr(0), v.ptr(j));
      auto o3dk = orient3d(a.ptr(0), b.ptr(0), c.ptr(0), v.ptr(k));
      // both on one side, both on the other side, or both *in* the plane all preclude a single intersection
      if ((o3dj < 0 && o3dk < 0) || (o3dj > 0 && o3dk > 0) || (o3dj == 0 && o3dk == 0)) ok = false;
      if (brille::approx::scalar(o3dj, 0.) && brille::approx::scalar(o3dk, 0.)) ok = false;
    }
    bArray<T> at;
    if (ok) {
      auto plane_n = three_point_normal(a, b, c);
      auto line_u = v.view(k) - v.view(j);
      auto denominator = dot(plane_n, line_u);
      auto numerator = dot(plane_n, a - v.view(j));
      assert(denominator.abs().sum() > 0);
      at = v.view(j) + (numerator * line_u) / denominator;
//      auto scalar = (numerator / denominator).sum();
//      at = v.view(j) + scalar * line_u;
//      debug_update_if(orient3d(a.ptr(0), b.ptr(0), c.ptr(0), at.ptr(0)) != 0.,
//                      "The found intersection point ", at.to_string(0),
//                      " is off the plane, proportional to ", orient3d(a.ptr(0), b.ptr(0), c.ptr(0), at.ptr(0)));
    }
    return at;
  }

}
#endif