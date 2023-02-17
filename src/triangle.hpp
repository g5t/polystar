#pragma once
#include "types.hpp"
#include "polygon_network.hpp"

namespace polystar::triangle {
#include "triangle.h"
#include "triangle_api.h"

  template<class X, class T, template<class> class A, class W>
    polystar::polygon::Network<X,T,A> triangulate(const A<T> & v, const W & ws){
      context *ctx;
//      behavior *b;
      triangleio in;

      // fill the triangleio:
      in.pointlist = new double[20*v.size(0)]; // allocate too much memory Triangle decides to add points
      in.numberofpoints = static_cast<int>(v.size(0));
      in.pointattributelist = (double *) nullptr;
      in.pointmarkerlist = (int *) nullptr;
      in.numberofpointattributes = 0;
      for (ind_t i=0; i<v.size(0); ++i){
        in.pointlist[2*i] = v.val(i,0);
        in.pointlist[2*i+1] = v.val(i,1);
      }

      in.trianglelist = (int *) nullptr;
      in.triangleattributelist = (double *) nullptr;
      in.trianglearealist = (double *) nullptr;
      in.neighborlist = (int *) nullptr;
      in.numberoftriangles = 0;
      in.numberofcorners = 0;
      in.numberoftriangleattributes = 0;

      auto border = ws.border();
      auto holes = ws.wires();
      size_t hole_edge_count{0};
      for (const auto & h: holes) hole_edge_count += h.size();
      in.segmentlist = new int[20*(border.size() + hole_edge_count)];
      in.numberofsegments = static_cast<int>(border.size()+hole_edge_count);
      in.segmentmarkerlist = nullptr;

      size_t offset{0};
      for (size_t i=0; i<border.size(); ++i){
        auto e = border.edge(i);
        in.segmentlist[offset++] = e.first;
        in.segmentlist[offset++] = e.second;
      }
      for (const auto & hole: holes){
        for (size_t i=0; i<hole.size(); ++i){
          auto e = hole.edge(i);
          in.segmentlist[offset++] = e.first;
          in.segmentlist[offset++] = e.second;
        }
      }

      // In Triangle speak, a hole is any point *inside* of a hole. With its boundary already specified above.
      in.numberofholes = static_cast<int>(holes.size());
      in.holelist = new double[2*holes.size()];
      in.numberofregions = 0;
      in.regionlist = (double *) nullptr;
      // take the average of the hole vertices as the in-hole point?
      // NO! Non-convex holes can have their center of mass outside their boundary.
      offset = 0;
      for (const auto & hole: holes){
        auto h = hole.find_contained_point(v);
//        auto h = 0 * v.view(0);
//        for (const auto & i: hole) h += v.view(i);
//        h /= static_cast<double>(hole.size());
        in.holelist[offset++] = h.val(0, 0);
        in.holelist[offset++] = h.val(0, 1);
      }

      in.edgelist = (int *) nullptr;
      in.edgemarkerlist = (int *) nullptr;
      in.numberofedges = 0;


      ctx = triangle_context_create();
      char cmdline[] = "zpN";
      int status;
      status = triangle_context_options(ctx, cmdline);
      if (status != TRI_OK){
        std::cout << "Error setting triangle context\n";
      }
      status = triangle_mesh_create(ctx, &in);
      if (status != TRI_OK){
        std::cout << "Error triangulating\n";
      }
      status = triangle_mesh_copy(ctx, &in, 0 /*write edges*/, 0 /*write neighbors*/);
      if (status != TRI_OK){
        std::cout << "Error copying mesh\n";
      }

      // the in structure also contains output results
      std::vector<X> triangles;
      if (in.numberoftriangles < 0) {
        throw std::runtime_error("Found less than no triangles");
      }
      triangles.reserve(in.numberoftriangles);
      offset=0;
      auto poly_size = static_cast<size_t>(in.numberofcorners);
      for (int i=0; i<in.numberoftriangles; ++i){
        X tri;
        tri.reserve(in.numberofcorners);
        for (size_t j=0; j<poly_size; ++j) tri.push_back(in.trianglelist[offset++]);
        // protect against inserted extra vertices?
        //if (std::all_of(tri.begin(), tri.end(), [&](const auto n){return n < v.size(0);}))
        triangles.push_back(tri);
      }
      auto out_v = A<T>(in.numberofpoints, 2);
      for (int i=0; i<in.numberofpoints; ++i){
        out_v.val(i, 0) = in.pointlist[2*i];
        out_v.val(i, 1) = in.pointlist[2*i + 1];
      }
      auto network = polystar::polygon::Network<X,T,A>(triangles, out_v);

      // cleanup
      triangle_context_destroy(ctx);
      // delete the memory we allocated
      delete[] in.pointlist;
      delete[] in.segmentlist;
      delete[] in.holelist;
      // plus free what Triangle (may have) allocated
      if (in.pointattributelist) free(in.pointattributelist);
      if (in.pointmarkerlist) free(in.pointmarkerlist);
      if (in.trianglelist) free(in.trianglelist);
      if (in.triangleattributelist) free(in.triangleattributelist);
      if (in.trianglearealist) free(in.trianglearealist);
      if (in.neighborlist) free(in.neighborlist);
      if (in.segmentmarkerlist) free(in.segmentmarkerlist);
      if (in.regionlist) free(in.regionlist);
      if (in.edgelist) free(in.edgelist);
      if (in.edgemarkerlist) free(in.edgemarkerlist);

    return network;
    }
}