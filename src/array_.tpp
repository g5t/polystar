/* This file is part of brille.

Copyright © 2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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

/*!
\file
\author Greg Tucker
\brief Implementations of template functions for Array2 and LatVec objects
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS // opeator overloading macros
// In order to overload operations between bArrays and LatVecs unambiguously
// these definitions can only be made after all types have been defined:

// 'pure' bArray<T> [+-*/] 'pure' bArray<R>
// (or derived classes with only extra static properties or methods)
#define ARRAY_ARRAY_OP(X) template<class T, class R, template<class> class A, class S = std::common_type_t<T,R>>\
std::enable_if_t<bareArrays<T,A,R,A>, A<S>>\
operator X (const A<T>& a, const A<R>& b){\
  auto itr = a.broadcastItr(b);\
  A<S> out(itr.shape());\
  for (auto [ox, ax, bx]: itr) out[ox] = a[ax] X b[bx];\
  return out;\
}
ARRAY_ARRAY_OP(+)
ARRAY_ARRAY_OP(-)
ARRAY_ARRAY_OP(*)
ARRAY_ARRAY_OP(/)
#undef ARRAY_ARRAY_OP

// any derived bArray<T> class with a copy constructor [+-*/] scalar:
#define ARRAY_SCALAR_OP(X,Y) template<class T, class R, template<class> class A, class S = std::common_type_t<T,R>>\
std::enable_if_t<isArray<T,A>, A<S>>\
operator X (const A<T>& a, const R b){\
  A<S> out = A<S>(a).decouple(); \
  out Y b;\
  return out;\
}
ARRAY_SCALAR_OP(+,+=)
ARRAY_SCALAR_OP(-,-=)
ARRAY_SCALAR_OP(*,*=)
ARRAY_SCALAR_OP(/,/=)
#undef ARRAY_SCALAR_OP

// scalar [+-*] any derived bArray<T> class with a copy constructor:
#define SCALAR_ARRAY_OP(X,Y) template<class T, class R, template<class> class A, class S = std::common_type_t<T,R>>\
std::enable_if_t<isArray<T,A>, A<S>>\
operator X (const R b, const A<T>& a){\
  A<S> out = A<S>(a).decouple();\
  out Y b;\
  return out;\
}
SCALAR_ARRAY_OP(+,+=)
SCALAR_ARRAY_OP(-,-=)
SCALAR_ARRAY_OP(*,*=)
#undef SCALAR_ARRAY_OP

// any derived bArray<T> with a copy constructor [+-*/] std::array<R,N>
// broadcasts the std::array but only if N matches the size of the last dimension of the bArray
#define ARRAY_STDA_OP(X,Y) template<class T, template<class> class A, class R, class S = std::common_type_t<T,R>, size_t Nel>\
std::enable_if_t<isArray<T,A>, A<S>>\
operator X (const A<T>& a, const std::array<R,Nel>& b){\
  assert(a.shape().back() == Nel);\
  A<S> out = A<S>(a).decouple();\
  auto sh = a.shape();\
  sh.back()=0;/*fix the last dimension of the iterator*/\
  for (auto x: out.subItr(sh)) for (size_t i=0; i<Nel; ++i){\
    x.back() = static_cast<brille::ind_t>(i);\
    out[x] Y b[i];\
  }\
  return out;\
}
ARRAY_STDA_OP(+,+=)
ARRAY_STDA_OP(-,-=)
ARRAY_STDA_OP(*,*=)
ARRAY_STDA_OP(/,/=)
#undef ARRAY_STDA_OP

#endif // operator overloading macros

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*==============================================================================
  The broadcasting of cross is intimately dependent on the implementation of
  the Array class. If Array2 is used then pop() is not a method of shape_t!
  So this variant must be used instead:
==============================================================================*/
  // cross (Array × Array)
  template<class T, class R, template<class> class L>
  std::enable_if_t<bareArrays<T,L,R,L>, L<double>>
  cross(const L<T>& a, const L<R>& b) {
    using namespace brille::utils;
    assert( a.size(1) == 3 && b.size(1)==3 );
    assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
    brille::ind_t aN=a.size(0), bN=b.size(0);
    assert( 1u==aN || 1u==bN || aN==bN );
    brille::ind_t oO = (1u == aN) ? bN : aN;
    bArray<double> oarray(oO, 3u); // row-ordered contiguous
    if (1u == aN || 1u == bN){
      if (1u == aN){
        for (brille::ind_t j=0; j<bN; ++j) vector_cross<double,T,R,3>(oarray.ptr(j), a.ptr(0), b.ptr(j));
      } else {
        for (brille::ind_t j=0; j<aN; ++j) vector_cross<double,T,R,3>(oarray.ptr(j), a.ptr(j), b.ptr(0));
      }
    } else {
      for (brille::ind_t j=0; j<aN; ++j) vector_cross<double,T,R,3>(oarray.ptr(j), a.ptr(j), b.ptr(j));
    }
    return oarray;
  }
#else
  /*! \brief Find the cross product of two Array2

  \param first The first Array2, with shape `(Na, 3)`
  \param second The second Array2, with shape `(Nb, 3)`

  Both arrays must hold the same number of 3-vectors, `Na = Nb`, or either `Na`
  or `Nb` must be one.
  Singleton vectors are broadcast across the other array.

  \return An array of the same data type as input with shape `(N,3)`
          where `N = max(Na, Nb)`
  */
  template<class T, class R, template<class> class A>
  Array2<double> cross(const A<T>& first, const A<R>& second);
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// dot (Array2 ⋅ Array2)
template<class T, class R, template<class> class A>
std::enable_if_t<bareArrays<T,A,R,A>, A<double>>
dot(const A<T>& a, const A<R>& b) {
  using namespace brille::utils;
  assert( a.size(1) == b.size(1) );
  assert( a.is_row_ordered() && b.is_row_ordered() && a.is_contiguous() && b.is_contiguous() );
  brille::ind_t aN=a.size(0), bN=b.size(0), d=a.size(1);
  assert( 1u==aN || 1u==bN || aN==bN );
  brille::ind_t oO = (1u == aN) ? bN : aN;
  bArray<double> oarray(oO, 1u);
  if (1u==aN || 1u==bN) {
    if (1u==aN){
      for (brille::ind_t i=0; i<bN; ++i) oarray.val(i,0) = vector_dot<double,T,R>(d, a.ptr(0), b.ptr(i));
    } else {
      for (brille::ind_t i=0; i<aN; ++i) oarray.val(i,0) = vector_dot<double,T,R>(d, a.ptr(i), b.ptr(0));
    }
  } else {
    for (brille::ind_t i=0; i<aN; ++i) oarray.val(i,0) = vector_dot<double,T,R>(d, a.ptr(i), b.ptr(i));
  }
  return oarray;
}
// #else
  /*! \brief Find the dot product of two Array2

  \param first The first Array2, with shape `(Na, M)`
  \param second The second Array2, with shape `(Nb, M)`

  Both arrays must have the same number of elements per vector, `M`, and either
  the same number of vectors, `Na = Nb`, or either `Na` or `Nb` must be one.
  Singleton vectors are broadcast across the other array.

  \return An Array2 with shape `(N,1)` where `N = max(Na, Nb)`
  */
  /* template<class T, class R, template<class> class A> */
  /* Array2<double> dot(const A<T>& first, const A<R>& second); */
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // [bArray] norm
  template<class T, template<class> class L>
  std::enable_if_t<isBareArray<T,L>, L<double>>
  norm(const L<T> &a){
    L<double> out = dot(a,a);
    for (auto& x : out.valItr()) x = std::sqrt(x);
    return out;
  }
#else
  /*! \brief Find the length of each vector in an Array2

  \param array The Array2 to find the norm(s) of, with shape `(N,M)`
  \return a `(N,1)` Array2 with the vector norms for each vector in `array`
  */
  template<class T, template<class> class A> Array2<double> norm(const A<T>& array);
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // [bArray],[bArray] cat
  template<class T, class R, template<class> class L, class S = std::common_type_t<T,R>>
  std::enable_if_t<bareArrays<T,L,R,L>, L<S>>
  cat(const brille::ind_t dim, const L<T>& a, const L<R>& b){
    return L<S>(a.decouple()).append(dim, b);
  }
#else
  /*! \brief Concatenate two Array2 or LatVec arrays along the specified direction

  \param dim The dimension to concatenate along
  \param a The first Array2
  \param b The second Array2
  \return The concatenation of `a` and `b` along the dimension `dim`
  */
  template<class T, class R, template<class> class A, class S=std::common_type_t<T,R>>
  A<S> cat(const ind_t dim, const A<T>& a, const A<R>& b);
#endif
/*! \brief variadic Array2 concatenation

\param dim The dimension to concatenate along (probably only 0 makes sense)
\param a The first Array2 or LatVec
\param b The second Array2 or LatVec
\param args Any number of additional Array2 or LatVecs
\return All arrays concatenated in order along the indicated dimension
*/
template<class T, class R, template<class> class A, class S=std::common_type_t<T,R>, class... Args>
std::enable_if_t<bothArrays<T,A,R,A>, A<S>>
cat(const brille::ind_t dim, const A<T>& a, const A<R>& b, Args... args){
  return cat(dim,cat(dim,a,b),args...);
}


#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // -('pure' brille:Array)
  template<class T, template<class> class L>
  std::enable_if_t<(isBareArray<T,L>&&!std::is_unsigned_v<T>), L<T>>
  operator-(const L<T>& a){
    return -1*a;
  }
#else
  /*! \brief Find the inverse of a bare Array2 or LatVec

  \param array Eitehr a bare Array2 or one of its LatVec subclasses
  \return -array
  \note This operator is not overloaded for arrays holding unsigned integers
  */
  template<class T, template<class> class A>
  A<T> operator-(const A<T>& array);
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // sum(bArray)
  template<class T, template<class> class A>
  std::enable_if_t<isBareArray<T,A>, A<T>>
  sum(const A<T>& a, typename A<T>::ind_t dim){
    return a.sum(dim);
  }
#else
  /*! \brief Sum over one of the dimensions of a bare Array2 or LatVec

  \param array Either a bare Array2 or one of its LatVec subclasses
  \param dim The dimension to sum over (and reduce to singleton)
  \return (1,M) or (N,1) array
  */
  template<class T, template<class> class A>
  A<T> sum(const A<T>& array, typename A<T>::ind_t dim);
#endif


#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // abs(bArray)
  template<class T, template<class> class L>
  std::enable_if_t<isBareArray<T,L>, L<T>>
  abs(const L<T>& a){
    L<T> out(a.shape());
    for (auto x : a.subItr()) out[x] = brille::utils::magnitude(a[x]);
    return out;
  }
#else
  /*! \brief Find the elementwise absolute value of a bare Array2 or LatVec

  \param array Either a bare Array or one of its LatVec subclasses
  \return |array|
  */
  template<class T, template<class> class L> L<T> abs(const L<T>& array);
#endif
