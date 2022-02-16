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

#ifndef BRILLE_ARRAY__HPP_
#define BRILLE_ARRAY__HPP_
#include <typeinfo> // for std::bad_cast
#include <exception>
/*! \file
    \brief Forward declarations for Array* classes which can be constructed
           from each other.
*/

namespace brille {
  // declare both Array and Array2 so that they can each contain conversion
  // constructors for the other.
  template<class T> class Array;
  template<class T> class Array2;

  // declare their itterators here
  template<class T> class ArrayIt;
  template<class T> class Array2It;
}

#include "array2.hpp"
#include "array.hpp"

namespace brille{
//! An alias used while deciding between Array and Array2 for data storage
template<class T>
using bArray = Array2<T>;


/*! \brief Templated struct to help differentiate between Array2 and its derived
           classes LQVec and LDVec

Since the derived classes add Lattice information some functions behave
differently for 'raw' Array2 versus LQVec or LDVec objects.
This struct is used in template arguments to cause substitution failure-based
differentiation of otherwise identical function signatures.
*/
template<class... T> struct ArrayTraits{
  static constexpr bool array = false;  //!< is this an Array2, LQVec, or LDVec
  static constexpr bool latvec = false; //!< is this an LQVec or LDVec
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T> struct ArrayTraits<bArray<T>>{
  static constexpr bool array = true;
  static constexpr bool latvec = false;
};
#endif

/*! Easier access to ArrayTraits for template substitution

true for Array2, LQVec, LDVec
*/
template<class T, template<class> class A>
inline constexpr bool isArray = ArrayTraits<A<T>>::array;
/*! Easier access to ArrayTraits for template substitution

true for LQVec, LDVec
*/
template<class T, template<class> class A>
inline constexpr bool isLatVec = ArrayTraits<A<T>>::latvec;
/*! Easier access to ArrayTraits for template substitution

true for Array2
*/
template<class T, template<class> class A>
inline constexpr bool isBareArray = isArray<T,A> && !isLatVec<T,A>;

/*! Easier access to ArrayTraits for double template substitution

true for any combination of two Array2, LQVec, LDVec objects
*/
template<class T, template<class> class A, class R, template<class> class B>
inline constexpr bool bothArrays = isArray<T,A> && isArray<R,B>;
/*! Easier access to ArrayTraits for double template substitution

true for two Array2 objects
*/
template<class T, template<class> class A, class R, template<class> class B>
inline constexpr bool bareArrays = isBareArray<T,A> && isBareArray<R,B>;
/*! Easier access to ArrayTraits for double template substitution

true for any combination of two LQVec, LDVec objects
*/
template<class T, template<class> class A, class R, template<class> class B>
inline constexpr bool bothLatVecs = isLatVec<T,A> && isLatVec<R,B>;

//! Check if both vectors are the same
template<class T>
bool equal_shapes(const std::vector<T>& a, const std::vector<T>&b){
  bool ok = a.size() == b.size();
  if (ok) ok = std::equal(a.begin(), a.end(), b.begin());
  return ok;
}
//! Check if both arrays are the same
template<class T, size_t Nel>
bool equal_shapes(const std::array<T,Nel>& a, const std::array<T,Nel>&b){
  return std::equal(a.begin(), a.end(), b.begin());
}


#include "array_.tpp"

} // end namespace

#endif
