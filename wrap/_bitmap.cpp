#include <pybind11/pybind11.h>
#include "_bitmap.hpp"
void wrap_bitmap(pybind11::module &m){
  using namespace polystar;
  using namespace polystar::bitmap;
  namespace py = pybind11;

  py::class_<Bitmap<double>> d(m, "BitmapD");
  define_bitmap_inits<double>(d);
  define_bitmap<double>(d);

  py::class_<Bitmap<float>> f(m, "BitmapF");
  define_bitmap_inits<float>(f);
  define_bitmap<float>(f);

  py::class_<Bitmap<int>> i(m, "BitmapI");
  define_bitmap_inits<int>(i);
  define_bitmap<int>(i);

}