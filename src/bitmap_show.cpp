#include <locale>
#include <codecvt>
#include <sstream>

#include "bitmap_show.hpp"
#include "bitmap.hpp"

using namespace polystar;
using namespace polystar::bitmap;

std::wstring polystar::bitmap::s2ws(const std::string& str){
  using convert_typeX = std::codecvt_utf8<wchar_t>;
  std::wstring_convert<convert_typeX, wchar_t> converterX;
  return converterX.from_bytes(str);
}

std::string polystar::bitmap::ws2s(const std::wstring& wstr) {
  using convert_typeX = std::codecvt_utf8<wchar_t>;
  std::wstring_convert<convert_typeX, wchar_t> converterX;
  return converterX.to_bytes(wstr);
}

std::array<wchar_t, 6> polystar::bitmap::border_chars(int choice){
  // return {top-left, top-right, bottom-left, bottom-right, horizontal, vertical} characters
  switch (choice) {
    case 1: return {{0x2554, 0x2557, 0x255A, 0x255D, 0x2550, 0x2551}}; // double lines
    case 2: return {{0x256D, 0x256E, 0x2570, 0x256F, 0x2500, 0x2502}}; // light arced lines
    case 3: return {{0x250F, 0x2513, 0x2517, 0x251B, 0x2501, 0x2503}}; // bold lines
    default: return {{0x250C, 0x2510, 0x2514, 0x2518, 0x2500, 0x2502}}; // light lines
  }
}