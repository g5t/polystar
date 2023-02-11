#ifndef POLYSTAR_BITMAP_SHOW_HPP
#define POLYSTAR_BITMAP_SHOW_HPP
#include <string>
#include <array>
#include <vector>
#include "bitmap.hpp"

namespace polystar::bitmap {
  std::wstring s2ws(const std::string& str);
  std::string ws2s(const std::wstring& wstr);
  std::array<wchar_t, 6> border_chars(int choice);

  template<class T>
  char16_t braille_char(T lt, T lm, T ll, T lb, T rt , T rm, T rl, T rb){
    char16_t b = 0x2800;
    if (lt) b += 0x0001;
    if (lm) b += 0x0002;
    if (ll) b += 0x0004;
    if (lb) b += 0x0040;
    if (rt) b += 0x0008;
    if (rm) b += 0x0010;
    if (rl) b += 0x0020;
    if (rb) b += 0x0080;
    return b;
  }

  template<class T>
  std::string braille_map(const Array2<T>& map) {
    // Complete with box-drawing character borders https://en.wikipedia.org/wiki/Box-drawing_character
    auto ny = map.size(0);
    auto nx = map.size(1);
    std::wstringstream ss;
    auto bc = border_chars(2);
    ss << bc[0];
    for (ind_t x=0; x<nx; x+=2) ss << bc[4];
    ss << bc[1] << "\n";
    for (ind_t y=0; y < ny; y+=4){
      ss << bc[5];
      for (ind_t x=0; x < nx; x+=2){
        auto tl =                        map.val(y, x);
        auto tr = x+1 < nx             ? map.val(y,   x+1) : T(0);
        auto ml = y+1 < ny             ? map.val(y+1, x  ) : T(0);
        auto mr = y+1 < ny && x+1 < nx ? map.val(y+1, x+1) : T(0);
        auto ll = y+2 < ny             ? map.val(y+2, x  ) : T(0);
        auto lr = y+2 < ny && x+1 < nx ? map.val(y+2, x+1) : T(0);
        auto bl = y+3 < ny             ? map.val(y+3, x  ) : T(0);
        auto br = y+3 < ny && x+1 < nx ? map.val(y+3, x+1) : T(0);
        ss << static_cast<wchar_t>(braille_char(tl, ml, ll, bl, tr, mr, lr, br));
      }
      ss << bc[5] << "\n";
    }
    ss << bc[2];
    for (ind_t x=0; x<nx; x+=2) ss << bc[4];
    ss << bc[3] << "\n";
    auto wstr = ss.str();
    return ws2s(wstr);
  }

  template<class T>
  char16_t block_char(T lt, T lb, T rt, T rb){
    auto code{(lt ? 0b1000 : 0b0000) + (rt ? 0b0100: 0b0000) + (lb ? 0b0010 : 0b0000) + (rb ? 0b0001: 0b0000)};
    switch (code){
      case 0b1111: return 0x2588;
      case 0b1110: return 0x259B;
      case 0b1101: return 0x259C;
      case 0b1100: return 0x2580;
      case 0b1011: return 0x2599;
      case 0b1010: return 0x258C;
      case 0b1001: return 0x259A;
      case 0b1000: return 0x2598;
      case 0b0111: return 0x259F;
      case 0b0101: return 0x2590;
      case 0b0100: return 0x259D;
      case 0b0011: return 0x2584;
      case 0b0010: return 0x2596;
      case 0b0001: return 0x2597;
      case 0b0000: return 0x2591;
      default: return 0x2591;
    }
  }

  template<class T>
  std::string block_map(const Array2<T>& map, bool border=false) {
    // Complete with box-drawing character borders https://en.wikipedia.org/wiki/Box-drawing_character
    auto ny = map.size(0);
    auto nx = map.size(1);
    std::wstringstream ss;
    if (border) {
      ss << wchar_t(0x2554);
      for (ind_t x = 0; x < nx; x += 2) ss << wchar_t(0x2550);
      ss << wchar_t(0x2557) << "\n";
    }
    for (ind_t y=0; y < ny; y+=2){
      if (border) ss << wchar_t(0x2551);
      for (ind_t x=0; x < nx; x+=2){
        auto tl =                        map.val(y, x);
        auto tr = x+1 < nx             ? map.val(y,   x+1) : T(0);
        auto ml = y+1 < ny             ? map.val(y+1, x  ) : T(0);
        auto mr = y+1 < ny && x+1 < nx ? map.val(y+1, x+1) : T(0);
        ss << static_cast<wchar_t>(block_char(tl, ml, tr, mr));
      }
      if (border) ss << wchar_t(0x2551);
      ss << "\n";
    }
    if (border) {
      ss << wchar_t(0x255A);
      for (ind_t x = 0; x < nx; x += 2) ss << wchar_t(0x2550);
      ss << wchar_t(0x255D) << "\n";
    }
    auto wstr = ss.str();
    return ws2s(wstr);
  }

}

#endif