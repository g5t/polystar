#include "process_id.hpp"

#if defined(_MSC_VER) || defined(__MINGW32__)
#include <process.h>
int polystar::process_id() {return _getpid();}
#else
#include <unistd.h>
int polystar::process_id() {return static_cast<int>(getpid());}
#endif

std::string polystar::pid_filename(std::string base, std::string ext){
  return base + std::to_string(polystar::process_id()) + ext;
}