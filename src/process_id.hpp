#ifndef POLYSTAR_PROCESS_ID_HPP_
#define POLYSTAR_PROCESS_ID_HPP_

#include <string>

namespace polystar{
    int process_id();
    std::string pid_filename(std::string base, std::string ext);
}

#endif