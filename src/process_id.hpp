#ifndef BRILLE_PROCESS_ID_HPP_
#define BRILLE_PROCESS_ID_HPP_

#include <string>

namespace brille{
    int process_id();
    std::string pid_filename(std::string base, std::string ext);
}

#endif