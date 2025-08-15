#pragma once

#define LOG(msg) \
    std::cerr << "[LOG] " << msg << " (en " << __FILE__ << ", " << __func__ << ", línea " << __LINE__ << ")\n"
