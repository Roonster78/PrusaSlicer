#ifndef TRIANGULATEWALL_HPP
#define TRIANGULATEWALL_HPP

#include "libslic3r/Polygon.hpp"

namespace Slic3r {

std::vector<size_t> triangulate_wall(
    const Polygon &          lower,
    const Polygon &          upper,
    double                   lower_z_mm,
    double                   upper_z_mm,
    double                   offset_diff_mm,
    std::function<bool(int)> statusfn = [](int) { return false; });

}

#endif // TRIANGULATEWALL_HPP
