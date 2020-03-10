#ifndef TRIANGULATEWALL_HPP
#define TRIANGULATEWALL_HPP

#include <Polygon.hpp>
#include <MTUtils.hpp>

namespace Slic3r {

inline const Point &getpt(const Polygon &p, size_t idx) { return p.points[idx]; }
inline size_t num_points(const Polygon &p) { return p.points.size(); }

struct Tri {
    size_t idx = 0, nextidx = 0, endidx = 0;
    bool started = false;
};

struct LowerTri : public Tri {};
struct UpperTri : public Tri {};

// This will be the flip switch to toggle between upper and lower triangle
// creation mode
enum class Proceed {
    UPPER, // A segment from the upper polygon and one vertex from the lower
    LOWER  // A segment from the lower polygon and one vertex from the upper
};

class Triangulator {
    const coord_t offset_diff_sq;
    const coord_t zdiff_sq;
    
public:
    LowerTri lower; UpperTri upper;
    Proceed dir = Proceed::UPPER;
     
    coord_t score(const Polygon &lpoly, const Polygon &upoly)
    {
        coord_t a = 0, b = 0;
        switch (dir) {
        case Proceed::UPPER:
            a = sq_distance(getpt(upoly, upper.idx), getpt(lpoly, lower.idx));
            b = sq_distance(getpt(upoly, upper.nextidx), getpt(lpoly, lower.idx));
            break;
        case Proceed::LOWER:
            a = sq_distance(getpt(upoly, upper.idx), getpt(lpoly, lower.idx));
            b = sq_distance(getpt(upoly, upper.idx), getpt(lpoly, lower.nextidx));
            break;
        }
        
        a = offset_diff_sq - a + zdiff_sq;
        b = offset_diff_sq - b + zdiff_sq;
        
        return (std::abs(a) + std::abs(b)) / 2;
    }

    template<class StatusFn, class BackInserterIt>
    void emplace_indices(size_t         upper_offs,
                         size_t         lower_offs,
                         BackInserterIt it)
    {
        using VT = typename std::iterator_traits<BackInserterIt>::value_type;
        std::array<VT, 3> tri;
        
        switch (dir) {
        case Proceed::UPPER:
            tri = {VT(upper_offs + upper.nextidx), VT(lower_offs + lower.idx),
                   VT(upper_offs + upper.idx)};
            break;
        case Proceed::LOWER:
            tri = {VT{upper_offs + upper.idx}, VT{lower_offs + lower.nextidx},
                   VT{lower_offs + lower.idx}};
            break;
        }
        std::copy(tri.begin(), tri.end(), it);
    }

    explicit Triangulator(coord_t odiff, coord_t zdiff)
        : offset_diff_sq{odiff * odiff}, zdiff_sq{zdiff * zdiff}
    {}
    
//    template<class StatusFn, class BackInserterIt>
//    void run(size_t upper_offs, size_t lower_offs, BackInserterIt it, 
//             const Polygon &lpoly, const Polygon &upoly)
//    {
        
//    }
};

template<class StatusFn>
std::vector<size_t> triangulate_walls(const Polygon &lower,
                                      const Polygon &upper,
                                      double         lower_z_mm,
                                      double         upper_z_mm,
                                      double         offset_difference_mm,
                                      StatusFn &&    status)
{
    if(num_points(upper) < 3 || num_points(lower) < 3) return {};
    
    std::vector<size_t> ret;

    Triangulator tri{scaled(offset_difference_mm), scaled(upper_z_mm - lower_z_mm)};

    coord_t current_fit = 0, prev_fit = 0;
    
    // The offset for indices from the lower polygon
    size_t offs = num_points(upper);
    
    // Create pointing indices into vertex arrays.
    tri.upper.idx     = 0;
    tri.lower.idx     = offs;
    tri.upper.nextidx = 1;
    tri.lower.nextidx = offs + 1;

    // Mark the current vertex iterator positions. If the iterators return to
    // the same position, the loop can be terminated.
    tri.upper.endidx = tri.upper.idx;
    tri.lower.endidx = tri.lower.idx;
    
    // If the Z levels are flipped, or the offset difference is negative, we
    // will interpret that as the triangles normals should be inverted.
    bool inverted = upper_z_mm < lower_z_mm || offset_difference_mm < 0;
    
    // We need to find the closest point on lower polygon to the first point on
    // the upper polygon. These will be our starting points.
    double distmin = std::numeric_limits<double>::max();
    for(size_t l = tri.lower.idx; l < num_points(lower); ++l) {
//        thr();
//        double d = sq_distance(rpts[l], rpts[uidx]);
//        if(d < distmin) { lidx = l; distmin = d; }
    }
    
    return ret;
}

}

#endif // TRIANGULATEWALL_HPP
