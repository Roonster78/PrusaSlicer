#ifndef TRIANGULATEWALL_HPP
#define TRIANGULATEWALL_HPP

#include <type_traits>
#include <array>
#include <vector>

namespace Slic3r {

namespace __trw {

template<class T> struct remove_cvref
{
    using type = std::remove_cv_t<std::remove_reference_t<T>>;
};

template<class T> using remove_cvref_t = typename remove_cvref<T>::type;

template<class Pt> struct _PointTraits {
    using CoordType = typename Pt::CoordType;
    static const constexpr size_t DIM = Pt::DIM;
    template<size_t dim> static const CoordType get(const Pt &p);
};

template<class T> using PointTraits = _PointTraits<remove_cvref_t<T>>;
template<class T> using TCoord = typename PointTraits<T>::CoordType;

template<class PSet> struct _PointSetTraits {
    using PointType   = typename PSet::PointType;
    using PointType3D = typename PSet::PointType3D;

    static const PointType &get(const PSet &p, size_t idx) { return p.points[idx]; }
    static size_t      size(const PSet &p) { return p.points.size(); }
    static PointType3D to_3d(const PointType &p, TCoord<PointType3D> z)
    {
        return {PointTraits<PointType>::get<0>(p), PointTraits<PointType>::get<1>(p), z};
    }
};

template<class T> using PointSetTraits = _PointSetTraits<remove_cvref_t<T>>;
template<class PSet> using TPoint = typename PointSetTraits<PSet>::PointType;

template<class PSet> const TPoint<PSet> &getpt(const PSet &p, size_t idx)
{
    return PointSetTraits<PSet>::get(p, idx);
}

template<class PSet> size_t num_points(const PSet &p)
{
    return PointSetTraits<PSet>::size(p);
}

struct Ring {
    size_t idx = 0, nextidx = 1, ringsize = 0;
    explicit Ring(size_t size): ringsize{size - 1} {}
    void inc() {
        if (nextidx > 0) nextidx++;
        if (nextidx == ringsize) nextidx = 0;
        idx ++;
        if (idx == ringsize) idx = 0;
    }
    bool finised() const { return nextidx == idx; }
};

// This will be the flip switch to toggle between upper and lower triangle
// creation mode
enum class Proceed {
    UPPER, // A segment from the upper polygon and one vertex from the lower
    LOWER  // A segment from the lower polygon and one vertex from the upper
};


// The concept of the algorithm is relatively simple. It will try to find
// the closest vertices from the upper and the lower polygon and use those
// as starting points. Then it will create the triangles sequentially using
// an edge from the upper polygon and a vertex from the lower or vice versa,
// depending on the resulting triangle's quality.
// The quality is measured by a scalar value. So far it looks like it is
// enough to derive it from the slope of the triangle's two edges connecting
// the upper and the lower part. A reference slope is calculated from the
// height and the offset difference.
template<class PSet> class Triangulator {
    const TCoord<PSet> offset_diff_sq;
    const TCoord<PSet> zdiff_sq;
    
public:
    Ring lower, upper;
    Proceed dir = Proceed::UPPER;
    
    void turn() { dir = dir == Proceed::UPPER ? Proceed::LOWER : Proceed::UPPER; }
    void proceed() { dir == Proceed::UPPER ? upper.inc() : lower.inc(); }
    bool finished() const { return lower.finised() && upper.finised(); }
    bool dir_finished() const { return dir == Proceed::UPPER ? upper.finised() : lower.finised(); }
    
    // Every triangle of the wall has two edges connecting the upper plate with
    // the lower plate. From the length of these two edges and the zdiff we
    // can calculate the momentary squared offset distance at a particular
    // position on the wall. The average of the differences from the reference
    // (squared) offset distance will give us the driving fitness value.
    //
    // This function measures the quality of the triangle derived from an upper
    // vertex and two lower vertices (or lower vertex and two upper vertices)
    // depending on the current 'dir'
    TCoord<PSet> score(const PSet &lpoly, const PSet &upoly)
    {
        TCoord<PSet> a = 0, b = 0;
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
    
    // Store the triangle derivable from the current ring positions and dir
    template<class BackInserterIt>
    void emplace_indices(BackInserterIt it)
    {
        using VT = typename std::iterator_traits<BackInserterIt>::value_type;
        std::array<VT, 3> tri;
        
        switch (dir) {
        case Proceed::UPPER:
            tri = {VT(upper.nextidx), VT(upper.ringsize + lower.idx), VT(upper.idx)};
            break;
        case Proceed::LOWER:
            tri = {VT{upper.idx}, VT{upper.ringsize + lower.nextidx}, VT{lower.idx}};
            break;
        }
        
        std::copy(tri.begin(), tri.end(), it);
    }

    explicit Triangulator(size_t  num_lower,
                          size_t  num_upper,
                          TCoord<PSet> odiff,
                          TCoord<PSet> zdiff)
        : offset_diff_sq{odiff * odiff}
        , zdiff_sq{zdiff * zdiff}
        , lower{num_lower}
        , upper{num_upper}
    {}
};

} // namespace __trw

// This function will return a triangulation of a sheet connecting an upper
// and a lower plate given as input polygons. It will not triangulate the
// plates themselves only the sheet. The caller has to specify the lower and
// upper z levels in world coordinates as well as the offset difference
// between the sheets. If the lower_z_mm is higher than upper_z_mm or the
// offset difference is negative, the resulting triangle orientation will be
// reversed.
//
// IMPORTANT: This is not a universal triangulation algorithm. It assumes
// that the lower and upper polygons are offsetted versions of the same
// original polygon. In general, it assumes that one of the polygons is
// completely inside the other. The offset difference is the reference
// distance from the inner polygon's perimeter to the outer polygon's
// perimeter. The real distance will be variable as the clipper offset has
// different strategies (rounding, etc...). This algorithm should have
// O(2n + 3m) complexity where n is the number of upper vertices and m is the
// number of lower vertices.
template<class PSet, class StatusFn>
std::vector<size_t> triangulate_walls(const PSet &        lower,
                                      const PSet &        upper,
                                      __trw::TCoord<PSet> lower_z,
                                      __trw::TCoord<PSet> upper_z,
                                      __trw::TCoord<PSet> offset_difference,
                                      StatusFn &&         status)
{
    using namespace __trw;
    using Coord = TCoord<PSet>;
    
    if(num_points(upper) < 3 || num_points(lower) < 3) return {};
    
    std::vector<size_t> ret;
    
    Triangulator t{offset_difference, upper_z - lower_z};
    
    if (status(std::max(t.lower.idx, t.upper.idx))) return {};
    
    // We need to find the closest point on lower polygon to the first point on
    // the upper polygon. These will be our starting points.
    auto distmin = std::numeric_limits<Coord>::max();
    auto refpt = PointSetTraits<PSet>::to_3d(getpt(upper, 0), upper_z);    
    for(size_t l = t.lower.idx; l < num_points(lower); ++l) {
        auto lowerpt = PointSetTraits<PSet>::to_3d(getpt(lower, l), lower_z);
        Coord d = sq_distance(lowerpt, refpt);
        if(d < distmin) { t.lower.idx = l; distmin = d; }
    }
    
    Coord current_fit = 0, prev_fit = 0;    
    while (!t.finished()) {
        prev_fit = current_fit;
        
        if (t.dir_finished() || (current_fit = t.score(lower, upper)) > prev_fit)
            t.turn();
        else {
            t.emplace_indices(std::back_inserter(ret));
            t.proceed();
        }
    }
    
    // If the Z levels are flipped, or the offset difference is negative, we
    // will interpret that as the triangles normals should be inverted.
    if (upper_z < lower_z || offset_difference < 0) {
        size_t num_triangles = ret.size() / 3;
        for (size_t i = 0; i < num_triangles; ++i)
            std::swap(ret[3*i], ret[3*i + 2]);
    }
    
    return ret;
}

} // namespace Slic3r

#endif // TRIANGULATEWALL_HPP
