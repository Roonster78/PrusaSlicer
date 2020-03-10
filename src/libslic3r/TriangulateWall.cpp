#include "TriangulateWall.hpp"
#include "TriangulateWallImpl.hpp"
#include "MTUtils.hpp"

namespace Slic3r {

namespace __trw {

template<> struct _PointTraits<Point> {
    using CoordType = coord_t;
    static const constexpr size_t DIM = 2;
    template<size_t dim> static coord_t get(const Point &p) { return p(dim); }
};

template<> struct _PointTraits<Vec3d> {
    using CoordType = double;
    static const constexpr size_t DIM = 3;
    template<size_t dim> static double get(const Vec3d &p) { return p(dim); }
};

template<> struct _PointSetTraits<Polygon> {
    using PointType   = Point;
    using PointType3D = Vec3d;
    
    static const Point &get(const Polygon &p, size_t idx) { return p.points[idx]; }
    static size_t      size(const Polygon &p) { return p.points.size(); }
    static Vec3d to_3d(const Point &p, double z) { return Slic3r::to_3d(unscaled(p), z); }
};

} // namespace __trw

std::vector<size_t> triangulate_wall(
    const Polygon &          lower,
    const Polygon &          upper,
    double                   lower_z_mm,
    double                   upper_z_mm,
    double                   offset_diff_mm,
    std::function<bool(int)> statusfn)
{
    return _triangulate_walls(lower, upper, scaled(lower_z_mm),
                              scaled(upper_z_mm), scaled(offset_diff_mm),
                              statusfn);
}

}
