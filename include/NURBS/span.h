#ifndef SPAN_H
#define SPAN_H

#include "declarations.h"
#include "utils.h"
#include <memory>

namespace NURBS {

class Span
{
public:
    ~Span() = default;
    Span(Eigen::Ref<Eigen::MatrixX3d> wpoints_,
         Eigen::Ref<Eigen::ArrayXd> knot_v,
         double start, double end, uint p);

    Eigen::Ref<Eigen::MatrixX3d> wpoints_;
    Eigen::Ref<Eigen::ArrayXd> knots;

    double start_t_, end_t_;
    Eigen::RowVectorXd cached_wbf_;
    Eigen::MatrixXd cached_vbf_;
    bool contains(double t) const;
    void update();
    void updateControlPoints();
    Eigen::MatrixXd getBasisFunction() const;
    PointVector polyline() const;
    Point valueAt(double u) const;
    Eigen::MatrixXd basis_function_;
    uint p_;
    /// Reset all privately cached data
    void resetCache();
private:
    mutable std::optional<PointVector> cached_polyline_;
};
}

#endif // SPAN_H
