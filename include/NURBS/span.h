#ifndef SPAN_H
#define SPAN_H

#include "declarations.h"

namespace NURBS {

class Span
{
public:
    ~Span() = default;
    Span(Eigen::Ref<Eigen::MatrixX3d> wpoints,
         Eigen::Ref<Eigen::ArrayXd> knot_v,
         double start, double end, uint p);

    Eigen::Ref<Eigen::MatrixX3d> wpoints;
    Eigen::Ref<Eigen::ArrayXd> knots;

    double start_t, end_t;
    uint p_;
    Eigen::MatrixXd basis_function_;
    Eigen::RowVectorXd cached_w_bf;
    Eigen::MatrixXd cached_v_bf;
    bool contains(double t) const;
    void update();
    void updateControlPoints();
    Eigen::MatrixXd getBasisFunction() const;
};
}

#endif // SPAN_H
