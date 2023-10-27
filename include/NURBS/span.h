#ifndef SPAN_H
#define SPAN_H

#include "declarations.h"
#include "utils.h"
#include <memory>

#include <unsupported/Eigen/FFT>

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
    Point derivativeAt(int n, double u) const;
    Point derivativeAt(double u) const;
    Eigen::MatrixXd basis_function_;
    uint p_;
    /// Reset all privately cached data
    void resetCache();
    double length(double t) const;
    double length() const;


private:
    mutable std::optional<PointVector> cached_polyline_;
    mutable std::optional<Eigen::VectorXd> cached_chebyshev_coeffs_; /*!  If generated, stores chebyshev coefficients
                                                                          for calculating the length of the curve */
    mutable std::optional<double> cached_length_;

    static double evaluate_chebyshev(double t, const Eigen::VectorXd& coeff);
};
}

#endif // SPAN_H
