#include "NURBS/span.h"


using namespace NURBS;

///// Curve::Span

Span::Span(Eigen::Ref<Eigen::MatrixX3d> wpoints,
                  Eigen::Ref<Eigen::ArrayXd> knot_v, double start, double end, uint p) :
    wpoints_(wpoints), knots(knot_v), start_t_(start), end_t_(end), p_(p)
{
    update();
}

bool Span::contains(double t) const
{
    if (start_t_==end_t_)
        return false;
    else
        return (t >= start_t_) && (t < end_t_);
}

Eigen::MatrixXd Span::getBasisFunction() const
{
    return basis_function_;
}

void Span::update()
{
    // update start and end
    start_t_ = knots(p_-1);
    end_t_ = knots(p_);

    // generate basis function
    Eigen::MatrixXd m(1, 1); m<<1;

    if (start_t_ != end_t_) {
        static const int i = p_-1;
        for (int k=2; k<=p_+1; k++)
        {
            Eigen::MatrixXd m1(k, k-1),
                    m2 = Eigen::MatrixXd::Zero(k-1, k),
                    m3(k, k-1),
                    m4 = Eigen::MatrixXd::Zero(k-1, k);

            m1 << m, Eigen::MatrixXd::Zero(1, k-1);
            m3 << Eigen::MatrixXd::Zero(1, k-1), m;

            Eigen::ArrayXd
                    d0 = Eigen::ArrayXd::Constant(k-1, start_t_),
                    d1 = Eigen::ArrayXd::Constant(k-1, end_t_ - start_t_),
                    ddwn = Eigen::ArrayXd::Zero(k-1);

            ddwn = knots.segment(i+1, k-1) - knots.segment(i-k+2, k-1);
            d0 -= knots.segment(i-k+2, k-1);

            d0 /= ddwn;
            d1 /= ddwn;

            m2.diagonal() = 1 - d0;
            m2.diagonal(1) = d0;

            m4.diagonal() = -d1;
            m4.diagonal(1) = d1;

            m = (m1*m2)+(m3*m4);
        }
    } else {
        m = Eigen::MatrixXd::Zero(p_+1, p_+1);
    }

    basis_function_ = m;

    updateControlPoints();
}

void Span::updateControlPoints()
{
    // generate w_bf, v_bf
    Eigen::MatrixX3d test = wpoints_;
    cached_vbf_ = basis_function_ * wpoints_.leftCols<2>();
    cached_wbf_ = basis_function_ * wpoints_.col(2);
}

PointVector Span::polyline() const {
    if (!cached_polyline_)
    {
        cached_polyline_ = PointVector();
        for(double u = 0.0; u < 1.0 + 0.005; u+=0.02) {
            cached_polyline_->emplace_back(valueAt(u));
        }
    }
    return *cached_polyline_;
}

Point Span::valueAt(double u) const {
    Eigen::RowVectorXd pw = _powSeries(u, p_);
    return (pw * cached_vbf_) / pw.dot(cached_wbf_);
}

void Span::resetCache() {
    cached_polyline_.reset();
}
