#include "NURBS/span.h"


using namespace NURBS;

///// Curve::Span

Span::Span(Eigen::Ref<Eigen::MatrixX3d> wpoints,
                  Eigen::Ref<Eigen::ArrayXd> knot_v, double start, double end, uint p) :
    wpoints(wpoints), knots(knot_v), start_t(start), end_t(end), p_(p)
{
    update();
}

bool Span::contains(double t) const
{
    if (start_t==end_t)
        return false;
    else
        return (t >= start_t) && (t < end_t);
}

Eigen::MatrixXd Span::getBasisFunction() const
{
    return basis_function_;
}

void Span::update()
{
    // update start and end
    start_t = knots(p_-1);
    end_t = knots(p_);

    // generate basis function
    Eigen::MatrixXd m(1, 1); m<<1;

    if (start_t != end_t) {
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
                    d0 = Eigen::ArrayXd::Constant(k-1, start_t),
                    d1 = Eigen::ArrayXd::Constant(k-1, end_t - start_t),
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
    Eigen::MatrixX3d test = wpoints;
    cached_v_bf = basis_function_ * wpoints.leftCols<2>();
    cached_w_bf = basis_function_ * wpoints.col(2);
}


