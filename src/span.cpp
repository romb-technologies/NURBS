#include "NURBS/span.h"

using namespace NURBS;

///// Curve::Span

Span::Span(Eigen::Ref<Eigen::MatrixX3d> wpoints, Eigen::Ref<Eigen::ArrayXd> knot_v, double start, double end, uint p)
    : wpoints_(wpoints), knots(knot_v), start_t_(start), end_t_(end), p_(p)
{
  update();
}

bool Span::contains(double t) const
{
  if (start_t_ == end_t_)
    return false;
  else
    return (t >= start_t_) && (t < end_t_);
}

Eigen::MatrixXd Span::getBasisFunction() const { return basis_function_; }

void Span::update()
{
  // update start and end
  start_t_ = knots(p_ - 1);
  end_t_ = knots(p_);

  // generate basis function
  Eigen::MatrixXd m(1, 1);
  m << 1;

  if (start_t_ != end_t_)
  {
    int i = p_ - 1;
    for (int k = 2; k <= p_ + 1; k++)
    {
      Eigen::MatrixXd m1(k, k - 1), m2 = Eigen::MatrixXd::Zero(k - 1, k), m3(k, k - 1),
                                    m4 = Eigen::MatrixXd::Zero(k - 1, k);

      m1 << m, Eigen::MatrixXd::Zero(1, k - 1);
      m3 << Eigen::MatrixXd::Zero(1, k - 1), m;

      Eigen::ArrayXd d0 = Eigen::ArrayXd::Constant(k - 1, start_t_),
                     d1 = Eigen::ArrayXd::Constant(k - 1, end_t_ - start_t_), ddwn = Eigen::ArrayXd::Zero(k - 1);

      ddwn = knots.segment(i + 1, k - 1) - knots.segment(i - k + 2, k - 1);
      d0 -= knots.segment(i - k + 2, k - 1);

      d0 /= ddwn;
      d1 /= ddwn;

      m2.diagonal() = 1 - d0;
      m2.diagonal(1) = d0;

      m4.diagonal() = -d1;
      m4.diagonal(1) = d1;

      m = (m1 * m2) + (m3 * m4);
    }
  }
  else
  {
    m = Eigen::MatrixXd::Zero(p_ + 1, p_ + 1);
  }

  basis_function_ = m;

  updateControlPoints();
}

void Span::updateControlPoints()
{
  // generate w_bf, v_bf
  cached_vbf_ = basis_function_ * wpoints_.leftCols<2>();
  cached_wbf_ = basis_function_ * wpoints_.col(2);
}

PointVector Span::polyline() const
{
  if (!cached_polyline_)
  {
    cached_polyline_ = PointVector();
    for (double u = 0.0; u < 1.0 + 0.01; u += 0.02)
    {
      cached_polyline_->emplace_back(valueAt(u));
    }
  }
  return *cached_polyline_;
}

Point Span::valueAt(double u) const
{
  Eigen::RowVectorXd pw = _powSeries(u, p_);
  return (pw * cached_vbf_) / pw.dot(cached_wbf_);
}

void Span::resetCache()
{
  cached_polyline_.reset();
  start_t_ = knots(p_ - 1);
  end_t_ = knots(p_);
  updateControlPoints();
}

Point Span::derivativeAt(int n, double u) const
{
  // temporary solution

  Eigen::RowVectorXd pw = _powSeries(u, p_);
  Eigen::RowVectorXd pwd1 = _powSeriesDerivative(u, p_, 1);
  Eigen::RowVectorXd pwd2 = _powSeriesDerivative(u, p_, 2);
  Eigen::RowVectorXd pwd3 = _powSeriesDerivative(u, p_, 3);

  Eigen::MatrixX2d r = cached_vbf_;
  Eigen::VectorXd s = cached_wbf_;

  Eigen::RowVector2d ru = pw * r;
  double su = pw.dot(s);

  // derivatives of 1/S(u) in point t
  double d1su = -pwd1.dot(s) / pow(su, 2);

  double d2su = -pwd2.dot(s) / pow(su, 2) + 2 * pow(pwd1.dot(s), 2) / pow(su, 3);

  double d3su = -(pwd3.dot(s) / pow(su, 2)) + 4 * (pwd1.dot(s) * pwd2.dot(s) / pow(su, 3)) -
                6 * (pow(pwd1.dot(s), 3) / pow(su, 4));

  switch (n)
  {
  case 1:
    return ru * d1su + (pwd1 * r) / su;
  case 2:
    return ru * d2su + 2 * (pwd1 * r) * d1su + (pwd2 * r) / su;
  case 3:
    return ru * d3su + 3 * (pwd1 * r) * d2su + 3 * (pwd2 * r) * d1su + (pwd3 * r) / su;
  default:
    return valueAt(u);
  }
}

Point Span::derivativeAt(double u) const { return derivativeAt(1, u); }

double Span::length(double t) const
{
  if (t < 0.0 || t > 1.0)
    throw std::logic_error{"Length can only be calculated for t within [0.0, 1.0] range."};

  if (!cached_chebyshev_coeffs_)
  {
    constexpr unsigned START_LOG_N = 10;
    unsigned log_n = START_LOG_N - 1;
    unsigned n = _exp2(START_LOG_N - 1);

    Eigen::VectorXd derivative_cache(2);
    auto updateDerivativeCache = [this, &derivative_cache](double n) {
      //            auto start = std::chrono::steady_clock::now();
      derivative_cache.conservativeResize(n + 1);
      derivative_cache.segment(n / 2 + 1, n / 2) =
          ((1 + Eigen::cos(Eigen::ArrayXd::LinSpaced(n / 2, 1, n - 1) * M_PI / n)) / 2).unaryExpr([this](double t) {
            return derivativeAt(t).norm();
          });
    };

    derivative_cache.head(2) << derivativeAt(1.0).norm(), derivativeAt(0.0).norm();
    for (unsigned k = 2; k <= n; k *= 2)
      updateDerivativeCache(k);

    Eigen::VectorXd chebyshev;
    Eigen::FFT<double> fft;
    Eigen::VectorXcd fft_out;
    do
    {
      n *= 2;
      log_n++;
      updateDerivativeCache(n);

      unsigned N = 2 * n;
      Eigen::VectorXd coeff(N);
      coeff(0) = derivative_cache(0);
      coeff(n) = derivative_cache(1);

      for (unsigned k = 1; k <= log_n; k++)
      {
        auto lin_spaced = Eigen::ArrayXi::LinSpaced(_exp2(k - 1), 0, _exp2(k - 1) - 1);
        auto index_c = _exp2(log_n + 1 - (k + 1)) + lin_spaced * _exp2(log_n + 1 - k);
        auto index_dc = _exp2(k - 1) + 1 + lin_spaced;
        // TODO: make use of slicing & indexing in Eigen3.4
        // coeff(index_c) = coeff(N - index_c) = derivative_cache(index_dc) / n;
        for (unsigned i = 0; i < lin_spaced.size(); i++)
          coeff(index_c(i)) = coeff(N - index_c(i)) = derivative_cache(index_dc(i)) / n;
      }

      fft.fwd(fft_out, coeff);
      chebyshev = (fft_out.real().head(n - 1) - fft_out.real().segment(2, n - 1)).array() /
                  Eigen::ArrayXd::LinSpaced(n - 1, 4, 4 * (n - 1));
    } while (std::fabs(chebyshev.tail<1>()[0]) > _epsilon * 1e-2);

    unsigned cut = 0;
    while (std::fabs(chebyshev(cut)) > _epsilon * 1e-2)
      cut++;
    cached_chebyshev_coeffs_ = Eigen::VectorXd(cut + 1);
    *cached_chebyshev_coeffs_ << 0, chebyshev.head(cut);
    (*cached_chebyshev_coeffs_)(0) = -evaluate_chebyshev(0, *cached_chebyshev_coeffs_);
  }
  return evaluate_chebyshev(t, *cached_chebyshev_coeffs_);
}

double Span::length() const { return length(1.0); }

double Span::evaluate_chebyshev(double t, const Eigen::VectorXd& coeff)
{
  t = 2 * t - 1;
  double tn{t}, tn_1{1}, res{coeff(0) + coeff(1) * t};
  for (unsigned k = 2; k < coeff.size(); k++)
  {
    std::swap(tn_1, tn);
    tn = 2 * t * tn_1 - tn;
    res += coeff(k) * tn;
  }
  return res;
}

Eigen::RowVectorXd Span::getCachedWBF() const { return cached_wbf_; }

Eigen::MatrixXd Span::getCachedVBF() const { return cached_vbf_; }
