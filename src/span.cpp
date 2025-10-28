#include "NURBS/span.h"

using namespace NURBS;

///// Curve::Span

Span::Span(Eigen::Ref<Eigen::MatrixX3d> wpoints, Eigen::Ref<Eigen::ArrayXd> knot_v, uint p)
    : wpoints_(wpoints), knots_(knot_v), p_(p)
{
  update();
}

Eigen::Ref<Eigen::MatrixX3d> Span::wpoints() const { return wpoints_; }

Eigen::Ref<Eigen::ArrayXd> Span::knots() const { return knots_; }

Eigen::MatrixXd Span::basisFunction() const { return basis_function_; }

Eigen::RowVectorXd Span::cachedWBF() const { return cached_wbf_; }

Eigen::MatrixXd Span::cachedVBF() const { return cached_vbf_; }

bool Span::contains(double t) const { return end() - start() > _epsilon && t >= start() && t < end(); }

void Span::update()
{
  if (end() - start() <= _epsilon) // start_t_ == end_t_
  {
    basis_function_ = Eigen::MatrixXd::Zero(p_ + 1, p_ + 1);
    updateControlPoints();
    return;
  }

  // generate basis function
  basis_function_.resize(1, 1);
  basis_function_ << 1;

  for (int k = 2, i = p_ - 1; k <= p_ + 1; k++)
  {

    Eigen::MatrixXd m1(k, k - 1), m2(k - 1, k), m3(k, k - 1), m4(k - 1, k);
    m1 << basis_function_, Eigen::MatrixXd::Zero(1, k - 1);
    m3 << Eigen::MatrixXd::Zero(1, k - 1), basis_function_;
    m2.setZero(), m4.setZero();

    Eigen::ArrayXd d0(k - 1), d1(k - 1), ddwn(k - 1);
    d0.setConstant(start()), d1.setConstant(end() - start());

    ddwn = knots_.segment(i + 1, k - 1) - knots_.segment(i - k + 2, k - 1);
    d0 -= knots_.segment(i - k + 2, k - 1);

    d0 /= ddwn, d1 /= ddwn;
    m2.diagonal() = 1 - d0;
    m2.diagonal(1) = d0;
    m4.diagonal() = -d1;
    m4.diagonal(1) = d1;
    basis_function_ = (m1 * m2) + (m3 * m4);
  }

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
      cached_polyline_->emplace_back(valueAt(u));
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
  updateControlPoints();
}

Point Span::derivativeAt(int n, double u) const
{
  if (p_ < n)
    return Point(0, 0);

  // Derivatives of t power series
  std::vector<Eigen::RowVectorXd> dt;
  for (int i = 0; i <= n; i++)
    dt.emplace_back(_powSeriesDerivative(u, p_, i));

  const auto& V = cached_vbf_;
  const auto& W = cached_wbf_;

  // g = 1/W
  std::vector<double> dg;
  dg.emplace_back(1 / dt[0].dot(W));

  // Generalized derivative of 1/W
  for (int i = 1; i <= n; i++)
  {
    double dg_temp = 0.0;
    for (int j = 1; j <= i; j++)
      dg_temp += _binomial(i, j) * dt[j].dot(W) * dg[i - j];
    dg.emplace_back(-dg[0] * dg_temp);
  }

  // Leibniz product rule
  Point deriv{0, 0};
  for (int i = 0; i <= n; i++)
  {
    deriv += _binomial(n, i) * (dt[n - i] * V) * dg[i];
  }
  return deriv;
}

Point Span::derivativeAt(double u) const { return derivativeAt(1, u); }

double Span::length(double t) const
{
  if (t > 1.0 || t < 0.0)
    throw std::logic_error{"Length can only be calculated for t within [0.0, 1.0] range."};

  auto evaluateChebyshev = [](double t, const Eigen::VectorXd& coeff) {
    t = 2 * t - 1;
    double tn{t}, tn_1{1}, res{coeff(0) + coeff(1) * t};
    for (unsigned k = 2; k < coeff.size(); k++)
    {
      std::swap(tn_1, tn);
      tn = 2 * t * tn_1 - tn;
      res += coeff(k) * tn;
    }
    return res;
  };

  if (!cached_chebyshev_coeffs_)
  {
    constexpr unsigned START_LOG_N = 10;
    unsigned log_n = START_LOG_N - 1;
    // todo: počet s manje točaka (kasnije)
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
    } while (std::fabs(chebyshev.tail<1>()(0)) > _epsilon * 1e-2);

    unsigned cut = 0;
    while (std::fabs(chebyshev(cut)) > _epsilon * 1e-2)
      cut++;

    cached_chebyshev_coeffs_ = Eigen::VectorXd(cut + 1);
    *cached_chebyshev_coeffs_ << 0, chebyshev.head(cut);
    (*cached_chebyshev_coeffs_)(0) = -evaluateChebyshev(0, *cached_chebyshev_coeffs_);
  }
  return evaluateChebyshev(t, *cached_chebyshev_coeffs_);
}

double Span::length() const { return length(1.0); }

void Span::reassign(Eigen::Ref<Eigen::MatrixX3d> wpoints, Eigen::Ref<Eigen::ArrayXd> knots)
{
  new (&wpoints_) Eigen::Ref<Eigen::MatrixX3d>{wpoints};
  new (&knots_) Eigen::Ref<Eigen::VectorXd>{knots};
}
