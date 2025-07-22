#include "NURBS/span.h"

using namespace NURBS;

///// Curve::Span

Span::Span(Eigen::Ref<Eigen::MatrixX3d> wpoints, Eigen::Ref<Eigen::ArrayXd> knot_v, uint p) : p_(p)
{
  update(knot_v, wpoints);
}

Span::Span(Eigen::Ref<const Eigen::MatrixXd> basis_func, Eigen::Ref<const Eigen::MatrixXd> vbf,
           Eigen::Ref<const Eigen::RowVectorXd> wbf, double start, double end)
    : basis_function_(basis_func), cached_wbf_(wbf), cached_vbf_(vbf), p_(basis_func.rows() - 1), start_(start),
      end_(end)
{
}

Eigen::MatrixXd Span::basisFunction() const { return basis_function_; }

Eigen::RowVectorXd Span::cachedWBF() const { return cached_wbf_; }

Eigen::MatrixXd Span::cachedVBF() const { return cached_vbf_; }

bool Span::contains(double t) const { return std::fabs(start_ - end_) > _epsilon && t >= start() && t < end(); }

void Span::update(Eigen::Ref<Eigen::ArrayXd> knots, Eigen::Ref<Eigen::MatrixX3d> wpoints)
{
  start_ = knots(p_ - 1);
  end_ = knots(p_);

  // generate basis function
  Eigen::MatrixXd m(1, 1);
  m << 1;

  if (std::fabs(start() - end()) > _epsilon) // start_t_ != end_t_
  {
    int i = p_ - 1;
    for (int k = 2; k <= p_ + 1; k++)
    {
      Eigen::MatrixXd m1(k, k - 1), m2 = Eigen::MatrixXd::Zero(k - 1, k), m3(k, k - 1),
                                    m4 = Eigen::MatrixXd::Zero(k - 1, k);

      m1 << m, Eigen::MatrixXd::Zero(1, k - 1);
      m3 << Eigen::MatrixXd::Zero(1, k - 1), m;

      Eigen::ArrayXd d0 = Eigen::ArrayXd::Constant(k - 1, start()),
                     d1 = Eigen::ArrayXd::Constant(k - 1, end() - start()), ddwn = Eigen::ArrayXd::Zero(k - 1);

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

  updateControlPoints(wpoints);
}

void Span::updateControlPoints(Eigen::Ref<Eigen::MatrixX3d> wpoints)
{
  // generate w_bf, v_bf
  cached_vbf_ = basis_function_ * wpoints.leftCols<2>();
  cached_wbf_ = basis_function_ * wpoints.col(2);
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
  cached_length_.reset();
  cached_chebyshev_coeffs_.reset();
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
  double d1su = -pwd1.dot(s) / _pow(su, 2);

  double d2su = -pwd2.dot(s) / _pow(su, 2) + 2 * _pow(pwd1.dot(s), 2) / _pow(su, 3);

  double d3su = -(pwd3.dot(s) / _pow(su, 2)) + 4 * (pwd1.dot(s) * pwd2.dot(s) / _pow(su, 3)) -
                6 * (_pow(pwd1.dot(s), 3) / _pow(su, 4));

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

  auto evaluate_chebyshev = [](double t, const Eigen::VectorXd& coeff) {
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
    (*cached_chebyshev_coeffs_)(0) = -evaluate_chebyshev(0, *cached_chebyshev_coeffs_);
  }
  return evaluate_chebyshev(t, *cached_chebyshev_coeffs_);
}

double Span::length() const { return length(1.0); }

std::vector<double> Span::extrema() const
{
  std::vector<double> extr;
  Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver;

  // d/du R(u)
  Eigen::MatrixX2d p1 = Eigen::MatrixXd::Zero(p_ + 1, 2);
  p1.topRows(p_) = (cachedVBF().array().colwise() * _powSeriesDerivative(1, p_, 1).transpose().array()).bottomRows(p_);

  // d/du S(u)
  Eigen::RowVectorXd pb = Eigen::VectorXd::Zero(p_ + 1);
  pb.head(p_) = (cachedWBF().array() * _powSeriesDerivative(1, p_, 1).array()).tail(p_);

  Eigen::MatrixX2d poly(2 * p_ + 1, 2);
  poly.col(0) = -_multiplyPolynomials(pb, cachedVBF().col(0)) + _multiplyPolynomials(p1.col(0), cachedWBF());
  poly.col(1) = -_multiplyPolynomials(pb, cachedVBF().col(1)) + _multiplyPolynomials(p1.col(1), cachedWBF());

  auto trimmed_x = _trimZeroes(poly.col(0));
  auto trimmed_y = _trimZeroes(poly.col(1));

  _PolynomialRoots roots(trimmed_x.size() + trimmed_y.size());
  if (trimmed_x.size() > 1)
  {
    poly_solver.compute(trimmed_x);
    poly_solver.realRoots(roots);
  }
  if (trimmed_y.size() > 1)
  {
    poly_solver.compute(trimmed_y);
    poly_solver.realRoots(roots);
  }
  for (int j = 0; j < trimmed_x.size() + trimmed_y.size(); j++)
    if (roots[j] >= 0.0 && roots[j] <= 1.0)
      extr.emplace_back(roots[j]);

  return extr;
}

BoundingBox Span::boundingBox() const
{
  if (!cached_bounding_box_)
  {
    auto ex = extrema();
    Eigen::MatrixXd extremes(ex.size() + 2, 2);
    for (unsigned k = 0; k < ex.size(); k++)
      extremes.row(k) = valueAt(ex[k]);

    extremes.row(extremes.rows() - 1) = valueAt(0.0);
    extremes.row(extremes.rows() - 2) = valueAt(1.0);

    cached_bounding_box_ = BoundingBox(Point(extremes.col(0).minCoeff(), extremes.col(1).minCoeff()),
                                       Point(extremes.col(0).maxCoeff(), extremes.col(1).maxCoeff()));
  }
  return cached_bounding_box_.value();
}

BoundingBox Span::fastBoundingBox(const Eigen::MatrixXd& vbf, const Eigen::RowVectorXd& wbf,
                                  const Eigen::MatrixXd& inverse_basis_function)
{
  Eigen::MatrixXd control_points_implicit = inverse_basis_function * vbf;
  Eigen::VectorXd weights_implicit = wbf * inverse_basis_function.transpose();
  control_points_implicit = control_points_implicit.array().colwise() / weights_implicit.array();

  return BoundingBox(Point(control_points_implicit.col(0).minCoeff(), control_points_implicit.col(1).minCoeff()),
                     Point(control_points_implicit.col(0).maxCoeff(), control_points_implicit.col(1).maxCoeff()));
}

PointVector Span::intersections(const Span& other) const
{
  PointVector intersections;
  const double epsilon_length = _epsilon * (length() + other.length()) / 2;

  if (!boundingBox().intersects(other.boundingBox()))
  {
    return intersections;
  }

  const unsigned max_intersections = p_ * other.p_;
  intersections.reserve(max_intersections);

  Point start_this = valueAt(0.0), start_other = other.valueAt(0.0);

  std::vector<SplitPair_> subcurve_pairs;

  // Splitting matrices based on the curves' degrees
  auto splittingCoeffs = [](unsigned p) {
    Eigen::MatrixXd zL = Eigen::MatrixXd::Zero(p + 1, p + 1);
    Eigen::MatrixXd zR = Eigen::MatrixXd::Zero(p + 1, p + 1);
    zL.diagonal() = _powSeries(0.5, p);
    zR(0, 0) = 1;
    for (int i = 1; i < p + 1; i++)
    {
      zR.col(i) = zR.col(i - 1);
      zR.col(i).tail(p) += zR.col(i - 1).head(p);
    }
    for (int i = 0; i < p + 1; i++)
    {
      zR.diagonal(i) *= _pow(0.5, i);
      zR.row(i) *= _pow(0.5, i);
    }
    return std::pair<Eigen::MatrixXd, Eigen::MatrixXd>{zL, zR};
  };

  auto [zL_a, zR_a] = splittingCoeffs(p_);
  auto [zL_b, zR_b] = splittingCoeffs(other.p_);

  // Self-intersections
  if (this == &other)
  {
    auto pair_init = SplitPair_(zL_a * cachedVBF(), zL_a * cachedWBF().transpose(), zR_a * cachedVBF(),
                                zR_a * cachedWBF().transpose());
    subcurve_pairs.emplace_back(pair_init);
    start_other = other.valueAt(0.5);
  }
  else
    subcurve_pairs.emplace_back(cachedVBF(), cachedWBF(), other.cachedVBF(), other.cachedWBF());

  // Check if intersection already exists, if not then add it
  auto addIntersection = [&intersections, start_this, start_other, epsilon_length](const Point& new_point) {
    if (std::none_of(
            intersections.begin(), intersections.end(),
            [&new_point, epsilon_length](const Point& point) { return (point - new_point).norm() < epsilon_length; }) &&
        ((new_point - start_this).norm() >= epsilon_length) && ((new_point - start_other).norm() >= epsilon_length))
    {
      intersections.push_back(new_point);
    }
  };

  // Cached basis function inverse for bounding box checks
  Eigen::MatrixXd inv_a = basisFunction().inverse();
  Eigen::MatrixXd inv_b = other.basisFunction().inverse();

  while (!subcurve_pairs.empty() && intersections.size() <= max_intersections)
  {
    SplitPair_ pair = std::move(subcurve_pairs.back());
    subcurve_pairs.pop_back();

    BoundingBox bbox1 = fastBoundingBox(pair.vbf_a, pair.wbf_a, inv_a);
    BoundingBox bbox2 = fastBoundingBox(pair.vbf_b, pair.wbf_b, inv_b);

    if (!bbox1.intersects(bbox2))
      continue;
    else if (bbox1.diagonal().norm() < _epsilon)
      addIntersection(bbox1.center());
    else if (bbox2.diagonal().norm() < _epsilon)
      addIntersection(bbox2.center());
    else
    {
      auto pair_a = SplitPair_(zL_a * pair.vbf_a, zL_a * pair.wbf_a.transpose(), zR_a * pair.vbf_a,
                               zR_a * pair.wbf_a.transpose());
      auto pair_b = SplitPair_(zL_b * pair.vbf_b, zL_b * pair.wbf_b.transpose(), zR_b * pair.vbf_b,
                               zR_b * pair.wbf_b.transpose());

      subcurve_pairs.emplace_back(pair_a.vbf_a, pair_a.wbf_a, pair_b.vbf_a, pair_b.wbf_a);
      subcurve_pairs.emplace_back(pair_a.vbf_a, pair_a.wbf_a, pair_b.vbf_b, pair_b.wbf_b);
      subcurve_pairs.emplace_back(pair_a.vbf_b, pair_a.wbf_b, pair_b.vbf_a, pair_b.wbf_a);
      subcurve_pairs.emplace_back(pair_a.vbf_b, pair_a.wbf_b, pair_b.vbf_b, pair_b.wbf_b);
    }
  }

  return intersections;
}
