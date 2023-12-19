#ifndef UTILS_H
#define UTILS_H

#include "declarations.h"

namespace NURBS
{

struct _PolynomialRoots : public std::vector<double>
{
  explicit _PolynomialRoots(unsigned reserve) { std::vector<double>::reserve(reserve); }
  void clear() {}          // no-op so that PolynomialSolver::RealRoots() doesn't clear it
  void push_back(double t) // only allow valid roots
  {
    if (t >= 0 && t <= 1)
      std::vector<double>::push_back(t);
  }
};

inline unsigned _exp2(unsigned exp) { return 1 << exp; }

inline double _pow(double base, unsigned exp)
{
  double result = exp & 1 ? base : 1;
  while (exp >>= 1)
    if (base *= base, exp & 1)
      result *= base;
  return result;
}

inline Eigen::RowVectorXd _powSeries(double base, unsigned exp)
{
  Eigen::RowVectorXd power_series(exp + 1);
  power_series(0) = 1;
  for (unsigned k = 1; k <= exp; k++)
    power_series(k) = power_series(k - 1) * base;
  return power_series;
}

inline Eigen::RowVectorXd _powSeriesDerivative(double base, unsigned exp, unsigned drv)
{
  Eigen::RowVectorXd power_series = Eigen::RowVectorXd::Ones(exp + 1);

  for (uint i = 0; i < drv; i++)
    for (uint j = 0; j <= exp; j++)
      power_series(j) *= j - i;

  for (unsigned k = drv; k <= exp; k++)
    power_series(k) *= _pow(base, k - drv);

  return power_series;
}

inline Eigen::VectorXd _trimZeroes(const Eigen::VectorXd& vec)
{
  auto idx = vec.size();
  while (idx && std::fabs(vec(idx - 1)) < _epsilon)
    --idx;
  return vec.head(idx);
}

inline Eigen::VectorXd _multiplyPolynomials(const Eigen::VectorXd& poly1, const Eigen::VectorXd& poly2)
{
  Eigen::VectorXd result = Eigen::VectorXd::Zero(poly1.size() + poly2.size() - 1);
  for (int i = 0; i < poly1.size(); i++)
    for (int j = 0; j < poly2.size(); j++)
      result(i + j) += poly1(i) * poly2(j);
  return result;
}

inline Eigen::MatrixX2d _pointVectorToMatrix (PointVector pv)
{
    Eigen::MatrixX2d out(pv.size(), 2);
    for (unsigned k = 0; k < pv.size(); k++)
        out.row(k) = pv[k];
    return out;
}

} // namespace NURBS

#endif
