#include "NURBS/utils.h"

namespace NURBS
{

unsigned _exp2(unsigned exp) { return 1 << exp; }

double _pow(double base, unsigned exp)
{
  double result = exp & 1 ? base : 1;
  while (exp >>= 1)
  {
    base *= base;
    if (exp & 1)
      result *= base;
  }
  return result;
}

Eigen::RowVectorXd _powSeries(double base, unsigned exp)
{
  Eigen::RowVectorXd power_series(exp + 1);
  power_series(0) = 1;
  for (unsigned k = 1; k <= exp; k++)
    power_series(k) = power_series(k - 1) * base;
  return power_series;
}

Eigen::RowVectorXd _powSeriesDerivative(double base, unsigned exp, unsigned drv)
{
  Eigen::RowVectorXd power_series = Eigen::RowVectorXd::Ones(exp + 1);
  for (uint i = 0; i < drv; i++)
  {
    for (uint j = 0; j <= exp; j++)
    {
      power_series(j) *= j - i;
    }
  }
  for (unsigned k = drv; k <= exp; k++)
    power_series(k) *= pow(base, k - drv);
  return power_series;
}

Eigen::VectorXd _trimZeroes(const Eigen::VectorXd& vec)
{
  auto idx = vec.size();
  while (idx && std::abs(vec(idx - 1)) < _epsilon)
    --idx;
  return vec.head(idx);
}

Eigen::VectorXd _multiplyPolynomials(const Eigen::VectorXd& poly1, const Eigen::VectorXd& poly2)
{
  Eigen::VectorXd result = Eigen::VectorXd::Zero(poly1.size() + poly2.size() - 1);
  for (int i = 0; i < poly1.size(); i++)
    for (int j = 0; j < poly2.size(); j++)
      result[i + j] += poly1[i] * poly2[j];

  return result;
}

double _dist(Point first, Point second)
{
  return std::sqrt(std::pow(first.x() - second.x(), 2) + std::pow(first.y() - second.y(), 2));
}

double _distSquared(Point first, Point second)
{
  return std::pow(first.x() - second.x(), 2) + std::pow(first.y() - second.y(), 2);
}

} // namespace NURBS
