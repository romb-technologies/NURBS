#ifndef UTILS_H
#define UTILS_H

#include "declarations.h"

namespace NURBS
{

unsigned _exp2(unsigned exp);
double _pow(double base, unsigned exp);
Eigen::RowVectorXd _powSeries(double base, unsigned exp);
Eigen::RowVectorXd _powSeriesDerivative(double base, unsigned exp, unsigned drv);
Eigen::VectorXd _trimZeroes(const Eigen::VectorXd& vec);
Eigen::VectorXd _multiplyPolynomials(const Eigen::VectorXd& poly1, const Eigen::VectorXd& poly2);
double _dist(Point first, Point second);
double _distSquared(Point first, Point second);
} // namespace NURBS

#endif
