#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <Eigen/Dense>
#include <vector>

namespace NURBS {

class Curve;

/*!
 * \brief Point in xy plane
 */
using Point = Eigen::Vector2d;

/*!
 * \brief A vector of Points
 */
using PointVector = std::vector<Point>;

/*!
 * \brief A Vector in xy plane
 */
using Vector = Eigen::Vector2d;

/*!
 * \brief Bounding box class
 */
using BoundingBox = Eigen::AlignedBox2d;

/*!
 * \brief Precision for numerical methods
 */
const double _epsilon = std::sqrt(std::numeric_limits<double>::epsilon());

}

#endif // DECLARATIONS_H
