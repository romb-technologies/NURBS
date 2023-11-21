#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <Eigen/Dense>
#include <vector>

/*!
 * Nominal namespace containing class pre-definitions and typedefs
 */
namespace NURBS
{
/*!
 * \brief A NURBS curve class
 *
 * A class for storing and using any-order NURBS curve.
 * It uses private and static caching for storing often accessed data.
 * Private caching is used for data concerning individual curves, while
 * static caching is used for common data (coefficient matrices)
 */
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

} // namespace NURBS

#endif // DECLARATIONS_H
