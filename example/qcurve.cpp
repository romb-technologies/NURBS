#include "qcurve.h"

#include <QPainter>
#include <QPen>

void qCurve::setDraw_control_points(bool value) { draw_control_points = value; }

void qCurve::setDraw_curvature_radious(bool value) { draw_curvature_radius = value; }

void qCurve::setDraw_knots(bool value) { draw_knots = value; }

bool qCurve::getDraw_control_points() const { return draw_control_points; }

bool qCurve::getDraw_curvature_radious() const { return draw_curvature_radius; }

bool qCurve::getLocked() const
{
    return locked;
}

void qCurve::setLocked(bool value)
{
    locked = value;
}

int qCurve::type() const { return QGraphicsItem::UserType + 1; }

void qCurve::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
{
  Q_UNUSED(option)
  Q_UNUSED(widget)

  setFlag(GraphicsItemFlag::ItemIsSelectable, true);

  painter->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform, true);

  QPen pen;
  pen.setStyle(isSelected() ? Qt::DashDotLine : Qt::SolidLine);
  pen.setColor(getLocked() ? Qt::red : Qt::black);
  painter->setPen(pen);
  QPainterPath curve;
  auto poly = polyline();
  curve.moveTo(poly[0].x(), poly[0].y());
  for (uint k = 1; k < poly.size(); k++)
    curve.lineTo(poly[k].x(), poly[k].y());
  painter->drawPath(curve);

  if (draw_control_points)
  {
    const int dot_size = 6;
    painter->setBrush(QBrush(Qt::blue, Qt::SolidPattern));
    NURBS::PointVector points = controlPoints();
    for (uint k = 1; k < points.size(); k++)
    {
      painter->setPen(Qt::blue);
      painter->drawEllipse(QRectF(points[k - 1].x() - dot_size / 2, points[k - 1].y() - dot_size / 2, dot_size, dot_size));
      painter->setPen(QPen(QBrush(Qt::gray), 1, Qt::DotLine));
      painter->drawLine(QLineF(points[k - 1].x(), points[k - 1].y(), points[k].x(), points[k].y()));
    }
    painter->setPen(Qt::blue);
    painter->drawEllipse(QRectF(points.back().x() - dot_size / 2, points.back().y() - dot_size / 2, dot_size, dot_size));
  }

  if (draw_knots) {
      const int dot_size = 4;
      painter->setBrush(QBrush(Qt::magenta, Qt::SolidPattern));
      Eigen::VectorXd knots = Curve::knotVector().matrix();
      for (uint k=order(); k<controlPoints().size()+1; k++) {
          painter->setPen(Qt::magenta);
          painter->drawEllipse(QRectF(valueAt(knots(k))(0) - dot_size / 2, valueAt(knots(k))(1) - dot_size / 2, dot_size, dot_size));
      }
  }

  if (draw_curvature_radius)
  {
    for (double t = 1.0 / 100; t <= 1.0; t += 1.0 / 200)
    {
      painter->setPen(QColor(abs(255 * (0.5 - t)), (int)(255 * t), (int)(255 * (1 - t))));
      auto p = valueAt(t);
      auto tangent = tangentAt(t);
      NURBS::Point normal(-tangent.y(), tangent.x());
      double kappa = curvatureAt(t);
      auto n1 = p + normal * kappa * 100;
      auto n2 = p - normal * kappa * 100;
      painter->drawLine(QLineF(n1.x(), n1.y(), n2.x(), n2.y()));
    }
  }
}


QRectF qCurve::boundingRect() const
{
  auto bbox = boundingBox();
  return QRectF(QPointF(bbox.min().x(), bbox.min().y()), QPointF(bbox.max().x(), bbox.max().y()));
}


std::vector<double> qCurve::knotVector()
{
  std::vector<double> out(Curve::knotVector().rows());
  Eigen::Map<Eigen::ArrayXd>(out.data(), out.size()) = Curve::knotVector();
  return out;
}

std::vector<double> qCurve::weights()
{
  std::vector<double> out(Curve::weights().rows());
  Eigen::Map<Eigen::VectorXd>(out.data(), out.size()) = Curve::weights();
  return out;
}

void qCurve::setKnot(int idx, double value)
{
  if (knot(idx) != value) {
      Curve::setKnot(idx, value);
      emit knotChanged(idx, knot(idx));
      emit curveChanged();
    }
}

void qCurve::setWeight(int idx, double value)
{
  if (weight(idx) != value) {
      Curve::setWeight(idx, value);
      emit weightChanged(idx, value);
      emit curveChanged();
    }
}
