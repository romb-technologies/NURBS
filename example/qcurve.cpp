#include "qcurve.h"

#include <QPainter>
#include <QPen>

void qCurve::setDraw_control_points(bool value) { draw_control_points = value; }

void qCurve::setDraw_curvature_radious(bool value) { draw_curvature_radious = value; }

bool qCurve::getDraw_control_points() const { return draw_control_points; }

bool qCurve::getDraw_curvature_radious() const { return draw_curvature_radious; }

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
      painter->setBrush(QBrush(Qt::green, Qt::SolidPattern));
      Eigen::VectorXd knots = knotVector().matrix();
      for (uint k=order(); k<controlPoints().size()+1; k++) {
          painter->setPen(Qt::green);
//          painter->drawEllipse(QRectF(valueAt(knots(k))(0) - dot_size / 2, valueAt(knots(k))(1) - dot_size / 2, dot_size, dot_size));
      }
  }
}


QRectF qCurve::boundingRect() const
{
  auto bbox = boundingBox();
  return QRectF(QPointF(bbox.min().x(), bbox.min().y()), QPointF(bbox.max().x(), bbox.max().y()));
}

