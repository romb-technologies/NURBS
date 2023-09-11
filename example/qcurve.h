#ifndef QCURVE_H
#define QCURVE_H

#include <QGraphicsItem>

#include "NURBS/declarations.h"
#include "NURBS/nurbs.h"

class qCurve : public QGraphicsItem, public NURBS::Curve
{
private:
  bool draw_control_points = true;
  bool draw_curvature_radius = false;
  bool draw_knots = true;
  bool locked = false;

public:
  qCurve(const Eigen::MatrixX2d& points) : QGraphicsItem(), NURBS::Curve(points) {}
  qCurve(const NURBS::Curve& curve) : QGraphicsItem(), NURBS::Curve(curve) {}
  qCurve(NURBS::Curve&& curve) : QGraphicsItem(), NURBS::Curve(curve) {}

  int type() const Q_DECL_OVERRIDE;
  void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) Q_DECL_OVERRIDE;
  QRectF boundingRect() const Q_DECL_OVERRIDE;

  void prepareGeometryChange() { QGraphicsItem::prepareGeometryChange(); }
  void setDraw_control_points(bool value);
  void setDraw_curvature_radious(bool value);
  void setDraw_knots(bool value);
  bool getDraw_control_points() const;
  bool getDraw_curvature_radious() const;
  bool getDraw_knots() const;
  std::shared_ptr<NURBS::Curve> getSharedPtr();
  bool getLocked() const;
  void setLocked(bool value);
};

#endif // QCURVE_H
