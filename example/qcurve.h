#ifndef QCURVE_H
#define QCURVE_H

#include <QGraphicsObject>

#include "NURBS/declarations.h"
#include "NURBS/nurbs.h"

class qCurve : public QGraphicsObject, public NURBS::Curve
{
  Q_OBJECT
private:
  bool draw_control_points = true;
  bool draw_curvature_radius = false;
  bool draw_knots = true;
  bool locked = false;

public:
  qCurve(const Eigen::MatrixX2d& points) : QGraphicsObject(), NURBS::Curve(points) {}
  qCurve(const NURBS::Curve& curve) : QGraphicsObject(), NURBS::Curve(curve) {}
  qCurve(NURBS::Curve&& curve) : QGraphicsObject(), NURBS::Curve(curve) {}

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
  std::vector<double> knotVector();
  std::vector<double> weights();

public slots:
  void setKnot(int idx, double value);
  void setWeight(int idx, double value);
  void elevateOrder(uint t);

signals:
  void knotChanged(int idx, double value);
  void weightChanged(int idx, double value);
  void curveChanged();
};

#endif // QCURVE_H
