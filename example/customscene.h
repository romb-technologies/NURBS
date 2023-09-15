#ifndef CUSTOMSCENE_H
#define CUSTOMSCENE_H

#include <QGraphicsScene>

#include "qcurve.h"
#include "qgraphicsviewzoom.h"

class CustomScene : public QGraphicsScene
{
  Q_OBJECT
public:
  bool draw_box_ = true;
  bool draw_inter_ = false;
  bool show_projection = false;
private:
  QGraphicsEllipseItem* dot;
  QGraphicsTextItem* number_display;
  QGraphicsLineItem* line;
  QGraphicsLineItem* tan;
  QGraphicsEllipseItem* byLength;

  std::pair<qCurve*, double> t_to_update;
  bool update_cp = false;
  std::pair<QGraphicsItem*, uint> cp_to_update;
  bool update_weights = false;
  double current_weight;

protected:
  void drawForeground(QPainter* painter, const QRectF& rect) Q_DECL_OVERRIDE;
  void mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;
  void mouseMoveEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;
  void mouseReleaseEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;

  void projectPointOntoCurve(qCurve* curve, NURBS::Point p);
  qCurve *getClosestCurve(NURBS::Point p, bool selected=false, double max_dist = 10);

signals:
  void cursorMove(QPointF pos);
  void pointToT(double t);

public slots:
  void selectItem(qCurve* c);
};

#endif // CUSTOMSCENE_H
