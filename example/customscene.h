#ifndef CUSTOMSCENE_H
#define CUSTOMSCENE_H

#include <QGraphicsScene>

#include "qcurve.h"
#include "qgraphicsviewzoom.h"

class CustomScene : public QGraphicsScene
{
private:
  QGraphicsEllipseItem* dot;
  QGraphicsTextItem* number_display;
  QMap<QGraphicsItem*, QGraphicsLineItem*> line;
  QMap<QGraphicsItem*, QGraphicsLineItem*> tan;
  QMap<QGraphicsItem*, QGraphicsEllipseItem*> byLength;
  bool draw_box_ = true;
  bool draw_inter_ = false;
  bool show_projection = false;
  std::pair<qCurve*, double> t_to_update;
  bool update_cp = false;
  std::pair<QGraphicsItem*, uint> cp_to_update;
  bool update_weights = false;
  double current_weight;

protected:
  void keyPressEvent(QKeyEvent* keyEvent) Q_DECL_OVERRIDE;

  void drawForeground(QPainter* painter, const QRectF& rect) Q_DECL_OVERRIDE;
  void mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;
  void mouseMoveEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;
  void mouseReleaseEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;
};

#endif // CUSTOMSCENE_H
