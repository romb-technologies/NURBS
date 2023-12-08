#include "customscene.h"

#include <QGraphicsSceneMouseEvent>
#include <QKeyEvent>
#include <QMessageBox>

#define is_curve (curve->type() == QGraphicsItem::UserType + 1)
#define is_curve_item(x) (x->type() == QGraphicsItem::UserType + 1)

#define c_curve (static_cast<qCurve*>(curve))
#define c_curve_item(x) (static_cast<qCurve*>(x))

void CustomScene::drawForeground(QPainter* painter, const QRectF& rect)
{
  Q_UNUSED(rect)

  if (draw_box_)
  {
    QPen pen;
    pen.setStyle(Qt::DotLine);
    pen.setColor(Qt::blue);
    painter->setPen(pen);
    for (auto&& curve : selectedItems())
    {
      NURBS::BoundingBox bbox;
      if (is_curve)
        bbox = c_curve->boundingBox();
      painter->drawRect(bbox.min().x(), bbox.min().y(), bbox.max().x() - bbox.min().x(),
                        bbox.max().y() - bbox.min().y());
    }
  }
  if (draw_inter_)
  {
    painter->setPen(Qt::red);
    painter->setBrush(QBrush(Qt::red, Qt::SolidPattern));
    for (int k = 0; k < items().size(); k++)
      for (int i = k; i < items().size(); i++)
      {
        if (is_curve_item(items()[k]) && is_curve_item(items()[i]))
        {
          NURBS::PointVector inter;
          inter = static_cast<qCurve*>(items()[i])->intersections(*static_cast<qCurve*>(items()[k]));

          for (auto& dot : inter)
            painter->drawEllipse(QPointF(dot.x(), dot.y()), 3, 3);
        }
      }
  }
}

void CustomScene::mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  const int sensitivity = 5;
  NURBS::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());

  if (show_projection || update_cp || update_weights)
    return;

  // point projection
  if (mouseEvent->button() == Qt::RightButton)
  {
    for (auto&& curve : selectedItems())
    {
      if (is_curve)
      {
        dot = addEllipse(QRectF(QPointF(p.x(), p.y()), QSizeF(6, 6)), QPen(Qt::yellow),
                         QBrush(Qt::red, Qt::SolidPattern));
        line = addLine(0, 0, 0, 0, QPen(Qt::red));
        tan = addLine(0, 0, 0, 0, QPen(Qt::blue));
        projectPointOntoCurve(c_curve, p);
      }
    }
    show_projection = true;
  }
  if (mouseEvent->button() == Qt::LeftButton)
  {
    for (auto&& curve : items())
    {
      if (is_curve)
      {
        auto pv = c_curve->controlPoints();
        for (uint k = 0; k < pv.size(); k++)
          // if control point is pressed
          if ((pv[k] - p).norm() < sensitivity && c_curve->getDraw_control_points())
          {
            selectItem(c_curve);

            number_display = addText(QString("(%1, %2)").arg(QString::number(p[0]), QString::number(p[1])));
            number_display->setPos(p[0], p[1]);
            cp_to_update = std::make_pair(curve, k);

            update_cp = true;
            break;
          }
      }
      if (update_cp)
        break;
    }
    if (!update_cp)
    {
      if (qCurve* c = getClosestCurve(p))
        selectItem(c);
      else
      {
        clearSelection();
        update();
      }
    }
  }
  if (mouseEvent->button() == Qt::MiddleButton)
  {
    for (auto&& curve : items())
    {
      if (!is_curve)
        continue;

      else
      {
        if (is_curve)
        {
          auto pv = c_curve->controlPoints();
          for (uint k = 0; k < pv.size(); k++)
            if ((pv[k] - p).norm() < sensitivity && c_curve->getDraw_control_points())
            {
              number_display = addText(QString::number(c_curve->weight(k)));
              number_display->setPos(p[0], p[1]);

              update_weights = true;
              current_weight = c_curve->weight(k);
              cp_to_update = std::make_pair(curve, k);
              break;
            }
        }
        if (update_weights)
          break;
      }
    }
  }
}

void CustomScene::mouseReleaseEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  if (mouseEvent->button() == Qt::RightButton)
  {
    if (show_projection)
    {
      for (auto&& curve : selectedItems())
        if (is_curve)
        {
          removeItem(dot);
          removeItem(line);
          removeItem(tan);
        }
      show_projection = false;
    }
  }
  if (mouseEvent->button() == Qt::LeftButton)
  {
    if (update_cp)
      removeItem(number_display);
    update_cp = false;
  }
  if (mouseEvent->button() == Qt::MiddleButton)
  {
    if (update_weights)
      removeItem(number_display);
    update_weights = false;
  }
}

void CustomScene::mouseMoveEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  NURBS::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());

  if (show_projection)
  {
    for (auto&& curve : selectedItems())
    {
      if (is_curve)
      {
        projectPointOntoCurve(c_curve, p);
      }
    }
  }

  // move control points
  else if (update_cp)
  {
    auto curve = cp_to_update.first;
    if (is_curve)
    {
      c_curve->prepareGeometryChange();
      c_curve->setControlPoint(cp_to_update.second, p);
    }
    number_display->setPos(p[0], p[1]);
    number_display->setPlainText(QString("(%1, %2)").arg(QString::number(p[0]), QString::number(p[1])));
    emit cursorMove(mouseEvent->scenePos());
    update();
  }
  else if (update_weights)
  {
    auto curve = cp_to_update.first;
    if (is_curve)
    {
      NURBS::Point pos = c_curve->controlPoint(cp_to_update.second);
      c_curve->prepareGeometryChange();
      c_curve->setWeight(cp_to_update.second, current_weight + (p[0] - pos[0]) * 0.005);
    }
    number_display->setPlainText(QString::number(c_curve->weight(cp_to_update.second)));
  }
  else
  {
    emit cursorMove(mouseEvent->scenePos());
  }
  QGraphicsScene::mouseMoveEvent(mouseEvent);
}

void CustomScene::projectPointOntoCurve(qCurve* curve, NURBS::Point p)
{
  dot->setRect(QRectF(QPointF(p.x() - 3, p.y() - 3), QSizeF(6, 6)));
  double t1 = curve->projectPoint(p);
  auto p1 = curve->valueAt(t1);
  auto tan1 = curve->tangentAt(t1);
  line->setLine(QLineF(QPointF(p.x(), p.y()), QPointF(p1.x(), p1.y())));
  tan->setLine(QLineF(QPointF(p1.x(), p1.y()) - 500 * QPointF(tan1.x(), tan1.y()),
                      QPointF(p1.x(), p1.y()) + 500 * QPointF(tan1.x(), tan1.y())));
  emit pointToT(t1);
}

qCurve* CustomScene::getClosestCurve(NURBS::Point p, bool selected, double max_dist)
{
  double dist = max_dist;
  qCurve* out = nullptr;
  for (auto&& curve : selected ? selectedItems() : items())
  {
    if (is_curve)
    {
      double t = c_curve->projectPoint(p);
      auto pt = c_curve->valueAt(t);
      if (double d = (pt - p).norm() < dist)
      {
        dist = d;
        out = c_curve;
      }
    }
  }
  return out;
}

void CustomScene::selectItem(qCurve* c)
{
  blockSignals(true);
  clearSelection();
  blockSignals(false);
  c->setSelected(true);
  update();
}
