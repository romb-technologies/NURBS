#include "customscene.h"

#include <QGraphicsSceneMouseEvent>
#include <QKeyEvent>
#include <QMessageBox>

#define is_curve (curve->type() == QGraphicsItem::UserType + 1)
#define is_curve_item(x) (x->type() == QGraphicsItem::UserType + 1)
#define is_poly (curve->type() == QGraphicsItem::UserType + 2)

#define c_curve (static_cast<qCurve*>(curve))
#define c_curve_item(x) (static_cast<qCurve*>(x))

void CustomScene::keyPressEvent(QKeyEvent* keyEvent)
{
  if (keyEvent->key() == 72) // key H
  {
    QMessageBox::information(nullptr, "Help", "\
Mouse controls:\n\
Right click - project mouse pointer on all curves\n\
Ctrl + Scroll - zoom in/out\n\
Ctrl + Left click - select curves\n\
Left click - deselect curves\n\
Left click + drag control point - manipulate control point\n\
Left click + drag curve - manipulate curve (only for 2nd and 3rd order)\n\
\n\
Keyboard shortcuts:\n\
H - display help\n\
B - toggle bounding box display\n\
I - toggle intesections display\n\
C - toggle curvature display (of selected curves)\n\
P - toggle control points display (of selected curves)\n\
Key Up - raise the order (of selected curves)\n\
Key Down - lower the order (of selected curves)\n\
Key + - join multiple curves into polycurve\n\
Delete - delete curve/polycurve");
  }

  if (keyEvent->key() == 80) // key P
  {
    for (auto&& curve : selectedItems())
      if (is_curve)
        c_curve->setDraw_control_points(!c_curve->getDraw_control_points());
    update();
  }
}

void CustomScene::mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  const int sensitivity = 5;
  NURBS::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());

  // point projection
  if (mouseEvent->button() == Qt::RightButton)
  {
    dot = addEllipse(QRectF(QPointF(p.x(), p.y()), QSizeF(6, 6)), QPen(Qt::yellow), QBrush(Qt::red, Qt::SolidPattern));
    for (auto&& curve : items())
    {
      if (is_curve)
      {
        auto t1 = c_curve->projectPoint(p);
        auto p1 = c_curve->valueAt(t1);
        auto tan1 = c_curve->tangentAt(t1);
        line.insert(curve, addLine(QLineF(QPointF(p.x(), p.y()), QPointF(p1.x(), p1.y())), QPen(Qt::red)));
        tan.insert(curve, addLine(QLineF(QPointF(p1.x(), p1.y()) - 150 * QPointF(tan1.x(), tan1.y()),
                                         QPointF(p1.x(), p1.y()) + 150 * QPointF(tan1.x(), tan1.y())),
                                  QPen(Qt::blue)));
//        auto t2 = c_curve->iterateByLength(t1, 50);
//        auto a = c_curve->valueAt(t2);
//        byLength.insert(curve, addEllipse(QRectF(QPointF(a.x() - 3, a.y() - 3), QSizeF(6, 6)), QPen(Qt::yellow), QBrush(Qt::red, Qt::SolidPattern)));
      }
    }
    show_projection = true;
  }
  if (mouseEvent->button() == Qt::LeftButton)
  {
      for (auto&& curve : items())
      {
        if (!is_curve && !is_poly)
          continue;

        if (mouseEvent->modifiers().testFlag(Qt::ControlModifier))
        {
          if (is_curve)
          {
            double t = c_curve->projectPoint(p);
            auto pt = c_curve->valueAt(t);
            if ((pt - p).norm() < 10)
              curve->setSelected(true);
//            update();
          }
        }
        else
        {
          for (auto&& item : selectedItems())
            item->setSelected(false);
          if (is_curve)
          {
            auto pv = c_curve->controlPoints();
            for (uint k = 0; k < pv.size(); k++)
              if ((pv[k] - p).norm() < sensitivity && c_curve->getDraw_control_points())
              {
                number_display = addText(
                            QString("(%1, %2)")
                            .arg(QString::number(p[0]), QString::number(p[1]))
                        );
                number_display->setPos(p[0], p[1]);

                update_cp = true;
                cp_to_update = std::make_pair(curve, k);
              }
          }
          if (update_cp)
            break;
        }
      }
  }
  if (mouseEvent->button() == Qt::MiddleButton) {
      for (auto&& curve : items())
      {
        if (!is_curve && !is_poly)
          continue;

        else
        {
          if (is_curve)
          {
            auto pv = c_curve->controlPoints();
            for (uint k = 0; k < pv.size(); k++)
              if ((pv[k] - p).norm() < sensitivity && c_curve->getDraw_control_points())
              {
                number_display = addText(
                            QString::number(c_curve->weight(k))
                        );
                number_display->setPos(p[0], p[1]);

                update_weights = true;
                last_weight = c_curve->weight(k);
                cp_to_update = std::make_pair(curve, k);
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
        removeItem(dot);
        for (auto&& curve : items())
          if (is_curve || is_poly)
          {
            removeItem(line[curve]);
            removeItem(tan[curve]);
            //removeItem(byLength[curve]);
          }
        line.clear();
        tan.clear();
        //byLength.clear();
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

  // project point
  if (show_projection)
  {
    dot->setRect(QRectF(QPointF(p.x() - 3, p.y() - 3), QSizeF(6, 6)));
    for (auto&& curve : items())
    {
      if (is_curve)
      {
        auto t1 = c_curve->projectPoint(p);
        auto p_c = c_curve->valueAt(t1);
        auto p1 = c_curve->valueAt(t1);
        auto tan1 = c_curve->tangentAt(t1);
        line[curve]->setLine(QLineF(QPointF(p.x(), p.y()), QPointF(p_c.x(), p_c.y())));
        tan[curve]->setLine(QLineF(QPointF(p1.x(), p1.y()) - 500 * QPointF(tan1.x(), tan1.y()),
                                   QPointF(p1.x(), p1.y()) + 500 * QPointF(tan1.x(), tan1.y())));
//        auto t2 = c_curve->iterateByLength(t1, 50);
//        auto a = c_curve->valueAt(t2);
//        byLength[curve]->setRect(QRectF(QPointF(a.x() - 3, a.y() - 3), QSizeF(6, 6)));
      }
    }
  }

  // move control points
  if (update_cp)
  {
    auto curve = cp_to_update.first;
    if (is_curve)
    {
      c_curve->prepareGeometryChange();
      c_curve->setControlPoint(cp_to_update.second, p);
    }
    number_display->setPos(p[0], p[1]);
    number_display->setPlainText(QString("(%1, %2)")
                                 .arg(QString::number(p[0]), QString::number(p[1]))
            );
    update();
  }
  if (update_weights) {
      auto curve = cp_to_update.first;
      if (is_curve)
      {
        NURBS::Point pos = c_curve->controlPoint(cp_to_update.second);
        c_curve->prepareGeometryChange();
        c_curve->setWeight(last_weight + (p[0] - pos[0])*0.005, cp_to_update.second);
      }
      number_display->setPlainText(
                  QString::number(c_curve->weight(cp_to_update.second))
              );
      update();
  }
  QGraphicsScene::mouseMoveEvent(mouseEvent);
}
