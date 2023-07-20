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
  if (mouseEvent->button() == Qt::LeftButton)
  {
      for (auto&& curve : items())
      {
        if (!is_curve)
          continue;

        if (mouseEvent->modifiers().testFlag(Qt::ControlModifier))
        {
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
                update_cp = true;
                cp_to_update = std::make_pair(curve, k);
              }
          }
          if (update_cp)
            break;
        }
      }
  }
}


void CustomScene::mouseReleaseEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  if (mouseEvent->button() == Qt::LeftButton)
  {
    update_cp = false;
    update_curvature = false;
  }
}

void CustomScene::mouseMoveEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  NURBS::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());
  if (update_cp)
  {
    auto curve = cp_to_update.first;
    if (is_curve)
    {
      c_curve->prepareGeometryChange();
      c_curve->setControlPoint(cp_to_update.second, p);
    }
    update();
  }
  QGraphicsScene::mouseMoveEvent(mouseEvent);
}
