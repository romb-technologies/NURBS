#ifndef QGRAPHICSVIEWZOOM_H
#define QGRAPHICSVIEWZOOM_H

#include <QDebug>
#include <QGraphicsView>
#include <QObject>

class QGraphicsViewZoom : public QObject
{
  Q_OBJECT
public:
  QGraphicsViewZoom(QGraphicsView* view);
  void gentle_zoom(double factor);
  void set_modifiers(Qt::KeyboardModifiers modifiers);
  void set_zoom_factor_base(double value);

private:
  QGraphicsView* m_view;
  Qt::KeyboardModifiers m_modifiers;
  double m_zoom_factor_base;
  QPointF target_scene_pos, target_viewport_pos;
  bool eventFilter(QObject* object, QEvent* event) Q_DECL_OVERRIDE;

signals:
  void zoomed();
};

#endif // QGRAPHICSVIEWZOOM_H
