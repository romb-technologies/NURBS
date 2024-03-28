#include "mainwindow.h"
#include "ui_mainwindow.h"

void MainWindow::on_actionDisplay_bounding_box_toggled(bool arg1)
{
  scene->draw_box_ = arg1;
  scene->update();
}

void MainWindow::on_actionDisplay_control_points_toggled(bool arg1)
{
  if (!activeCurve)
    return;
  activeCurve->setDraw_control_points(arg1);
  scene->update();
}

void MainWindow::on_actionDisplay_curvature_toggled(bool arg1)
{
  if (!activeCurve)
    return;
  activeCurve->setDraw_curvature_radious(arg1);
  scene->update();
}

void MainWindow::on_actionDisplay_intersections_toggled(bool arg1)
{
  scene->draw_inter_ = arg1;
  scene->update();
}

void MainWindow::on_actionDisplay_knots_toggled(bool arg1)
{
  if (!activeCurve)
    return;
  activeCurve->setDraw_knots(arg1);
  scene->update();
}

void MainWindow::on_actionDelete_triggered()
{
    removeActiveCurveFromScene();
  scene->update();
}

void MainWindow::on_actionConvert_to_Beziers_triggered()
{
  if (!activeCurve)
    return;
  static int counter = 1;
  std::vector<NURBS::Curve> beziers = activeCurve->piecewiseBezier();
  removeActiveCurveFromScene();
  for (NURBS::Curve b : beziers)
    addCurveToScene(b, QString("Segment %1").arg(QString::number(counter++)));
  scene->update();
}

void MainWindow::on_actionRaise_order_triggered()
{
  if (!activeCurve)
    return;
  activeCurve->elevateOrder();
  ui->knotVector->setData(activeCurve->knotVector(), 0.0, 1.0);
  scene->update();
}

void MainWindow::on_actionLower_order_triggered()
{
  if (!activeCurve)
    return;
  activeCurve->lowerOrder();
  scene->update();
}


void MainWindow::on_actionInsert_knot_triggered()
{

}
