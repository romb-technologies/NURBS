#include "mainwindow.h"

void MainWindow::on_actionDisplay_bounding_box_toggled(bool arg1)
{
    scene->draw_box_ = arg1;
    scene->update();
}

void MainWindow::on_actionDisplay_control_points_toggled(bool arg1)
{
    activeCurve->setDraw_control_points(arg1);
    scene->update();
}

void MainWindow::on_actionDisplay_curvature_toggled(bool arg1)
{
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
    activeCurve->setDraw_knots(arg1);
    scene->update();
}

void MainWindow::on_actionDelete_triggered()
{
    removeActiveCurveFromScene();
    scene->update();
}
