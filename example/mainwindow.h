#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDoubleSpinBox>
#include <QListWidgetItem>

#include "qcustomplot.h"
#include "customscene.h"

namespace Ui
{
class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  explicit MainWindow(QWidget* parent = nullptr);
  ~MainWindow();

private slots:
  void on_actionDisplay_bounding_box_toggled(bool arg1);
  void on_actionDisplay_control_points_toggled(bool arg1);
  void on_actionDisplay_curvature_toggled(bool arg1);
  void on_actionDisplay_intersections_toggled(bool arg1);
  void on_actionDisplay_knots_toggled(bool arg1);
  void on_actionDelete_triggered();

  void on_actionConvert_to_Beziers_triggered();

private:
  Ui::MainWindow* ui;
  CustomScene* scene;

  std::vector<QDoubleSpinBox*> knotVectorField;
  std::vector<QDoubleSpinBox*> weightField;

  qCurve* activeCurve;

  QDoubleSpinBox* makeSpinBox(double min, double max, double step, std::vector<QDoubleSpinBox*>& addTo);
  void displayKnotVector(qCurve *curve);
  void displayWeights(qCurve *curve);
  void onUpdateKnot();
  void onUpdateWeight();
  qCurve *addCurveToScene(NURBS::Curve c);
  qCurve *addCurveToScene(NURBS::Curve c, QString name);
  void removeActiveCurveFromScene();
  void graphKnotVector(qCurve *curve);
};

#endif // MAINWINDOW_H
