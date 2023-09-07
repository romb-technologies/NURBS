#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QDebug>
#include <chrono>
#include <iostream>
#include "curvelistwidgetitem.h"
#include <QStatusBar>

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
    , scene(new CustomScene)
{
  ui->setupUi(this);

  ui->graphicsView->setScene(scene);

  new QGraphicsViewZoom(ui->graphicsView);

  Eigen::MatrixX2d cp1, cp2;
  cp1.resize(4, 2);
  cp2.resize(5, 2);
  cp1 << 84, 162,
      246, 30,
      48, 236,
      180, 110;

  cp2 << 180, 110,
      175, 160,
      60, 48,
      164, 165,
      124, 134;

  NURBS::Curve curve1(cp1 * 5);
  NURBS::Curve curve2(cp2 * 5);

  curve1.appendPoint({600, 1000});
  auto split = curve2.splitCurve(0.4);
  auto beziers = curve2.piecewiseBezier();

  addCurveToScene(curve1);
  addCurveToScene(curve2);
//  addCurveToScene(split.first);
//  addCurveToScene(split.second);
//  for (NURBS::Curve b: beziers) {
//      addCurveToScene(b);
//  }

  // init knot vector and weight display
  connect(ui->curveList,
          &QListWidget::currentItemChanged,
          this,
          [this](QListWidgetItem *current) {
            activeCurve = static_cast<CurveListWidgetItem*>(current)->curve;
            displayKnotVector(activeCurve);
            displayWeights(activeCurve);
  });
  displayKnotVector(activeCurve);
  displayWeights(activeCurve);

  ui->graphicsView->centerOn(scene->itemsBoundingRect().center());
}

qCurve* MainWindow::addCurveToScene(NURBS::Curve c) {
    qCurve* qc = new qCurve(c);
    scene->addItem(qc);
    new CurveListWidgetItem(qc, ui->curveList);
    activeCurve = qc;
    return qc;
}

QDoubleSpinBox* MainWindow::makeSpinBox(double min, double max, double step, std::vector<QDoubleSpinBox*>& addTo) {
    QDoubleSpinBox* box = new QDoubleSpinBox();
    box->setMinimum(min);
    box->setMaximum(max);
    box->setSingleStep(step);
    addTo.emplace_back(box);
    return box;
}


void MainWindow::displayKnotVector(qCurve *curve) {
    //Clear previous knot vector
    for (int i=0; i<knotVectorField.size(); i++) {
        delete knotVectorField[i];
    }
    knotVectorField.clear();

    for (int i=0; i<curve->knotVector().rows(); i++){
        QDoubleSpinBox* box = makeSpinBox(0.0, 1.0, 0.01, knotVectorField);
        ui->knotVectorLayout->insertWidget(
                    ui->knotVectorLayout->count() - 1,
                    box);
        box->setValue(curve->knotVector()(i));
        connect(box,
                QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this,
                [this, i, curve](double d){
                        curve->setKnot(i, d);
                        qobject_cast<QDoubleSpinBox*>(sender())->setValue(curve->knot(i));
                        scene->update();
                });
    }
}

void MainWindow::displayWeights(qCurve* curve) {
    // Clear previous weights
    for (int i=0; i<weightField.size(); i++) {
        delete weightField[i];
    }
    weightField.clear();

    for (int i=0; i<curve->weights().rows(); i++){
        QDoubleSpinBox* box = makeSpinBox(-100.0, 100.0, 0.01, weightField);
        ui->weightLayout->insertWidget(
                    ui->weightLayout->count() - 1,
                    box);
        box->setValue(curve->weight(i));
        connect(box,
                QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this,
                [this, i, curve](double d){
                        curve->setWeight(d, i);
                        scene->update();
                });
    }
}

MainWindow::~MainWindow() { delete ui; }
