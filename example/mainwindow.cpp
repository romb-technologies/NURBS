#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QDebug>
#include <chrono>
#include <iostream>
#include "curvelistwidgetitem.h"

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

  qCurve* curve1 = new qCurve(cp1 * 5);
  qCurve* curve2 = new qCurve(cp2 * 5);

  scene->addItem(curve1);
  scene->addItem(curve2);

  activeCurve = curve2;

  new CurveListWidgetItem(curve1, ui->curveList);
  new CurveListWidgetItem(curve2, ui->curveList);

  connect(ui->curveList,
          &QListWidget::currentItemChanged,
          this,
          [this](QListWidgetItem *current) {
            activeCurve = static_cast<CurveListWidgetItem*>(current)->curve;
            displayKnotVector(activeCurve);
            displayWeights(activeCurve);
  });

  // knot vector and weights
  displayKnotVector(activeCurve);
  displayWeights(activeCurve);


  // testiranje
  qCurve c = qCurve(cp1*5);

  // cached basis fucntions
  auto start = std::chrono::steady_clock::now();
  for(int k = 0; k < 1000;k++)
      {
      volatile auto asd = c.valueAt(0.5);
  }
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "cached basis functions: " << elapsed_seconds.count() << "\n";

  // cached basis fucntions (again)
  for(int k = 0; k < 1000;k++)
      {
      volatile auto asd = c.valueAt3(0.5);
  }
  end = std::chrono::steady_clock::now();
  elapsed_seconds = end - start;
  std::cout << "de boor bla: " << elapsed_seconds.count() << "\n";

  //de boor method
  start = std::chrono::steady_clock::now();
  for(int k = 0; k < 1000;k++)
      {
      volatile auto asd = c.valueAt2(0.5);
  }
  end = std::chrono::steady_clock::now();
  elapsed_seconds = end - start;
  std::cout << "de Boor method: " << elapsed_seconds.count() << "\n";


  ui->graphicsView->centerOn(scene->itemsBoundingRect().center());
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
    // Clear previous knot vector
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
