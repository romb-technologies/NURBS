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
  curve2.insertKnot(0.25, 3);
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
            if (current) {
              if (activeCurve) {
                  activeCurve->disconnect(ui->knotVector);
                  activeCurve->disconnect(ui->weights);
                }
              activeCurve = static_cast<CurveListWidgetItem*>(current)->curve;

              ui->knotVector->setData(activeCurve->knotVector(), 0.0, 1.0);
              connect(ui->knotVector,
                      &qVectorField::valueChanged,
                      activeCurve,
                      &qCurve::setKnot);
              connect(activeCurve,
                      &qCurve::knotChanged,
                      ui->knotVector,
                      &qVectorField::setValue);

              ui->weights->setData(activeCurve->weights(), 0.0, 100.0);
              connect(ui->weights,
                      &qVectorField::valueChanged,
                      activeCurve,
                      &qCurve::setWeight);
              connect(activeCurve,
                      &qCurve::weightChanged,
                      ui->weights,
                      &qVectorField::setValue);

              ui->infoText->setText(QString("Order: %1 | Points: %2").arg(
                                      QString::number(activeCurve->order()),
                                      QString::number(activeCurve->controlPoints().size()))
                                    );

              ui->customPlot->clearGraphs();
              int N = activeCurve->controlPoints().size();
              for (int pt=0; pt<N; pt++) {
                  ui->customPlot->addGraph()->setPen(QPen(QColor(255.0*pt/(N-1), 0, 255.0*(1-pt/(N-1)))));
                }
              graphKnotVector(activeCurve);
              connect(activeCurve,
                      &qCurve::knotChanged,
                      this,
                      [this]() {
                  graphKnotVector(activeCurve);
                });
              }
            else {
                activeCurve = nullptr;
              }
  });

  ui->customPlot->xAxis->setRange(0, 1);
  ui->customPlot->yAxis->setRange(0, 1);

  ui->curveList->setCurrentRow(0);
  ui->graphicsView->scene()->setItemIndexMethod(QGraphicsScene::NoIndex);
  ui->graphicsView->centerOn(scene->itemsBoundingRect().center());
}

qCurve* MainWindow::addCurveToScene(NURBS::Curve c)
{
    static uint counter = 1;
    qCurve *qc = addCurveToScene(c, QString("Curve %1").arg(QString::number(counter++)));
    return qc;
}

qCurve* MainWindow::addCurveToScene(NURBS::Curve c, QString name)
{
    qCurve* qc = new qCurve(c);
    scene->addItem(qc);
    new CurveListWidgetItem(qc,
                            name,
                            ui->curveList);
    connect(qc,
            &qCurve::curveChanged,
            this,
            [this]() {
        scene->update();
      });
    activeCurve = qc;
    return qc;
}

void MainWindow::removeActiveCurveFromScene()
{
    if (!activeCurve)
      return;
    scene->removeItem(activeCurve);
    delete activeCurve;
    activeCurve = nullptr;

    scene->selectedItems().clear();

    ui->knotVector->clear();
    ui->weights->clear();

    ui->customPlot->clearGraphs();

    ui->infoText->setText("");
    ui->curve_name->setText("");

    QListWidgetItem *it = ui->curveList->takeItem(ui->curveList->currentRow());
    delete it;

}


void MainWindow::graphKnotVector(qCurve *curve)
{
    std::vector<QVector<double>> x, y;

    int m = curve->order() + 1;
    int N = curve->controlPoints().size();
    double detail = 1000;

    for (int pt=0; pt<N; pt++) {
        x.emplace_back(QVector<double>(detail+1));
        y.emplace_back(QVector<double>(detail+1));
      }

    for (int i=0; i<=detail; i++) {
        Eigen::VectorXd bf = curve->getBasisFunctionsAt(i/detail);

        int sp = curve->getKnotSpanIndex(i/detail)-m+1;
        for (int pt=0; pt<m; pt++) {
            x[sp+pt][i] = i/detail;
            y[sp+pt][i] = bf(pt);
        }
      }
    for (int pt=0; pt<N; pt++)
      ui->customPlot->graph(pt)->setData(x[pt], y[pt]);

    ui->customPlot->replot();
}

MainWindow::~MainWindow() { delete ui; }

