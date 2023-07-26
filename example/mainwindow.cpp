#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QDebug>
#include <chrono>
#include <iostream>

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

  scene->addItem(new qCurve(cp1 * 5));
  scene->addItem(new qCurve(cp2 * 5));

  qCurve c = qCurve(cp1*5);

  // cached basis fucntions
  auto start = std::chrono::steady_clock::now();
//  for (double t=0.0; t<=1.0; t+=0.0001) {
//      c.valueAt(t);
//  }
  for(int k = 0; k < 1000;k++)
      {
      volatile auto asd = c.valueAt(0.5);
  }
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "cached basis functions: " << elapsed_seconds.count() << "\n";

  // cached basis fucntions (again)
  start = std::chrono::steady_clock::now();
//  for (double t=0.0; t<=1.0; t+=0.0001) {
//      c.valueAt(t);
//  }
  for(int k = 0; k < 1000;k++)
      {
      volatile auto asd = c.valueAt3(0.5);
  }
  end = std::chrono::steady_clock::now();
  elapsed_seconds = end - start;
  std::cout << "de boor bla: " << elapsed_seconds.count() << "\n";

  //de boor method
  start = std::chrono::steady_clock::now();
//  for (double t=0.0; t<=1.0; t+=0.0001) {
//      c.valueAt2(t);
//  }
  for(int k = 0; k < 1000;k++)
      {
      volatile auto asd = c.valueAt2(0.5);
  }
  end = std::chrono::steady_clock::now();
  elapsed_seconds = end - start;
  std::cout << "de Boor method: " << elapsed_seconds.count() << "\n";


  ui->graphicsView->centerOn(scene->itemsBoundingRect().center());
}

MainWindow::~MainWindow() { delete ui; }
