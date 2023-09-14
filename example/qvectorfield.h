#ifndef QVECTORFIELD_H
#define QVECTORFIELD_H

#include <QDoubleSpinBox>
#include <QHBoxLayout>
#include <QSpacerItem>

class qVectorField : public QWidget
{
  Q_OBJECT
public:
  qVectorField(QWidget *parent = nullptr);
  virtual ~qVectorField() {}
  void clear();
  void setData(const std::vector<double> data, double min, double max);
signals:
  void valueChanged(int idx, double value);
public slots:
  void setValue(int idx, double value);
private:
  std::vector<QDoubleSpinBox*> boxes;
  QDoubleSpinBox* makeSpinBox(double min, double max, double step);
};

#endif // QVECTORFIELD_H
