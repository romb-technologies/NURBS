#ifndef QVECTORFIELD_H
#define QVECTORFIELD_H

#include <QDoubleSpinBox>
#include <QHBoxLayout>
#include <QSpacerItem>

class qVectorField : public QWidget
{
  Q_OBJECT
public:
  qVectorField(QWidget* parent = nullptr);
  virtual ~qVectorField() {}
signals:
  void valueChanged(int idx, double value);
public slots:
  void setValue(int idx, double value);
  void clear();
  void setData(const std::vector<double> data, double min, double max);

private:
  std::vector<QDoubleSpinBox*> boxes;
  QDoubleSpinBox* makeSpinBox(double min, double max, double step);
  QHBoxLayout* h_layout;
};

#endif // QVECTORFIELD_H
