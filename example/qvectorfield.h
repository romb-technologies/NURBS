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
  virtual ~qVectorField() { clear(); }
  void clear();
signals:
  void valueChanged(double d, int i);
public slots:
  void setValue(double d, int i);
private:
  std::vector<QDoubleSpinBox*> boxes;
  QDoubleSpinBox* makeSpinBox(double min, double max, double step);
};

#endif // QVECTORFIELD_H
