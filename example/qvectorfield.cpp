#include "qvectorfield.h"

qVectorField::qVectorField(QWidget* parent)
{
  h_layout = new QHBoxLayout(this);
  setLayout(h_layout);
  h_layout->setMargin(0);
  h_layout->setSpacing(0);
  h_layout->addStretch();
}

void qVectorField::setData(const std::vector<double> data, double min, double max)
{
  clear();
  for (int i = 0; i < data.size(); i++)
  {
    QDoubleSpinBox* box = makeSpinBox(min, max, 0.01);
    box->setValue(data[i]);

    // pass valueChanged to parent
    connect(box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            [this, i](double d) { emit valueChanged(i, d); });
  }
}

void qVectorField::setValue(int idx, double value)
{
  if (boxes[idx]->value() != value)
  {
    boxes[idx]->setValue(value);
    emit valueChanged(idx, value);
  }
}

QDoubleSpinBox* qVectorField::makeSpinBox(double min, double max, double step)
{
  QDoubleSpinBox* box = new QDoubleSpinBox();
  box->setMinimum(min);
  box->setMaximum(max);
  box->setSingleStep(step);
  boxes.emplace_back(box);
  h_layout->insertWidget(h_layout->count() - 1, box);
  return box;
}

void qVectorField::clear()
{
  for (int i = 0; i < boxes.size(); i++)
    delete boxes[i];

  boxes.clear();
  disconnect();
}
