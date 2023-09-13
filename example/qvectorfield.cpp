#include "qvectorfield.h"

qVectorField::qVectorField(QWidget *parent)
{
  QHBoxLayout *layout = new QHBoxLayout(this);
  layout->setMargin(0);
  layout->setSpacing(0);
  for (int i=0; i<5; i++){
      QDoubleSpinBox* box = makeSpinBox(0.0, 1.0, 0.01);
      layout->addWidget(box);
      connect(box, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
              this, [this, i](double d) {
                  emit valueChanged(d, i);
              });
  }
  layout->addStretch();
}

void qVectorField::setValue(double d, int i) {
  if (boxes[i]->value() != d) {
      boxes[i]->setValue(d);
      emit valueChanged(d, i);
    }
}

QDoubleSpinBox* qVectorField::makeSpinBox(double min, double max, double step) {
    QDoubleSpinBox* box = new QDoubleSpinBox();
    box->setMinimum(min);
    box->setMaximum(max);
    box->setSingleStep(step);
    boxes.emplace_back(box);
    return box;
}

void qVectorField::clear() {
      for (int i=0; i<boxes.size(); i++) {
          delete boxes[i];
      }
      boxes.clear();
}
