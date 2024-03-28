#include "qinsertknotdialog.h"

QInsertKnotDialog::QInsertKnotDialog(QWidget *parent) : QDialog(parent)
{
    setModal(true);
    layout = new QFormLayout();
    knot = new QDoubleSpinBox();
    multi = new QSpinBox();
    layout.addRow(tr("Knot"), knot);
    layout.addRow(tr("Multiplicity"), multi);
}
