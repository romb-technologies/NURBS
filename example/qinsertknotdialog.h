#ifndef QINSERTKNOTDIALOG_H
#define QINSERTKNOTDIALOG_H

#include <QDialog>
#include <QObject>
#include <QFormLayout>
#include <QDoubleSpinBox>
#include <QSpinBox>


class QInsertKnotDialog : public QDialog
{
    Q_OBJECT
public:
    QInsertKnotDialog(QWidget *parent = nullptr);
private:
    QFormLayout layout;
    QDoubleSpinBox knot;
    QSpinBox multi;
};

#endif // QINSERTKNOTDIALOG_H
