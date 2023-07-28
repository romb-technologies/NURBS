#ifndef CURVELISTWIDGETITEM_H
#define CURVELISTWIDGETITEM_H

#include <QWidget>
#include <QListWidgetItem>
#include <qcurve.h>

class CurveListWidgetItem: public QListWidgetItem
{
public:
    CurveListWidgetItem(qCurve *c, QListWidget* parent);
    qCurve *curve;
};

#endif // CURVELISTWIDGETITEM_H
