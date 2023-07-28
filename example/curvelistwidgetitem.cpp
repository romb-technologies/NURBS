#include "curvelistwidgetitem.h"

CurveListWidgetItem::CurveListWidgetItem(qCurve *c, QListWidget *parent):
    QListWidgetItem("Curve", parent),
    curve(c) {}
