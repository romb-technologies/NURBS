#include "curvelistwidgetitem.h"

CurveListWidgetItem::CurveListWidgetItem(qCurve* c, QString name, QListWidget* parent)
    : QListWidgetItem(name, parent), curve(c)
{
}
