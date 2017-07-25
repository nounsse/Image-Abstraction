#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include <image.h>
#include <misc.h>
#include <pnmfile.h>

#include <QImage>
#include "disjoint-set.h"

class Segmentation
{
public:
    Segmentation(const QImage & input_image);
    ~Segmentation();
    QImage segment(float sigma, float c, int min_size);
    QImage removeRegionUnder( const QPoint &clickedPoint );
    QImage getResult( );
protected:
    image<rgb> *input;
    image<rgb> *segmentation;
    Universe *u ;
    edge *edges;
    QImage result;
    int num;
    std::map<int, QColor> component_colors;
    std::map<int, bool> removed_regions;

};

#endif // SEGMENTATION_H
