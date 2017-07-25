
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "segment-image.h"

#include "Segmentation.h"
#include <QColor>

Segmentation::Segmentation(const QImage & input_image )
{

    int width = input_image.width();
    int height = input_image.height();

    input = new image<rgb> ( width, height );

    // smooth each color channel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {

            QColor color =input_image.pixel( x, y );
            imRef(input, x, y).r = (uchar)color.red();
            imRef(input, x, y).g = (uchar)color.green();
            imRef(input, x, y).b = (uchar)color.blue();
        }
    }

    result = QImage(input_image);
}

Segmentation::~Segmentation(){
    delete [] edges;
    delete [] u;
}

QImage Segmentation::segment(float sigma, float c, int min_size) {
    int width = input->width();
    int height = input->height();

    image<float> *r = new image<float>(width, height);
    image<float> *g = new image<float>(width, height);
    image<float> *b = new image<float>(width, height);

    // smooth each color channel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            imRef(r, x, y) = imRef(input, x, y).r;
            imRef(g, x, y) = imRef(input, x, y).g;
            imRef(b, x, y) = imRef(input, x, y).b;
        }
    }
    image<float> *smooth_r = smooth(r, sigma);
    image<float> *smooth_g = smooth(g, sigma);
    image<float> *smooth_b = smooth(b, sigma);
    delete r;
    delete g;
    delete b;

    // build graph
    edges = new edge[width*height*4];
    num = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (x < width-1) {
                edges[num].a = y * width + x;
                edges[num].b = y * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
                num++;
            }

            if (y < height-1) {
                edges[num].a = y * width + x;
                edges[num].b = (y+1) * width + x;
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x, y+1);
                num++;
            }

            if ((x < width-1) && (y < height-1)) {
                edges[num].a = y * width + x;
                edges[num].b = (y+1) * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y+1);
                num++;
            }

            if ((x < width-1) && (y > 0)) {
                edges[num].a = y * width + x;
                edges[num].b = (y-1) * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y-1);
                num++;
            }
        }
    }
    delete smooth_r;
    delete smooth_g;
    delete smooth_b;

    // segment
    u = segment_graph(width*height, num, edges, c);


    // post process small components
    for (int i = 0; i < num; i++) {
        int a = u->find(edges[i].a);
        int b = u->find(edges[i].b);
        if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
            u->join(a, b);
    }

    int num_ccs = u->num_sets();

    segmentation = new image<rgb>(width, height);


    // pick random colors for each component
    std::map<int, float> average_r;
    std::map<int, float> average_g;
    std::map<int, float> average_b;

    std::map<int, int> nb_pixels;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int comp = u->find(y * width + x);

            average_r[comp] += imRef(input, x, y).r;
            average_g[comp] += imRef(input, x, y).g;
            average_b[comp] += imRef(input, x, y).b;
            nb_pixels[comp]++;
        }
    }

    // pick random colors for each component
    rgb *colors = new rgb[width*height];
    for (int i = 0; i < width*height; i++)
        colors[i] = random_rgb();

    std::cout << "Num ccs "<< num_ccs << " and " << nb_pixels.size() << std::endl;

    for( std::map<int, int>::iterator it = nb_pixels.begin() ; it != nb_pixels.end() ; it++ ){
        int comp = it->first;
        component_colors[comp] = QColor(int(float(average_r[comp])/nb_pixels[comp]), int(float(average_g[comp])/nb_pixels[comp]), int(float(average_b[comp])/nb_pixels[comp]));
        removed_regions[comp] = false;
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int comp = u->find(y * width + x);

            imRef(segmentation, x, y) = colors[comp];
            QColor color (colors[comp].r, colors[comp].g, colors[comp].b);
            result.setPixel(x, y , qRgb(color.red(), color.green(), color.blue()));
        }
    }



    return result;
}

/*
QImage Segmentation::removeRegionUnder( const QPoint &A, const QPoint &B ){

    int width = input->width();
    int height = input->height();

    if( clickedPoint.x()< width && clickedPoint.y()<height ){

        int region_id = u->find(clickedPoint.y() * width + clickedPoint.x());
        removed_regions[region_id] = true;

        std::map<int, bool> isolated_region;
        for(std::map<int, bool>::iterator it = removed_regions.begin(); it != removed_regions.end() ; it++){
            isolated_region[it->first] = true;
        }

        // post process small components
        for (int i = 0; i < num; i++) {
            int a = u->find(edges[i].a);
            int b = u->find(edges[i].b);
            if ( a != b ){
                if( !removed_regions[a] ) isolated_region [b] = false;
                if( !removed_regions[b] ) isolated_region [a] = false;
            }
        }

        for(std::map<int, bool>::iterator it = isolated_region.begin(); it != isolated_region.end() ; it++){
            if( it->second )
                removed_regions[it->first] = true;
        }

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int comp = u->find(y * width + x);

                //QColor color;
                if( comp == region_id || isolated_region [comp]){
                    //color = component_colors[comp];
                    result.setPixel(x, y , qRgb(255, 255, 255));
                }

            }
        }
        //  result.setPixel(clickedPoint.x(), clickedPoint.y(), qRgb(255, 0, 0));

    }
    return result;
}
*/
QImage Segmentation::removeRegionUnder( const QPoint &clickedPoint ){

    int width = input->width();
    int height = input->height();

    if( clickedPoint.x()< width && clickedPoint.y()<height ){

        int region_id = u->find(clickedPoint.y() * width + clickedPoint.x());
        removed_regions[region_id] = true;

        std::map<int, bool> isolated_region;
        for(std::map<int, bool>::iterator it = removed_regions.begin(); it != removed_regions.end() ; it++){
            isolated_region[it->first] = true;
        }

        // post process small components
        for (int i = 0; i < num; i++) {
            int a = u->find(edges[i].a);
            int b = u->find(edges[i].b);
            if ( a != b ){
                if( !removed_regions[a] ) isolated_region [b] = false;
                if( !removed_regions[b] ) isolated_region [a] = false;
            }
        }

        for(std::map<int, bool>::iterator it = isolated_region.begin(); it != isolated_region.end() ; it++){
            if( it->second )
                removed_regions[it->first] = true;
        }

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int comp = u->find(y * width + x);

                //QColor color;
                if( comp == region_id || isolated_region [comp]){
                    //color = component_colors[comp];
                    result.setPixel(x, y , qRgb(255, 255, 255));
                }

            }
        }
        //  result.setPixel(clickedPoint.x(), clickedPoint.y(), qRgb(255, 0, 0));

    }
    return result;
}

QImage Segmentation::getResult( ){

    int width = input->width();
    int height = input->height();



    // pick random colors for each component
    float average_r =0.;
    float average_g =0.;
    float average_b =0.;

    int nb_pixels = 0;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int comp = u->find(y * width + x);
            if( removed_regions[comp] ){
                average_r += imRef(input, x, y).r;
                average_g += imRef(input, x, y).g;
                average_b += imRef(input, x, y).b;
                nb_pixels++;
            }
        }
    }

    average_r = int(average_r/nb_pixels);
    average_g = int(average_g/nb_pixels);
    average_b = int(average_b/nb_pixels);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int comp = u->find(y * width + x);
            if( ! removed_regions[comp] ){

                rgb color = imRef(input, x, y);
                result.setPixel(x, y , qRgb(color.r, color.g, color.b));
            } else
                result.setPixel(x, y , qRgb(average_r, average_g, average_b));
        }
    }
    return result;
}
