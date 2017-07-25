#ifndef _TREE_OF_SHAPES_H_
#define _TREE_OF_SHAPES_H_

#ifdef __cplusplus
extern "C" {
#endif
#include "mw3.h"

typedef struct Info{
    int index;
    int label;
    int show;
    Flist Mu;
    Flist Sigma;
    int *label_pixel;
    int label_max;
    float *attribute;
    float lambda1;
    float lambda2;
    float *boundingbox;
    float oren;
    float x0;
    float y0;
    float xShift;
    float yShift;
    float rotation;
    float contrast;
    float r;
    float g;
    float b;
    float DistMin;
    int n_cc;
    int child;
}Info;

typedef struct TOSParameters {
    int order; //rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2"
    int model; //synthesis model: orignal shape: m=0; ellipse: m=1; rectangle: m=2; random: m=3;
    float alpha; //alpha for transparent",
    int ns; //  "scale ratio order for color filtering",
    float threshold; //  "threshold for color filtering",
    int smodel; //shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
    float shift; //add a random shiftS to each shape, shiftS = shift*rand()",
    float theta; //add a random rotation to each shape, thetaS = theta*rand()",
    int mpixel; //minimal area (in pixel) for FLST",
    int maxarea; //large shape",
    float kappa; //compactness parameter of the attribute filtering on the orignal image",
    int relief; //     "add relief effects, if relief =1",
    float reliefOrientation; //    "relief orentation, in degree",
    float reliefHeight; //     "relief height",
    int blur; //     "add blur effects, if blur =1",
    int median; //kernel size for median filter",
    int kerSize; //kernel size for Gaussian blur",
    float kerStd; //std for the gaussian kernel",
    int color_sketch; //     "compute the sketch: filter shapes based on contrast if=1",
    float eps; // "-log10(max number of false alarms)",

} TOSParameters;

typedef struct DictionaryParameters {
    int equal; //"scaling shape with equal aspect ratio or not, if equal=1, scaling shape with equal aspect ratio, otherwise scaling    shape with different aspect ratio on x-axis and y-axis",
    int mcolor; //   "color model: mcolor=1, use the color of tranferred image, otherwise use the color of the orignal image",
    int randS; //     "selection model: randS=0, randomly select shapes;
    //              randS=1, select shapes according to elongation, compactness and scale;
    //             randS=2, select shapes according to elongation, compactness, scale and color",
    float paS2I; //   "parameter for transfer: areaOFshape / areaOFimage, if parameter > paS2I, shape is removed",
    float paC2S; //   "parameter for transfer: areaOFconnectcomponent / areaOFshape, if parameter < paC2S, shape is removed",
    float paS2P; //   "parameter for transfer: areaOFshape / areaOFshapeparent, if parameter > paS2P, shape is removed",
    float kappaDict; // "compactness parameter of the attribute filtering on the transferred image",
} DictionaryParameters;

TOSParameters getDefaultTOSParameters();
TOSParameters getWaterColorTOSParameters ();
TOSParameters getShapeShakingTOSParameters ();
TOSParameters getAbstractionTOSParameters ();
TOSParameters getStyleTransferTOSParameters();
TOSParameters getShapeSmoothingTOSParameters();

DictionaryParameters getDefaultDictionaryParameters ();
TOSParameters getDictionaryTOSParameters ();

void shapeInitialize(Shapes pTree);

#ifdef __cplusplus
}
#endif

#endif /* !_TREE_OF_SHAPES_H_ */
