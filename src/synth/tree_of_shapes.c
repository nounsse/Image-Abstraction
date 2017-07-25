#include "tree_of_shapes.h"
#include <limits.h>

/*====== Shape Initialization  ===*/
void shapeInitialize(Shapes pTree)
{
    int i,j;
    Shape pShape;
    Info *ShapeInfo;
    float *boundingbox;
    float *attribute;
    int *label;

    flst_pixels(pTree);

    for(i = pTree->nb_shapes-1; i>= 0; i--)
    {
        pShape = pTree->the_shapes + i;

        if(pShape == NULL )
            continue;

        pShape->data_size = sizeof(struct Info);
        ShapeInfo = (struct Info *) malloc(sizeof(struct Info));

        attribute = (float *) malloc(5*sizeof(float));

        for(j = 0; j< 5; j++)
            attribute[j] = 0.0;

        ShapeInfo->attribute = attribute;

        boundingbox = (float *) malloc(4*sizeof(float));

        for(j = 0; j< 4; j++)
            boundingbox[j] = 0.0;

        ShapeInfo->boundingbox = boundingbox;
        ShapeInfo->index = i;

        label = (int *) malloc(pShape->area*sizeof(int));

        for(j = 0; j< pShape->area; j++)
            label[j] = 0;

        ShapeInfo->label_pixel  = label;
        ShapeInfo->label_max = 0;

        ShapeInfo->show= 1;
        ShapeInfo->lambda1 = 0.0;
        ShapeInfo->lambda2 = 0.0;
        ShapeInfo->xShift =0.0;
        ShapeInfo->yShift =0.0;
        ShapeInfo->rotation =0.0;
        ShapeInfo->contrast= 0.0;
        ShapeInfo->oren = 0.0;
        ShapeInfo->x0 = 0.0;
        ShapeInfo->y0 = 0.0;
        ShapeInfo->DistMin = 1.0;
        ShapeInfo->r = 0.0;
        ShapeInfo->g = 0.0;
        ShapeInfo->b = 0.0;
        ShapeInfo->n_cc = 0;

        ShapeInfo->Mu = mw_new_flist();
        ShapeInfo->Sigma = mw_new_flist();

        ShapeInfo->child = 0;

        pShape->data = (void *)ShapeInfo;
    }
}

TOSParameters getDefaultTOSParameters (){

    TOSParameters tosParam;
    tosParam.order=0; //rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2"
    tosParam.model=2; //synthesis model: rectangle: m=0; ellipse: m=1; orignal shape: m=2
    tosParam.alpha=0.; //alpha for transparent",
    tosParam.ns=3; //  "scale ratio order for color filtering",
    tosParam.threshold=0.6; //  "threshold for color filtering",
    tosParam.smodel=0; //shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
    tosParam.shift=2.0; //add a random shiftS to each shape, shiftS = shift*rand()",
    tosParam.theta=0.0; //add a random rotation to each shape, thetaS = theta*rand()",
    tosParam.mpixel=5; //minimal area (in pixel) for FLST",
    tosParam.maxarea=INT_MAX; //large shape",
    tosParam.kappa=0.; //compactness parameter of the attribute filtering on the orignal image",
    tosParam.relief=0; //     "add relief effects, if relief =1",
    tosParam.reliefOrientation=45; //    "relief orentation, in degree",
    tosParam.reliefHeight=3; //     "relief height",
    tosParam.blur=0;  //     "add blur effects, if blur =1",
    tosParam.median=13; //kernel size for median filter",
    tosParam.kerSize=3; //kernel size for Gaussian blur",
    tosParam.kerStd=0.8; //std for the gaussian kernel",
    tosParam.color_sketch=0; // "compute the sketch: filter shapes based on contrast if=1",
    tosParam.eps=0.0; // "-log10(max number of false alarms)",
    return tosParam;
}


TOSParameters getWaterColorTOSParameters (){

    TOSParameters tosParam;
    tosParam.order=0; //rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2"
    tosParam.model=0; //synthesis model: rectangle: m=0; ellipse: m=1; orignal shape: m=2
    tosParam.alpha=0.; //alpha for transparent",
    tosParam.ns=3; //  "scale ratio order for color filtering",
    tosParam.threshold=0.6; //  "threshold for color filtering",
    tosParam.smodel=0; //shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
    tosParam.shift=10.0; //add a random shiftS to each shape, shiftS = shift*rand()",
    tosParam.theta=0.0; //add a random rotation to each shape, thetaS = theta*rand()",
    tosParam.mpixel=5; //minimal area (in pixel) for FLST",
    tosParam.maxarea=INT_MAX; //large shape",
    tosParam.kappa=0.; //compactness parameter of the attribute filtering on the orignal image",
    tosParam.relief=0; //     "add relief effects, if relief =1",
    tosParam.reliefOrientation=45; //    "relief orentation, in degree",
    tosParam.reliefHeight=3; //     "relief height",
    tosParam.blur=1;  //     "add blur effects, if blur =1",
    tosParam.median=13; //kernel size for median filter",
    tosParam.kerSize=3; //kernel size for Gaussian blur",
    tosParam.kerStd=0.8; //std for the gaussian kernel",
    tosParam.color_sketch=0; // "compute the sketch: filter shapes based on contrast if=1",
    tosParam.eps=0.0; // "-log10(max number of false alarms)",
    return tosParam;
}

TOSParameters getShapeShakingTOSParameters (){

    TOSParameters tosParam;
    tosParam.order=0; //rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2"
    tosParam.model=0; //synthesis model: rectangle: m=0; ellipse: m=1; orignal shape: m=2
    tosParam.alpha=0.; //alpha for transparent",
    tosParam.ns=3; //  "scale ratio order for color filtering",
    tosParam.threshold=0.6; //  "threshold for color filtering",
    tosParam.smodel=0; //shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
    tosParam.shift=10.0; //add a random shiftS to each shape, shiftS = shift*rand()",
    tosParam.theta=0.0; //add a random rotation to each shape, thetaS = theta*rand()",
    tosParam.mpixel=5; //minimal area (in pixel) for FLST",
    tosParam.maxarea=INT_MAX; //large shape",
    tosParam.kappa=0.; //compactness parameter of the attribute filtering on the orignal image",
    tosParam.relief=0; //     "add relief effects, if relief =1",
    tosParam.reliefOrientation=45; //    "relief orentation, in degree",
    tosParam.reliefHeight=3; //     "relief height",
    tosParam.blur=0;  //     "add blur effects, if blur =1",
    tosParam.median=3; //kernel size for median filter",
    tosParam.kerSize=3; //kernel size for Gaussian blur",
    tosParam.kerStd=0.4; //std for the gaussian kernel",
    tosParam.color_sketch=0; // "compute the sketch: filter shapes based on contrast if=1",
    tosParam.eps=0.0; // "-log10(max number of false alarms)",
    return tosParam;
}

TOSParameters getAbstractionTOSParameters (){

    TOSParameters tosParam;
    tosParam.order=1; //rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2"
    tosParam.model=2; //synthesis model: rectangle: m=0; ellipse: m=1; orignal shape: m=2
    tosParam.alpha=0.; //alpha for transparent",
    tosParam.ns=3; //  "scale ratio order for color filtering",
    tosParam.threshold=0.5; //  "threshold for color filtering",
    tosParam.smodel=0; //shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
    tosParam.shift=0.0; //add a random shiftS to each shape, shiftS = shift*rand()",
    tosParam.theta=0.0; //add a random rotation to each shape, thetaS = theta*rand()",
    tosParam.mpixel=20; //minimal area (in pixel) for FLST",
    tosParam.maxarea=INT_MAX; //large shape",
    tosParam.kappa=0.; //compactness parameter of the attribute filtering on the orignal image",
    tosParam.relief=0; //     "add relief effects, if relief =1",
    tosParam.reliefOrientation=45; //    "relief orentation, in degree",
    tosParam.reliefHeight=3; //     "relief height",
    tosParam.blur=0;  //     "add blur effects, if blur =1",
    tosParam.median=3; //kernel size for median filter",
    tosParam.kerSize=3; //kernel size for Gaussian blur",
    tosParam.kerStd=0.5; //std for the gaussian kernel",
    tosParam.color_sketch=1; // "compute the sketch: filter shapes based on contrast if=1",
    tosParam.eps=0.0; // "-log10(max number of false alarms)",
    return tosParam;
}

TOSParameters getStyleTransferTOSParameters (){

    TOSParameters tosParam;
    tosParam.order=1; //rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2"
    tosParam.model=4; //synthesis model: rectangle: m=0; ellipse: m=1; orignal shape: m=2
    tosParam.alpha=0.2; //alpha for transparent",
    tosParam.ns=3; //  "scale ratio order for color filtering",
    tosParam.threshold=0.; //  "threshold for color filtering",
    tosParam.smodel=0; //shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
    tosParam.shift=0.; //add a random shiftS to each shape, shiftS = shift*rand()",
    tosParam.theta=0.0; //add a random rotation to each shape, thetaS = theta*rand()",
    tosParam.mpixel=5; //minimal area (in pixel) for FLST",
    tosParam.maxarea=INT_MAX; //large shape",
    tosParam.kappa=0.; //compactness parameter of the attribute filtering on the orignal image",
    tosParam.relief=0; //     "add relief effects, if relief =1",
    tosParam.reliefOrientation=45; //    "relief orentation, in degree",
    tosParam.reliefHeight=3; //     "relief height",
    tosParam.blur=1;  //     "add blur effects, if blur =1",
    tosParam.median=3; //kernel size for median filter",
    tosParam.kerSize=3; //kernel size for Gaussian blur",
    tosParam.kerStd=0.4; //std for the gaussian kernel",
    tosParam.color_sketch=0; // "compute the sketch: filter shapes based on contrast if=1",
    tosParam.eps=0.0; // "-log10(max number of false alarms)",
    return tosParam;
}

TOSParameters getDictionaryTOSParameters (){

    TOSParameters tosParam;
    tosParam.order=0; //rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2"
    tosParam.model=2; //synthesis model: rectangle: m=0; ellipse: m=1; orignal shape: m=2
    tosParam.alpha=0.; //alpha for transparent",
    tosParam.ns=3; //  "scale ratio order for color filtering",
    tosParam.threshold=0.6; //  "threshold for color filtering",
    tosParam.smodel=0; //shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
    tosParam.shift=0.0; //add a random shiftS to each shape, shiftS = shift*rand()",
    tosParam.theta=0.0; //add a random rotation to each shape, thetaS = theta*rand()",
    tosParam.mpixel=20; //minimal area (in pixel) for FLST",
    tosParam.maxarea=INT_MAX; //large shape",
    tosParam.kappa=0.; //compactness parameter of the attribute filtering on the orignal image",
    tosParam.relief=0; //     "add relief effects, if relief =1",
    tosParam.reliefOrientation=45; //    "relief orentation, in degree",
    tosParam.reliefHeight=3; //     "relief height",
    tosParam.blur=0;  //     "add blur effects, if blur =1",
    tosParam.median=13; //kernel size for median filter",
    tosParam.kerSize=3; //kernel size for Gaussian blur",
    tosParam.kerStd=0.8; //std for the gaussian kernel",
    tosParam.color_sketch=0; // "compute the sketch: filter shapes based on contrast if=1",
    tosParam.eps=0.0; // "-log10(max number of false alarms)",
    return tosParam;
}

TOSParameters getShapeSmoothingTOSParameters (){

    TOSParameters tosParam;
    tosParam.order=10; //rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2"
    tosParam.model=0; //synthesis model: rectangle: m=0; ellipse: m=1; orignal shape: m=2
    tosParam.alpha=0.; //alpha for transparent",
    tosParam.ns=3; //  "scale ratio order for color filtering",
    tosParam.threshold=0.8; //  "threshold for color filtering",
    tosParam.smodel=0; //shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
    tosParam.shift=0.; //add a random shiftS to each shape, shiftS = shift*rand()",
    tosParam.theta=0.0; //add a random rotation to each shape, thetaS = theta*rand()",
    tosParam.mpixel=200; //minimal area (in pixel) for FLST",
    tosParam.maxarea=INT_MAX; //large shape",
    tosParam.kappa=0.; //compactness parameter of the attribute filtering on the orignal image",
    tosParam.relief=0; //     "add relief effects, if relief =1",
    tosParam.reliefOrientation=45; //    "relief orentation, in degree",
    tosParam.reliefHeight=3; //     "relief height",
    tosParam.blur=0;  //     "add blur effects, if blur =1",
    tosParam.median=3; //kernel size for median filter",
    tosParam.kerSize=3; //kernel size for Gaussian blur",
    tosParam.kerStd=0.4; //std for the gaussian kernel",
    tosParam.color_sketch=0; // "compute the sketch: filter shapes based on contrast if=1",
    tosParam.eps=0.0; // "-log10(max number of false alarms)",
    return tosParam;
}


DictionaryParameters getDefaultDictionaryParameters (){

    DictionaryParameters dictionaryParam;
    dictionaryParam.equal=1; //"scaling shape with equal aspect ratio or not, if equal=1, scaling shape with equal aspect ratio, otherwise scaling    shape with different aspect ratio on x-axis and y-axis",
    dictionaryParam.mcolor=1; //   "color model: mcolor=1, use the color of tranferred image, otherwise use the color of the orignal image",
    dictionaryParam.randS=2; //     "selection model: randS=0, randomly select shapes;
    //              randS=1, select shapes according to elongation, compactness and scale;
    //             randS=2, select shapes according to elongation, compactness, scale and color",
    dictionaryParam.paS2I=1.00; //   "parameter for transfer: areaOFshape / areaOFimage, if parameter > paS2I, shape is removed",
    dictionaryParam.paC2S=0.00; //   "parameter for transfer: areaOFconnectcomponent / areaOFshape, if parameter < paC2S, shape is removed",
    dictionaryParam.paS2P=1.00; //   "parameter for transfer: areaOFshape / areaOFshapeparent, if parameter > paS2P, shape is removed",
    dictionaryParam.kappaDict=0; // "compactness parameter of the attribute filtering on the transferred image",
    return dictionaryParam;
}
