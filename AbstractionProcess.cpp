#include "AbstractionProcess.h"

#include "tree_of_shapes.h"
#include <iostream>
#include <sys/time.h>
#include <QColor>

extern void flst();
#define PI 3.1415926
#define _MIN(x, y) ( (x)<(y) ? (x) : (y) )
#define _MAX(x, y) ( (x)>(y) ? (x) : (y) )

AbstractionProcess::AbstractionProcess( std::string fileNameIn )
{

    /**************************************************/
    /*************   READ INPUT IMAGE   ***************/
    /**************************************************/

    _imgin = cfimageread(fileNameIn.c_str());

    /**************************************************/
    /*************   SET DEFAULT INPUT OPTIONS   **************/
    /**************************************************/

    // init(_imgin);
    _image_loaded = true;
    _tosParameters = getDefaultTOSParameters();
    _tree_computed = false;
}

Cfimage AbstractionProcess::cfimages_from_qimage( const QImage &input_image  ){

    int nx = input_image.width(),ny = input_image.height();
    Cfimage out = mw_change_cfimage(NULL,ny,nx);

    float * red = (float *)malloc(nx*ny*sizeof(float));
    float * green = (float *)malloc(nx*ny*sizeof(float));
    float * blue = (float *)malloc(nx*ny*sizeof(float));

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            QColor color = input_image.pixel( x, y );
            int comp = y * nx + x;
            red[comp]= color.red();
            green[comp]= color.green();
            blue[comp]= color.blue();
        }
    }
    free(out->red);
    free(out->green);
    free(out->blue);
    out->red = red;
    out->green = green;
    out->blue = blue;


    return out;

}

QImage AbstractionProcess::qimages_from_cfimage( Cfimage input_img  ){
    QImage result_image( QSize(input_img->ncol, input_img->nrow), QImage::Format_RGB32 );

    for( int j= 0; j< input_img->nrow; j++)
        for( int i= 0; i< input_img->ncol; i++)
        {
            int comp = j*input_img->ncol + i;

            QColor color (input_img->red[comp], input_img->green[comp], input_img->blue[comp]);
            result_image.setPixel(i, j , qRgb(color.red(), color.green(), color.blue()));
        }

    return result_image;
}

QImage AbstractionProcess::qimages_from_ccimage( Ccimage input_img  ){
    QImage result_image( QSize(input_img->ncol, input_img->nrow), QImage::Format_RGB32 );

    for( int j= 0; j< input_img->nrow; j++)
        for( int i= 0; i< input_img->ncol; i++)
        {
            int comp = j*input_img->ncol + i;

            QColor color (input_img->red[comp], input_img->green[comp], input_img->blue[comp]);
            result_image.setPixel(i, j , qRgb(color.red(), color.green(), color.blue()));
        }

    return result_image;
}

AbstractionProcess::AbstractionProcess( const QImage &imageIn ){


    /**************************************************/
    /*************   READ INPUT IMAGE   ***************/
    /**************************************************/

    _imgin = cfimages_from_qimage(imageIn);

    /**************************************************/
    /*************   SET DEFAULT INPUT OPTIONS   **************/
    /**************************************************/

    //init(_imgin);
    _image_loaded = true;
    _tosParameters = getDefaultTOSParameters();
    _tree_computed = false;

    _treeOfShapes = new TreeOfShapes( cfimages_from_qimage(imageIn) );
}

void AbstractionProcess::addDictionnary(TreeOfShapes * dictionary ){
    /**************************************************/
    /*************   READ INPUT DICT   ***************/
    /**************************************************/

    _dictionnary = dictionary;


}

void AbstractionProcess::init(Cfimage inputImg, Shapes &pTree){

    Shape pShape;
    Fimage imgIntensity;

    if  ( ((imgIntensity = mw_new_fimage()) == NULL) ||
          (mw_alloc_fimage(imgIntensity,inputImg->nrow,inputImg->ncol) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");

    if((pTree = mw_new_shapes()) == NULL)
        mwerror(FATAL, 1,
                "fgrain --> Not enough memory to allocate the tree of shapes");
    if  ( ((_NormOfDu = mw_new_fimage()) == NULL) ||
          (mw_alloc_fimage(_NormOfDu,inputImg->nrow,inputImg->ncol) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");
    /*=======================================*/
    /*==    Compute Intensity imag        ===*/
    /*=======================================*/
    for( int i= 0; i< inputImg->ncol; i++)
        for( int j= 0; j< inputImg->nrow; j++)
        {
            //imgIntensity->gray[j*imgin->ncol + i] = (imgin->blue[j*imgin->ncol + i] * 255.0 );
            imgIntensity->gray[j*inputImg->ncol + i] = (int)(inputImg->blue[j*inputImg->ncol + i]
                    + inputImg->red[j*inputImg->ncol + i]
                    + inputImg->green[j*inputImg->ncol + i])/3;
        }
    float fzero = 0.; int nsize = 3;
    fderiv(imgIntensity,NULL,NULL,NULL,NULL,NULL,NULL,_NormOfDu,NULL,&fzero,&nsize);

    /*=======================================*/
    /*=== Compute FLST on Intensity image ===*/
    /*=======================================*/
    int minArea = 5;
    flst(&minArea, imgIntensity, pTree);

    /*=======================================*/
    /*==    Initialization                ===*/
    /*=======================================*/
    shapeInitialize(pTree);

    /*=======================================*/
    /*==   Assign color to each shape     ===*/
    /*=======================================*/
    // assign color to each shape
    Shape * ppShapeOfPixel = pTree->smallest_shape;

    for(int i= 0; i< pTree->ncol*pTree->nrow; i++){
        pShape = ppShapeOfPixel[i];

        ((Info*)(pShape->data))->r += inputImg->red[i];
        ((Info*)(pShape->data))->g += inputImg->green[i];
        ((Info*)(pShape->data))->b += inputImg->blue [i];
        ((Info*)(pShape->data))->n_cc++;
    }

    for(int i= 0; i< pTree->nb_shapes; i++){
        pShape = pTree->the_shapes + i;

        ((Info*)(pShape->data))->r /= (float)((Info*)(pShape->data))->n_cc;
        ((Info*)(pShape->data))->g /= (float)((Info*)(pShape->data))->n_cc;
        ((Info*)(pShape->data))->b /= (float)((Info*)(pShape->data))->n_cc;
        //printf("%f, %f\n", ((Info*)(pShape->data))->b, pShape->value/255.0);
    }

    mw_delete_fimage(imgIntensity);
}




AbstractionProcess::~AbstractionProcess(){
    if(_tree_computed){
        mw_delete_shapes(_pTree);
    }
}

Cfimage AbstractionProcess::cfimageread(const char* name)
{

    QImage image(name);
    return cfimages_from_qimage(image);
}


/* This removes the shapes from the tree associated to pFloatImageInput
that are too small (threshold *pMinArea). As a consequence all the remaining
shapes of pFloatImageOutput are of area larger or equal than *pMinArea */

void AbstractionProcess::mw_fgrain_side(int *pMinArea, Fimage pFloatImageInput, Fimage pFloatImageOutput, int sideflag)
{
    int i;
    Shapes pTree;

    if(mw_change_fimage(pFloatImageOutput, pFloatImageInput->nrow,
                        pFloatImageInput->ncol) == NULL)
        mwerror(FATAL, 1,
                "fgrain --> Not enough memory to allocate the output image");
    if((pTree = mw_new_shapes()) == NULL)
        mwerror(FATAL, 1,
                "fgrain --> Not enough memory to allocate the tree of shapes");

    /* Compute the Level Sets Transform of the input image */
    flst(NULL, pFloatImageInput, pTree);

    /* Kill too small grains.
     Bound i>0 because it is forbidden to delete the root, at index 0 */
    for(i = pTree->nb_shapes-1; i > 0; i--)
        if(pTree->the_shapes[i].area < *pMinArea && ( (sideflag >0) ^ pTree->the_shapes[i].inferior_type ))
            pTree->the_shapes[i].removed = (char)1;

    /* Reconstruct in pFloatImageOutput the modified tree of shapes */
    flst_reconstruct(pTree, pFloatImageOutput);

    mw_delete_shapes(pTree);
}


/*in and out must be allocated*/
void AbstractionProcess::fgrain_side(int MinArea, float *in, int nx, int ny, float *out, int sideflag) {
    int i;
    Fimage mwin = mw_change_fimage(NULL,ny,nx);
    Fimage mwout = mw_new_fimage();
    for(i=0;i<nx*ny;i++) mwin->gray[i] = in[i];

    mw_fgrain_side( &MinArea, mwin, mwout, sideflag);

    for(i=0;i<nx*ny;i++) out[i] = mwout->gray[i];
    mw_delete_fimage(mwin);
    mw_delete_fimage(mwout);
}


/*========= Compute the orientation and Elongation of shapes ========*/
void AbstractionProcess::shape_orilam(Shape pShape, float *out_ori, float *out_e, float *out_k)
{
    float size;
    float a11, a20, a02, x0, y0, sumx, sumy, lambda1, lambda2;
    int i,j;

    size = (float)pShape->area;

    sumx = 0;
    sumy = 0;
    for(i = 0; i< size; i++)
    {
        sumx += (double)((pShape->pixels+i)->x);
        sumy += (double)((pShape->pixels+i)->y);
    }
    x0 = sumx /(float)size;
    y0 = sumy /(float)size;

    a11 = 0;
    a20 = 0;
    a02 = 0;
    for(i = 0; i< size; i++)
    {
        a11 += (float) ((pShape->pixels+i)->x - x0)*((pShape->pixels+i)->y - y0);
        a20 += (float) ((pShape->pixels+i)->x - x0)*((pShape->pixels+i)->x - x0);
        a02 += (float) ((pShape->pixels+i)->y - y0)*((pShape->pixels+i)->y - y0);
    }
    a11 = a11 /(float)pow(size,1);
    a20 = a20 /(float)pow(size,1)+ 1.0/12.0;
    a02 = a02 /(float)pow(size,1)+ 1.0/12.0;

    *out_ori = (0.5*atan2(2*a11,(a20-a02)) + PI/2)/PI;

    lambda1 =0.5*( a02 + a20 + sqrt((a20-a02)*(a20-a02) + 4*a11*a11));
    lambda2 =0.5*( a02 + a20 - sqrt((a20-a02)*(a20-a02) + 4*a11*a11));

    *out_e = lambda2/lambda1;

    *out_k = ((float) pShape->area)/(sqrt(lambda2*lambda1)*4*PI);
}


/*========= Compute the orientation and engienvalues of shapes ========*/
void AbstractionProcess::shape_orilam(Shape pShape, float *out_ori, float *out_e, float *out_k, float *pX0, float *pY0)
{
    float size;
    float a11, a20, a02, x0, y0, sumx, sumy, lambda1, lambda2;
    int i,j;

    size = (float)pShape->area;

    sumx = 0;
    sumy = 0;
    for(i = 0; i< size; i++)
    {
        sumx += (double)((pShape->pixels+i)->x);
        sumy += (double)((pShape->pixels+i)->y);
    }
    x0 = sumx /(float)size;
    y0 = sumy /(float)size;

    a11 = 0;
    a20 = 0;
    a02 = 0;
    for(i = 0; i< size; i++)
    {
        a11 += (float) ((pShape->pixels+i)->x - x0)*((pShape->pixels+i)->y - y0);
        a20 += (float) ((pShape->pixels+i)->x - x0)*((pShape->pixels+i)->x - x0);
        a02 += (float) ((pShape->pixels+i)->y - y0)*((pShape->pixels+i)->y - y0);
    }
    a11 = a11 /(float)pow(size,1);
    a20 = a20 /(float)pow(size,1)+ 1.0/12.0;
    a02 = a02 /(float)pow(size,1)+ 1.0/12.0;

    *out_ori = (0.5*atan2(2*a11,(a20-a02)));

    lambda1 =0.5*( a02 + a20 + sqrt((a20-a02)*(a20-a02) + 4*a11*a11));
    lambda2 =0.5*( a02 + a20 - sqrt((a20-a02)*(a20-a02) + 4*a11*a11));

    *out_e = lambda1;
    *out_k = lambda2;

    *pX0 = x0;
    *pY0 = y0;
}

/* sort two shapes according to their scales*/
void AbstractionProcess::Order(Fsignal t2b_index,
                               int *p, int *q)
{
    int temp;
    Shape pShape1, pShape2;
    pShape1 = _pTree->the_shapes + (int) t2b_index->values[*p];
    pShape2 = _pTree->the_shapes + (int) t2b_index->values[*q];

    if( pShape1->area < pShape2->area)
    {
        temp = t2b_index->values[*p];

        t2b_index->values[*p] =  t2b_index->values[*q];
        t2b_index->values[*q] =  temp;
    }
}
/*====== Indexing the mn-order parent of the pShape ===*/
Shape AbstractionProcess::m_order_parent(Shape pShape,
                                         int *mn,
                                         bool dict)
{
    int t;
    Shape pShapeTemp;
    pShapeTemp = pShape;
    for(t=0; t<(*mn); t++)
    {
        if(pShapeTemp->parent == NULL)
            break;

        if( dict ){
            for(pShapeTemp=pShape->parent; ((Info*)(pShapeTemp->data))->show != 1 && pShapeTemp !=NULL; pShapeTemp=pShapeTemp->parent)
            {;
            }
        } else {
            pShapeTemp = pShapeTemp->parent;
        }
    }

    return pShapeTemp;
}

/* sampling the leaf node of the tree according to
   given distribution */
void AbstractionProcess::random_leaf(Fsignal leaf,
                                     Fsignal t2b_index,
                                     int *k_ind)
{
    Shape pShape, pShapeTemp;
    int i, size, select_i;
    float pb_sum, temp;
    Fsignal pb=0;

    size = leaf->size;

    if(size <= 0)
        return;

    srand ( time(NULL) );
    /*   for(i=0; i<size; i++) */
    /*     printf("leaf: %f, %f \n", leaf->values[2*i +0], leaf->values[2*i +1]); */

    select_i = size-1;

    if(size > 1) {
        pb_sum = 0.0;
        pb = mw_change_fsignal(pb,size);
        mw_clear_fsignal(pb,0.0);

        for(i=0; i< size; i++){
            /*// sampling by scale distribution; */
            //pb->values[i] = 1/(pow((float)leaf->values[2*i +1], -1.0));

            /*// sampling by uniform distribution; */
            pb->values[i] = 1;
            pb_sum += pb->values[i];
        }

        for(i=0; i< size; i++){
            pb->values[i] /= pb_sum;
        }

        for(i=1; i< size; i++){
            pb->values[i] += pb->values[i-1];
        }

        temp = ((float)rand())/RAND_MAX;

        for(i= 0; i< size; i++){
            if( temp <= pb->values[i]){
                select_i = i;
                break;
            }
        }
    }

    t2b_index->values[*k_ind] = leaf->values[2*select_i +0];
    //printf("select_i: %d %d %f %f\n", *k_ind, select_i, leaf->values[2*select_i +0], t2b_index->values[*k_ind]);

    pShape = _pTree->the_shapes + (int) leaf->values[2*select_i +0];
    pShapeTemp = pShape->parent;
    if(pShapeTemp != NULL){
        ((Info*)(pShapeTemp->data))->child --;
        if(((Info*)(pShapeTemp->data))->child ==0){
            leaf->values[2*select_i +0] = (float) ((Info*)(pShapeTemp->data))->index;
            leaf->values[2*select_i +1] = (float) pShapeTemp->area;
        }
        else{
            for(i = select_i+1; i< size; i++){
                leaf->values[2*(i-1) +0] = leaf->values[2*(i) +0];
                leaf->values[2*(i-1) +1] = leaf->values[2*(i) +1];
            }
            size = size - 1;
            mw_change_fsignal(leaf,size);
        }
    }
    else{
        size = 0;
        mw_change_fsignal(leaf,size);
    }
}

/* visit the tree in random order */
void AbstractionProcess::random_tree_order(Fsignal t2b_index)
{
    int i,j, a,b, nleaf, k_ind;
    Shape pShape, pShapeTemp;
    Fsignal leaf=0;

    leaf = mw_change_fsignal(leaf, 2*_pTree->nb_shapes);
    mw_clear_fsignal(leaf,0.0);

    mw_change_fsignal(t2b_index, _pTree->nb_shapes);
    mw_clear_fsignal(t2b_index,0.0);

    nleaf = 0;
    for(i=0; i< _pTree->nb_shapes; i++) {
        pShape = _pTree->the_shapes + i;
        ((Info*)(pShape->data))->index = i;
        ((Info*)(pShape->data))->child = 0;

        for(pShapeTemp=pShape->child; pShapeTemp!=NULL;
            pShapeTemp=pShapeTemp->next_sibling){
            ((Info*)(pShape->data))->child += 1;
        }

        if( ((Info*)(pShape->data))->child == 0){
            leaf->values[2*nleaf +0] = (float) ((Info*)(pShape->data))->index;
            leaf->values[2*nleaf +1] = (float)  pShape->area;
            nleaf++;
        }

        //printf("index: %d, %d\n", i, pShape->area);
    }

    mw_change_fsignal(leaf, nleaf);
    k_ind = 0;

    while(k_ind < _pTree->nb_shapes){
        random_leaf(leaf, t2b_index, &k_ind);
        //printf("%d, select: %f\n",pTree->nb_shapes,  t2b_index->values[k_ind]);
        k_ind++;
    }

    mw_change_fsignal(leaf, _pTree->nb_shapes);
    mw_clear_fsignal(leaf,0.0);
    for(i=0; i< _pTree->nb_shapes; i++){
        leaf->values[i] = t2b_index->values[i];
    }
    for(i=0; i< _pTree->nb_shapes; i++){
        t2b_index->values[i] =  leaf->values[_pTree->nb_shapes-1-i];
        //printf("%d, select: %f\n",i, t2b_index->values[i]);
    }

    mw_delete_fsignal(leaf);
}

/*Indexing the tree by the Breadth-first order*/
void AbstractionProcess::top2bottom_index_tree(Fsignal t2b_index)
{
    int queArr[_pTree->nb_shapes];
    int i, qInd;
    int *head, *rear;
    Shape pShape, pShapeTemp;

    mw_change_fsignal(t2b_index, _pTree->nb_shapes);

    for(i=0; i< _pTree->nb_shapes; i++) {
        pShape = _pTree->the_shapes + i;
        ((Info*)(pShape->data))->index = i;
        queArr[i] = -i;
    }

    //printf("---top 1----!\n");

    head = queArr;
    rear = queArr;
    queArr[0] = 0;
    rear++;
    t2b_index->values[0] = 0;

    i=0; qInd =1;

    while(head != rear){
        t2b_index->values[i++] = (float) *head;
        pShape = _pTree->the_shapes + *head;
        //printf("%d %d \n", *head, *rear);
        head++;

        for(pShapeTemp = pShape->child; pShapeTemp != NULL;
            pShapeTemp = pShapeTemp->next_sibling){
            queArr[qInd++] = (((Info*)(pShapeTemp->data))->index);
            rear++;
        }
    }
}

/*===== Compute the mean contrast of the curve l =====*/
float AbstractionProcess::mean_contrast(Shape pShape)
{
    double per;
    float mu,meanmu,x,y,ox,oy;
    int i,ix,iy;

    Flist pBoundary = NULL;
    pBoundary = mw_change_flist(pBoundary, 4*pShape->area+1, 0, 2);
    flst_boundary(_pTree, pShape, pBoundary);

    per = 0.;
    meanmu = 0.0;
    /*   meanmu = FLT_MAX; */

    for (i=0; i<pBoundary->size;i++) {

        x = pBoundary->values[i*2];
        y = pBoundary->values[i*2+1];

        if (i>0) per += sqrt((double)(x-ox)*(x-ox)+(y-oy)*(y-oy));
        ox = x; oy = y;

        ix = (int)rint((double)x)-1;
        iy = (int)rint((double)y)-1;
        if (ix>=0 && iy>=0 && ix<_NormOfDu->ncol && iy<_NormOfDu->nrow) {
            mu = _NormOfDu->gray[_NormOfDu->ncol*iy+ix];
            meanmu += mu;
            /*       if (mu<meanmu) meanmu=mu; */
        }
    }
    meanmu /= pBoundary->size;
    if (meanmu == FLT_MAX) meanmu = 0.;

    free(pBoundary->data);
    mw_delete_flist(pBoundary);
    return(meanmu);
}

void AbstractionProcess::adaptive_shift_shape(float *shift,
                                              float *theta)
{
    int i,j;
    Shape pShape;
    float SHIFT, THETA, tempx, tempy, tempt, tempShiftx, tempShifty;
    float a,b,phi, CONTRAST;
    SHIFT = *shift;
    THETA = *theta;

    srand( time(NULL) );
    for(i = _pTree->nb_shapes-1; i>= 0; i--)
    {
        pShape = _pTree->the_shapes + i;

        if(pShape == NULL )
            continue;

        if(i==0)
        {
            ((Info*)(pShape->data))->xShift = 0;
            ((Info*)(pShape->data))->yShift = 0;
            ((Info*)(pShape->data))->rotation = 0;
            continue;
        }

        if( (rand()%10) >5 )
            tempx = -((float)rand())/RAND_MAX;
        else
            tempx = ((float)rand())/RAND_MAX;

        if( (rand()%10) >5 )
            tempy = -((float)rand())/RAND_MAX;
        else
            tempy = ((float)rand())/RAND_MAX;

        if( (rand()%10) >5 )
            tempt = -((float)rand())/RAND_MAX;
        else
            tempt = ((float)rand())/RAND_MAX;

        a = 2.0 * sqrt(((Info*)(pShape->data))->lambda1);
        b = 2.0 * sqrt(((Info*)(pShape->data))->lambda2);
        phi = ((Info*)(pShape->data))->oren;

        tempShiftx = tempx*SHIFT;
        tempShifty = tempy*SHIFT*pow((b/a),2);

        ((Info*)(pShape->data))->xShift = tempShiftx*cos(phi) + tempShifty*sin(phi);
        ((Info*)(pShape->data))->yShift = tempShifty*cos(phi) - tempShiftx*sin(phi);

        ((Info*)(pShape->data))->rotation = tempt*THETA*PI*(b/a);

        CONTRAST = mean_contrast(pShape);

        ((Info*)(pShape->data))->xShift *= 1/pow(CONTRAST, 0.5);
        ((Info*)(pShape->data))->yShift *= 1/pow(CONTRAST, 0.5);
        ((Info*)(pShape->data))->rotation *= 1/pow(CONTRAST, 0.5);

    }
}


void AbstractionProcess::random_shift_shape(float *shift)
{
    int i,j;
    Shape pShape;
    float SHIFT, tempx, tempy;
    SHIFT = *shift;

    srand( time(NULL) );
    for(i = _pTree->nb_shapes-1; i>= 0; i--)
    {
        pShape = _pTree->the_shapes + i;

        if(pShape == NULL )
            continue;
        if(i==0)
        {
            ((Info*)(pShape->data))->xShift = 0;
            ((Info*)(pShape->data))->yShift = 0;
            ((Info*)(pShape->data))->rotation = 0;
            continue;
        }

        if( (rand()%10) >5 )
            tempx = -((float)rand())/RAND_MAX;
        else
            tempx = ((float)rand())/RAND_MAX;

        if( (rand()%10) >5 )
            tempy = -((float)rand())/RAND_MAX;
        else
            tempy = ((float)rand())/RAND_MAX;

        ((Info*)(pShape->data))->xShift = tempx*SHIFT;
        ((Info*)(pShape->data))->yShift = tempy*SHIFT;

        /*       printf("temp: %f, shift: %f \n", ((Info*)(pShape->data))->xShift,  */
        /*       ((Info*)(pShape->data))->yShift);*/
    }
}

/*===== Compute the mean contrast of the curve l =====*/
float AbstractionProcess::min_contrast(Shape pShape)
{
    double per;
    float mu,meanmu,x,y,ox,oy;
    int i,ix,iy;

    Flist pBoundary = NULL;
    pBoundary = mw_change_flist(pBoundary, 4*pShape->area+1, 0, 2);
    flst_boundary(_pTree, pShape, pBoundary);

    per = 0.;
    meanmu = FLT_MAX;

    for(i=0; i<pBoundary->size;i++)
    {

        x = pBoundary->values[i*2];
        y = pBoundary->values[i*2+1];

        if (i>0) per += sqrt((double)(x-ox)*(x-ox)+(y-oy)*(y-oy));
        ox = x; oy = y;

        ix = (int)rint((double)x)-1;
        iy = (int)rint((double)y)-1;
        if (ix>=0 && iy>=0 && ix<_NormOfDu->ncol && iy<_NormOfDu->nrow)
        {
            mu =_NormOfDu->gray[_NormOfDu->ncol*iy+ix];
            if (mu<meanmu) meanmu=mu;
        }
    }
    if (meanmu == FLT_MAX) meanmu = 0.;

    free(pBoundary->data);
    mw_delete_flist(pBoundary);
    return(meanmu);
}


/*====== Compute shape attribute ===*/
void AbstractionProcess::compute_shape_attribute()
{
    int i;
    float oren, lamb1, lamb2, x0, y0;
    Shape pShape;

    for(i = _pTree->nb_shapes-1; i>=0; i--)  {
        pShape = _pTree->the_shapes + i;

        shape_orilam(pShape, &oren, &lamb1, &lamb2, &x0, &y0);

        ((Info*)(pShape->data))->lambda1 = lamb1;
        ((Info*)(pShape->data))->lambda2 = lamb2;
        ((Info*)(pShape->data))->oren = oren;

        ((Info*)(pShape->data))->x0 = x0;
        ((Info*)(pShape->data))->y0 = y0;

        pShape->removed = 0;

        if(i != 0)
            ((Info*)(pShape->data))->contrast = fabs(min_contrast(pShape));
    }
}

/*================ Compute the perimeter of one shape ================== */
float AbstractionProcess::peri_shape(Shape pShape)
{
    float fPerimeter;
    int nr, nc, size, i;
    Flist pBoundary = NULL;


    nr = _pTree->nrow;
    nc = _pTree->ncol;
    size = nr*nc;

    if (pShape->area >= size)
        return((float) 2*(nr+nc));

    if (!pShape)
    {
        printf("shape is empty\n");
    }

    pBoundary = mw_change_flist(pBoundary, 4*pShape->area+1, 0, 2);

    flst_boundary(_pTree, pShape, pBoundary);
    fPerimeter = (float)pBoundary->size;

    free(pBoundary->data);
    mw_delete_flist(pBoundary);

    return (fPerimeter);
}


/*====== Compute shape attribute ===*/
void AbstractionProcess::compute_shape_attribute(int *ns)
{
    int i, j, kn, nb, nn, nsize, num;
    float oren, elg, kap, stdtemp, meantemp, mmin, mmax;
    int *pKN;
    Shape pShape, pShapeTemp, pShapeCh;

    nn = *ns; nsize = 3;

    for(i = _pTree->nb_shapes-1; i>=0; i--)  {
        pShape = _pTree->the_shapes + i;

        shape_orilam(pShape, &oren, &elg, &kap);
        pShapeTemp = m_order_parent(pShape, &nn);
        ((Info*)(pShape->data))->attribute[0] = ((float) pShape->area)/((float) pShapeTemp->area);

        if( ((Info*)(pShape->data))->attribute[0] <= ((float) pShape->area)/((float) pShapeTemp->area)){
            ((Info*)(pShape->data))->attribute[0] = ((float) pShape->area)/((float) pShapeTemp->area);
            //std::cout << ((Info*)(pShape->data))->attribute[0] << " it happens 2222222 ????  "<<((float) pShape->area)/((float) pShapeTemp->area)<< std::endl;

        }
        if( ((Info*)(pShapeTemp->data))->attribute[0] <= ((float) pShape->area)/((float) pShapeTemp->area)){
            ((Info*)(pShapeTemp->data))->attribute[0] = ((float) pShape->area)/((float) pShapeTemp->area);
            //std::cout <<((Info*)(pShapeTemp->data))->attribute[0] << "it happens 2222222 ????  "<<((float) pShape->area)/((float) pShapeTemp->area )<< std::endl;
        }

        ((Info*)(pShape->data))->attribute[1] = kap;
        ((Info*)(pShape->data))->attribute[2] = elg;
        ((Info*)(pShape->data))->attribute[3] = oren;
    }
}

/*=== Synthesis a shape using a rectangle with the same parameters (second-order moments)===*/
void AbstractionProcess::synshapeRect(Shape pShape,
                                      Ccimage imgsyn,
                                      float *alpha,
                                      int *relief,
                                      float *reliefOrentation, float *reliefHeight)
{
    int xi, yi;
    float a, b, x0temp, y0temp, top, right, left, bottom, bLimit, ALPHA;
    float phi, xi_e, yi_e;
    float xShift, yShift, theta, tR, tG, tB, tr, tg, tb, TR, TG, TB;

    ALPHA = *alpha;
    x0temp = (((Info*)(pShape->data))->x0);
    y0temp = (((Info*)(pShape->data))->y0);

    xShift = (((Info*)(pShape->data))->xShift);
    yShift = (((Info*)(pShape->data))->yShift);
    theta  = (((Info*)(pShape->data))->rotation);

    x0temp += xShift;
    y0temp += yShift;

    a = (sqrt(3.0 * ((Info*)(pShape->data))->lambda1));
    b = (sqrt(3.0 * ((Info*)(pShape->data))->lambda2));
    phi = ((Info*)(pShape->data))->oren;

    bLimit = sqrt(2.0)*a;

    left   = _MAX(0, x0temp - bLimit);
    right  = _MIN(_pTree->ncol -1, x0temp + bLimit);
    top    = _MAX(0, y0temp - bLimit);
    bottom = _MIN(_pTree->nrow - 1, y0temp + bLimit);

    TR  = ((Info*)(pShape->data))->r;
    TG  = ((Info*)(pShape->data))->g;
    TB  = ((Info*)(pShape->data))->b;


    /*=== Synthesis ===*/
    if(*relief == 1 && pShape->area > 10)
    {
        float shLambda, shTR, shTG, shTB;
        int xsh, ysh, shiftsh;
        if(pShape->area > 10)
            shiftsh = *reliefHeight;
        //   shiftsh = 3*(1 - 1/((Info*)(pShape->data))->contrast);
        else
            shiftsh = (*reliefHeight)*( (float) pShape->area /10.0);
        //       shiftsh = pow(pShape->area, 0.25);
        //     shiftsh = pow(pShape->area, 0.25)*(1 - 1/((Info*)(pShape->data))->contrast);

        //     printf("%f %d \n", ((Info*)(pShape->data))->contrast, shiftsh);
        shLambda = 0.3;
        for(xi = ceil(left); xi <= right; xi++)
            for(yi = ceil(top); yi <= bottom; yi++)
            {

                xi_e = ((float)xi - x0temp)*cos(phi+theta) + ((float)yi - y0temp)*sin(phi+theta);
                yi_e = ((float)yi - y0temp)*cos(phi+theta) - ((float)xi - x0temp)*sin(phi+theta);

                if( xi_e >= -a && xi_e <= +a && yi_e >= -b && yi_e <= +b )
                {
                    tr = TR* shLambda;
                    tg = TG* shLambda;
                    tb = TB* shLambda;


                    xsh = xi + shiftsh*cos( PI*(*reliefOrentation)/180.0 );
                    ysh = yi - shiftsh*sin( PI*(*reliefOrentation)/180.0 );
                    xsh = _MAX(0, xsh);
                    xsh = _MIN(imgsyn->ncol - 1, xsh);
                    ysh = _MAX(0, ysh);
                    ysh = _MIN(imgsyn->nrow - 1, ysh);

                    tR = ((float) imgsyn->red[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*tr;
                    imgsyn->red[ysh*imgsyn->ncol + xsh]   = (int) tR;/*(int)rint((double) tR); */

                    tG = ((float) imgsyn->green[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*tg;
                    imgsyn->green[ysh*imgsyn->ncol + xsh] = (int) tG; /*(int)rint((double) tG); */

                    tB = ((float) imgsyn->blue[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*tb;
                    imgsyn->blue[ysh*imgsyn->ncol + xsh]  = (int) tB; /*(int)rint((double) tB);*/
                }
            }
    }

    for( xi= ceil(left); xi<= right; xi++)
        for( yi= ceil(top); yi<= bottom; yi++)
        {
            xi_e = ((float)xi - x0temp)*cos(phi+theta) + ((float)yi - y0temp)*sin(phi+theta);
            yi_e = ((float)yi - y0temp)*cos(phi+theta) - ((float)xi - x0temp)*sin(phi+theta);

            if( xi_e >= -a && xi_e <= +a && yi_e >= -b && yi_e <= +b )
            {
                tR = ((float) imgsyn->red[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->r;
                imgsyn->red[yi*_pTree->ncol + xi] = (int) tR; /*rint((double) tR);*/

                tG = ((float) imgsyn->green[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->g;
                imgsyn->green[yi*_pTree->ncol + xi] = (int) tG; /*(int)rint((double) tG); */

                tB = ((float) imgsyn->blue[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->b;
                imgsyn->blue[yi*_pTree->ncol + xi] = (int) tB; /*(int)rint((double) tB); */
            }
        }
}

/*=== Synthesis a shape using a ellipse with the same parameters (second-order moments) ===*/
void AbstractionProcess::synshapeEllipse(Shape pShape,
                                         Ccimage imgsyn,
                                         float *alpha,
                                         int *relief,
                                         float *reliefOrentation, float *reliefHeight)
{
    int xi, yi;
    float a, b, x0temp, y0temp, top, right, left, bottom, ALPHA;
    float phi, xi_e, yi_e;
    float xShift, yShift, theta, tR, tG, tB, TR, TG, TB;

    ALPHA = *alpha;
    x0temp = (((Info*)(pShape->data))->x0);
    y0temp = (((Info*)(pShape->data))->y0);

    xShift = (((Info*)(pShape->data))->xShift);
    yShift = (((Info*)(pShape->data))->yShift);
    theta  = (((Info*)(pShape->data))->rotation);

    x0temp += xShift;
    y0temp += yShift;

    a = 2.0 * sqrt(((Info*)(pShape->data))->lambda1);
    b = 2.0 * sqrt(((Info*)(pShape->data))->lambda2);
    phi = ((Info*)(pShape->data))->oren;

    left   = _MAX(0, x0temp - a);
    right  = _MIN(_pTree->ncol -1, x0temp + a);
    top    = _MAX(0, y0temp - a);
    bottom = _MIN(_pTree->nrow - 1, y0temp + a);

    TR  = ((Info*)(pShape->data))->r;
    TG  = ((Info*)(pShape->data))->g;
    TB  = ((Info*)(pShape->data))->b;

    /*=== Synthesis ===*/
    if(*relief == 1 && pShape->area > 10)
    {
        float shLambda, shTR, shTG, shTB;
        int xsh, ysh, shiftsh;
        if(pShape->area > 10)
            shiftsh = *reliefHeight;
        else
            shiftsh = (*reliefHeight)*( (float) pShape->area /10.0);

        shLambda = 0.3;
        for( xi= ceil(left); xi<= right; xi++)
            for( yi= ceil(top); yi<= bottom; yi++)
            {
                xi_e = ((float)xi - x0temp)*cos(phi+theta) + ((float)yi - y0temp)*sin(phi+theta);
                yi_e = ((float)yi - y0temp)*cos(phi+theta) - ((float)xi - x0temp)*sin(phi+theta);

                if( xi_e*xi_e/(a*a) + yi_e*yi_e/(b*b) <= 1 ){

                    if(xi<0 || xi>= imgsyn->ncol ||
                            yi<0 || yi>= imgsyn->nrow )
                        continue;

                    shTR = TR* shLambda;
                    shTG = TG* shLambda;
                    shTB = TB* shLambda;


                    xsh = xi + shiftsh*cos( PI*(*reliefOrentation)/180.0 );
                    ysh = yi - shiftsh*sin( PI*(*reliefOrentation)/180.0 );
                    xsh = _MAX(0, xsh);
                    xsh = _MIN(imgsyn->ncol - 1, xsh);
                    ysh = _MAX(0, ysh);
                    ysh = _MIN(imgsyn->nrow - 1, ysh);

                    tR = ((float) imgsyn->red[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTR;
                    imgsyn->red[ysh*imgsyn->ncol + xsh]   = (int) tR;/*(int)rint((double) tR); */

                    tG = ((float) imgsyn->green[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTG;
                    imgsyn->green[ysh*imgsyn->ncol + xsh] = (int) tG; /*(int)rint((double) tG); */

                    tB = ((float) imgsyn->blue[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTB;
                    imgsyn->blue[ysh*imgsyn->ncol + xsh]  = (int) tB; /*(int)rint((double) tB);*/
                }
            }
    }

    for( xi= ceil(left); xi<= right; xi++)
        for( yi= ceil(top); yi<= bottom; yi++)
        {
            xi_e = ((float)xi - x0temp)*cos(phi+theta) + ((float)yi - y0temp)*sin(phi+theta);
            yi_e = ((float)yi - y0temp)*cos(phi+theta) - ((float)xi - x0temp)*sin(phi+theta);

            if( xi_e*xi_e/(a*a) + yi_e*yi_e/(b*b) <= 1 ){

                tR = ((float) imgsyn->red[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->r;
                imgsyn->red[yi*_pTree->ncol + xi] = (int) tR; /*(int)rint((double) tR); */

                tG = ((float) imgsyn->green[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->g;
                imgsyn->green[yi*_pTree->ncol + xi] = (int) tG; /*(int)rint((double) tG); */

                tB = ((float) imgsyn->blue[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->b;
                imgsyn->blue[yi*_pTree->ncol + xi] = (int) tB; /*(int)rint((double) tB); */
            }
        }
}

/*=== Synthesis by Shape Shaking ===*/
/*=== Before the Shaking, smooth the shape with a gaussian kernel or a median filter ===*/
void AbstractionProcess::synshapeOriginal(Shape pShape,
                                          Ccimage imgsyn,
                                          Cimage imgShapeLabelSyn,
                                          Fimage imgShapeBlurSyn,
                                          Fsignal gaussKernel,
                                          int *median,
                                          float *alpha,
                                          int *relief,
                                          float *reliefOrentation, float *reliefHeight)
{
    int i, xi, yi, x, y, iKer, jKer, KerSize, MedSize, xKer, yKer, numMedain;
    float xr, yr, x0temp, y0temp, ALPHA, BETA, a, b;
    float xShift, yShift, theta, tR, tG, tB, tr, tg, tb, TR, TG, TB;
    float minX, maxX, minY, maxY;
    float top, right, left, bottom;

    ALPHA = *alpha;
    x0temp = (((Info*)(pShape->data))->x0);
    y0temp = (((Info*)(pShape->data))->y0);

    xShift = (((Info*)(pShape->data))->xShift);
    yShift = (((Info*)(pShape->data))->yShift);
    theta  = (((Info*)(pShape->data))->rotation);

    TR  = ((Info*)(pShape->data))->r;
    TG  = ((Info*)(pShape->data))->g;
    TB  = ((Info*)(pShape->data))->b;

    right = 0.0;  bottom= 0.0;
    left = (float)(imgShapeLabelSyn->ncol-1);
    top  = (float)(imgShapeLabelSyn->nrow-1);
    for(i=0; i< pShape->area; i++)
    {
        x = (pShape->pixels+i)->x;
        y = (pShape->pixels+i)->y;

        xr = (x - x0temp)*cos(theta) + (y - y0temp)*sin(theta);
        yr = (y - y0temp)*cos(theta) - (x - x0temp)*sin(theta);

        xi = floor(xShift + x0temp + xr);
        yi = floor(yShift + y0temp + yr);

        if(xi<0 || xi>= imgShapeLabelSyn->ncol ||
                yi<0 || yi>= imgShapeLabelSyn->nrow )
            continue;

        imgShapeLabelSyn->gray[yi*imgShapeLabelSyn->ncol + xi] = 1;

        left   = _MIN(xi, left);
        top    = _MIN(yi, top);
        right  = _MAX(xi, right);
        bottom = _MAX(yi, bottom);
    }

    /*=== Median Filter ===*/
    MedSize = (int)((*median)/2.0);
    for(x = left; x <= right; x++)
        for(y = top; y <= bottom; y++)
        {
            numMedain = 0;
            for(iKer = - MedSize; iKer <= MedSize; iKer++)
                for(jKer = - MedSize; jKer <= MedSize; jKer++)
                {
                    xKer = x + iKer;
                    yKer = y + jKer;

                    if(xKer<0 || xKer>= imgShapeLabelSyn->ncol ||
                            yKer<0 || yKer>= imgShapeLabelSyn->nrow )
                        continue;

                    imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] +=
                            (float)(imgShapeLabelSyn->gray[yKer*imgShapeLabelSyn->ncol + xKer]);
                    numMedain++;
                }
            if( imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] < ((float) numMedain)/2.0 )
                imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] = 0.0;
            else
                imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] = 1.0;
        }

    for(x = left; x <= right; x++)
        for(y = top; y <= bottom; y++)
        {
            imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] =
                    (int) imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x];

            imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] = 0.0;
        }

    /*=== Add Gaussian Bulr ===*/
    KerSize = (int) ( sqrt( (double) gaussKernel->size) /2.0 );
    for(x = left; x <= right; x++)
        for(y = top; y <= bottom; y++)
        {
            for(iKer = -KerSize; iKer <= KerSize; iKer++)
                for(jKer = -KerSize; jKer <= KerSize; jKer++)
                {
                    xKer = x + iKer;
                    yKer = y + jKer;

                    if(xKer<0 || xKer>= imgShapeLabelSyn->ncol ||
                            yKer<0 || yKer>= imgShapeLabelSyn->nrow )
                        continue;

                    imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] +=
                            gaussKernel->values[(iKer + KerSize)*KerSize + (jKer + KerSize)]*
                            (float)(imgShapeLabelSyn->gray[yKer*imgShapeLabelSyn->ncol + xKer]);
                }
        }

    /*=== Synthesis ===*/
    if(*relief == 1 && pShape->area > 10)
    {
        float shLambda, shTR, shTG, shTB;
        int xsh, ysh, shiftsh;
        if(pShape->area > 10)
            shiftsh = *reliefHeight;
        //   shiftsh = 3*(1 - 1/((Info*)(pShape->data))->contrast);
        else
            shiftsh = (*reliefHeight)*( (float) pShape->area /10.0);
        //       shiftsh = pow(pShape->area, 0.25);
        //     shiftsh = pow(pShape->area, 0.25)*(1 - 1/((Info*)(pShape->data))->contrast);

        //     printf("%f %d \n", ((Info*)(pShape->data))->contrast, shiftsh);
        shLambda = 0.3;
        for(x = ceil(left); x <= right; x++)
            for(y = ceil(top); y <= bottom; y++)
            {
                if(imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] == 0)
                    continue;

                BETA = imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x];

                shTR = TR* shLambda;
                shTG = TG* shLambda;
                shTB = TB* shLambda;


                xsh = x + shiftsh*cos( PI*(*reliefOrentation)/180.0 );
                ysh = y - shiftsh*sin( PI*(*reliefOrentation)/180.0 );
                xsh = _MAX(0, xsh);
                xsh = _MIN(imgsyn->ncol - 1, xsh);
                ysh = _MAX(0, ysh);
                ysh = _MIN(imgsyn->nrow - 1, ysh);

                tr = ((float) imgsyn->red[ysh*imgsyn->ncol + xsh])*(1-BETA)   + BETA*shTR;
                tg = ((float) imgsyn->green[ysh*imgsyn->ncol + xsh])*(1-BETA) + BETA*shTG;
                tb = ((float) imgsyn->blue[ysh*imgsyn->ncol + xsh])*(1-BETA)  + BETA*shTB;

                tR = ((float) imgsyn->red[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*tr;
                imgsyn->red[ysh*imgsyn->ncol + xsh]   = (int) tR;/*(int)rint((double) tR); */

                tG = ((float) imgsyn->green[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*tg;
                imgsyn->green[ysh*imgsyn->ncol + xsh] = (int) tG; /*(int)rint((double) tG); */

                tB = ((float) imgsyn->blue[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*tb;
                imgsyn->blue[ysh*imgsyn->ncol + xsh]  = (int) tB; /*(int)rint((double) tB);*/

            }
    }

    for(x = left; x <= right; x++)
        for(y = top; y <= bottom; y++)
        {
            if(imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] == 0)
                continue;

            BETA = imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x];

            tr = ((float) imgsyn->red[y*imgsyn->ncol + x])*(1-BETA)   + BETA*TR;
            tg = ((float) imgsyn->green[y*imgsyn->ncol + x])*(1-BETA) + BETA*TG;
            tb = ((float) imgsyn->blue[y*imgsyn->ncol + x])*(1-BETA)  + BETA*TB;

            tR = ((float) imgsyn->red[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tr;
            imgsyn->red[y*imgsyn->ncol + x]   = (int) tR;/*(int)rint((double) tR); */

            tG = ((float) imgsyn->green[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tg;
            imgsyn->green[y*imgsyn->ncol + x] = (int) tG; /*(int)rint((double) tG); */

            tB = ((float) imgsyn->blue[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tb;
            imgsyn->blue[y*imgsyn->ncol + x]  = (int) tB; /*(int)rint((double) tB);*/

            imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x]  = 0.0;
            imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] = 0;
        }
}

void AbstractionProcess::synshapeOriginal( Shape pShape,
                                           Ccimage imgsyn,
                                           float *alpha)
{
    int i, xi, yi;
    float x, y, xr, yr, x0temp, y0temp, ALPHA;
    float xShift, yShift, theta, tR, tG, tB;

    ALPHA = *alpha;
    x0temp = (((Info*)(pShape->data))->x0);
    y0temp = (((Info*)(pShape->data))->y0);

    xShift = (((Info*)(pShape->data))->xShift);
    yShift = (((Info*)(pShape->data))->yShift);
    theta  = (((Info*)(pShape->data))->rotation);

    for(i=0; i< pShape->area; i++){

        x = (float)((pShape->pixels+i)->x);
        y = (float)((pShape->pixels+i)->y);

        xr = (x - x0temp)*cos(theta) + (y - y0temp)*sin(theta);
        yr = (y - y0temp)*cos(theta) - (x - x0temp)*sin(theta);

        xi = floor(xShift + x0temp + xr);
        yi = floor(yShift + y0temp + yr);

        if(xi>= 0 && xi< _pTree->ncol &&
                yi>= 0 && yi< _pTree->nrow){
            tR = ((float) imgsyn->red[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->r;
            imgsyn->red[yi*_pTree->ncol + xi] = (int)rint((double) tR);

            tG = ((float) imgsyn->green[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->g;
            imgsyn->green[yi*_pTree->ncol + xi] = (int)rint((double) tG);

            tB = ((float) imgsyn->blue[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->b;
            imgsyn->blue[yi*_pTree->ncol + xi] = (int)rint((double) tB);
        }
    }
}

/*Index the tree in coast-to-fine order*/
/*Buuble sorting of integer array */
void AbstractionProcess::Bubble( Fsignal t2b_index)
{
#if 1
    int i,j, a, b;
    Shape pShape;

    mw_change_fsignal(t2b_index, _pTree->nb_shapes);

    for(i=0; i< _pTree->nb_shapes; i++)
    {
        pShape = _pTree->the_shapes + i;
        ((Info*)(pShape->data))->index = i;
        t2b_index->values[i] = i;
    }

    for (i=0; i< _pTree->nb_shapes; i++)
        for (j = _pTree->nb_shapes-1; i<j; j--)
        {
            a = j-1; b = j;
            Order(t2b_index, &a, &b);
        }
#else
    int i,j, a, b;
    Shape pShape;
    mw_change_fsignal(t2b_index, pTree->nb_shapes);

    for(i=0; i< pTree->nb_shapes; i++)
    {
        pShape = pTree->the_shapes + i;
        ((Info*)(pShape->data))->index = i;
        t2b_index->values[i] = i;
    }

    for (i = 0; i< pTree->nb_shapes-1; i++)
        for (j = i+1; j< pTree->nb_shapes; j++)
        {
            a = i; b = j;
            Order_v4(pTree, t2b_index, &a, &b);
        }

    for(i=0; i< pTree->nb_shapes; i++)
    {
        pShape = pTree->the_shapes + (int) t2b_index->values[i];
        printf("%d \n", pShape->area);
    }
#endif
}

/*=== Generate a 1D Gaussian kernel ===*/
Fsignal AbstractionProcess::sgauss(float *std,
                                   Fsignal out,
                                   int *size)
{
    int i,n;
    double sum, v;

    n = *size;
    v = (0.5*(double)(n-1)) / (double)(*std);
    v = 0.5*v*v/log(10.);

    out = mw_change_fsignal(out, n);
    if (!out) mwerror(FATAL,1,"Not enough memory.");

    out->shift = -0.5*(float)(n-1);

    if (n==1)
    {
        out->values[0]=1.0;
    }
    else
    {
        /*=== store Gaussian signal ===*/
        for(i=(n+1)/2; i--; )
        {
            v = ((double)i+(double)out->shift)/(double)(*std);
            out->values[i] = out->values[n-1-i] = (float)exp(-0.5*v*v);
        }
        /*=== normalize to get unit mass ===*/
        for (sum=0.0,i=n; i--; )
            sum += (double)out->values[i];
        for (i=n; i--; )
        {
            out->values[i] /= (float)sum;
        }
    }

    return(out);
}

/*=== Generate a 2D Gaussian kernel ===*/
Fsignal AbstractionProcess::Sgauss(float *std, Fsignal out, int *size)
{
    Fsignal sgaussX, sgaussY;
    int i, j, n;
    float sum;
    n = *size;

    if  ( ((sgaussX = mw_new_fsignal()) == NULL) ||
          (mw_alloc_fsignal(sgaussX, n) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");
    if  ( ((sgaussY = mw_new_fsignal()) == NULL) ||
          (mw_alloc_fsignal(sgaussY, n) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");

    sgaussX = sgauss(std, sgaussX, size);
    sgaussY = sgauss(std, sgaussY, size);

    out = mw_change_fsignal(out, n*n);
    for (i = 0; i<n; i++)
        for (j = i; j<n; j++)
        {
            out->values[j*n + i] = sgaussX->values[i] * sgaussY->values[j];
            out->values[i*n + j] = out->values[j*n + i];
        }
    sum = 0.0;
    for (i = 0; i<n*n; i++)
        sum += out->values[i];
    for (i = 0; i<n*n; i++)
        out->values[i] /= sum;

    mw_delete_fsignal(sgaussX);
    mw_delete_fsignal(sgaussY);

    return(out);
}

void AbstractionProcess::get_shapes_truearea(Shape s, Shape root,
                                             int *truearea)
{
    int index;
    Shape t;

    index = s-root;
    truearea[index] = s->area;
    t = mw_get_first_child_shape(s);
    while(t){
        get_shapes_truearea(t,root,truearea);
        truearea[index] -= t->area;
        t = mw_get_next_sibling_shape(t);
    }
}


void AbstractionProcess::filter_shapes( Fimage sketch,
                                        Cfimage out,
                                        char *local,
                                        float *eps)
{
    std::cout <<"AbstractionProcess::filter_shapes()::begin"<< std::endl;
    Fimage Fv;
    float lgeo,step,hstep,tcleanup,farea, std;
    float *gray,*red,*green,*blue;
    int prec,visit,nrow,ncol,i,j,indexshape,*truearea;
    char all;
    Shapes tree;
    Shape s,t;

    tree = mw_new_shapes();

    Fv = mw_change_fimage(NULL,_imgin->nrow,_imgin->ncol);
    sketch = mw_change_fimage(NULL,_imgin->nrow,_imgin->ncol);

    if(!(Fv && sketch)) mwerror(FATAL,1,"Not enough memory.\n");


    /*   eps = 0.; lgeo = 10.; step = 0.00392; prec = 2; hstep = 0.0000392; */
    /* eps = 15.;  */
    lgeo = 10.; step = 1.; prec = 2; hstep = 0.01; std=0.5;
    tcleanup = 1.; all=(char) 1; visit = 100;

    nrow = _imgin->nrow;
    ncol = _imgin->ncol;
    for(i=0;i<ncol*nrow;i++)
        Fv->gray[i] = (_imgin->red[i]+_imgin->green[i]+_imgin->blue[i])/3.;

    /*   if(boundaries) */
    /*     mw_copy_flists(ll_boundaries2(Fv,&eps,NULL,&step,&prec,&std,&hstep,NULL,&visit,local,NULL,tree),boundaries); */
    /* /\*     mw_copy_flists(lll_bdha_multiscale(Fv,&eps,nscale,NULL,&lgeo,NULL,&step, *\/ */
    /* /\* 			     &prec,&hstep,NULL,&visit,local,NULL, *\/ */
    /* /\* 			     NULL,tree,NULL,&tcleanup),boundaries); *\/ */
    /*   else    */
    ll_boundaries2(Fv,eps,NULL,&step,&prec,&std,&hstep,NULL,&visit,local,NULL,tree);

    //        lll_bdha_multiscale(Fv,&eps,nscale,NULL,&lgeo,NULL,&step,
    //              &prec,&hstep,NULL,&visit,local,NULL,NULL,
    //             tree,NULL,&tcleanup);

    /*   ll_boundaries2(_imgin,eps,tree,step,prec,std,hstep,all,visit,loc,image_out,keep_tree); */
    /*   lll_bdha_multiscale(_imgin,eps,nscale,reghisto,lgeo,tree,step, */
    /* 			   prec,hstep,all,visit,loc,gc_only,image_out, */
    /* 		      keep_tree,cleanup,tcleanup); */

    /* compute recursively integral of
     gray level
     saturation
     cos and sin of hue */
    red = (float*)calloc(tree->nb_shapes,sizeof(float));
    green = (float*)calloc(tree->nb_shapes,sizeof(float));
    blue = (float*)calloc(tree->nb_shapes,sizeof(float));
    gray = (float*)calloc(tree->nb_shapes,sizeof(float));
    if(!(gray && red  && green && blue))
        mwerror(FATAL,1,"Not enough memory.\n");


    /* integrate grey level saturation and hue */
    for(i=0;i<ncol;i++)
        for(j=0;j<nrow;j++){
            s = mw_get_smallest_shape(tree,i,j);
            indexshape = s-tree->the_shapes;
            gray[indexshape] += Fv->gray[i+j*ncol];
            red[indexshape] += _imgin->red[i+j*ncol];
            green[indexshape] += _imgin->green[i+j*ncol];
            blue[indexshape] += _imgin->blue[i+j*ncol];

        }

    /* recursively compute area of meaningful shapes when holes are removed */
    truearea = (int*)malloc(sizeof(int)*tree->nb_shapes);
    if(!truearea)     mwerror(FATAL,1,"Not enough memory.\n");
    get_shapes_truearea(tree->the_shapes,tree->the_shapes,truearea);

    for(i=0;i<tree->nb_shapes;i++){
        s = tree->the_shapes+i;
        farea = (float) (s->area);
        t = mw_get_first_child_shape(s);
        while(t){
            farea -= (float)t->area;
            t = mw_get_next_sibling_shape(t);
        }

        gray[i] /= (float)truearea[i];
        red[i] /= (float)truearea[i];
        green[i] /= (float)truearea[i];
        blue[i] /= (float)truearea[i];
    }

    /* replace gray saturation and hue by average */
    out = mw_change_cfimage(out,_imgin->nrow,_imgin->ncol);
    for(i=0;i<ncol;i++)
        for(j=0;j<nrow;j++){
            s = mw_get_smallest_shape(tree,i,j);
            indexshape = s-tree->the_shapes;
            sketch->gray[i+j*ncol] = gray[indexshape];
            out->red[i+j*ncol]= red[indexshape];
            out->green[i+j*ncol]= green[indexshape];
            out->blue[i+j*ncol]= blue[indexshape];
        }

    free(gray);
    free(red);
    free(green);
    free(blue);
    free(truearea);
    std::cout <<"AbstractionProcess::filter_shapes()::end"<< std::endl;
}


/*=== Filtering the image ===*/
void AbstractionProcess::filter_image(int *ns,
                                      float *threshold,
                                      int *mpixel){
    //@Declare variables here.==========================
    int i, kl ,j, rmn, lableTemp, nn;
    float thre;
    float  R,G,B,H,S,L, CONTR;
    float elong, elong_pre, kappa, kappa_pre, oren, oren_pre, sca, sca_pre, Dist;
    Shape pShape;

    thre = *threshold;
    nn = *ns;

    compute_shape_attribute(&nn);

    //filtering the image
    for(i = 0; i<=_pTree->nb_shapes-1; i++)  {
        pShape = _pTree->the_shapes + i;
        if(pShape->parent == NULL)
            continue;

        CONTR = sqrt(pow((((Info*)(pShape->data))->r - ((Info*)(pShape->parent->data))->r), 2.0) +
                     pow((((Info*)(pShape->data))->g - ((Info*)(pShape->parent->data))->g), 2.0) +
                     pow((((Info*)(pShape->data))->b - ((Info*)(pShape->parent->data))->b), 2.0));

        elong = ((Info*)(pShape->data))->attribute[2];
        kappa = ((Info*)(pShape->data))->attribute[1];
        oren  = ((Info*)(pShape->data))->attribute[3];
        sca   = (float) pShape->area;

        elong_pre = ((Info*)(pShape->parent->data))->attribute[2];
        kappa_pre = ((Info*)(pShape->parent->data))->attribute[1];
        oren_pre  = ((Info*)(pShape->parent->data))->attribute[3];
        sca_pre   = (float) pShape->parent->area;

        Dist = sqrt((elong - elong_pre)*(elong - elong_pre) +
                    (kappa - kappa_pre)*(kappa - kappa_pre) +
                    (oren - oren_pre)*(oren - oren_pre)/(PI*PI) +
                    (1 - _MIN(sca_pre/sca, sca/sca_pre))*(1 - _MIN(sca_pre/sca, sca/sca_pre)));
        Dist /= 4;

        if(pShape->area <= *mpixel
                || (((Info*)(pShape->data))->attribute[0])*CONTR<= thre
                || Dist*CONTR < 0.
                /*        ||( fabs(((Info*)(pShape->data))->attribute[0] - ((Info*)(pShape->parent->data))->attribute[0]) >= thre */
                /* 	   && fabs(((Info*)(pShape->data))->attribute[1] - ((Info*)(pShape->parent->data))->attribute[1]) >= thre */
                /* 	   && fabs(((Info*)(pShape->data))->attribute[2] - ((Info*)(pShape->parent->data))->attribute[2]) >= thre */
                /* 	   && fabs(((Info*)(pShape->data))->attribute[3] - ((Info*)(pShape->parent->data))->attribute[3]) >= thre */
                /* 	   && fabs(((Info*)(pShape->data))->attribute[4] - ((Info*)(pShape->parent->data))->attribute[4]) >= thre */
                /* 	   )//thre   */
                //|| DisVectL2(((Info*)(pShape->data))->attribute, ((Info*)(pShape->parent->data))->attribute, 4) >= thre
                /* 1/(1 + thre*((float)peri_shape(_pTree, pShape))/(float)pShape->area) */)
        {
            lableTemp = 1;
            pShape->removed = 1;
            kl++;
        }
        if(i ==0 )
            pShape->removed = 0;
    }
}

int AbstractionProcess::random_number(int *M)
{
    int i, size, select_i;
    float pb_sum, temp;
    Fsignal pb=0;

    size = *M;

    //   srand ( time(NULL) );
    /*   for(i=0; i<size; i++) */
    /*     printf("leaf: %f, %f \n", leaf->values[2*i +0], leaf->values[2*i +1]); */

    select_i = size-1;

    pb_sum = 0.0;
    pb = mw_change_fsignal(pb,size);
    mw_clear_fsignal(pb,0.0);

    for(i=0; i< size; i++)
    {
        /*// sampling by scale distribution; */
        //pb->values[i] = 1/(pow((float)leaf->values[2*i +1], -1.0));

        /*// sampling by uniform distribution; */
        pb->values[i] = 1;
        pb_sum += pb->values[i];
    }

    for(i=0; i< size; i++)
        pb->values[i] /= pb_sum;

    for(i=1; i< size; i++)
        pb->values[i] += pb->values[i-1];

    temp = ((float)rand())/RAND_MAX;

    for(i= 0; i< size; i++)
    {
        if( temp <= pb->values[i])
        {
            select_i = i;
            break;
        }
    }

    return select_i;
}

/*=================================================*/
/*=====  Select Shape according to the      =======*/
/*=====  distance of its attributes         =======*/
/*=================================================*/
/*=================================================*/
/*=== randS=0, randomly select shapes;          ===*/
/*=== randS=1, select shapes according to       ===*/
/*===        elongation, compactness and scale; ===*/
/*=== randS=2, select shapes according to       ===*/
/*===  elongation, compactness, scale and color ===*/
/*=================================================*/
Shape AbstractionProcess::selectShapeDict(Shapes pTreeDict,
                                          Shape pShape,
                                          float *paDict,
                                          int *randS)
{
    Shape pShapeDict, pShapeTemp;
    float lambda1, lambda2;

    float elong, kappa, elongDict, kappaDict, Dist, minDist;
    float sca, scaDict, pa;
    int i, index, mn, temp;

    lambda1 = ((Info*)(pShape->data))->lambda1;
    lambda2 = ((Info*)(pShape->data))->lambda2;
    elong = lambda2 / lambda1;
    kappa = ((float) pShape->area)/(sqrt(lambda2*lambda1)*4*PI);
    sca = ((float) pShape->area);

    index = 1; minDist = 10000.0;
    if(*randS == 0)
    {
        temp = pTreeDict->nb_shapes -1;
        index = random_number(&temp) + 1;
    }
    else
    {
        for(i= 1; i<pTreeDict->nb_shapes; i++)
        {
            pShapeDict = pTreeDict->the_shapes + i;

            mn=3;
            pShapeTemp =  m_order_parent(pShapeDict, &mn, true);
            pa = ((float) pShapeDict->area)/((float) pShapeTemp->area);

            if(pa < *paDict)
                continue;

            lambda1 = ((Info*)(pShapeDict->data))->lambda1;
            lambda2 = ((Info*)(pShapeDict->data))->lambda2;
            elongDict = lambda2 / lambda1;
            kappaDict = ((float) pShapeDict->area)/(sqrt(lambda2*lambda1)*4*PI);
            scaDict = ((float) pShapeDict->area);

            Dist = pow((elong - elongDict), 2.0) +
                    pow((kappa - kappaDict), 2.0) +
                    pow((1 - _MIN(sca/scaDict, scaDict/sca)), 2.0);

            if(*randS == 2)
            {
                Dist +=  ( pow((1 - _MIN(((Info*)(pShape->data))->r/((Info*)(pShapeDict->data))->r,
                                         ((Info*)(pShapeDict->data))->r/((Info*)(pShape->data))->r)), 2.0) +
                           pow((1 - _MIN(((Info*)(pShape->data))->g/((Info*)(pShapeDict->data))->g,
                                         ((Info*)(pShapeDict->data))->g/((Info*)(pShape->data))->g)), 2.0) +
                           pow((1 - _MIN(((Info*)(pShape->data))->b/((Info*)(pShapeDict->data))->b,
                                         ((Info*)(pShapeDict->data))->b/((Info*)(pShape->data))->b)), 2.0) )/3.0;
            }

            //         Dist =  sqrt(Dist);

            if(minDist > Dist)
            {
                minDist = Dist;
                index = i;
            }
        }
    }


    pShapeDict = pTreeDict->the_shapes + index;
    //    printf("index: %d \n", index);

    return pShapeDict;
}

QImage AbstractionProcess::render(TOSParameters tosParameters, bool &tree_recomputed, DictionaryParameters dictionaryParameters, TreeOfShapes * dictionnary){


    Ccimage imgsyn =NULL;
    if( tosParameters.model == 4 )
        imgsyn = _treeOfShapes->render( tosParameters, tree_recomputed, dictionnary, dictionaryParameters );
    else
        imgsyn = _treeOfShapes->render( tosParameters, tree_recomputed );
    QImage result_image( QSize(imgsyn->ncol, imgsyn->nrow), QImage::Format_RGB32 );

    _tosParameters = tosParameters;

    for( int j= 0; j< imgsyn->nrow; j++)
        for( int i= 0; i< imgsyn->ncol; i++)
        {
            int comp = j*imgsyn->ncol + i;

            QColor color (imgsyn->red[comp], imgsyn->green[comp], imgsyn->blue[comp]);
            result_image.setPixel(i, j , qRgb(color.red(), color.green(), color.blue()));
        }

    if( imgsyn != NULL )
        mw_delete_ccimage(imgsyn);
    return result_image;
}

void AbstractionProcess::save_shapes( QString folder_name, TOSParameters tosParameters, DictionaryParameters dictionaryParameters, TreeOfShapes * dictionnary){
    Ccimage imgsyn =NULL;

    bool tree_recomputed;
    imgsyn = _treeOfShapes->render( tosParameters, tree_recomputed, dictionnary, dictionaryParameters, true, folder_name );

    _tosParameters = tosParameters;

    if( imgsyn != NULL )
        mw_delete_ccimage(imgsyn);

}

std::vector<QImage> AbstractionProcess::render_shape_by_shape(TOSParameters tosParameters, DictionaryParameters dictionaryParameters, TreeOfShapes * dictionnary){


    std::vector<QImage> result_images;
    if( tosParameters.model == 4 )
        result_images = _treeOfShapes->render_shape_by_shape( tosParameters, dictionnary, dictionaryParameters );
    else
        result_images = _treeOfShapes->render_shape_by_shape( tosParameters );

    return result_images;
}



QImage AbstractionProcess::run (int process, char* dictionnary_name)
{
#if 0
    QImage result_image;

    /**************************************************/
    /*****  TIME AND CARRY OUT GRAIN RENDERING   ******/
    /**************************************************/

    struct timeval start, end;
    gettimeofday(&start, NULL);
    std::cout << "AbstractionProcess::run ()::mode "<< process << std::endl;

    QImage image_dict ( dictionnary_name );

    enum Abstraction_Mode mode =(enum Abstraction_Mode)process;

    if( mode == SHOW_TREE ){
        std::cout << "AbstractionProcess::run ()::SHOW_TREE" << std::endl;
        Cfimage out = mw_change_cfimage(NULL,_imgin->nrow,_imgin->ncol);
        showtree3D_matlab(_imgin, out);
        result_image = qimages_from_cfimage(out);

    } else if( mode == COLOR_SKETCH ) {
        std::cout << "AbstractionProcess::run ()::COLOR_SKETCH"<< std::endl;
#if 1
        float eps=10.; // "-log10(max number of false alarms)",
        char local=NULL ; //"local boundaries, default NULL",
        int nscale=2; // "number of scales for meaningful boundaries",
        Flists boundaries=mw_new_flists(); // "save meaningful boundaries",

        Cfimage out = mw_change_cfimage(NULL,_imgin->nrow,_imgin->ncol);
        color_sketch(_imgin,out,
                     boundaries,
                     &nscale,
                     &local,
                     &eps);

        result_image = qimages_from_cfimage(out);
#else


        Fimage in = fimageread(image_name);
        Fimage image_out = mw_change_fimage(NULL,in->nrow,in->ncol);
        printf("ll_boundaries2\n");
        fflush(stdout);
        float eps=0.; // "-log10(max number of false alarms), default 0.",
        float step=1.; // "quantization step, default 1.",
        int prec=2; //  "sampling precision for flst_bilinear, default 2",
        float std=0.5; //  "standard dev. for preliminary convolution, default 0.5",
        char all=NULL; //       "keep all meaningful level lines (not only maximal ones)",
        float hstep=0.01; // "step for contrast histogram, default 0.01",
        Shapes tree =NULL; //          "use a precomputed bilinear tree, default NULL",
        int visit=100; //  "maximal number of visits for a boundary, default 100",
        char loc= NULL; //                  "force local research",
        Shapes keep_tree = mw_new_shapes(); // "keep meaningful tree",

        Flists boundaries = ll_boundaries2(in,&eps,tree,
                                           &step,&prec,&std,&hstep,
                                           &all,&visit,&loc,image_out,
                                           keep_tree);
        printf("ll_boundaries2 done\n");
        fflush(stdout);
        fimagewrite(dictionnary_name, image_out);

#endif
    } else if( mode == FILTER_COLOR ){
        std::cout << "AbstractionProcess::run ()::FILTER_COLOR"<< std::endl;
        Ccimage out = mw_change_ccimage(NULL,_imgin->nrow,_imgin->ncol);
        // fgrain_side(min_area, in->gray, in->ncol, in->nrow, out->gray, sideflag);

        int ns=3;   // "scale ratio order",
        float alpha=.6; // "threshold",

        attribute_filter_color(&ns, &alpha, _imgin, out);

        result_image = qimages_from_ccimage(out);

    }else if( mode == SYNTEXTURE_COLOR ){
        std::cout << "AbstractionProcess::run ()::SYNTEXTURE_COLOR"<< std::endl;
        int model = 0;
        Cfimage out = mw_change_cfimage(NULL,_imgin->nrow,_imgin->ncol);
        syntexturecolor(_imgin, _pTree, &model, out);

        result_image = qimages_from_cfimage(out);
    } else if( mode == SYNTEXTURE_COLOR_WA ){
        std::cout << "AbstractionProcess::run ()::SYNTEXTURE_COLOR_WA"<< std::endl;

        int order=1; //rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2"
        int model=2; //synthesis model: rectangle: m=0; ellipse: m=1; orignal shape: m=2
        float alpha=0.; //alpha for transparent",
        int smodel=1; //shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
        float shift=5.0; //add a random shiftS to each shape, shiftS = shift*rand()",
        float theta=0.0; //add a random rotation to each shape, thetaS = theta*rand()",
        int mpixel=5; //minimal area (in pixel) for FLST",
        float kappa=0.; //compactness parameter of the attribute filtering on the orignal image",
        int median=13; //kernel size for median filter",
        int kerSize=3; //kernel size for Gaussian blur",
        float KerStd=0.8; //std for the gaussian kernel",

        Ccimage out = mw_change_ccimage(NULL,_imgin->nrow,_imgin->ncol);
        syntexturecolor_wa(&order,
                           &model,
                           &alpha,
                           &smodel,
                           &shift, &theta,
                           &mpixel,
                           &kappa,
                           &median,&kerSize,
                           &KerStd,
                           _imgin,
                           out);
        result_image = qimages_from_ccimage(out);
    } else if( mode == SYNTEXTURE_COLOR_V2 ){
        std::cout << "AbstractionProcess::run ()::SYNTEXTURE_COLOR_V2"<< std::endl;
        int order=1; //"rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2",
        int model=0; // "synthesis model: rectangle: m=0; ellipse: m=1;",
        float alpha=0.5; // "alpha for transparent",
        int smodel=1; // "shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
        float shift=0.; // "add a random shiftS to each shape, shiftS = shift*rand()",
        float theta=0.; // "add a random rotation to each shape, thetaS = theta*rand()",

        Ccimage imgsyn = mw_change_ccimage(NULL,_imgin->nrow,_imgin->ncol);
        syntexturecolor_v2(&order,
                           &model,
                           &alpha,
                           &smodel,
                           &shift,
                           &theta,
                           _imgin,
                           imgsyn);
        result_image = qimages_from_ccimage(imgsyn);
    } else if( mode == SYNTEXTURE_COLOR_V3 ){
        std::cout << "AbstractionProcess::run ()::SYNTEXTURE_COLOR_V3"<< std::endl;
        int order=1; //  "rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2",
        int model=3; //       "synthesis model: rectangle: m=0; ellipse: m=1; orignal shape: m=2",
        float alpha=0.; //     "alpha for transparent",
        int smodel=1; //     "shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
        float shift=0.0; //    "add a random shiftS to each shape, shiftS = shift*rand()",
        float theta=0.0; //     "add a random rotation to each shape, thetaS = theta*rand()",
        int mpixel=10; //    "minimal area (in pixel) for FLST",
        int relief=1; //     "add relief effects, if relief =1",
        float reliefOrentation=45; //    "relief orentation, in degree",
        float reliefHeight=3; //     "relief height",
        float kappa=0.; //    "compactness parameter of the attribute filtering on the orignal image",
        int median=3; //    "kernel size for median filter",
        int kerSize=3; //  "kernel size for Gaussian blur",
        float KerStd=0.4; //  "std for the gaussian kernel",

        Ccimage imgsyn = mw_change_ccimage(NULL,_imgin->nrow,_imgin->ncol);
        syntexturecolor_v3(&order,
                           &model,
                           &alpha,
                           &smodel,
                           &shift, &theta,
                           &mpixel, &relief,
                           &reliefOrentation, &reliefHeight, &kappa,
                           &median, &kerSize,
                           &KerStd,
                           _imgin,
                           imgsyn);
        result_image = qimages_from_ccimage(imgsyn);
    } else if( mode == SYNTEXTURE_COLOR_V4 ){
        std::cout << "AbstractionProcess::run ()::SYNTEXTURE_COLOR_V4"<< std::endl;
        int order=0; //  "rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2",
        int model=2; //   "synthesis model: rectangle: m=0; ellipse: m=1; orignal shape: m=2",
        float alpha=0.; //   "alpha for transparent",
        int smodel=1; //   "shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
        float shift=0.0; //  "add a random shiftS to each shape, shiftS = shift*rand()",
        float theta=0.0; //   "add a random rotation to each shape, thetaS = theta*rand()",
        int mpixel=5; //    "minimal area (in pixel) for FLST",
        int relief=0; //   "add relief effects, if relief =1",
        float reliefOrentation=45; //   "relief orentation, in degree",
        float reliefHeight=3; //   "relief height",
        float kappa=0.; //   "compactness parameter of the attribute filtering on the orignal image",
        int median=3; //   "kernel size for median filter",
        int kerSize=3; //  "kernel size for Gaussian blur",
        float KerStd=0.4; // "std for the gaussian kernel",

        Ccimage imgsyn = mw_change_ccimage(NULL,_imgin->nrow,_imgin->ncol);
        syntexturecolor_v4(&order,
                           &model,
                           &alpha,
                           &smodel,
                           &shift, &theta,
                           &mpixel, &relief,
                           &reliefOrentation, &reliefHeight, &kappa,
                           &median, &kerSize,
                           &KerStd,
                           _imgin,
                           imgsyn);
        result_image = qimages_from_ccimage(imgsyn);
    } else if (mode == SYNTEXTURE_COLOR_MULT ) {
        std::cout << "AbstractionProcess::run ()::SYNTEXTURE_COLOR_MULT"<< std::endl;
        int order=1; // "rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2",
        int model=1; //  "synthesis model: rectangle: m=0; ellipse: m=1;",
        float alpha=0.0; //  "alpha for transparent",
        float threshold=0.9; //  "compactness threshold for approximating shapes",
        int smallest=10; //   "smallest cc for approximating shapes",
        int smodel=1; //  "shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
        float shift=0.0; //  "add a random shiftS to each shape, shiftS = shift*rand()",
        float theta=0.0; //  "add a random rotation to each shape, thetaS = theta*rand()",

        Ccimage imgsyn = mw_change_ccimage(NULL,_imgin->nrow,_imgin->ncol);
        syntexturecolor_mult(&order,
                             &model,
                             &alpha,
                             &threshold,
                             &smallest,
                             &smodel,
                             &shift,
                             &theta,
                             _imgin,
                             imgsyn);
        result_image = qimages_from_ccimage(imgsyn);
    } else if( mode == SYNTEXTURE_COLOR_DICT ){

        std::cout << "AbstractionProcess::run ()::SYNTEXTURE_COLOR_DICT"<< std::endl;
        int order=1; //"rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2",
        float alpha=0.5; //"alpha for transparent",
        int equal=0; //"scaling shape with equal aspect ratio or not, if equal=1, scaling shape with equal aspect ratio, otherwise scaling    shape with different aspect ratio on x-axis and y-axis",
        int mcolor=0; //"color model: mcolor=1, use the color of tranferred image, otherwise use the color of the orignal image",
        float kappa=0; //"compactness parameter of the attribute filtering on the orignal image",
        float kappaDict=0; //"compactness parameter of the attribute filtering on the transferred image",
        int median=0; //"select shape smoothing methods: M=0: smooth using gaussian kernel, otherwise smooth with median filter",
        int kerSize=5; //"kernel size for shape smoothing methods",
        float KerStd=0.6; //"std for the gaussian kernel",

        Cfimage imgDict = cfimages_from_qimage(image_dict);
        Ccimage imgsyn = mw_change_ccimage(NULL,_imgin->nrow,_imgin->ncol);
        syntexturecolor_dict(&order,
                             &alpha,
                             &equal, &mcolor,
                             &kappa, &kappaDict,
                             &median, &kerSize,
                             &KerStd,
                             _imgin,
                             imgDict,
                             imgsyn);
        result_image = qimages_from_ccimage(imgsyn);
    } else if (mode == SYNTEXTURE_COLOR_DICT2 ){
        std::cout << "AbstractionProcess::run ()::SYNTEXTURE_COLOR_DICT2"<< std::endl;
        int order=1; //    "rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2",
        float alpha=0.; //  "alpha for transparent",
        int equal=1; //     "scaling shape with equal aspect ratio or not, if equal=1, scaling shape with equal aspect ratio, otherwise scaling    shape with different aspect ratio on x-axis and y-axis",
        int mcolor=1; //   "color model: mcolor=1, use the color of tranferred image, otherwise use the color of the orignal image",
        int randS=2; //     "selection model: randS=0, randomly select shapes;
        //              randS=1, select shapes according to elongation, compactness and scale;
        //             randS=2, select shapes according to elongation, compactness, scale and color",
        int mpixel=5; //    "minimal area (in pixel) for FLST",
        float paS2I=1.00; //   "parameter for transfer: areaOFshape / areaOFimage, if parameter > paS2I, shape is removed",
        float paC2S=0.00; //   "parameter for transfer: areaOFconnectcomponent / areaOFshape, if parameter < paC2S, shape is removed",
        float paS2P=1.00; //   "parameter for transfer: areaOFshape / areaOFshapeparent, if parameter > paS2P, shape is removed",
        int relief=0; //    "add relief effects, if relief =1",
        float reliefOrentation=45; //    "relief orentation, in degree",
        float reliefHeight=3; //    "relief height",
        float kappa=0; //      "compactness parameter of the attribute filtering on the orignal image",
        float kappaDict=0; // "compactness parameter of the attribute filtering on the transferred image",
        int median=3; //    "kernel size for median filter",
        int kerSize=3; //   "kernel size for Gaussian blur",
        float KerStd=0.4; // "std for the gaussian kernel",

        Cfimage imgDict = cfimages_from_qimage(image_dict);
        Ccimage imgsyn = mw_change_ccimage(NULL,_imgin->nrow,_imgin->ncol);
        syntexturecolor_dict2(&order,
                              &alpha,
                              &equal, &mcolor, &randS, &mpixel,
                              &paS2I, &paC2S, &paS2P,
                              &relief,
                              &reliefOrentation, &reliefHeight, &kappa, &kappaDict,
                              &median, &kerSize,
                              &KerStd,
                              _imgin,
                              imgDict,
                              imgsyn);
        result_image = qimages_from_ccimage(imgsyn);
    } else if( mode == SYNTEXTURE_COLOR_DICT2_OUTPUT ){
        std::cout << "AbstractionProcess::run ()::SYNTEXTURE_COLOR_DICT2_OUTPUT"<< std::endl;
        //SEGMENTATION FAULT !!!!

        int order=1; // "rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2",
        float alpha=0.; // "alpha for transparent",
        int equal=1; // "scaling shape with equal aspect ratio or not, if equal=1, scaling shape with equal aspect ratio, otherwise scaling    shape with different aspect ratio on x-axis and y-axis",
        int mcolor=0; // "color model: mcolor=1, use the color of tranferred image, otherwise use the color of the orignal image",
        int randS=2; // "selection model: randS=0, randomly select shapes;
        //                              randS=1, select shapes according to elongation, compactness and scale;
        //                              randS=2, select shapes according to elongation, compactness, scale and color",
        int mpixel=5; //  "minimal area (in pixel) for FLST",
        float paS2I=1.00; // "parameter for transfer: areaOFshape / areaOFimage, if parameter > paS2I, shape is removed",
        float paC2S=0.00; //  "parameter for transfer: areaOFconnectcomponent / areaOFshape, if parameter < paC2S, shape is removed",
        float paS2P=1.00; // "parameter for transfer: areaOFshape / areaOFshapeparent, if parameter > paS2P, shape is removed",
        int relief=0; //  "add relief effects, if relief =1",
        float reliefOrentation=45; // "relief orentation, in degree",
        float reliefHeight=3; //  "relief height",
        float kappa=0; //  "compactness parameter of the attribute filtering on the orignal image",
        float kappaDict=0; // "compactness parameter of the attribute filtering on the transferred image",
        int median=3; //   "kernel size for median filter",
        int kerSize=3; // "kernel size for Gaussian blur",
        float KerStd=0.4; // "std for the gaussian kernel",

        Cfimage imgDict = cfimages_from_qimage(image_dict);
        Ccimage imgsyn = mw_change_ccimage(NULL,_imgin->nrow,_imgin->ncol);
        Ccmovie synmovie = mw_new_ccmovie();
        syntexturecolor_dict2_output(&order,
                                     &alpha,
                                     &equal, &mcolor, &randS, &mpixel,
                                     &paS2I, &paC2S, &paS2P,
                                     &relief,
                                     &reliefOrentation, &reliefHeight, &kappa, &kappaDict,
                                     &median, &kerSize,
                                     &KerStd,
                                     _imgin,
                                     imgDict,
                                     imgsyn,
                                     synmovie);
        result_image = qimages_from_ccimage(imgsyn);
    } else if ( mode == SYNTEXTURE_COLOR_DICT3 ){
        std::cout << "AbstractionProcess::run ()::SYNTEXTURE_COLOR_DICT3"<< std::endl;

        int order=1; // "rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2",
        float alpha=0.; // "alpha for transparent",
        int equal=1; // "scaling shape with equal aspect ratio or not, if equal=1, scaling shape with equal aspect ratio, otherwise                              scaling shape with different aspect ratio on x-axis and y-axis",
        int mcolor=0; // "color model: mcolor=1, use the color of tranferred image, otherwise use the color of the orignal image",
        int smodel=1; //  "shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
        float shift=0.0; //  "add a random shiftS to each shape, shiftS = shift*rand()",
        int randS=2; //   "selection model: randS=0, randomly select shapes;
        //                 randS=1, select shapes according to elongation, compactness and scale;
        //                 randS=2, select shapes according to elongation, compactness, scale and color",
        int mpixel=20; //  "minimal area (in pixel) for FLST",
        float paS2I=1.00; //  "parameter for transfer: areaOFshape / areaOFimage, if parameter > paS2I, shape is removed",
        float paC2S=0.00; //  "parameter for transfer: areaOFconnectcomponent / areaOFshape, if parameter < paC2S, shape is removed",
        float paS2P=1.00; // "parameter for transfer: areaOFshape / areaOFshapeparent, if parameter > paS2P, shape is removed",
        float kappa=0.0; //    "compactness parameter of the attribute filtering on the orignal image",
        float kappaDict=0.0; //"compactness parameter of the attribute filtering on the transferred image",
        int median=3; //   "kernel size for median filter",
        int kerSize=3; // "kernel size for Gaussian blur",
        float KerStd=0.4; // "std for the gaussian kernel",


        Cfimage imgDict = cfimages_from_qimage(image_dict);
        Ccimage imgsyn = mw_change_ccimage(NULL,_imgin->nrow,_imgin->ncol);

        syntexturecolor_dict3(&order, &alpha, &equal, &mcolor, &smodel,
                              &shift, &randS, &mpixel, &paS2I, &paC2S, &paS2P, &kappa,
                              &kappaDict, &median, &kerSize, &KerStd, _imgin, imgDict, imgsyn);

        result_image = qimages_from_ccimage(imgsyn);
    } else if (mode == SYNTEXTURE_COLOR_TT ) {
        std::cout << "AbstractionProcess::run ()::SYNTEXTURE_COLOR_TT"<< std::endl;

        int order=1; // "rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2",
        float alpha=0.; // "alpha for transparent",
        int equal=1; //   "scaling shape with equal aspect ratio or not, if equal=1, scaling shape with equal aspect ratio, otherwise                              scaling shape with different aspect ratio on x-axis and y-axis",
        int mcolor=0; //  "color model: mcolor=1, use the color of tranferred image, otherwise use the color of the orignal image",
        int smodel=1; //   "shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
        float shift=0.0; //   "add a random shiftS to each shape, shiftS = shift*rand()",
        int randS=2; //    "selection model: randS=0, randomly select shapes;
        //                                       randS=1, select shapes according to elongation, compactness and scale;
        //                                       randS=2, select shapes according to elongation, compactness, scale and color",
        int mpixel=10; //  "minimal area (in pixel) for FLST",
        float paS2I=1.00; //  "parameter for transfer: areaOFshape / areaOFimage, if parameter > paS2I, shape is removed",
        float paC2S=0.00; // "parameter for transfer: areaOFconnectcomponent / areaOFshape, if parameter < paC2S, shape is removed",
        float paS2P=1.00; //  "parameter for transfer: areaOFshape / areaOFshapeparent, if parameter > paS2P, shape is removed",
        float kappa=0.0; //    "compactness parameter of the attribute filtering on the orignal image",
        float kappaDict=0.0; //"compactness parameter of the attribute filtering on the transferred image",
        int median=3; //  "kernel size for median filter",
        int kerSize=3; // "kernel size for Gaussian blur",
        float KerStd=0.4; // "std for the gaussian kernel",

        Cfimage imgDict = cfimages_from_qimage(image_dict);
        Ccimage imgsyn = mw_change_ccimage(NULL,_imgin->nrow,_imgin->ncol);
        syntexturecolor_tt(&order, &alpha, &equal, &mcolor, &smodel,
                           &shift, &randS, &mpixel, &paS2I,
                           & paC2S,& paS2P, &kappa,& kappaDict,
                           &median, &kerSize, &KerStd,
                           _imgin, imgDict, imgsyn);
        result_image = qimages_from_ccimage(imgsyn);
    }


    gettimeofday(&end, NULL);
    double elapsedTime = (end.tv_sec  - start.tv_sec) +
            (end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "time elapsed : " << elapsedTime <<" seconds"<< std::endl;
    std::cout << "***************************" << std::endl << std::endl << std::endl;

    return result_image;
#endif
}

