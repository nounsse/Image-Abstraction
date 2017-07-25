#include "TreeOfShapes.h"

#include "tree_of_shapes.h"
#include <iostream>
#include <sys/time.h>
#include <QColor>
#include <queue>
#include <algorithm>

extern void flst();
#define PI 3.1415926
#define _MIN(x, y) ( (x)<(y) ? (x) : (y) )
#define _MAX(x, y) ( (x)>(y) ? (x) : (y) )

int TreeOfShapes::_tree_count = 0;

TreeOfShapes::TreeOfShapes( Cfimage imageIn ){
    /**************************************************/
    /*************   READ INPUT IMAGE   ***************/
    /**************************************************/

    _imgin = imageIn;

    /**************************************************/
    /*************   SET DEFAULT INPUT OPTIONS   **************/
    /**************************************************/
    _tosParameters = getDefaultTOSParameters();
    _dictionaryParameters = getDefaultDictionaryParameters();

    _tree_computed = false;
    _tree_recomputed = false;

    _tree_id = _tree_count++;

    _large_to_small_index = NULL;
    _large_to_small_index_computed = false;
    _texture_image_loaded = false;
    _use_kdtree = false;
}

TreeOfShapes::TreeOfShapes( Cfimage imageIn, Cfimage texture_image ){


    /**************************************************/
    /*************   READ INPUT IMAGE   ***************/
    /**************************************************/

    _imgin = imageIn;
    _texture_image = texture_image;
    if( imageIn->nrow != texture_image->nrow || imageIn->ncol != texture_image->ncol )
        mwerror(FATAL,1,"TreeOfShapes()::Texture image of different size than input image.\n");

    _texture_image_loaded = true;
    /**************************************************/
    /*************   SET DEFAULT INPUT OPTIONS   **************/
    /**************************************************/
    _tosParameters = getDefaultTOSParameters();
    _dictionaryParameters = getDefaultDictionaryParameters();

    _tree_computed = false;
    _tree_recomputed = false;

    _tree_id = _tree_count++;

    _large_to_small_index = NULL;
    _large_to_small_index_computed = false;
    _use_kdtree = false;

}

void TreeOfShapes::init(Cfimage inputImg, Shapes &pTree){

    std::cout << "TreeOfShapes::init(Cfimage inputImg, Shapes &pTree)" << std::endl;
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

    /*=======================================*/
    /*====   Compute shape attribute     === */
    /*=======================================*/

    std::cout << "Compute shape attributes" << std::endl;
    compute_shape_attribute();
    //synthesize the image
    printf("---1---- Synthesize the image ..........\n");
}




TreeOfShapes::~TreeOfShapes(){
    if( _tree_computed ){
        mw_delete_shapes(_pTree);
        mw_delete_fimage(_NormOfDu);
        mw_delete_cfimage(_imgin);
        if( _large_to_small_index_computed )
            mw_delete_fsignal(_large_to_small_index);
        for (std::map<int, Fsignal>::iterator it = _dictionary_selections.begin(); it !=  _dictionary_selections.end(); ++it){
            mw_delete_fsignal( it->second );
        }

        if( _texture_image_loaded )
            mw_delete_cfimage(_texture_image);
    }
}

void TreeOfShapes::getTreeInfo(std::vector<QPoint> &positions, std::vector<QColor> &colors,
                               std::vector< std::vector< std::pair<int, int> > > &pixels, std::vector<int> & heights, std::vector<std::pair<int, int> > &edges){

    int queArr[_pTree->nb_shapes];
    int i, qInd;
    int *head, *rear;
    Shape pShape, pShapeTemp;

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
    i=0; qInd =1;

    edges.clear();
    while(head != rear){
        int parent_id = *head;
        pShape = _pTree->the_shapes + *head;
        //printf("%d %d \n", *head, *rear);
        head++;

        for(pShapeTemp = pShape->child; pShapeTemp != NULL;
            pShapeTemp = pShapeTemp->next_sibling){
            queArr[qInd++] = (((Info*)(pShapeTemp->data))->index);
            edges.push_back(std::make_pair(parent_id, (((Info*)(pShapeTemp->data))->index)));
            std::pair<int, int > pair = std::make_pair(parent_id, (((Info*)(pShapeTemp->data))->index));
            if( pair.first >= _pTree->nb_shapes || pair.second >= _pTree->nb_shapes )
                std::cout <<pair.first << " , " << pair.second << std::endl;
            rear++;
        }
    }

    positions.clear();
    heights.clear();
    colors.clear();
    pixels.clear();

    positions.resize(_pTree->nb_shapes);
    heights.resize(_pTree->nb_shapes);
    colors.resize(_pTree->nb_shapes);
    pixels.resize(_pTree->nb_shapes);
    int height;
    for(int j= 0; j< _pTree->nb_shapes; j++){
        pShape = _pTree->the_shapes + j;

        height = 0;
        for(pShapeTemp = pShape; pShapeTemp != NULL; pShapeTemp = pShapeTemp->parent)
            height++;

        heights[j] = height;
        positions[j] = QPoint((((Info*)(pShape->data))->x0),(((Info*)(pShape->data))->y0));
        colors[j] = QColor((((Info*)(pShape->data))->r), (((Info*)(pShape->data))->g), (((Info*)(pShape->data))->b));

        if( j==0 ){
            pixels[j].push_back(std::make_pair(_pTree->ncol, _pTree->nrow));
        } else {
            pixels[j].resize(pShape->area);
            for(i=0; i< pShape->area; i++){
                pixels[j][i]=std::make_pair(((pShape->pixels+i)->x), ((pShape->pixels+i)->y));
            }
        }
    }

}


/* This removes the shapes from the tree associated to pFloatImageInput
that are too small (threshold *pMinArea). As a consequence all the remaining
shapes of pFloatImageOutput are of area larger or equal than *pMinArea */

void TreeOfShapes::mw_fgrain_side(int *pMinArea, int sideflag)
{
    int i;
    /* Kill too small grains.
     Bound i>0 because it is forbidden to delete the root, at index 0 */
    for(i = _pTree->nb_shapes-1; i > 0; i--)
        if(_pTree->the_shapes[i].area < *pMinArea && ( (sideflag >0) ^ _pTree->the_shapes[i].inferior_type ))
            _pTree->the_shapes[i].removed = (char)1;

    /* Reconstruct in pFloatImageOutput the modified tree of shapes */
    //    flst_reconstruct(pTree, pFloatImageOutput);

}


/*========= Compute the orientation and Elongation of shapes ========*/
void TreeOfShapes::shape_orilam(Shape pShape, float *out_ori, float *out_e, float *out_k)
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
void TreeOfShapes::shape_orilam(Shape pShape, float *out_ori, float *out_e, float *out_k, float *pX0, float *pY0)
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
void TreeOfShapes::Order(Fsignal t2b_index,
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
Shape TreeOfShapes::m_order_parent(Shape pShape,
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
void TreeOfShapes::random_leaf(Fsignal leaf,
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

    mw_delete_fsignal(pb);
}

/* visit the tree in random order */
void TreeOfShapes::random_tree_order(Fsignal t2b_index)
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
void TreeOfShapes::top2bottom_index_tree(Fsignal t2b_index)
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
float TreeOfShapes::mean_contrast(Shape pShape)
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

void TreeOfShapes::adaptive_shift_shape(float *shift,
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


void TreeOfShapes::random_shift_shape(float *shift, float * theta)
{
    int i,j;
    Shape pShape;
    float SHIFT, THETA, tempx, tempy, temprotation;
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

//        if( (rand()%10) >5 )
//            temprotation = -((float)rand())/RAND_MAX;
//        else
//            temprotation = ((float)rand())/RAND_MAX;

        ((Info*)(pShape->data))->xShift = tempx*SHIFT;
        ((Info*)(pShape->data))->yShift = tempy*SHIFT;
 //       ((Info*)(pShape->data))->rotation = temprotation*THETA;
        /*       printf("temp: %f, shift: %f \n", ((Info*)(pShape->data))->xShift,  */
        /*       ((Info*)(pShape->data))->yShift);*/
    }
}

/*===== Compute the mean contrast of the curve l =====*/
float TreeOfShapes::min_contrast(Shape pShape)
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
/*========================================================*/
/*====== Compute boundingbox of each shape on the tree ===*/
/*====== contains two functions:                       ===*/
/*====== shape_boundingbox(pTree, pShape),  and        ===*/
/*====== tree_boundingbox(pTree, pShape)               ===*/
/*========================================================*/

/*====== Compute boundingbox of a shape===*/
/* Note: the setp 'flst_pixels(pTree)' should be done before this function; */
void TreeOfShapes::shape_boundingbox(Shape pShape)
{
    float xmin, xmax, ymin, ymax, theta, x0temp, y0temp;
    float x, y, xr, yr, sx, sy;
    int size, ncol, nrow, i;

    ncol = _pTree->ncol;
    nrow = _pTree->nrow;

    size = pShape->area;
    theta = ((Info*)(pShape->data))->oren;
    x0temp = ((Info*)(pShape->data))->x0;
    y0temp = ((Info*)(pShape->data))->y0;

    xmin = ncol;  xmax = -ncol;
    ymin = nrow;  ymax = -nrow;
    sx = 1.0; sy = 1.0;
    for(i = 0; i< size; i++)
    {
        x = (float)((pShape->pixels+i)->x - x0temp);
        y = (float)((pShape->pixels+i)->y - y0temp);

        //scaling_rotation_dict2(&x, &y, &sx, &sy, &theta, &xr, &yr);
        xr = (int)( x*cos(theta) + y*sin(theta));
        yr = (int)( y*cos(theta) - x*sin(theta));

        xmin = _MIN(xr, xmin);
        ymin = _MIN(yr, ymin);
        xmax = _MAX(xr, xmax);
        ymax = _MAX(yr, ymax);
    }

    ((Info*)(pShape->data))->boundingbox[0] = (int) xmin;
    ((Info*)(pShape->data))->boundingbox[1] = (int) xmax;
    ((Info*)(pShape->data))->boundingbox[2] = (int) ymin;
    ((Info*)(pShape->data))->boundingbox[3] = (int) ymax;
}

/*====== Compute boundingbox of each shape on the tree ===*/
void TreeOfShapes::tree_boundingbox()
{
    Shape pShape;
    int i;

    for(i = _pTree->nb_shapes-1; i>=0; i--)
    {
        pShape = _pTree->the_shapes + i;
        shape_boundingbox(pShape);
    }
}

/*====== Compute shape attribute ===*/
void TreeOfShapes::compute_shape_attribute()
{
    int i;
    float oren, lamb1, lamb2, x0, y0;
    Shape pShape;

    _average_r = 0.;
    _average_g = 0.;
    _average_b = 0.;
    _maxArea = 0;
    for(i = _pTree->nb_shapes-1; i>=0; i--)  {
        pShape = _pTree->the_shapes + i;

        shape_orilam(pShape, &oren, &lamb1, &lamb2, &x0, &y0);

        ((Info*)(pShape->data))->lambda1 = lamb1;
        ((Info*)(pShape->data))->lambda2 = lamb2;
        ((Info*)(pShape->data))->oren = oren;

        ((Info*)(pShape->data))->x0 = x0;
        ((Info*)(pShape->data))->y0 = y0;

        pShape->removed = 0;

        _average_r += ((Info*)(pShape->data))->r ;
        _average_g += ((Info*)(pShape->data))->g ;
        _average_b += ((Info*)(pShape->data))->b ;

        if(i != 0){
            ((Info*)(pShape->data))->contrast = fabs(min_contrast(pShape));
            _maxArea = std::max( pShape->area, _maxArea );
        }
    }

    _average_r /= _pTree->nb_shapes;
    _average_g /= _pTree->nb_shapes;
    _average_b /=  _pTree->nb_shapes;


}

/*================ Compute the perimeter of one shape ================== */
float TreeOfShapes::peri_shape(Shape pShape)
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
void TreeOfShapes::compute_shape_attribute(int *ns)
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

/*=== Synthesis by Shape Shaking ===*/
/*=== Before the Shaking, smooth the shape with a gaussian kernel or a median filter ===*/
void TreeOfShapes::synshapeRect(Shape pShape,
                                Ccimage imgsyn,
                                Cimage imgShapeLabelSyn,
                                Fimage imgShapeBlurSyn,
                                Fsignal gaussKernel,
                                int *median,
                                float *alpha,
                                int *relief,
                                float *reliefOrentation, float *reliefHeight)
{

    int xi, yi, x, y, iKer, jKer, KerSize, MedSize, xKer, yKer, numMedain;
    float bLimit, xr, yr, ALPHA, BETA, a, b, x0temp, y0temp, top, right, left, bottom;
    float phi, xi_e, yi_e;
    float xShift, yShift, theta, tR, tG, tB, TR, TG, TB, tr, tg, tb;

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
    for( xi= ceil(left); xi<= right; xi++)
        for( yi= ceil(top); yi<= bottom; yi++)
        {
            xi_e = ((float)xi - x0temp)*cos(phi+theta) + ((float)yi - y0temp)*sin(phi+theta);
            yi_e = ((float)yi - y0temp)*cos(phi+theta) - ((float)xi - x0temp)*sin(phi+theta);

            if( xi_e >= -a && xi_e <= +a && yi_e >= -b && yi_e <= +b )
            {

                if(xi<0 || xi>= imgsyn->ncol ||
                        yi<0 || yi>= imgsyn->nrow )
                    continue;

                imgShapeLabelSyn->gray[yi*imgShapeLabelSyn->ncol + xi] = 1;

                left   = _MIN(xi, left);
                top    = _MIN(yi, top);
                right  = _MAX(xi, right);
                bottom = _MAX(yi, bottom);
            }
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


/*=== Synthesis a shape using a rectangle with the same parameters (second-order moments)===*/
void TreeOfShapes::synshapeRect(Shape pShape,
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
/*=== Synthesis by Shape Shaking ===*/
/*=== Before the Shaking, smooth the shape with a gaussian kernel or a median filter ===*/
void TreeOfShapes::synshapeEllipse(Shape pShape,
                                   Ccimage imgsyn,
                                   Cimage imgShapeLabelSyn,
                                   Fimage imgShapeBlurSyn,
                                   Fsignal gaussKernel,
                                   int *median,
                                   float *alpha,
                                   int *relief,
                                   float *reliefOrentation, float *reliefHeight)
{

    int xi, yi, x, y, iKer, jKer, KerSize, MedSize, xKer, yKer, numMedain;
    float xr, yr, ALPHA, BETA, a, b, x0temp, y0temp, top, right, left, bottom;
    float phi, xi_e, yi_e;
    float xShift, yShift, theta, tR, tG, tB, TR, TG, TB, tr, tg, tb;

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
    for( xi= ceil(left); xi<= right; xi++)
        for( yi= ceil(top); yi<= bottom; yi++)
        {
            xi_e = ((float)xi - x0temp)*cos(phi+theta) + ((float)yi - y0temp)*sin(phi+theta);
            yi_e = ((float)yi - y0temp)*cos(phi+theta) - ((float)xi - x0temp)*sin(phi+theta);

            if( xi_e*xi_e/(a*a) + yi_e*yi_e/(b*b) <= 1 ){

                if(xi<0 || xi>= imgsyn->ncol ||
                        yi<0 || yi>= imgsyn->nrow )
                    continue;

                imgShapeLabelSyn->gray[yi*imgShapeLabelSyn->ncol + xi] = 1;

                left   = _MIN(xi, left);
                top    = _MIN(yi, top);
                right  = _MAX(xi, right);
                bottom = _MAX(yi, bottom);
            }
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

/*=== Synthesis by Shape Shaking ===*/
/*=== Before the Shaking, smooth the shape with a gaussian kernel or a median filter ===*/
void TreeOfShapes::synshapeCircle(Shape pShape,
                                  Ccimage imgsyn,
                                  Cimage imgShapeLabelSyn,
                                  Fimage imgShapeBlurSyn,
                                  Fsignal gaussKernel,
                                  int *median,
                                  float *alpha,
                                  int *relief,
                                  float *reliefOrentation, float *reliefHeight)
{

    int xi, yi, x, y, iKer, jKer, KerSize, MedSize, xKer, yKer, numMedain;
    float xr, yr, ALPHA, BETA, a, b, x0temp, y0temp, top, right, left, bottom;
    float phi, xi_e, yi_e;
    float xShift, yShift, theta, tR, tG, tB, TR, TG, TB, tr, tg, tb;

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

    left   = _MAX(0, x0temp - b);
    right  = _MIN(_pTree->ncol -1, x0temp + b);
    top    = _MAX(0, y0temp - b);
    bottom = _MIN(_pTree->nrow - 1, y0temp + b);

    TR  = ((Info*)(pShape->data))->r;
    TG  = ((Info*)(pShape->data))->g;
    TB  = ((Info*)(pShape->data))->b;
    for( xi= ceil(left); xi<= right; xi++)
        for( yi= ceil(top); yi<= bottom; yi++)
        {
            xi_e = ((float)xi - x0temp)*cos(phi+theta) + ((float)yi - y0temp)*sin(phi+theta);
            yi_e = ((float)yi - y0temp)*cos(phi+theta) - ((float)xi - x0temp)*sin(phi+theta);

            if( xi_e*xi_e/(b*b) + yi_e*yi_e/(b*b) <= 1 ){

                if(xi<0 || xi>= imgsyn->ncol ||
                        yi<0 || yi>= imgsyn->nrow )
                    continue;

                imgShapeLabelSyn->gray[yi*imgShapeLabelSyn->ncol + xi] = 1;

                left   = _MIN(xi, left);
                top    = _MIN(yi, top);
                right  = _MAX(xi, right);
                bottom = _MAX(yi, bottom);
            }
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

/*=== Synthesis a shape using a ellipse with the same parameters (second-order moments) ===*/
void TreeOfShapes::synshapeCircle(Shape pShape,
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

    left   = _MAX(0, x0temp - b);
    right  = _MIN(_pTree->ncol -1, x0temp + b);
    top    = _MAX(0, y0temp - a);
    bottom = _MIN(_pTree->nrow - 1, y0temp + b);

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

                if( xi_e*xi_e/(b*b) + yi_e*yi_e/(b*b) <= 1 ){

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

            if( xi_e*xi_e/(b*b) + yi_e*yi_e/(b*b) <= 1 ){

                tR = ((float) imgsyn->red[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->r;
                imgsyn->red[yi*_pTree->ncol + xi] = (int) tR; /*(int)rint((double) tR); */

                tG = ((float) imgsyn->green[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->g;
                imgsyn->green[yi*_pTree->ncol + xi] = (int) tG; /*(int)rint((double) tG); */

                tB = ((float) imgsyn->blue[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->b;
                imgsyn->blue[yi*_pTree->ncol + xi] = (int) tB; /*(int)rint((double) tB); */
            }
        }
}


/*=== Synthesis a shape using a ellipse with the same parameters (second-order moments) ===*/
void TreeOfShapes::synshapeEllipse(Shape pShape,
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


/*========================================================*/
/*====== Synthesize a shape by a given shape from      ===*/
/*====== the dictionary.                               ===*/
/*=== Before the Shaking, smooth the shape with a      ===*/
/*=== gaussian kernel or a median filter               ===*/
/*========================================================*/
void TreeOfShapes::synShapeDict(Shape pShapeDict, Shape pShape,
                                Ccimage imgsyn,
                                Cfimage imgDict, Cfimage imgShapeColorSyn,
                                Cimage imgShapeLabel, Cimage imgShapeLabelSyn,
                                float *alpha,
                                int *equal, int *mcolor, int *relief,
                                float *reliefOrentation, float *reliefHeight)
{
#if 1
    int i, x, y,  xi, yi;
    float xr, yr, xt, yt, x0temp, y0temp, x0tempDict, y0tempDict, x0tempTransformed, y0tempTransformed;
    float tR, tG, tB, ALPHA, TR, TG, TB;
    float xBound_l, xBound_r, yBound_t, yBound_b;
    float SCALEx, SCALEy, theta, thetaDict, rotation, xShift, yShift;
    float top, right, left, bottom;
    float a, b, bLimit, tempx;
    ALPHA = *alpha;

    xBound_l = (int) ((Info*)(pShapeDict->data))->boundingbox[0];
    xBound_r = (int) ((Info*)(pShapeDict->data))->boundingbox[1];
    yBound_t = (int) ((Info*)(pShapeDict->data))->boundingbox[2];
    yBound_b = (int) ((Info*)(pShapeDict->data))->boundingbox[3];

    /*=== begin=== Label pShapeDict ========*/
    for(i=0; i < pShapeDict->area; i++)
    {
        x = (pShapeDict->pixels+i)->x;
        y = (pShapeDict->pixels+i)->y;
        imgShapeLabel->gray[y*imgShapeLabel->ncol + x] = 1;
    }

    srand( time(NULL) );
    if(*mcolor == 1)
    {
        tempx = (rand()%10);
        TR = ((Info*)(pShapeDict->data))->r + tempx*0;
        TG = ((Info*)(pShapeDict->data))->g + tempx*0;
        TB = ((Info*)(pShapeDict->data))->b + tempx*0;
    }

    else
    {
        TR  = ((Info*)(pShape->data))->r + tempx*0;
        TG  = ((Info*)(pShape->data))->g + tempx*0;
        TB  = ((Info*)(pShape->data))->b + tempx*0;
    }

    if(*mcolor == 2){
        TR  = ((Info*)(pShapeDict->data))->r;
        TG  = ((Info*)(pShapeDict->data))->g;
        TB  = ((Info*)(pShapeDict->data))->b;
    }

    if(*equal == 1)
    {
        SCALEx = sqrt( (double) pShape->area / (double) pShapeDict->area );
        SCALEy = SCALEx;
    }
    else
    {
        SCALEx =  sqrt( ((Info*)(pShape->data))->lambda1) /  sqrt(((Info*)(pShapeDict->data))->lambda1 );
        SCALEy =  sqrt( ((Info*)(pShape->data))->lambda2) /  sqrt(((Info*)(pShapeDict->data))->lambda2 );
    }

    a = SCALEx * (xBound_r - xBound_l);
    b = SCALEy * (yBound_b - yBound_t);

    bLimit = (1.0*sqrt(a*a + b*b)/2.0);

    x0temp = ((Info*)(pShape->data))->x0;
    y0temp = ((Info*)(pShape->data))->y0;

    x0tempDict =  ((Info*)(pShapeDict->data))->x0;
    y0tempDict =  ((Info*)(pShapeDict->data))->y0;

    theta  = ((Info*)(pShape->data))->oren;
    thetaDict = ((Info*)(pShapeDict->data))->oren;

    xShift = (((Info*)(pShape->data))->xShift);
    yShift = (((Info*)(pShape->data))->yShift);

    rotation  = (((Info*)(pShape->data))->rotation);

    left   = _MAX(0, x0temp - bLimit);
    right  = _MIN(imgsyn->ncol - 1, x0temp + bLimit );
    top    = _MAX(0, y0temp - bLimit);
    bottom = _MIN(imgsyn->nrow - 1, y0temp + bLimit);


    /*=== Transformation ===*/
    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++)
        {
            xt = (x - x0temp);
            yt = (y - y0temp);

            /*====  Rotate to the major axis for scaling    ====*/
            xi = ( (1/SCALEx)* (xt*cos(theta) + yt*sin(theta)));
            yi = ( (1/SCALEy)* (yt*cos(theta) - xt*sin(theta)));

            /*====  Inverse-transformation    ====*/
            xr = ( (xi*cos(thetaDict) - yi*sin(thetaDict)) + x0tempDict);
            yr = ( (yi*cos(thetaDict) + xi*sin(thetaDict)) + y0tempDict);

            if(xr <0 || xr >=imgShapeLabel->ncol ||
                    yr <0 || yr >=imgShapeLabel->nrow)
                continue;

            if( imgShapeLabel->gray[(int)(yr)*imgShapeLabel->ncol + (int)(xr)] == 1 )
            {
                imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] = 1;

                imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x]   =  imgDict->red[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x] =  imgDict->green[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x]  =  imgDict->blue[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
            }
        }

    x0tempTransformed = (x0tempDict - x0temp)*cos(theta) + (y0tempDict - y0temp)*sin(theta) +x0temp;
    y0tempTransformed = (y0tempDict - y0temp)*cos(theta) - (x0tempDict - x0temp)*sin(theta) +y0temp;

    if(*relief == 1 && pShape->area > 30)
    {
        float shLambda, shTR, shTG, shTB;
        int xsh, ysh, shiftsh;
        if(pShape->area > 30)
            shiftsh = (*reliefHeight);
        else
            shiftsh = (*reliefHeight)*( (float) pShape->area /30.0);;

        shLambda = 0.3;
        for(x = ceil(left); x <= right; x++)
            for(y = ceil(top); y <= bottom; y++)
            {
                if(imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] == 0)
                    continue;

                if(*mcolor == 1){
                    TR  = imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x];
                    TG  = imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x];
                    TB  = imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x];
                }

                xr = (x - x0tempTransformed)*cos(rotation) + (y - y0tempTransformed)*sin(rotation);
                yr = (y - y0tempTransformed)*cos(rotation) - (x - x0tempTransformed)*sin(rotation);

                xi = floor(xShift + x0tempTransformed + xr);
                yi = floor(yShift + y0tempTransformed + yr);

                shTR = TR* shLambda;
                shTG = TG* shLambda;
                shTB = TB* shLambda;


                xsh = xi + shiftsh*cos( PI*(*reliefOrentation)/180.0 );
                ysh = yi - shiftsh*sin( PI*(*reliefOrentation)/180.0 );
                xsh = _MAX(0, xsh);
                xsh = _MIN(imgsyn->ncol - 1, xsh);
                ysh = _MAX(0, ysh);
                ysh = _MIN(imgsyn->nrow - 1, ysh);

                if(xsh>= 0 && xsh< _pTree->ncol &&
                        ysh>= 0 && ysh< _pTree->nrow){
                    tR = ((float) imgsyn->red[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTR;
                    imgsyn->red[ysh*imgsyn->ncol + xsh]   = (int) tR;/*(int)rint((double) tR); */

                    tG = ((float) imgsyn->green[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTG;
                    imgsyn->green[ysh*imgsyn->ncol + xsh] = (int) tG; /*(int)rint((double) tG); */

                    tB = ((float) imgsyn->blue[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTB;
                    imgsyn->blue[ysh*imgsyn->ncol + xsh]  = (int) tB; /*(int)rint((double) tB);*/
                }

            }
    }

    TR  = ((Info*)(pShape->data))->r + tempx*0;
    TG  = ((Info*)(pShape->data))->g + tempx*0;
    TB  = ((Info*)(pShape->data))->b + tempx*0;

    if(*mcolor == 2){
        TR  = ((Info*)(pShapeDict->data))->r;
        TG  = ((Info*)(pShapeDict->data))->g;
        TB  = ((Info*)(pShapeDict->data))->b;
    }
    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++)
        {
            if(imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] == 0)
                continue;

            if(*mcolor == 1){
                TR  = imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x];
                TG  = imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x];
                TB  = imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x];
            }


            xr = (x - x0tempTransformed)*cos(rotation) + (y - y0tempTransformed)*sin(rotation);
            yr = (y - y0tempTransformed)*cos(rotation) - (x - x0tempTransformed)*sin(rotation);

            xi = floor(xShift + x0tempTransformed + xr);
            yi = floor(yShift + y0tempTransformed + yr);

            if(xi>= 0 && xi< _pTree->ncol &&
                    yi>= 0 && yi< _pTree->nrow){
                tR = ((float) imgsyn->red[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*TR;
                imgsyn->red[yi*_pTree->ncol + xi] = (int)rint((double) tR);

                tG = ((float) imgsyn->green[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*TG;
                imgsyn->green[yi*_pTree->ncol + xi] = (int)rint((double) tG);

                tB = ((float) imgsyn->blue[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*TB;
                imgsyn->blue[yi*_pTree->ncol + xi] = (int)rint((double) tB);
            }
        }
#else

    int i, x, y;
    float xi, yi, xr, yr, xt, yt, x0temp, y0temp, x0tempDict, y0tempDict;
    float tR, tG, tB, ALPHA, TR, TG, TB;
    float xBound_l, xBound_r, yBound_t, yBound_b;
    float SCALEx, SCALEy, theta, thetaDict;
    float top, right, left, bottom;
    float a, b, bLimit, tempx;
    ALPHA = *alpha;

    xBound_l = (int) ((Info*)(pShapeDict->data))->boundingbox[0];
    xBound_r = (int) ((Info*)(pShapeDict->data))->boundingbox[1];
    yBound_t = (int) ((Info*)(pShapeDict->data))->boundingbox[2];
    yBound_b = (int) ((Info*)(pShapeDict->data))->boundingbox[3];

    /*=== begin=== Label pShapeDict ========*/
    for(i=0; i < pShapeDict->area; i++)
    {
        x = (pShapeDict->pixels+i)->x;
        y = (pShapeDict->pixels+i)->y;
        imgShapeLabel->gray[y*imgShapeLabel->ncol + x] = 1;
    }

    srand( time(NULL) );
    if(*mcolor == 1)
    {
        tempx = (rand()%10);
        TR = ((Info*)(pShapeDict->data))->r + tempx*0;
        TG = ((Info*)(pShapeDict->data))->g + tempx*0;
        TB = ((Info*)(pShapeDict->data))->b + tempx*0;
    }

    else
    {
        TR  = ((Info*)(pShape->data))->r + tempx*0;
        TG  = ((Info*)(pShape->data))->g + tempx*0;
        TB  = ((Info*)(pShape->data))->b + tempx*0;
    }

    if(*equal == 1)
    {
        SCALEx = sqrt( (double) pShape->area / (double) pShapeDict->area );
        SCALEy = SCALEx;
    }
    else
    {
        SCALEx =  sqrt( ((Info*)(pShape->data))->lambda1) /  sqrt(((Info*)(pShapeDict->data))->lambda1 );
        SCALEy =  sqrt( ((Info*)(pShape->data))->lambda2) /  sqrt(((Info*)(pShapeDict->data))->lambda2 );
    }

    a = SCALEx * (xBound_r - xBound_l);
    b = SCALEy * (yBound_b - yBound_t);

    bLimit = (1.0*sqrt(a*a + b*b)/2.0);

    x0temp = ((Info*)(pShape->data))->x0;
    y0temp = ((Info*)(pShape->data))->y0;
    x0tempDict =  ((Info*)(pShapeDict->data))->x0;
    y0tempDict =  ((Info*)(pShapeDict->data))->y0;

    theta  = ((Info*)(pShape->data))->oren;
    thetaDict = ((Info*)(pShapeDict->data))->oren;

    left   = _MAX(0, x0temp - bLimit);
    right  = _MIN(imgsyn->ncol - 1, x0temp + bLimit );
    top    = _MAX(0, y0temp - bLimit);
    bottom = _MIN(imgsyn->nrow - 1, y0temp + bLimit);


    /*=== Transformation ===*/
    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++)
        {
            xt = (x - x0temp);
            yt = (y - y0temp);

            /*====  Rotate to the major axis for scaling    ====*/
            xi = ( (1/SCALEx)* (xt*cos(theta) + yt*sin(theta)));
            yi = ( (1/SCALEy)* (yt*cos(theta) - xt*sin(theta)));

            /*====  Inverse-transformation    ====*/
            xr = ( (xi*cos(thetaDict) - yi*sin(thetaDict)) + x0tempDict);
            yr = ( (yi*cos(thetaDict) + xi*sin(thetaDict)) + y0tempDict);

            if(xr <0 || xr >=imgShapeLabel->ncol ||
                    yr <0 || yr >=imgShapeLabel->nrow)
                continue;

            if( imgShapeLabel->gray[(int)(yr)*imgShapeLabel->ncol + (int)(xr)] == 1 )
            {
                imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] = 1;

                imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x]   =  imgDict->red[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x] =  imgDict->green[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x]  =  imgDict->blue[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
            }
        }


    if(*relief == 1 && pShape->area > 30)
    {
        float shLambda, shTR, shTG, shTB;
        int xsh, ysh, shiftsh;
        if(pShape->area > 30)
            shiftsh = (*reliefHeight);
        else
            shiftsh = (*reliefHeight)*( (float) pShape->area /30.0);;

        shLambda = 0.3;
        for(x = ceil(left); x <= right; x++)
            for(y = ceil(top); y <= bottom; y++)
            {
                if(imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] == 0)
                    continue;

                if(*mcolor == 1){
                    TR  = imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x];
                    TG  = imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x];
                    TB  = imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x];
                }


                shTR = TR* shLambda;
                shTG = TG* shLambda;
                shTB = TB* shLambda;


                xsh = x + shiftsh*cos( PI*(*reliefOrentation)/180.0 );
                ysh = y - shiftsh*sin( PI*(*reliefOrentation)/180.0 );
                xsh = _MAX(0, xsh);
                xsh = _MIN(imgsyn->ncol - 1, xsh);
                ysh = _MAX(0, ysh);
                ysh = _MIN(imgsyn->nrow - 1, ysh);

                if(xsh>= 0 && xsh< _pTree->ncol &&
                        ysh>= 0 && ysh< _pTree->nrow){
                    tR = ((float) imgsyn->red[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTR;
                    imgsyn->red[ysh*imgsyn->ncol + xsh]   = (int) tR;/*(int)rint((double) tR); */

                    tG = ((float) imgsyn->green[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTG;
                    imgsyn->green[ysh*imgsyn->ncol + xsh] = (int) tG; /*(int)rint((double) tG); */

                    tB = ((float) imgsyn->blue[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTB;
                    imgsyn->blue[ysh*imgsyn->ncol + xsh]  = (int) tB; /*(int)rint((double) tB);*/
                }

            }
    }

    TR  = ((Info*)(pShape->data))->r + tempx*0;
    TG  = ((Info*)(pShape->data))->g + tempx*0;
    TB  = ((Info*)(pShape->data))->b + tempx*0;

    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++)
        {
            if(imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] == 0)
                continue;

            if(*mcolor == 1){
                TR  = imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x];
                TG  = imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x];
                TB  = imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x];
            }

            tR = ((float) imgsyn->red[y*_pTree->ncol + x])*ALPHA + (1-ALPHA)*TR;
            imgsyn->red[y*_pTree->ncol + x] = (int)rint((double) tR);

            tG = ((float) imgsyn->green[y*_pTree->ncol + x])*ALPHA + (1-ALPHA)*TG;
            imgsyn->green[y*_pTree->ncol + x] = (int)rint((double) tG);

            tB = ((float) imgsyn->blue[y*_pTree->ncol + x])*ALPHA + (1-ALPHA)*TB;
            imgsyn->blue[y*_pTree->ncol + x] = (int)rint((double) tB);
        }
#endif
}

void TreeOfShapes::synShapeDict(Shape pShapeDict, Shape pShape,
                                Ccimage imgsyn,
                                Cfimage imgDict, Cfimage imgShapeColorSyn,
                                Cimage imgShapeLabel, Cimage imgShapeLabelSyn,
                                Fimage imgShapeBlurSyn,
                                Fsignal gaussKernel,
                                int *median,
                                float *alpha,
                                int *equal, int *mcolor, int *relief,
                                float *reliefOrentation, float *reliefHeight)
{
#if 0
    int i, x, y,  xi, yi;
    float xr, yr, xt, yt, x0temp, y0temp, x0tempDict, y0tempDict, x0tempTransformed, y0tempTransformed;
    float tR, tG, tB, ALPHA, TR, TG, TB;
    float xBound_l, xBound_r, yBound_t, yBound_b;
    float SCALEx, SCALEy, theta, thetaDict, rotation, xShift, yShift;
    float top, right, left, bottom;
    float a, b, bLimit, tempx;
    ALPHA = *alpha;

    xBound_l = (int) ((Info*)(pShapeDict->data))->boundingbox[0];
    xBound_r = (int) ((Info*)(pShapeDict->data))->boundingbox[1];
    yBound_t = (int) ((Info*)(pShapeDict->data))->boundingbox[2];
    yBound_b = (int) ((Info*)(pShapeDict->data))->boundingbox[3];

    /*=== begin=== Label pShapeDict ========*/
    for(i=0; i < pShapeDict->area; i++)
    {
        x = (pShapeDict->pixels+i)->x;
        y = (pShapeDict->pixels+i)->y;
        imgShapeLabel->gray[y*imgShapeLabel->ncol + x] = 1;
    }

    srand( time(NULL) );
    if(*mcolor == 1)
    {
        tempx = (rand()%10);
        TR = ((Info*)(pShapeDict->data))->r + tempx*0;
        TG = ((Info*)(pShapeDict->data))->g + tempx*0;
        TB = ((Info*)(pShapeDict->data))->b + tempx*0;
    }

    else
    {
        TR  = ((Info*)(pShape->data))->r + tempx*0;
        TG  = ((Info*)(pShape->data))->g + tempx*0;
        TB  = ((Info*)(pShape->data))->b + tempx*0;
    }

    if(*equal == 1)
    {
        SCALEx = sqrt( (double) pShape->area / (double) pShapeDict->area );
        SCALEy = SCALEx;
    }
    else
    {
        SCALEx =  sqrt( ((Info*)(pShape->data))->lambda1) /  sqrt(((Info*)(pShapeDict->data))->lambda1 );
        SCALEy =  sqrt( ((Info*)(pShape->data))->lambda2) /  sqrt(((Info*)(pShapeDict->data))->lambda2 );
    }

    a = SCALEx * (xBound_r - xBound_l);
    b = SCALEy * (yBound_b - yBound_t);

    bLimit = (1.0*sqrt(a*a + b*b)/2.0);

    x0temp = ((Info*)(pShape->data))->x0;
    y0temp = ((Info*)(pShape->data))->y0;

    x0tempDict =  ((Info*)(pShapeDict->data))->x0;
    y0tempDict =  ((Info*)(pShapeDict->data))->y0;

    theta  = ((Info*)(pShape->data))->oren;
    thetaDict = ((Info*)(pShapeDict->data))->oren;

    xShift = (((Info*)(pShape->data))->xShift);
    yShift = (((Info*)(pShape->data))->yShift);

    rotation  = (((Info*)(pShape->data))->rotation);

    left   = _MAX(0, x0temp - bLimit);
    right  = _MIN(imgsyn->ncol - 1, x0temp + bLimit );
    top    = _MAX(0, y0temp - bLimit);
    bottom = _MIN(imgsyn->nrow - 1, y0temp + bLimit);


    /*=== Transformation ===*/
    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++)
        {
            xt = (x - x0temp);
            yt = (y - y0temp);

            /*====  Rotate to the major axis for scaling    ====*/
            xi = ( (1/SCALEx)* (xt*cos(theta) + yt*sin(theta)));
            yi = ( (1/SCALEy)* (yt*cos(theta) - xt*sin(theta)));

            /*====  Inverse-transformation    ====*/
            xr = ( (xi*cos(thetaDict) - yi*sin(thetaDict)) + x0tempDict);
            yr = ( (yi*cos(thetaDict) + xi*sin(thetaDict)) + y0tempDict);

            if(xr <0 || xr >=imgShapeLabel->ncol ||
                    yr <0 || yr >=imgShapeLabel->nrow)
                continue;

            if( imgShapeLabel->gray[(int)(yr)*imgShapeLabel->ncol + (int)(xr)] == 1 )
            {
                imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] = 1;

                imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x]   =  imgDict->red[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x] =  imgDict->green[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x]  =  imgDict->blue[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
            }
        }

    x0tempTransformed = (x0tempDict - x0temp)*cos(theta) + (y0tempDict - y0temp)*sin(theta) +x0temp;
    y0tempTransformed = (y0tempDict - y0temp)*cos(theta) - (x0tempDict - x0temp)*sin(theta) +y0temp;

    if(*relief == 1 && pShape->area > 30)
    {
        float shLambda, shTR, shTG, shTB;
        int xsh, ysh, shiftsh;
        if(pShape->area > 30)
            shiftsh = (*reliefHeight);
        else
            shiftsh = (*reliefHeight)*( (float) pShape->area /30.0);;

        shLambda = 0.3;
        for(x = ceil(left); x <= right; x++)
            for(y = ceil(top); y <= bottom; y++)
            {
                if(imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] == 0)
                    continue;

                if(*mcolor == 1){
                    TR  = imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x];
                    TG  = imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x];
                    TB  = imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x];
                }

                xr = (x - x0tempTransformed)*cos(rotation) + (y - y0tempTransformed)*sin(rotation);
                yr = (y - y0tempTransformed)*cos(rotation) - (x - x0tempTransformed)*sin(rotation);

                xi = floor(xShift + x0tempTransformed + xr);
                yi = floor(yShift + y0tempTransformed + yr);

                shTR = TR* shLambda;
                shTG = TG* shLambda;
                shTB = TB* shLambda;


                xsh = xi + shiftsh*cos( PI*(*reliefOrentation)/180.0 );
                ysh = yi - shiftsh*sin( PI*(*reliefOrentation)/180.0 );
                xsh = _MAX(0, xsh);
                xsh = _MIN(imgsyn->ncol - 1, xsh);
                ysh = _MAX(0, ysh);
                ysh = _MIN(imgsyn->nrow - 1, ysh);

                if(xsh>= 0 && xsh< _pTree->ncol &&
                        ysh>= 0 && ysh< _pTree->nrow){
                    tR = ((float) imgsyn->red[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTR;
                    imgsyn->red[ysh*imgsyn->ncol + xsh]   = (int) tR;/*(int)rint((double) tR); */

                    tG = ((float) imgsyn->green[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTG;
                    imgsyn->green[ysh*imgsyn->ncol + xsh] = (int) tG; /*(int)rint((double) tG); */

                    tB = ((float) imgsyn->blue[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTB;
                    imgsyn->blue[ysh*imgsyn->ncol + xsh]  = (int) tB; /*(int)rint((double) tB);*/
                }

            }
    }

    TR  = ((Info*)(pShape->data))->r + tempx*0;
    TG  = ((Info*)(pShape->data))->g + tempx*0;
    TB  = ((Info*)(pShape->data))->b + tempx*0;

    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++)
        {
            if(imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] == 0)
                continue;

            if(*mcolor == 1){
                TR  = imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x];
                TG  = imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x];
                TB  = imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x];
            }

            xr = (x - x0tempTransformed)*cos(rotation) + (y - y0tempTransformed)*sin(rotation);
            yr = (y - y0tempTransformed)*cos(rotation) - (x - x0tempTransformed)*sin(rotation);

            xi = floor(xShift + x0tempTransformed + xr);
            yi = floor(yShift + y0tempTransformed + yr);

            tR = ((float) imgsyn->red[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*TR;
            imgsyn->red[yi*_pTree->ncol + xi] = (int)rint((double) tR);

            tG = ((float) imgsyn->green[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*TG;
            imgsyn->green[yi*_pTree->ncol + xi] = (int)rint((double) tG);

            tB = ((float) imgsyn->blue[yi*_pTree->ncol + xi])*ALPHA + (1-ALPHA)*TB;
            imgsyn->blue[yi*_pTree->ncol + xi] = (int)rint((double) tB);
        }
#else
    int i, x, y, iKer, jKer, KerSize, MedSize, xKer, yKer, numMedain;
    float xi, yi, xr, yr, xt, yt, x0temp, y0temp, x0tempDict, y0tempDict;
    float tR, tG, tB, ALPHA, tr, tg, tb, BETA, TR, TG, TB;
    float xBound_l, xBound_r, yBound_t, yBound_b;
    float SCALEx, SCALEy, theta, thetaDict;
    float top, right, left, bottom;
    float a, b, bLimit, tempx;
    ALPHA = *alpha;

    xBound_l = (int) ((Info*)(pShapeDict->data))->boundingbox[0];
    xBound_r = (int) ((Info*)(pShapeDict->data))->boundingbox[1];
    yBound_t = (int) ((Info*)(pShapeDict->data))->boundingbox[2];
    yBound_b = (int) ((Info*)(pShapeDict->data))->boundingbox[3];

    /*=== begin=== Label pShapeDict ========*/
    for(i=0; i < pShapeDict->area; i++)
    {
        x = (pShapeDict->pixels+i)->x;
        y = (pShapeDict->pixels+i)->y;
        imgShapeLabel->gray[y*imgShapeLabel->ncol + x] = 1;
    }

    srand( time(NULL) );
    if(*mcolor == 1)
    {
        tempx = (rand()%10);
        TR = ((Info*)(pShapeDict->data))->r + tempx*0;
        TG = ((Info*)(pShapeDict->data))->g + tempx*0;
        TB = ((Info*)(pShapeDict->data))->b + tempx*0;
    }

    else
    {
        TR  = ((Info*)(pShape->data))->r + tempx*0;
        TG  = ((Info*)(pShape->data))->g + tempx*0;
        TB  = ((Info*)(pShape->data))->b + tempx*0;
    }

    if(*mcolor == 2)
    {
        TR = ((Info*)(pShapeDict->data))->r ;
        TG = ((Info*)(pShapeDict->data))->g;
        TB = ((Info*)(pShapeDict->data))->b;
    }

    if(*equal == 1)
    {
        SCALEx = sqrt( (double) pShape->area / (double) pShapeDict->area );
        SCALEy = SCALEx;
    }
    else
    {
        SCALEx =  sqrt( ((Info*)(pShape->data))->lambda1) /  sqrt(((Info*)(pShapeDict->data))->lambda1 );
        SCALEy =  sqrt( ((Info*)(pShape->data))->lambda2) /  sqrt(((Info*)(pShapeDict->data))->lambda2 );
    }

    a = SCALEx * (xBound_r - xBound_l);
    b = SCALEy * (yBound_b - yBound_t);

    bLimit = (1.0*sqrt(a*a + b*b)/2.0);

    x0temp = ((Info*)(pShape->data))->x0;
    y0temp = ((Info*)(pShape->data))->y0;
    x0tempDict =  ((Info*)(pShapeDict->data))->x0;
    y0tempDict =  ((Info*)(pShapeDict->data))->y0;

    theta  = ((Info*)(pShape->data))->oren;
    thetaDict = ((Info*)(pShapeDict->data))->oren;

    left   = _MAX(0, x0temp - bLimit);
    right  = _MIN(imgsyn->ncol - 1, x0temp + bLimit );
    top    = _MAX(0, y0temp - bLimit);
    bottom = _MIN(imgsyn->nrow - 1, y0temp + bLimit);

    /*=== Transformation ===*/
    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++)
        {
            xt = (x - x0temp);
            yt = (y - y0temp);

            /*====  Rotate to the major axis for scaling    ====*/
            xi = ( (1/SCALEx)* (xt*cos(theta) + yt*sin(theta)));
            yi = ( (1/SCALEy)* (yt*cos(theta) - xt*sin(theta)));

            /*====  Inverse-transformation    ====*/
            xr = ( (xi*cos(thetaDict) - yi*sin(thetaDict)) + x0tempDict);
            yr = ( (yi*cos(thetaDict) + xi*sin(thetaDict)) + y0tempDict);

            if(xr <0 || xr >=imgShapeLabel->ncol ||
                    yr <0 || yr >=imgShapeLabel->nrow)
                continue;

            if( imgShapeLabel->gray[(int)(yr)*imgShapeLabel->ncol + (int)(xr)] == 1 )
            {
                imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] = 1;

                imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x]   =  imgDict->red[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x] =  imgDict->green[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x]  =  imgDict->blue[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
            }
        }

    /*=== Median Filter ===*/
    MedSize = (int)((*median)/2.0);
    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++)
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
            {
                imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] = 0.0;
            }
            else
            {
                imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] = 1.0;
            }
        }

    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++)
        {
            if( imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] == 0  ){
                xt = (x - x0temp);
                yt = (y - y0temp);

                /*====  Rotate to the major axis for scaling    ====*/
                xi = ( (1/SCALEx)* (xt*cos(theta) + yt*sin(theta)));
                yi = ( (1/SCALEy)* (yt*cos(theta) - xt*sin(theta)));

                /*====  Inverse-transformation    ====*/
                xr = ( (xi*cos(thetaDict) - yi*sin(thetaDict)) + x0tempDict);
                yr = ( (yi*cos(thetaDict) + xi*sin(thetaDict)) + y0tempDict);


                imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x]   =  imgDict->red[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x] =  imgDict->green[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x]  =  imgDict->blue[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
            }
            (imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x]) =
                    (int) imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x];

            imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] = 0.0;
        }

    /*=== Add Gaussian Bulr ===*/
    KerSize = (int) ( sqrt( (double) gaussKernel->size) /2.0 );
    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++)
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

                    if( imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] == 0  ){
                        xt = (x - x0temp);
                        yt = (y - y0temp);

                        /*====  Rotate to the major axis for scaling    ====*/
                        xi = ( (1/SCALEx)* (xt*cos(theta) + yt*sin(theta)));
                        yi = ( (1/SCALEy)* (yt*cos(theta) - xt*sin(theta)));

                        /*====  Inverse-transformation    ====*/
                        xr = ( (xi*cos(thetaDict) - yi*sin(thetaDict)) + x0tempDict);
                        yr = ( (yi*cos(thetaDict) + xi*sin(thetaDict)) + y0tempDict);


                        imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x]   =  imgDict->red[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                        imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x] =  imgDict->green[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                        imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x]  =  imgDict->blue[(int)(yr)*imgShapeLabel->ncol + (int)(xr)];
                    }
                }

        }


    if(*relief == 1 && pShape->area > 30)
    {
        float shLambda, shTR, shTG, shTB;
        int xsh, ysh, shiftsh;
        if(pShape->area > 30)
            //       shiftsh = *reliefHeight;
            shiftsh = (*reliefHeight);
        else
            //       shiftsh = (*reliefHeight)*( (float) pShape->area /30.0);
            shiftsh = (*reliefHeight)*( (float) pShape->area /30.0);;
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

                if(*mcolor == 1){

                    TR  = imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x];
                    TG  = imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x];
                    TB  = imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x];
                }

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

    TR  = ((Info*)(pShape->data))->r + tempx*0;
    TG  = ((Info*)(pShape->data))->g + tempx*0;
    TB  = ((Info*)(pShape->data))->b + tempx*0;

    if(*mcolor == 2)
    {
        TR = ((Info*)(pShapeDict->data))->r ;
        TG = ((Info*)(pShapeDict->data))->g;
        TB = ((Info*)(pShapeDict->data))->b;
    }

    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++)
        {
            if(imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] == 0)
                continue;

            BETA = imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x];

            if(*mcolor == 1){
                TR  = imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x];
                TG  = imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x];
                TB  = imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x];
            }

            tr = ((float) imgsyn->red[y*imgsyn->ncol + x])*(1-BETA)   + BETA*TR;
            tg = ((float) imgsyn->green[y*imgsyn->ncol + x])*(1-BETA) + BETA*TG;
            tb = ((float) imgsyn->blue[y*imgsyn->ncol + x])*(1-BETA)  + BETA*TB;

            tR = ((float) imgsyn->red[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tr;
            imgsyn->red[y*imgsyn->ncol + x]   = (int) tR;/*(int)rint((double) tR); */

            tG = ((float) imgsyn->green[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tg;
            imgsyn->green[y*imgsyn->ncol + x] = (int) tG; /*(int)rint((double) tG); */

            tB = ((float) imgsyn->blue[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tb;
            imgsyn->blue[y*imgsyn->ncol + x]  = (int) tB; /*(int)rint((double) tB);*/

            //         printf("%f %f %f; %d %d %d\n", tR, tG, tB, imgsyn->red[y*imgsyn->ncol + x],
            //                imgsyn->green[y*imgsyn->ncol + x], imgsyn->blue[y*imgsyn->ncol + x]);

            imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x]   = 0.0;
            imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] = 0;
        }

    for(i=0; i < pShapeDict->area; i++)
    {
        x = (pShapeDict->pixels+i)->x;
        y = (pShapeDict->pixels+i)->y;
        imgShapeLabel->gray[y*imgShapeLabel->ncol + x] = 0;
    }
#endif
}

/*=== Synthesis by Shape Shaking ===*/
/*=== Before the Shaking, smooth the shape with a gaussian kernel or a median filter ===*/
void TreeOfShapes::synshapeOriginal(Shape pShape,
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

void TreeOfShapes::synshapeOriginal( Shape pShape,
                                     Ccimage imgsyn,
                                     float *alpha,
                                     int *relief,
                                     float *reliefOrentation, float *reliefHeight)
{
    int i, xi, yi;
    float x, y, xr, yr, x0temp, y0temp, ALPHA;
    float xShift, yShift, theta, tR, tG, tB, tb, TR, TG, TB;

    ALPHA = *alpha;
    x0temp = (((Info*)(pShape->data))->x0);
    y0temp = (((Info*)(pShape->data))->y0);

    xShift = (((Info*)(pShape->data))->xShift);
    yShift = (((Info*)(pShape->data))->yShift);
    theta  = (((Info*)(pShape->data))->rotation);


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
        for(i=0; i< pShape->area; i++){

            x = (float)((pShape->pixels+i)->x);
            y = (float)((pShape->pixels+i)->y);

            xr = (x - x0temp)*cos(theta) + (y - y0temp)*sin(theta);
            yr = (y - y0temp)*cos(theta) - (x - x0temp)*sin(theta);

            xi = floor(xShift + x0temp + xr);
            yi = floor(yShift + y0temp + yr);

            shTR = TR* shLambda;
            shTG = TG* shLambda;
            shTB = TB* shLambda;


            xsh = xi + shiftsh*cos( PI*(*reliefOrentation)/180.0 );
            ysh = yi - shiftsh*sin( PI*(*reliefOrentation)/180.0 );
            xsh = _MAX(0, xsh);
            xsh = _MIN(imgsyn->ncol - 1, xsh);
            ysh = _MAX(0, ysh);
            ysh = _MIN(imgsyn->nrow - 1, ysh);

            if(xsh>= 0 && xsh< _pTree->ncol &&
                    ysh>= 0 && ysh< _pTree->nrow){
                tR = ((float) imgsyn->red[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTR;
                imgsyn->red[ysh*imgsyn->ncol + xsh]   = (int) tR;/*(int)rint((double) tR); */

                tG = ((float) imgsyn->green[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTG;
                imgsyn->green[ysh*imgsyn->ncol + xsh] = (int) tG; /*(int)rint((double) tG); */

                tB = ((float) imgsyn->blue[ysh*imgsyn->ncol + xsh])*ALPHA + (1-ALPHA)*shTB;
                imgsyn->blue[ysh*imgsyn->ncol + xsh]  = (int) tB; /*(int)rint((double) tB);*/
            }
        }
    }

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


void TreeOfShapes::sortShapes(Fsignal t2b_index){


    struct timeval start, end;
    gettimeofday(&start, NULL);


    std::priority_queue < std::pair<int, int>, std::deque< std::pair<int, int> > , std::less<std::pair<int, int> > > AreaShapeIDQueue;
    Shape pShape;

    for(int i=0; i< _pTree->nb_shapes; i++) {
        pShape = _pTree->the_shapes + i;
        AreaShapeIDQueue.push( std::make_pair(pShape->area, i));
    }

    int i = 0;
    while(!AreaShapeIDQueue.empty()){

        std::pair<int, int> area_id_pair = AreaShapeIDQueue.top();
        AreaShapeIDQueue.pop();

        t2b_index->values[i++] = area_id_pair.second;
    }
    gettimeofday(&end, NULL);
    double elapsedTime = (end.tv_sec  - start.tv_sec) +
            (end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "priority sortShapes::time elapsed : " << elapsedTime <<" seconds"<< std::endl;
    std::cout << "***************************" << std::endl << std::endl << std::endl;
}

/*Index the tree in coast-to-fine order*/
/*Buuble sorting of integer array */
void TreeOfShapes::Bubble( Fsignal t2b_index)
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
Fsignal TreeOfShapes::sgauss(float *std,
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
Fsignal TreeOfShapes::Sgauss(float *std, Fsignal out, int *size)
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

void TreeOfShapes::get_shapes_truearea(Shape s, Shape root,
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


void TreeOfShapes::filter_shapes( Cfimage out,
                                  char *local,
                                  float *eps)
{
    std::cout <<"TreeOfShapes::filter_shapes()::begin"<< std::endl;
    Fimage Fv;
    float lgeo,step,hstep,tcleanup,farea, std;
    float *gray,*red,*green,*blue;
    int prec,visit,nrow,ncol,i,j,indexshape,*truearea;
    char all;
    Shapes tree;
    Shape s,t;

    tree = mw_new_shapes();

    Fv = mw_change_fimage(NULL,_imgin->nrow,_imgin->ncol);

    if(!Fv ) mwerror(FATAL,1,"Not enough memory.\n");


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
    ll_boundaries2(Fv,eps,NULL,&step,&prec,&std,&hstep,NULL,&visit,NULL,NULL,tree);

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
            out->red[i+j*ncol]= red[indexshape];
            out->green[i+j*ncol]= green[indexshape];
            out->blue[i+j*ncol]= blue[indexshape];
        }

    free(gray);
    free(red);
    free(green);
    free(blue);
    free(truearea);
    std::cout <<"TreeOfShapes::filter_shapes()::end"<< std::endl;
}


/*=== Filtering the image ===*/
void TreeOfShapes::filter_image(int *ns,
                                float *threshold,
                                int *mpixel,
                                int *maxpixel){
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

        if(pShape->area <= *mpixel || *maxpixel < pShape->area
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
        } else
            pShape->removed = 0;

        if(i ==0 )
            pShape->removed = 0;
    }
}

int TreeOfShapes::random_number(int *M)
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

    mw_delete_fsignal(pb);
    return select_i;
}

void TreeOfShapes::computeKdTree(float average_r, float average_g, float average_b ){
    _use_kdtree = false;

    if( _pTree->nb_shapes > 1000 ){
        std::cout <<"Building KD tree" << std::endl;
        _annTree = BasicANNkdTree(6);
        Shape pShape;


        std::vector<std::vector<float>> values;
        for(int i=1; i < _pTree->nb_shapes; i++)  {
            pShape = _pTree->the_shapes + i;

            float lambda1, lambda2;

            float elong, kappa, elongDict, kappaDict, Dist, minDist;
            float sca, scaDict, pa;
            int mn, temp;



            lambda1 = ((Info*)(pShape->data))->lambda1;
            lambda2 = ((Info*)(pShape->data))->lambda2;
            elongDict = lambda2 / lambda1;
            kappaDict = ((float) pShape->area)/(sqrt(lambda2*lambda1)*4*PI);
            scaDict = ((float) pShape->area);

            std::vector<float> value;
            value.push_back(elongDict);
            value.push_back(kappaDict);
            value.push_back(scaDict);
            value.push_back(((Info*)(pShape->data))->r/average_r);
            value.push_back(((Info*)(pShape->data))->g/average_g);
            value.push_back(((Info*)(pShape->data))->b/average_b);
            values.push_back(value);

            //        Dist = pow((elong - elongDict), 2.0) +
            //                pow((kappa - kappaDict), 2.0) +
            //                pow((1 - _MIN(sca/scaDict, scaDict/sca)), 2.0);


            //            Dist +=  ( pow((1 - _MIN(((Info*)(pShape->data))->r/((Info*)(pShapeDict->data))->r,
            //                                     ((Info*)(pShapeDict->data))->r/((Info*)(pShape->data))->r)), 2.0) +
            //                       pow((1 - _MIN(((Info*)(pShape->data))->g/((Info*)(pShapeDict->data))->g,
            //                                     ((Info*)(pShapeDict->data))->g/((Info*)(pShape->data))->g)), 2.0) +
            //                       pow((1 - _MIN(((Info*)(pShape->data))->b/((Info*)(pShapeDict->data))->b,
            //                                     ((Info*)(pShapeDict->data))->b/((Info*)(pShape->data))->b)), 2.0) )/3.0;


            //         Dist =  sqrt(Dist);
        }


        _annTree.build(values);

        _use_kdtree = true;
    }

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
Shape TreeOfShapes::selectShapeDict(Shape pShape,
                                    float *paDict,
                                    int *randS,
                                    int &index,
                                    float average_r,
                                    float average_g,
                                    float average_b)
{
    Shape pShapeDict, pShapeTemp;
    float lambda1, lambda2;

    float elong, kappa, elongDict, kappaDict, Dist, minDist;
    float sca, scaDict, pa;
    int i, mn, temp;

    lambda1 = ((Info*)(pShape->data))->lambda1;
    lambda2 = ((Info*)(pShape->data))->lambda2;
    elong = lambda2 / lambda1;
    kappa = ((float) pShape->area)/(sqrt(lambda2*lambda1)*4*PI);
    sca = ((float) pShape->area);


    index = 1; minDist = 10000.0;
    if(*randS == 0)
    {
        temp = _pTree->nb_shapes -1;
        index = random_number(&temp) + 1;
    }
    else
    {
        if (!_use_kdtree){
            for(i= 1; i<_pTree->nb_shapes; i++)
            {
                pShapeDict = _pTree->the_shapes + i;

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
        }else{
            std::vector<float> value;
            value.push_back(elong);
            value.push_back(kappa);
            value.push_back(sca);
            value.push_back(((Info*)(pShape->data))->r/average_r);
            value.push_back(((Info*)(pShape->data))->g/average_g);
            value.push_back(((Info*)(pShape->data))->b/average_b);

            int k = std::min(10, _pTree->nb_shapes);
            ANNidxArray neighbors = new ANNidx[ k ];
            ANNdistArray neighbors_sqr_dists= new ANNdist[ k ];
            _annTree.knearest(value, k, neighbors, neighbors_sqr_dists);

            for(i= 1; i<k; i++)
            {
                pShapeDict = _pTree->the_shapes + (int)neighbors[i];

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
                    index = (int)neighbors[i];
                }
            }

            //index = _annTree.nearest(value);
        }

    }


    pShapeDict = _pTree->the_shapes + index;
    //    printf("index: %d \n", index);

    return pShapeDict;
}

Shape TreeOfShapes::getShape(int index){
    if( index < _pTree->nb_shapes ){
        return _pTree->the_shapes + index;
    }
    return NULL;
}
void TreeOfShapes::compute_tree( TOSParameters tosParameters, bool dictionary ){
    struct timeval start, end;
    gettimeofday(&start, NULL);

    double elapsedTime = 0., current_time = 0.;
    std::cout <<"TreeOfShapes::Abstraction started"<< std::endl;
    if( !_tree_computed || _tosParameters.color_sketch != tosParameters.color_sketch ||
            ( _tosParameters.color_sketch == 1 & _tosParameters.eps != tosParameters.eps ) ){

        if( tosParameters.color_sketch == 1 ){

            char local=NULL ; //"local boundaries, default NULL",

            Cfimage out = mw_change_cfimage(NULL,_imgin->nrow,_imgin->ncol);
            filter_shapes(out,
                          &local,
                          &tosParameters.eps);
            std::cout << "color_sketch done" <<std::endl;
            init(out, _pTree);
            mw_delete_cfimage(out);

        } else if( _tosParameters.color_sketch == 1 & tosParameters.color_sketch == 0 || !_tree_computed){
            init(_imgin, _pTree);
        }
        _tree_computed = true;
        _tree_recomputed = true;
        _use_kdtree = false;
        if( _large_to_small_index_computed )
            mw_delete_fsignal(_large_to_small_index);

        for (std::map<int, Fsignal>::iterator it = _dictionary_selections.begin(); it !=  _dictionary_selections.end(); ++it){
            mw_delete_fsignal( it->second );
        }
        _dictionary_selections.clear();

        _large_to_small_index_computed = false;

        gettimeofday(&end, NULL);
        current_time = (end.tv_sec  - start.tv_sec) +
                (end.tv_usec - start.tv_usec) / 1.e6;
        std::cout << std::endl<<"TreeOfShapes::Tree of " << _pTree->nb_shapes << " shapes computed : " << current_time - elapsedTime <<" seconds"<< std::endl;
        elapsedTime = current_time;
    }

    if( dictionary ){
        /*=======================================*/
        /*====   Compute shape attribute     === */
        /*=======================================*/

        std::cout << "Compute shape attributes" << std::endl;
        compute_shape_attribute();
        tree_boundingbox();
    }
}


void TreeOfShapes::save_shapes( QString folder_name, bool average_color ){

    Fsignal t2b_index = NULL;
    Shape pShape;
    std::cout << "folder_name " << folder_name.toStdString() << std::endl;
    int x, y, i, j ;

    if  ( ((t2b_index = mw_new_fsignal()) == NULL) ||
          (mw_alloc_fsignal(t2b_index,_pTree->nb_shapes) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");

    top2bottom_index_tree(t2b_index);
    for(i= 1; i<_pTree->nb_shapes; i++)
    {
        QString image_name(folder_name);
        QImage result_image( QSize(_imgin->ncol, _imgin->nrow), QImage::Format_ARGB32);
        result_image.fill(QColor(255,255,255, 0.));
        pShape = _pTree->the_shapes + (int)t2b_index->values[i];

        for(j=0; j< pShape->area; j++){

            x = ((pShape->pixels+j)->x);
            y = ((pShape->pixels+j)->y);

            QColor color;
            if( average_color )
                color = QColor (((Info*)(pShape->data))->r, ((Info*)(pShape->data))->g, ((Info*)(pShape->data))->b);
            else
                color = QColor((float) _imgin->red[y*_pTree->ncol + x], (float) _imgin->green[y*_pTree->ncol + x], (float) _imgin->blue[y*_pTree->ncol + x]);


            result_image.setPixel(x, y, qRgb(color.red(), color.green(), color.blue()));
        }
        image_name.append(QString ( "/" ));
        QString num;
        num.setNum(i);
        num.append(QString(".png"));

        image_name.append(num);

        result_image.save(image_name);
    }

    mw_delete_fsignal(t2b_index);
    std::cout << "Saved" << std::endl;
}

std::vector<QImage> TreeOfShapes::render_shape_by_shape(TOSParameters tosParameters, TreeOfShapes *tosDictionary, DictionaryParameters dictionaryParameters ){

}

Ccimage TreeOfShapes::render(TOSParameters tosParameters, bool &tree_recomputed, TreeOfShapes *tosDictionary, DictionaryParameters dictionaryParameters, bool save_shapes, QString folder_name ){

    struct timeval start, end;
    gettimeofday(&start, NULL);

    double elapsedTime = 0., current_time = 0.;
    std::cout <<"TreeOfShapes::Abstraction started"<< std::endl;

    compute_tree(tosParameters, false);

    tree_recomputed = _tree_recomputed;

    printf("---0---- syntexturecolor ..........\n");

    fflush(stdout);
    //@Declare variables here.==========================
    int i,j, k, minArea, mn, nsize;
    Shape pShape, pShapeTemp, pShapeDict;

    float pa, fzero, ALPHA;
    Cimage imgShapeLabel;
    Fimage imgShapeBlur;
    Fsignal t2b_index = NULL, gaussKernel;
    ALPHA = 0.0;

    Ccimage imgsyn = mw_change_ccimage(imgsyn, _imgin->nrow, _imgin->ncol);
    if  ( ((imgShapeLabel = mw_new_cimage()) == NULL) ||
          (mw_alloc_cimage(imgShapeLabel, _imgin->nrow, _imgin->ncol) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");

    if  ( ((imgShapeBlur = mw_new_fimage()) == NULL) ||
          (mw_alloc_fimage(imgShapeBlur, _imgin->nrow, _imgin->ncol) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");

    imgsyn = mw_change_ccimage(imgsyn, _imgin->nrow, _imgin->ncol);
    imgShapeLabel = mw_change_cimage(imgShapeLabel, _imgin->nrow, _imgin->ncol);
    imgShapeBlur  = mw_change_fimage(imgShapeBlur, _imgin->nrow, _imgin->ncol);

    /*=======================================*/
    /*=== Compute FLST on Intensity image ===*/
    /*=======================================*/

    if  ( ((t2b_index = mw_new_fsignal()) == NULL) ||
          (mw_alloc_fsignal(t2b_index,_pTree->nb_shapes) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");

    for( i= 0; i< _pTree->ncol; i++)
        for( j= 0; j< _pTree->nrow; j++)
        {
            imgsyn->red[j*_pTree->ncol + i] = 1;
            imgsyn->green[j*_pTree->ncol + i] = 1;
            imgsyn->blue[j*_pTree->ncol + i] = 1;
        }


    /*=======================================*/
    /*==   Image filtering                 ==*/
    /*=======================================*/

    std::cout << "Image filtering" << std::endl;
    int max_area = tosParameters.maxarea;
    if( _tree_recomputed ) max_area = INT_MAX;
    filter_image(&tosParameters.ns,&tosParameters.threshold, &tosParameters.mpixel, &max_area);

    gettimeofday(&end, NULL);
    current_time = (end.tv_sec  - start.tv_sec) +
            (end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << std::endl<<"TreeOfShapes::Image filtered: " << current_time - elapsedTime<<" seconds"<< std::endl;
    elapsedTime = current_time;
    /*=======================================*/
    /*==   Select the rendering order *\/ ===*/
    /*=======================================*/
    std::cout << "Rendering order " << tosParameters.order <<std::endl;
    if(tosParameters.order == 0)
        top2bottom_index_tree(t2b_index);
    else if(tosParameters.order == 1){
        if( !_large_to_small_index_computed ){

            if  ( ((_large_to_small_index = mw_new_fsignal()) == NULL) ||
                  (mw_alloc_fsignal(_large_to_small_index,_pTree->nb_shapes) == NULL) )
                mwerror(FATAL,1,"Not enough memory.\n");
            sortShapes(_large_to_small_index);
            _large_to_small_index_computed = true;
        }
        mw_copy_fsignal_values(_large_to_small_index, t2b_index);

    } else if(tosParameters.order == 2)
        random_tree_order(t2b_index);
    else
        top2bottom_index_tree(t2b_index);

    gettimeofday(&end, NULL);
    current_time = (end.tv_sec  - start.tv_sec) +
            (end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << std::endl<<"TreeOfShapes::Shape sorted: " << current_time - elapsedTime <<" seconds"<< std::endl;
    elapsedTime = current_time;

    /*=======================================*/
    /*==  Add a random shift to each shape ==*/
    /*=======================================*/
    if(tosParameters.smodel == 0)
        random_shift_shape(&tosParameters.shift, &tosParameters.theta);
    else if(tosParameters.smodel == 1)
        adaptive_shift_shape(&tosParameters.shift, &tosParameters.theta);
    else
        adaptive_shift_shape(&tosParameters.shift, &tosParameters.theta);

    gettimeofday(&end, NULL);
    current_time = (end.tv_sec  - start.tv_sec) +
            (end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << std::endl<<"TreeOfShapes::Shaking computed: " << current_time - elapsedTime <<" seconds"<< std::endl;
    elapsedTime = current_time;

    /*=======================================*/
    /*==  Compute a Gaussian kernel        ==*/
    /*=======================================*/
    mw_clear_cimage(imgShapeLabel,0);
    mw_clear_fimage(imgShapeBlur,0.0);
    if  ( ((gaussKernel = mw_new_fsignal()) == NULL) ||
          (mw_alloc_fsignal(gaussKernel, tosParameters.kerSize*tosParameters.kerSize) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");
    gaussKernel = Sgauss(&tosParameters.kerStd, gaussKernel, &tosParameters.kerSize);

    gettimeofday(&end, NULL);
    current_time = (end.tv_sec  - start.tv_sec) +
            (end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << std::endl<<"TreeOfShapes::Kernel generated: " << current_time - elapsedTime<<" seconds"<< std::endl;
    elapsedTime = current_time;


    Fsignal dictionary_correspondance;
    if  ( ((dictionary_correspondance = mw_new_fsignal()) == NULL) ||
          (mw_alloc_fsignal(dictionary_correspondance,_pTree->nb_shapes) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");

    bool correspondance_computed = false ;
    //Check if correspondance computed for dictionary
    if( tosParameters.model == 4 ){
        std::map<int, Fsignal>::iterator it = _dictionary_selections.find( tosDictionary->getTreeId() );

        if( it != _dictionary_selections.end() ){
            mw_copy_fsignal_values(it->second, dictionary_correspondance);
            correspondance_computed = true;
        } else {
            if  ( ((_dictionary_selections[ tosDictionary->getTreeId() ] = mw_new_fsignal()) == NULL) ||
                  (mw_alloc_fsignal(_dictionary_selections[ tosDictionary->getTreeId() ],_pTree->nb_shapes) == NULL) )
                mwerror(FATAL,1,"Not enough memory.\n");
            mw_clear_fsignal(_dictionary_selections[ tosDictionary->getTreeId() ],-1.0);
            tosDictionary->computeKdTree(_average_r, _average_g, _average_b);
        }
    }

    /*=======================================*/
    /*==  Shape Shaking Filtering          ==*/
    /*=======================================*/
    for(i=0; i < _pTree->nb_shapes; i++)  {
        pShape = _pTree->the_shapes + (int)t2b_index->values[i];

        if((int)t2b_index->values[i] == 0 ) {
            float r=((Info*)(pShape->data))->r, g= ((Info*)(pShape->data))->g, b= ((Info*)(pShape->data))->b;

            if (tosParameters.model == 4 && (dictionaryParameters.mcolor == 1 || dictionaryParameters.mcolor ==2)){
                pShapeDict = tosDictionary->getShape(0);
                ((Info*)(pShape->data))->r = ((Info*)(pShapeDict->data))->r;
                ((Info*)(pShape->data))->g = ((Info*)(pShapeDict->data))->g;
                ((Info*)(pShape->data))->b = ((Info*)(pShapeDict->data))->b;
            }
            synshapeRect(pShape, imgsyn, &ALPHA, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
            if (tosParameters.model == 4 && dictionaryParameters.mcolor == 1 || dictionaryParameters.mcolor ==2){
                ((Info*)(pShape->data))->r = r;
                ((Info*)(pShape->data))->g = g;
                ((Info*)(pShape->data))->b = b;
            }
        } else if(pShape->removed != 1){

            /*=======================================*/
            /*====     Attribute filtering       === */
            /*=======================================*/
            mn=3;
            pShapeTemp =  m_order_parent(pShape, &mn);
            pa = ((float) pShape->area)/((float) pShapeTemp->area);
            //     pa = ((float) pShape->area)/((float) (pTree->ncol*pTree->nrow));

            if(pa < tosParameters.kappa)
                continue;

            if(save_shapes){

                for( int x= 0; x< _pTree->ncol; x++)
                    for( int y= 0; y< _pTree->nrow; y++)
                    {
                        imgsyn->red[y*_pTree->ncol + x] = 1;
                        imgsyn->green[y*_pTree->ncol + x] = 2;
                        imgsyn->blue[y*_pTree->ncol + x] = 3;
                    }
            }
            /*=======================================*/
            /*===  Select the rendering model    === */
            /*=======================================*/
            if(tosParameters.model == 0){
                if(tosParameters.blur == 1)
                    synshapeOriginal(pShape, imgsyn, imgShapeLabel, imgShapeBlur, gaussKernel, &tosParameters.median, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
                else
                    synshapeOriginal(pShape, imgsyn, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
            } else if(tosParameters.model == 1)
                if(tosParameters.blur == 1)
                    synshapeEllipse(pShape, imgsyn, imgShapeLabel, imgShapeBlur, gaussKernel, &tosParameters.median, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
                else
                    synshapeEllipse(pShape, imgsyn, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
            else if(tosParameters.model == 2){
                if(tosParameters.blur == 1)
                    synshapeRect(pShape, imgsyn, imgShapeLabel, imgShapeBlur, gaussKernel, &tosParameters.median, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
                else
                    synshapeRect(pShape, imgsyn, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
            } else if(tosParameters.model == 3){
                if(tosParameters.blur == 1)
                    synshapeCircle(pShape, imgsyn, imgShapeLabel, imgShapeBlur, gaussKernel, &tosParameters.median, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
                else
                    synshapeCircle(pShape, imgsyn, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
            } else if(tosParameters.model == 4){
                //Dictionary
                Cfimage imgDict = tosDictionary->getCfImage();
                Cfimage imgShapeColorSyn;
                Cimage imgShapeLabel, imgShapeLabelSyn;
                Fimage imgShapeBlurSyn;

                if  ( ((imgShapeColorSyn = mw_new_cfimage()) == NULL) ||
                      (mw_alloc_cfimage(imgShapeColorSyn, _imgin->nrow, _imgin->ncol) == NULL) )
                    mwerror(FATAL,1,"Not enough memory.\n");
                if  ( ((imgShapeLabel = mw_new_cimage()) == NULL) ||
                      (mw_alloc_cimage(imgShapeLabel, imgDict->nrow, imgDict->ncol) == NULL) )
                    mwerror(FATAL,1,"Not enough memory.\n");

                if  ( ((imgShapeLabelSyn = mw_new_cimage()) == NULL) ||
                      (mw_alloc_cimage(imgShapeLabelSyn, _imgin->nrow, _imgin->ncol) == NULL) )
                    mwerror(FATAL,1,"Not enough memory.\n");
                if  ( ((imgShapeBlurSyn = mw_new_fimage()) == NULL) ||
                      (mw_alloc_fimage(imgShapeBlurSyn, _imgin->nrow, _imgin->ncol) == NULL) )
                    mwerror(FATAL,1,"Not enough memory.\n");

                imgShapeColorSyn = mw_change_cfimage(imgShapeColorSyn, _imgin->nrow, _imgin->ncol);
                imgShapeLabel = mw_change_cimage(imgShapeLabel, imgDict->nrow, imgDict->ncol);
                imgShapeLabelSyn = mw_change_cimage(imgShapeLabelSyn, _imgin->nrow, _imgin->ncol);
                imgShapeBlurSyn  = mw_change_fimage(imgShapeBlurSyn, _imgin->nrow, _imgin->ncol);

                mw_clear_cimage(imgShapeLabel,0);
                mw_clear_cimage(imgShapeLabelSyn,0);
                mw_clear_fimage(imgShapeBlurSyn,0.0);

                if( correspondance_computed ){
                    int shape_id = (int)dictionary_correspondance->values[(int)t2b_index->values[i]];
                    if ( shape_id >= 0 )
                        pShapeDict = tosDictionary->getShape(shape_id);
                    if ( pShapeDict == NULL || shape_id < 0 ){

                        pShapeDict = tosDictionary->selectShapeDict(pShape, &dictionaryParameters.kappaDict, &dictionaryParameters.randS, shape_id, _average_r, _average_g, _average_b);
                        _dictionary_selections[ tosDictionary->getTreeId() ]->values[(int)t2b_index->values[i]] = shape_id;
                    }
                } else {
                    int shape_id;
                    pShapeDict = tosDictionary->selectShapeDict(pShape, &dictionaryParameters.kappaDict, &dictionaryParameters.randS, shape_id, _average_r, _average_g, _average_b);
                    dictionary_correspondance->values[(int)t2b_index->values[i]] = shape_id;
                }
                if( tosParameters.blur == 1 )
                    synShapeDict( pShapeDict, pShape, imgsyn, imgDict, imgShapeColorSyn, imgShapeLabel, imgShapeLabelSyn, imgShapeBlurSyn, gaussKernel, &tosParameters.median, &tosParameters.alpha,
                                  &dictionaryParameters.equal, &dictionaryParameters.mcolor,&tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
                else
                    synShapeDict( pShapeDict, pShape, imgsyn, imgDict, imgShapeColorSyn, imgShapeLabel, imgShapeLabelSyn, &tosParameters.alpha,
                                  &dictionaryParameters.equal, &dictionaryParameters.mcolor,&tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);

                mw_delete_cfimage(imgShapeColorSyn);
                mw_delete_cimage(imgShapeLabel);
                mw_delete_cimage(imgShapeLabelSyn);
                mw_delete_fimage(imgShapeBlurSyn);
            } else if(tosParameters.model == 5){
                mn = rand()%10;
                if(mn<=4)
                    synshapeEllipse(pShape, imgsyn, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
                else if(mn<=9)
                    synshapeRect(pShape, imgsyn, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
                else
                    synshapeOriginal(pShape, imgsyn, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
            } else if(tosParameters.model == 6){
                if( ((float)pShape->area) / ((float)(imgsyn->nrow*imgsyn->ncol)) > 0.2 ||
                        ((float) pShape->area)/(sqrt(((Info*)(pShape->data))->lambda1 * ((Info*)(pShape->data))->lambda2)*4*PI) <0.5)
                    synshapeOriginal(pShape, imgsyn, imgShapeLabel, imgShapeBlur, gaussKernel, &tosParameters.median, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
                else
                    synshapeEllipse(pShape, imgsyn, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);
            }
            else
                synshapeEllipse(pShape, imgsyn, &tosParameters.alpha, &tosParameters.relief, &tosParameters.reliefOrientation, &tosParameters.reliefHeight);

            if( save_shapes ){

                QString image_name(folder_name);
                QImage result_image( QSize(_imgin->ncol, _imgin->nrow), QImage::Format_ARGB32);
                result_image.fill(QColor(255,255,255, 0.));


                for( int y= 0; y< imgsyn->nrow; y++)
                    for( int x= 0; x< imgsyn->ncol; x++)
                    {
                        int comp = y*imgsyn->ncol + x;

                        QColor color (imgsyn->red[comp], imgsyn->green[comp], imgsyn->blue[comp]);

                        //color = QColor (_imgin->red[comp], _imgin->green[comp], _imgin->blue[comp]);

                        if( imgsyn->red[comp] != 1 || imgsyn->green[comp] != 2 || imgsyn->blue[comp] != 3)
                            result_image.setPixel(x, y , qRgb(color.red(), color.green(), color.blue()));
                    }
                QString num;
                num.setNum(i);
                num.append(QString(".png"));
                image_name.append(QString ( "/" ));
                image_name.append(num);

                result_image.save(image_name);
            }
        }
    }


    gettimeofday(&end, NULL);
    elapsedTime = (end.tv_sec  - start.tv_sec) +
            (end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "TreeOfShapes::time elapsed : " << elapsedTime <<" seconds"<< std::endl;
    std::cout << "***************************" << std::endl << std::endl << std::endl;

    mw_delete_fsignal(t2b_index);
    mw_delete_cimage(imgShapeLabel);
    mw_delete_fimage(imgShapeBlur);
    mw_delete_fsignal(gaussKernel);

    if( tosParameters.model == 4 && !correspondance_computed ){
        mw_copy_fsignal_values( dictionary_correspondance, _dictionary_selections[ tosDictionary->getTreeId() ]);
    }

    _tosParameters = tosParameters;

    _tree_recomputed = false;


    return imgsyn;
}
