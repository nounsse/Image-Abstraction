#ifndef TREEOFSHAPES_H
#define TREEOFSHAPES_H

#include <QImage>

#include "mw3.h"
#include "mw3-modules.h"
#include "stdio.h"
#include "stdlib.h"
#include "tree_of_shapes.h"
#include "jmkdtree.h"
#include <cfloat>
#include <map>
#include <QColor>


class TreeOfShapes
{
public:

    static int _tree_count;

    TreeOfShapes( Cfimage imageIn );
    TreeOfShapes( Cfimage imageIn, Cfimage texture_image );
    ~TreeOfShapes();

    Ccimage render(TOSParameters tosParameters, bool &tree_recomputed, TreeOfShapes *tosDictionary=NULL, DictionaryParameters dictionaryParameters=getDefaultDictionaryParameters(), bool save_shapes=false, QString folder_name=QString()  );
    void save_shapes( QString folder_name, bool average_color );
    std::vector<QImage> render_shape_by_shape(TOSParameters tosParameters, TreeOfShapes *tosDictionary=NULL, DictionaryParameters dictionaryParameters=getDefaultDictionaryParameters()  );
    void compute_tree( TOSParameters tosParameters, bool dictionary=false );
    void computeKdTree(float average_r, float average_g, float average_b );
    Cfimage getCfImage(){ if( _texture_image_loaded ) return _texture_image; else return _imgin; }
    Shape selectShapeDict(Shape pShape,
                          float *paDict,
                          int *randS,
                          int &index,
                          float average_r, float average_g, float average_b);

    Shape getShape(int index);
    int getTreeId(){ return _tree_id; }
    int getMaxArea(){ return _maxArea; }
    void getTreeInfo(std::vector<QPoint> &positions, std::vector<QColor> &colors,
                     std::vector< std::vector< std::pair<int, int> > > &pixels, std::vector<int> & heights, std::vector<std::pair<int, int> > &edges);
protected:
    bool _tree_computed;
    bool _texture_image_loaded;
    bool _tree_recomputed;
    bool _large_to_small_index_computed;
    bool _use_kdtree;
    BasicANNkdTree _annTree;
    int _tree_id;
    int _maxArea;

    Cfimage _imgin;
    Cfimage _texture_image;
    Shapes _pTree;
    Fimage _NormOfDu;
    Fsignal _large_to_small_index;
    std::map<int, Fsignal> _dictionary_selections;

    TOSParameters _tosParameters;
    DictionaryParameters _dictionaryParameters;

    float _average_r;
    float _average_g;
    float _average_b;

    void init(Cfimage inputImg, Shapes &pTree);

    /* This removes the shapes from the tree associated to pFloatImageInput
    that are too small (threshold *pMinArea). As a consequence all the remaining
    shapes of pFloatImageOutput are of area larger or equal than *pMinArea */

    void mw_fgrain_side(int *pMinArea, int sideflag);

    void sortShapes(Fsignal t2b_index);
    /*in and out must be allocated*/
    void fgrain_side(int MinArea, float *in, int nx, int ny, float *out, int sideflag);

    Shape m_order_parent(Shape pShape,
                         int *mn,
                         bool dict = false);
    void Order(Fsignal t2b_index,
               int *p, int *q);
    void Bubble( Fsignal t2b_index);

    void shape_orilam(Shape pShape,float *out_ori, float *out_e, float *out_k);
    void shape_orilam(Shape pShape, float *out_ori, float *out_e, float *out_k, float *pX0, float *pY0);
    void compute_shape_attribute();

    float min_contrast(Shape pShape);
    void synshapeCircle(Shape pShape,
                        Ccimage imgsyn,
                        float *alpha,
                        int *relief,
                        float *reliefOrentation, float *reliefHeight);
    void synshapeCircle(Shape pShape,
                        Ccimage imgsyn,
                        Cimage imgShapeLabelSyn,
                        Fimage imgShapeBlurSyn,
                        Fsignal gaussKernel,
                        int *median,
                        float *alpha,
                        int *relief,
                        float *reliefOrentation, float *reliefHeight);
    void synshapeEllipse(Shape pShape,
                         Ccimage imgsyn,
                         float *alpha,
                         int *relief,
                         float *reliefOrentation, float *reliefHeight);
    void synshapeEllipse(Shape pShape,
                         Ccimage imgsyn,
                         Cimage imgShapeLabelSyn,
                         Fimage imgShapeBlurSyn,
                         Fsignal gaussKernel,
                         int *median,
                         float *alpha,
                         int *relief,
                         float *reliefOrentation, float *reliefHeight);
    void synshapeRect(Shape pShape,
                      Ccimage imgsyn,
                      float *alpha,
                      int *relief,
                      float *reliefOrentation, float *reliefHeight);
    void synshapeRect(Shape pShape,
                      Ccimage imgsyn,
                      Cimage imgShapeLabelSyn,
                      Fimage imgShapeBlurSyn,
                      Fsignal gaussKernel,
                      int *median,
                      float *alpha,
                      int *relief,
                      float *reliefOrentation, float *reliefHeight);
    void synshapeOriginal(Shape pShape,
                          Ccimage imgsyn,
                          Cimage imgShapeLabelSyn,
                          Fimage imgShapeBlurSyn,
                          Fsignal gaussKernel,
                          int *median,
                          float *alpha,
                          int *relief,
                          float *reliefOrentation, float *reliefHeight);
    void synshapeOriginal( Shape pShape,
                           Ccimage imgsyn,
                           float *alpha,
                           int *relief,
                           float *reliefOrentation, float *reliefHeight);

    void top2bottom_index_tree(Fsignal t2b_index);
    void random_tree_order(Fsignal t2b_index);

    void random_leaf(Fsignal leaf,
                     Fsignal t2b_index,
                     int *k_ind);

    void random_shift_shape(float *shift, float * theta);
    void adaptive_shift_shape(float *shift,
                              float *theta);
    float mean_contrast(Shape pShape);

    Fsignal sgauss(float *std,
                   Fsignal out,
                   int *size);
    Fsignal Sgauss(float *std, Fsignal out, int *size);

    float peri_shape(Shape pShape);

    void compute_shape_attribute(int *ns);
    void filter_image(int *ns,
                      float *alpha,
                      int *mpixel,
                      int *maxpixel);
    void filter_shapes( Cfimage out,
                        char *local,
                        float *eps);
    void get_shapes_truearea(Shape s, Shape root,
                             int *truearea);

    int random_number(int *M);

    void shape_boundingbox(Shape pShape);
    void tree_boundingbox();

    void synShapeDict(Shape pShapeDict, Shape pShape,
                      Ccimage imgsyn,
                      Cfimage imgDict, Cfimage imgShapeColorSyn,
                      Cimage imgShapeLabel, Cimage imgShapeLabelSyn,
                      Fimage imgShapeBlurSyn,
                      Fsignal gaussKernel,
                      int *median,
                      float *alpha,
                      int *equal, int *mcolor, int *relief,
                      float *reliefOrentation, float *reliefHeight);
    void synShapeDict(Shape pShapeDict, Shape pShape,
                      Ccimage imgsyn,
                      Cfimage imgDict, Cfimage imgShapeColorSyn,
                      Cimage imgShapeLabel, Cimage imgShapeLabelSyn,
                      float *alpha,
                      int *equal, int *mcolor, int *relief,
                      float *reliefOrentation, float *reliefHeight);
};

#endif // TREEOFSHAPES_H
