/*--------------------------- MegaWave2 Module -----------------------------*/
/* mwcommand
  name = {fgrain_side};
  version = {"1.2"};
  author = {"Pascal Monasse, Frederic Guichard, G. Facciolo"};
  function = {"Grain filter of an image"};
  usage = {
    'a': [min_area=20]-> pMinArea   "Min area of grains we keep",
    image_in -> pFloatImageInput    "Input fimage",
    image_out <- pFloatImageOutput  "Output fimage"
    };
*/
/*----------------------------------------------------------------------
 v1.2 (04/2007): simplified header (LM)
 v1.3 (2013): portable version (GF)
 v1.4 (2014): fgrain_side: removes only upper (bright) or lower (dark) courves
----------------------------------------------------------------------*/

#ifndef ABSTRACTIONGRAINPROCESS_H
#define ABSTRACTIONGRAINPROCESS_H

#include <QImage>

#include "mw3.h"
#include "mw3-modules.h"
#include "stdio.h"
#include "stdlib.h"
#include "tree_of_shapes.h"
#include <cfloat>
#include "TreeOfShapes.h"

enum Abstraction_Mode { COLOR_SKETCH=0, FILTER_COLOR=1, SHOW_TREE=2, SYNTEXTURE_COLOR=3, SYNTEXTURE_COLOR_WA=4,
                        SYNTEXTURE_COLOR_DICT=5, SYNTEXTURE_COLOR_DICT2=6, SYNTEXTURE_COLOR_V2=7,
                        SYNTEXTURE_COLOR_V3=8, SYNTEXTURE_COLOR_V4=9,
                        SYNTEXTURE_COLOR_DICT2_OUTPUT=10, SYNTEXTURE_COLOR_DICT3=11, SYNTEXTURE_COLOR_MULT=12,
                        SYNTEXTURE_COLOR_TT=13};

class AbstractionProcess
{
public:
    AbstractionProcess(){ _image_loaded = false; }
    AbstractionProcess( std::string fileNameIn );
    AbstractionProcess( const QImage &imageIn );
    ~AbstractionProcess();

    QImage run (int process, char* dictionnary_name);
    QImage render(TOSParameters tosParameters, bool &tree_recomputed, DictionaryParameters dictionaryParameters = getDefaultDictionaryParameters(), TreeOfShapes * dictionnary=NULL);
    std::vector<QImage> render_shape_by_shape(TOSParameters tosParameters, DictionaryParameters dictionaryParameters, TreeOfShapes * dictionnary);
    void addDictionnary( TreeOfShapes * dictionary );

    void save_shapes( QString folder_name, bool average_color ){ _treeOfShapes->save_shapes(folder_name, average_color); }
    void save_shapes( QString folder_name, TOSParameters tosParameters, DictionaryParameters dictionaryParameters = getDefaultDictionaryParameters(), TreeOfShapes * dictionnary=NULL );
    TOSParameters getParameters(){return _tosParameters;}
    int getMaxArea(){return _treeOfShapes->getMaxArea();}
    void getTreeInfo(std::vector<QPoint> &positions, std::vector<QColor> &colors,
                     std::vector< std::vector< std::pair<int, int> > > &pixels, std::vector<int> & heights, std::vector<std::pair<int, int> > &edges){
        _treeOfShapes->getTreeInfo(positions, colors, pixels, heights, edges);
    }
protected:
    bool _tree_computed;
    Cfimage _imgin;
    bool _image_loaded;

    TreeOfShapes *_treeOfShapes;
    TreeOfShapes *_dictionnary;

    Shapes _pTree;

    Fimage   _NormOfDu;

    TOSParameters _tosParameters;

    void init(Cfimage inputImg, Shapes &pTree);

    Cfimage cfimageread(const char* name);

    Cfimage cfimages_from_qimage( const QImage &input_image  );

    QImage qimages_from_cfimage( Cfimage input_img  );

    QImage qimages_from_ccimage( Ccimage input_img  );

    /* This removes the shapes from the tree associated to pFloatImageInput
    that are too small (threshold *pMinArea). As a consequence all the remaining
    shapes of pFloatImageOutput are of area larger or equal than *pMinArea */

    void mw_fgrain_side(int *pMinArea, Fimage pFloatImageInput, Fimage pFloatImageOutput, int sideflag);


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

    void synshapeEllipse(Shape pShape,
                         Ccimage imgsyn,
                         float *alpha,
                         int *relief,
                         float *reliefOrentation, float *reliefHeight);
    void synshapeRect(Shape pShape,
                      Ccimage imgsyn,
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
                           float *alpha);

    void top2bottom_index_tree(Fsignal t2b_index);
    void random_tree_order(Fsignal t2b_index);

    void random_leaf(Fsignal leaf,
                     Fsignal t2b_index,
                     int *k_ind);

    void random_shift_shape(float *shift);
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
                      int *mpixel);
    void filter_shapes( Fimage sketch,
                        Cfimage out,
                        char *local,
                        float *eps);
    void get_shapes_truearea(Shape s, Shape root,
                             int *truearea);

    Shape selectShapeDict(Shapes pTreeDict,
                          Shape pShape,
                          float *paDict,
                          int *randS);
    int random_number(int *M);
};

#endif // ABSTRACTIONGRAINPROCESS_H
