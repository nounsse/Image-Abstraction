#ifndef OPTIONWIDGET_H
#define OPTIONWIDGET_H

#include <QDockWidget>
#include <QSlider>
#include <QSpinBox>
#include <QRadioButton>
#include <QCheckBox>
#include <QGroupBox>
#include <QListWidget>
#include <QPushButton>
#include <QGridLayout>
#include <QDialog>
#include <QComboBox>
#include "tree_of_shapes.h"

typedef struct SegParameters {
    float sigma; // sigma: to smooth the image.
    int c; /// c: constant for treshold function.
    int min_size; // min_size: minimum component size (enforced by post-processing stage).

} SegParameters;

class OptionWidget : public QDockWidget
{
    Q_OBJECT
public:
    OptionWidget(QWidget * parent );
    ~OptionWidget(){}


    int getProcessingMode(){ return _processingMode; }
    TOSParameters getTOSParameters(){ return _TOSParameters; }
    DictionaryParameters getDictionaryParameters(){ return _dictionaryParameters; }
    SegParameters getSegParameters(){ return _segParameters; }
    void addDictionary(const QString &dict_name);
    int getSelectedDictionary(){ return _dictionariesWidget->currentRow();}
    void updateFromParameters(TOSParameters tosParam);
    void setNbShapes(int shapeNb);
protected:

    QDialog * _segmentationDialog;

    QListWidget * _dictionariesWidget;

    DictionaryParameters _dictionaryParameters;
    SegParameters _segParameters;
    QDoubleSpinBox * _sigmaRSpinBox;
    QSlider * _sigmaRSlider;

    QDoubleSpinBox * _transparencySpinBox;
    QSlider * _transparencySlider;

    QSpinBox * _scaleRatioSpinBox;
    QSlider * _scaleRatioSlider;

    QDoubleSpinBox * _thresholdSpinBox;
    QSlider * _thresholdSlider;

    QDoubleSpinBox * _compactnessSpinBox;
    QSlider * _compactnessSlider;

    QSpinBox *_minAreaSpinBox;
    QSlider *_minAreaSlider;

    QSpinBox *_largeShapeSpinBox;
    QSlider *_largeShapeSlider;

    QGroupBox * _reliefGroupBox;
    QGroupBox * _blurGroupBox;

    QRadioButton * _colorSketchRadioButton;
    QSpinBox *_epsSpinBox;
    QSlider *_epsSlider;

    QDoubleSpinBox * _segSigmaSpinBox;
    QSlider * _segSigmaSlider;

    QSpinBox * _segConstantSpinBox;
    QSlider * _segConstantSlider;

    QSpinBox * _segMinSizeSpinBox;
    QSlider * _segMinSizeSlider;

    int _processingMode;

    int _effect_intensity;
    QSlider * _effectIntensitySlider;

    TOSParameters _TOSParameters;

    QDoubleSpinBox * _shiftSpinBox;
    QSlider * _shiftSlider;

    QDoubleSpinBox * _thetaSpinBox;
    QSlider * _thetaSlider;

    QDoubleSpinBox * _compactnessDictSpinBox;
    QSlider * _compactnessDictSlider;

    QDoubleSpinBox * _kerStdSpinBox;
    QSlider * _kerStdSlider;

    QSpinBox * _kerSizeSpinBox;
    QSlider * _kerSizeSlider;

    QSpinBox * _medianSizeSpinBox;
    QSlider * _medianSizeSlider;


    QSpinBox * _reliefHeightSpinBox;
    QSlider * _reliefHeightSlider;


    QSpinBox * _reliefOrientationSpinBox;
    QSlider * _reliefOrientationSlider;

    QGroupBox * _dictionaryParameterGroupBox;

    QPushButton * _segmentPushButton;
    QPushButton * _updateInputPushButton;

    QComboBox * _shakingTypeComboBox;
    QComboBox * _renderModelComboBox;
    QComboBox * _renderOrderComboBox;
    QComboBox * _modeComboBox;

    int _effect_intensity_nb;
    void setInLayout( QGridLayout*layout, char* text, QWidget *first, QWidget *second, int row );
private slots:
    void setTransparencySlider();
    void setTransparencySpinBox();
    void setScaleRatioSlider();
    void setScaleRatioSpinBox();
    void setThresholdSlider();
    void setThresholdSpinBox();
    void setSegSigmaSlider();
    void setSegSigmaSpinBox();
    void setSegMinSizeSlider();
    void setSegMinSizeSpinBox();
    void setSegConstantSlider();
    void setSegConstantSpinBox();
    void setCompactnessSlider();
    void setCompactnessSpinBox();
    void setShiftSlider();
    void setShiftSpinBox();
    void setThetaSlider();
    void setThetaSpinBox();
    void setMinAreaSlider();
    void setMinAreaSpinBox();
    void setLargeShapeSlider();
    void setLargeShapeSpinBox();
    void setMode( int mode );
    void setRenderModel( int renderModel );
    void setRenderOrder( int renderOrder );
    void setRelief( bool addRelief );
    void setBlur( bool addBlur );
    void setKerStdSlider();
    void setKerStdSpinBox();
    void setKerSizeSlider();
    void setKerSizeSpinBox();
    void setMedianSizeSlider();
    void setMedianSizeSpinBox();
    void setReliefOrientationSlider();
    void setReliefOrientationSpinBox();
    void setReliefHeightSlider();
    void setReliefHeightSpinBox();
    void setColorSketch( bool color_sketch );
    void setEpsSlider();
    void setEpsSpinBox();
    void setShakingType( int smodel );
    void setShapeScaling( int shapeScaling );
    void setSelectionModel( int selectionModel );
    void setColorModel( int colorModel );
    void setCompactnessDictSlider();
    void setCompactnessDictSpinBox();
    void setEffectIntensity( int effect_intensity );
};

#endif // OPTIONWIDGET_H
