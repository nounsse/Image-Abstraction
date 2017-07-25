#include "OptionWidget.h"

#include <QLabel>
#include <QToolBox>
#include <QFileInfo>
#include <iostream>


void OptionWidget::setInLayout( QGridLayout*layout, char* text, QWidget *first, QWidget *second, int row ){

    layout->addWidget(new QLabel(text), row,0);
    layout->addWidget(first, row,1);
    layout->addWidget(second, row,2);

}

OptionWidget::OptionWidget(QWidget * parent ):QDockWidget(parent)
{

    /**************************************************/
    /*************   SET DEFAULT INPUT OPTIONS   **************/
    /**************************************************/

    _TOSParameters = getAbstractionTOSParameters();
    _dictionaryParameters = getDefaultDictionaryParameters();

    _segParameters.c = 500;
    _segParameters.min_size = 50;
    _segParameters.sigma = 0.5;

    _effect_intensity = 5;
    _effect_intensity_nb = 10;
    QWidget * contents = new QWidget();

    QVBoxLayout * contentLayout = new QVBoxLayout(contents);

    _segmentationDialog = new QDialog();
    _segmentationDialog->setModal(true);

    QVBoxLayout * segLayout = new QVBoxLayout(_segmentationDialog);

    QGroupBox * segmentationParameterGroupBox = new QGroupBox("Shape segmentation" );
    segLayout->addWidget(segmentationParameterGroupBox );


    QVBoxLayout * segmentationVLayout = new QVBoxLayout( segmentationParameterGroupBox );
    QGridLayout * segmentationGroupBoxGLayout = new QGridLayout( );
    segmentationVLayout->addLayout(segmentationGroupBoxGLayout);

    _segSigmaSpinBox = new QDoubleSpinBox;
    _segSigmaSpinBox->setRange(0.0, 1.0);
    _segSigmaSpinBox->setSingleStep(0.1);
    _segSigmaSpinBox->setValue(_segParameters.sigma);

    _segSigmaSlider = new QSlider( Qt::Horizontal);
    _segSigmaSlider->setRange(0, 10);
    _segSigmaSlider->setSingleStep(1);
    _segSigmaSlider->setValue( int(_segParameters.sigma/0.1) );

    int count = 0;
    setInLayout(segmentationGroupBoxGLayout, "Sigma (smoothing)", _segSigmaSpinBox, _segSigmaSlider, count++ );

    connect(_segSigmaSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setSegSigmaSpinBox()));
    connect(_segSigmaSlider, SIGNAL(valueChanged(int)), this, SLOT(setSegSigmaSlider()));

    _segConstantSpinBox = new QSpinBox;
    _segConstantSpinBox->setRange(1, 2000);
    _segConstantSpinBox->setSingleStep(1);
    _segConstantSpinBox->setValue(_segParameters.c);

    _segConstantSlider = new QSlider( Qt::Horizontal);
    _segConstantSlider->setRange(1, 2000);
    _segConstantSlider->setSingleStep(1);
    _segConstantSlider->setValue(_segParameters.c);

    setInLayout(segmentationGroupBoxGLayout, "Threshold Function", _segConstantSpinBox, _segConstantSlider, count++ );

    connect(_segConstantSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setSegConstantSpinBox()));
    connect(_segConstantSlider, SIGNAL(valueChanged(int)), this, SLOT(setSegConstantSlider()));

    _segMinSizeSpinBox = new QSpinBox;
    _segMinSizeSpinBox->setRange(1, 1500);
    _segMinSizeSpinBox->setSingleStep(1);
    _segMinSizeSpinBox->setValue(_segParameters.min_size);

    _segMinSizeSlider = new QSlider( Qt::Horizontal);
    _segMinSizeSlider->setRange(0, 1500);
    _segMinSizeSlider->setSingleStep(1);
    _segMinSizeSlider->setValue( _segParameters.min_size );

    setInLayout(segmentationGroupBoxGLayout, "Min size", _segMinSizeSpinBox, _segMinSizeSlider, count++ );

    connect(_segMinSizeSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setSegMinSizeSpinBox()));
    connect(_segMinSizeSlider, SIGNAL(valueChanged(int)), this, SLOT(setSegMinSizeSlider()));

    _segmentPushButton = new QPushButton("Segment");
    connect(_segmentPushButton, SIGNAL(clicked()), parent, SLOT(segment()));
    connect(_segmentPushButton, SIGNAL(clicked()), _segmentationDialog, SLOT(close()));

    _updateInputPushButton = new QPushButton("Use as input");
    connect(_updateInputPushButton, SIGNAL(clicked()), parent, SLOT(updateFromSegmentation()));
    connect(_updateInputPushButton, SIGNAL(clicked()), _segmentationDialog, SLOT(close()));

    segmentationVLayout->addWidget(_segmentPushButton);


    QPushButton * segmentationPushButton = new QPushButton( "Segment Input image" );
    contentLayout->addWidget(segmentationPushButton);
    connect(segmentationPushButton, SIGNAL(clicked()), _segmentationDialog, SLOT(open()));
    contentLayout->addWidget(_updateInputPushButton);

    QGroupBox * outputImageParameterGroupBox = new QGroupBox("Image Abstraction");
    QVBoxLayout * outputVBoxLayout = new QVBoxLayout(outputImageParameterGroupBox);
    contentLayout->addWidget(outputImageParameterGroupBox);

    _modeComboBox = new QComboBox ();

    // modeComboBox->addItem ("COLOR_SKETCH");
    _modeComboBox->addItem ("Shape abstraction");
    _modeComboBox->addItem ("Watercolor effect");
    _modeComboBox->addItem ("Brush/Painting");
    _modeComboBox->addItem ("Shape smoothing");
    _modeComboBox->addItem ("Style transfer");


    _processingMode = 0;

    connect (_modeComboBox, SIGNAL (activated (int)),
             this, SLOT (setMode (int)));
    outputVBoxLayout->addWidget (_modeComboBox);
    outputVBoxLayout->addWidget (new QLabel("Effect intensity"));
    _effectIntensitySlider = new QSlider( Qt::Horizontal);
    _effectIntensitySlider->setRange(0, _effect_intensity_nb);
    _effectIntensitySlider->setSingleStep(1);
    _effectIntensitySlider->setValue(_effect_intensity);
    outputVBoxLayout->addWidget(_effectIntensitySlider);

    connect(_effectIntensitySlider, SIGNAL(valueChanged(int)), this, SLOT(setEffectIntensity(int)));


    _renderModelComboBox = new QComboBox ();

    _renderModelComboBox->addItem ("Original Shapes");
    _renderModelComboBox->addItem ("Ellipses");
    _renderModelComboBox->addItem ("Rectangles");
    _renderModelComboBox->addItem ("Circle");
    _renderModelComboBox->addItem ("From Dictionary");
    _renderModelComboBox->addItem ("Random");

    _renderModelComboBox->setCurrentIndex(_TOSParameters.model);

    connect (_renderModelComboBox, SIGNAL (activated (int)),
             this, SLOT (setRenderModel (int)));

    outputVBoxLayout->addWidget(new QLabel("Replace shapes with"));
    outputVBoxLayout->addWidget(_renderModelComboBox);


    _dictionaryParameterGroupBox = new QGroupBox("Dictionary options");
    outputVBoxLayout->addWidget(_dictionaryParameterGroupBox);


    _dictionariesWidget = new QListWidget(this);
    _dictionariesWidget->setVisible(false);

    outputVBoxLayout->addWidget(_dictionariesWidget);
    QGridLayout * dictionaryParamGroupBoxGLayout = new QGridLayout( _dictionaryParameterGroupBox );

    QComboBox * selectionModelComboBox = new QComboBox ();

    selectionModelComboBox->addItem ("Randomly");
    selectionModelComboBox->addItem ("Elong/Compactness/Scale");
    selectionModelComboBox->addItem ("Elong/Compactness/Scale/Color");

    selectionModelComboBox->setCurrentIndex(_dictionaryParameters.randS);

    connect (selectionModelComboBox, SIGNAL (activated (int)),
             this, SLOT (setSelectionModel(int)));

    count = 0;
    dictionaryParamGroupBoxGLayout->addWidget(new QLabel("Selection model"), count,0);
    dictionaryParamGroupBoxGLayout->addWidget(selectionModelComboBox, count++,1);


    QComboBox * colorModelComboBox = new QComboBox ();

    colorModelComboBox->addItem ("from input");
    colorModelComboBox->addItem ("from dictionary");
    colorModelComboBox->addItem ("average dictionary");

    colorModelComboBox->setCurrentIndex(_dictionaryParameters.mcolor);

    connect (colorModelComboBox, SIGNAL (activated (int)),
             this, SLOT (setColorModel(int)));
    dictionaryParamGroupBoxGLayout->addWidget(new QLabel("Use color"), count,0);
    dictionaryParamGroupBoxGLayout->addWidget(colorModelComboBox, count++,1);

    QComboBox * shapeScalingComboBox = new QComboBox ();

    shapeScalingComboBox->addItem ("Equal");
    shapeScalingComboBox->addItem ("Different");

    shapeScalingComboBox->setCurrentIndex(_dictionaryParameters.equal);

    connect (shapeScalingComboBox, SIGNAL (activated (int)),
             this, SLOT (setShapeScaling(int)));
    dictionaryParamGroupBoxGLayout->addWidget(new QLabel("Shape x/y scaling"), count,0);
    dictionaryParamGroupBoxGLayout->addWidget(shapeScalingComboBox, count++,1);

    _compactnessDictSpinBox = new QDoubleSpinBox;
    _compactnessDictSpinBox->setRange(0.0, 1.0);
    _compactnessDictSpinBox->setSingleStep(0.1);
    _compactnessDictSpinBox->setValue(_dictionaryParameters.kappaDict);

    _compactnessDictSlider = new QSlider( Qt::Horizontal);
    _compactnessDictSlider->setRange(0, 10);
    _compactnessDictSlider->setSingleStep(1);
    _compactnessDictSlider->setValue( int(_dictionaryParameters.kappaDict/0.1) );

    dictionaryParamGroupBoxGLayout->addWidget(new QLabel("Compactness"), count,0);
    dictionaryParamGroupBoxGLayout->addWidget(_compactnessDictSpinBox, count++,1);

    connect(_compactnessDictSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setCompactnessSpinBox()));
    connect(_compactnessDictSlider, SIGNAL(valueChanged(int)), this, SLOT(setCompactnessSlider()));

    /*
    modeComboBox->addItem ("COLOR_SKETCH");
    modeComboBox->addItem ("FILTER_COLOR");

    modeComboBox->addItem ("SHOW_TREE");

    modeComboBox->addItem ("SYNTEXTURE_COLOR");
    modeComboBox->addItem (" SYNTEXTURE_COLOR_WA");

    modeComboBox->addItem ("SYNTEXTURE_COLOR_DICT");
    modeComboBox->addItem ("SYNTEXTURE_COLOR_DICT2");

    modeComboBox->addItem ("SYNTEXTURE_COLOR_V2");
    modeComboBox->addItem ("SYNTEXTURE_COLOR_V3");
    modeComboBox->addItem ("SYNTEXTURE_COLOR_V4");

    modeComboBox->addItem ("SYNTEXTURE_COLOR_DICT2_OUTPUT");
    modeComboBox->addItem ("SYNTEXTURE_COLOR_DICT3");
    modeComboBox->addItem ("SYNTEXTURE_COLOR_MULT");

    modeComboBox->addItem ("SYNTEXTURE_COLOR_TT");

    _processingMode = 0;

    connect (modeComboBox, SIGNAL (activated (int)),
             this, SLOT (setMode (int)));
    contentLayout->addWidget (modeComboBox);
*/

    // outputImageParamGroupBoxGLayout->addWidget(new QLabel("Mode"), 0,1);
    // outputImageParamGroupBoxGLayout->addWidget(modeComboBox, 0,2);

    QGroupBox * parametersGroupBox = new QGroupBox("Image Abstraction Parameters");
    contentLayout->addWidget(parametersGroupBox);
    QVBoxLayout * paramVBoxLayout = new QVBoxLayout(parametersGroupBox);

    QToolBox * abstractionToolBox = new QToolBox();
    paramVBoxLayout->addWidget(abstractionToolBox);

    QGroupBox * filteringParameterGroupBox = new QGroupBox();
    abstractionToolBox->addItem(filteringParameterGroupBox, "Shape selection" );

    QGridLayout * filteringGroupBoxGLayout = new QGridLayout( filteringParameterGroupBox );

    _minAreaSpinBox = new QSpinBox;
    _minAreaSpinBox->setRange(1, 1500);
    _minAreaSpinBox->setSingleStep(1);
    _minAreaSpinBox->setValue(_TOSParameters.mpixel);

    _minAreaSlider = new QSlider( Qt::Horizontal);
    _minAreaSlider->setRange(1, 1500);
    _minAreaSlider->setSingleStep(1);
    _minAreaSlider->setValue(_TOSParameters.mpixel);

    count = 0;
    setInLayout(filteringGroupBoxGLayout, "Min area (pixel)", _minAreaSpinBox, _minAreaSlider, count++ );

    connect(_minAreaSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setMinAreaSpinBox()));
    connect(_minAreaSlider, SIGNAL(valueChanged(int)), this, SLOT(setMinAreaSlider()));

    _largeShapeSpinBox = new QSpinBox;
    _largeShapeSpinBox->setRange(1000, 1500);
    _largeShapeSpinBox->setSingleStep(1);
    _largeShapeSpinBox->setValue(_TOSParameters.maxarea);

    _largeShapeSlider = new QSlider( Qt::Horizontal);
    _largeShapeSlider->setRange(1000, 1500);
    _largeShapeSlider->setSingleStep(1);
    _largeShapeSlider->setValue(_TOSParameters.maxarea);

    setInLayout(filteringGroupBoxGLayout, "Max area", _largeShapeSpinBox, _largeShapeSlider, count++ );

    connect(_largeShapeSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setLargeShapeSpinBox()));
    connect(_largeShapeSlider, SIGNAL(valueChanged(int)), this, SLOT(setLargeShapeSlider()));

    _scaleRatioSpinBox = new QSpinBox;
    _scaleRatioSpinBox->setRange(0, 10);
    _scaleRatioSpinBox->setSingleStep(1);
    _scaleRatioSpinBox->setValue(_TOSParameters.ns);

    _scaleRatioSlider = new QSlider( Qt::Horizontal);
    _scaleRatioSlider->setRange(0, 10);
    _scaleRatioSlider->setSingleStep(1);
    _scaleRatioSlider->setValue( _TOSParameters.ns );

    setInLayout(filteringGroupBoxGLayout, "Scale Ratio Order", _scaleRatioSpinBox, _scaleRatioSlider, count++ );

    connect(_scaleRatioSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setScaleRatioSpinBox()));
    connect(_scaleRatioSlider, SIGNAL(valueChanged(int)), this, SLOT(setScaleRatioSlider()));


    _thresholdSpinBox = new QDoubleSpinBox;
    _thresholdSpinBox->setRange(0.0, 1.0);
    _thresholdSpinBox->setSingleStep(0.1);
    _thresholdSpinBox->setValue(_TOSParameters.threshold);

    _thresholdSlider = new QSlider( Qt::Horizontal);
    _thresholdSlider->setRange(0, 10);
    _thresholdSlider->setSingleStep(1);
    _thresholdSlider->setValue( int(_TOSParameters.threshold/0.1) );

    setInLayout(filteringGroupBoxGLayout, "Threshold", _thresholdSpinBox, _thresholdSlider, count++ );

    connect(_thresholdSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setThresholdSpinBox()));
    connect(_thresholdSlider, SIGNAL(valueChanged(int)), this, SLOT(setThresholdSlider()));

    _compactnessSpinBox = new QDoubleSpinBox;
    _compactnessSpinBox->setRange(0.0, 1.0);
    _compactnessSpinBox->setSingleStep(0.1);
    _compactnessSpinBox->setValue(_TOSParameters.kappa);

    _compactnessSlider = new QSlider( Qt::Horizontal);
    _compactnessSlider->setRange(0, 10);
    _compactnessSlider->setSingleStep(1);
    _compactnessSlider->setValue( int(_TOSParameters.kappa/0.1) );

    setInLayout(filteringGroupBoxGLayout, "Compactness", _compactnessSpinBox, _compactnessSlider, count++ );

    connect(_compactnessSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setCompactnessSpinBox()));
    connect(_compactnessSlider, SIGNAL(valueChanged(int)), this, SLOT(setCompactnessSlider()));

    _colorSketchRadioButton = new QRadioButton("Keep meaningful boundaries");
    if( _TOSParameters.color_sketch == 1)
        _colorSketchRadioButton->setChecked(true);
    else
        _colorSketchRadioButton->setChecked(false);

    connect(_colorSketchRadioButton, SIGNAL(toggled(bool)), this, SLOT(setColorSketch(bool)));

    filteringGroupBoxGLayout->addWidget(_colorSketchRadioButton, count++,0);

    _epsSpinBox = new QSpinBox;
    _epsSpinBox->setRange(-10, 10);
    _epsSpinBox->setSingleStep(1);
    _epsSpinBox->setValue(_TOSParameters.eps);

    _epsSlider = new QSlider( Qt::Horizontal);
    _epsSlider->setRange(-10, 10);
    _epsSlider->setSingleStep(1);
    _epsSlider->setValue( int(_TOSParameters.eps));

    setInLayout(filteringGroupBoxGLayout, "Eps", _epsSpinBox, _epsSlider, count++ );

    connect(_epsSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setEpsSpinBox()));
    connect(_epsSlider, SIGNAL(valueChanged(int)), this, SLOT(setEpsSlider()));

    QGroupBox * shakingParameterGroupBox = new QGroupBox();
    abstractionToolBox->addItem(shakingParameterGroupBox, "Shape displacement" );

    QGridLayout * shakingParamGroupBoxGLayout = new QGridLayout( shakingParameterGroupBox );

    _shiftSpinBox = new QDoubleSpinBox;
    _shiftSpinBox->setRange(0.0, 50.0);
    _shiftSpinBox->setSingleStep(0.5);
    _shiftSpinBox->setValue(_TOSParameters.shift);

    _shiftSlider = new QSlider( Qt::Horizontal);
    _shiftSlider->setRange(0, 500);
    _shiftSlider->setSingleStep(1);
    _shiftSlider->setValue( int(_TOSParameters.shift/0.1) );

    count = 0;
    setInLayout(shakingParamGroupBoxGLayout, "Shift", _shiftSpinBox, _shiftSlider, count++ );

    connect(_shiftSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setShiftSpinBox()));
    connect(_shiftSlider, SIGNAL(valueChanged(int)), this, SLOT(setShiftSlider()));

    _thetaSpinBox = new QDoubleSpinBox;
    _thetaSpinBox->setRange(0.0, 1.);
    _thetaSpinBox->setSingleStep(0.1);
    _thetaSpinBox->setValue(_TOSParameters.theta);

    _thetaSlider = new QSlider( Qt::Horizontal);
    _thetaSlider->setRange(0, 10);
    _thetaSlider->setSingleStep(1);
    _thetaSlider->setValue( int(_TOSParameters.theta/0.1) );

    setInLayout(shakingParamGroupBoxGLayout, "Rotation", _thetaSpinBox, _thetaSlider, count++ );

    connect(_thetaSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setThetaSpinBox()));
    connect(_thetaSlider, SIGNAL(valueChanged(int)), this, SLOT(setThetaSlider()));

    _shakingTypeComboBox = new QComboBox ();
    _shakingTypeComboBox->addItem ("Uniform");
    _shakingTypeComboBox->addItem ("Dominant");

    _shakingTypeComboBox->setCurrentIndex(_TOSParameters.smodel);

    connect (_shakingTypeComboBox, SIGNAL (activated (int)),
             this, SLOT (setShakingType (int)));

    setInLayout(shakingParamGroupBoxGLayout, "Shaking type", _shakingTypeComboBox, _shakingTypeComboBox, count++ );

    QGroupBox * renderingParameterGroupBox = new QGroupBox();
    abstractionToolBox->addItem(renderingParameterGroupBox, "Rendering");

    QVBoxLayout * renderingLayout = new QVBoxLayout( renderingParameterGroupBox );
    _renderOrderComboBox = new QComboBox ();

    _renderOrderComboBox->addItem ("Top -> down");
    _renderOrderComboBox->addItem ("Large -> small");
    _renderOrderComboBox->addItem ("Random");

    _renderOrderComboBox->setCurrentIndex(_TOSParameters.order);

    connect (_renderOrderComboBox, SIGNAL (activated (int)),
             this, SLOT (setRenderOrder (int)));
    renderingLayout->addWidget(new QLabel("Rendering order"));
    renderingLayout->addWidget(_renderOrderComboBox);

    QGridLayout * trLayout = new QGridLayout();
    renderingLayout->addLayout(trLayout);
    _transparencySpinBox = new QDoubleSpinBox;
    _transparencySpinBox->setRange(0.0, 1.0);
    _transparencySpinBox->setSingleStep(0.1);
    _transparencySpinBox->setValue(_TOSParameters.alpha);

    _transparencySlider = new QSlider( Qt::Horizontal);
    _transparencySlider->setRange(0, 10);
    _transparencySlider->setSingleStep(1);
    _transparencySlider->setValue( int(_TOSParameters.alpha/0.1) );

    setInLayout(trLayout, "Transparency", _transparencySpinBox, _transparencySlider, 0 );

    connect(_transparencySpinBox, SIGNAL(valueChanged(double)), this, SLOT(setTransparencySpinBox()));
    connect(_transparencySlider, SIGNAL(valueChanged(int)), this, SLOT(setTransparencySlider()));

    _blurGroupBox = new QGroupBox("Add Blur");
    _blurGroupBox->setCheckable(true);

    QGridLayout * blurLayout = new QGridLayout(_blurGroupBox);

    if( _TOSParameters.blur == 1)
        _blurGroupBox->setChecked(true);
    else
        _blurGroupBox->setChecked(false);

    connect(_blurGroupBox, SIGNAL(toggled(bool)), this, SLOT(setBlur(bool)));

    renderingLayout->addWidget(_blurGroupBox);

    _kerStdSpinBox = new QDoubleSpinBox;
    _kerStdSpinBox->setRange(0.0, 1.0);
    _kerStdSpinBox->setSingleStep(0.1);
    _kerStdSpinBox->setValue(_TOSParameters.kerStd);

    _kerStdSlider = new QSlider( Qt::Horizontal);
    _kerStdSlider->setRange(0, 10);
    _kerStdSlider->setSingleStep(1);
    _kerStdSlider->setValue( int(_TOSParameters.kerStd/0.1) );

    int blur_count = 0;
    setInLayout(blurLayout, "Gaussian kernel", _kerStdSpinBox, _kerStdSlider, blur_count++ );

    connect(_kerStdSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setKerStdSpinBox()));
    connect(_kerStdSlider, SIGNAL(valueChanged(int)), this, SLOT(setKerStdSlider()));

    _kerSizeSpinBox = new QSpinBox;
    _kerSizeSpinBox->setRange(0, 30);
    _kerSizeSpinBox->setSingleStep(1);
    _kerSizeSpinBox->setValue(_TOSParameters.kerSize);

    _kerSizeSlider = new QSlider( Qt::Horizontal);
    _kerSizeSlider->setRange(0, 30);
    _kerSizeSlider->setSingleStep(1);
    _kerSizeSlider->setValue( int(_TOSParameters.kerSize));

    // setInLayout(blurLayout, "Gaussian Blur size", _kerSizeSpinBox, _kerSizeSlider, blur_count++ );

    connect(_kerSizeSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setKerSizeSpinBox()));
    connect(_kerSizeSlider, SIGNAL(valueChanged(int)), this, SLOT(setKerSizeSlider()));

    _medianSizeSpinBox = new QSpinBox;
    _medianSizeSpinBox->setRange(0, 30);
    _medianSizeSpinBox->setSingleStep(1);
    _medianSizeSpinBox->setValue(_TOSParameters.median);

    _medianSizeSlider = new QSlider( Qt::Horizontal);
    _medianSizeSlider->setRange(0, 30);
    _medianSizeSlider->setSingleStep(1);
    _medianSizeSlider->setValue( int(_TOSParameters.median));

    setInLayout(blurLayout, "Median filter Size", _medianSizeSpinBox, _medianSizeSlider, blur_count++ );

    connect(_medianSizeSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setMedianSizeSpinBox()));
    connect(_medianSizeSlider, SIGNAL(valueChanged(int)), this, SLOT(setMedianSizeSlider()));

    _reliefGroupBox = new QGroupBox("Add Relief");
    _reliefGroupBox->setCheckable(true);

    QGridLayout * reliefLayout = new QGridLayout(_reliefGroupBox);

    if( _TOSParameters.relief == 1)
        _reliefGroupBox->setChecked(true);
    else
        _reliefGroupBox->setChecked(false);

    connect(_reliefGroupBox, SIGNAL(toggled(bool)), this, SLOT(setRelief(bool)));

    renderingLayout->addWidget(_reliefGroupBox);

    _reliefOrientationSpinBox = new QSpinBox;
    _reliefOrientationSpinBox->setRange(0, 365);
    _reliefOrientationSpinBox->setSingleStep(1);
    _reliefOrientationSpinBox->setValue(_TOSParameters.reliefOrientation);

    _reliefOrientationSlider = new QSlider( Qt::Horizontal);
    _reliefOrientationSlider->setRange(0, 365);
    _reliefOrientationSlider->setSingleStep(1);
    _reliefOrientationSlider->setValue( int(_TOSParameters.reliefOrientation));

    int relief_count = 0;
    setInLayout(reliefLayout, "Relief Orientation", _reliefOrientationSpinBox, _reliefOrientationSlider, relief_count++ );

    connect(_reliefOrientationSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setReliefOrientationSpinBox()));
    connect(_reliefOrientationSlider, SIGNAL(valueChanged(int)), this, SLOT(setReliefOrientationSlider()));

    _reliefHeightSpinBox = new QSpinBox;
    _reliefHeightSpinBox->setRange(1, 30);
    _reliefHeightSpinBox->setSingleStep(1);
    _reliefHeightSpinBox->setValue(_TOSParameters.reliefHeight);

    _reliefHeightSlider = new QSlider( Qt::Horizontal);
    _reliefHeightSlider->setRange(1, 30);
    _reliefHeightSlider->setSingleStep(1);
    _reliefHeightSlider->setValue( int(_TOSParameters.reliefHeight));

    setInLayout(reliefLayout, "Relief Height", _reliefHeightSpinBox, _reliefHeightSlider, relief_count++ );

    connect(_reliefHeightSpinBox, SIGNAL(valueChanged(int)), this, SLOT(setReliefHeightSpinBox()));
    connect(_reliefHeightSlider, SIGNAL(valueChanged(int)), this, SLOT(setReliefHeightSlider()));


    _dictionaryParameterGroupBox->setVisible(false);
    QPushButton * runPushButton = new QPushButton("Run");

    connect(runPushButton, SIGNAL(clicked()), parent, SLOT(run()));
    contentLayout->addWidget( runPushButton );



    contentLayout->addStretch(0);
    this->setWidget(contents);
    abstractionToolBox->setCurrentIndex(2);
}

void OptionWidget::updateFromParameters(TOSParameters tosParam){

    _TOSParameters = tosParam;
    _minAreaSpinBox->setValue(_TOSParameters.mpixel);
    _minAreaSlider->setValue(_TOSParameters.mpixel);
    _scaleRatioSpinBox->setValue(_TOSParameters.ns);
    _scaleRatioSlider->setValue( _TOSParameters.ns );
    _thresholdSpinBox->setValue(_TOSParameters.threshold);
    _thresholdSlider->setValue( int(_TOSParameters.threshold/0.1) );
    _compactnessSpinBox->setValue(_TOSParameters.kappa);
    _compactnessSlider->setValue( int(_TOSParameters.kappa/0.1) );
    if( _TOSParameters.color_sketch == 1)
        _colorSketchRadioButton->setChecked(true);
    else
        _colorSketchRadioButton->setChecked(false);
    _epsSpinBox->setValue(_TOSParameters.eps);
    _epsSlider->setValue( int(_TOSParameters.eps));
    _shiftSpinBox->setValue(_TOSParameters.shift);
    _shiftSlider->setValue( int(_TOSParameters.shift/0.1) );
    _thetaSpinBox->setValue(_TOSParameters.theta);
    _thetaSlider->setValue( int(_TOSParameters.theta/0.1) );
    _shakingTypeComboBox->setCurrentIndex(_TOSParameters.smodel);
    _renderModelComboBox->setCurrentIndex(_TOSParameters.model);
    setRenderModel(_TOSParameters.model);
    _transparencySpinBox->setValue(_TOSParameters.alpha);
    _transparencySlider->setValue( int(_TOSParameters.alpha/0.1) );
    if( _TOSParameters.relief == 1)
        _reliefGroupBox->setChecked(true);
    else
        _reliefGroupBox->setChecked(false);
    if( _TOSParameters.blur == 1)
        _blurGroupBox->setChecked(true);
    else
        _blurGroupBox->setChecked(false);
    _kerStdSpinBox->setValue(_TOSParameters.kerStd);
    _kerStdSlider->setValue( int(_TOSParameters.kerStd/0.1) );
    _medianSizeSpinBox->setValue(_TOSParameters.kerSize);
    _medianSizeSlider->setValue( int(_TOSParameters.kerSize));
    _medianSizeSpinBox->setValue(_TOSParameters.median);
    _medianSizeSlider->setValue( int(_TOSParameters.median));
    _renderOrderComboBox->setCurrentIndex(_TOSParameters.order);

}

void OptionWidget::addDictionary(const QString &dict_name){

    QFileInfo info (dict_name);
    new QListWidgetItem(QIcon(dict_name), info.baseName(), _dictionariesWidget);
}


void OptionWidget::setScaleRatioSlider(){
    _scaleRatioSpinBox->setValue( _scaleRatioSlider->value() );
    _TOSParameters.ns = _scaleRatioSpinBox->value();
}

void OptionWidget::setScaleRatioSpinBox(){
    _scaleRatioSlider->setValue( _scaleRatioSpinBox->value() );
    _TOSParameters.ns = _scaleRatioSpinBox->value();
}

void OptionWidget::setThresholdSlider(){
    _thresholdSpinBox->setValue( _thresholdSlider->value()*0.1 );
    _TOSParameters.threshold = _thresholdSpinBox->value();
}

void OptionWidget::setThresholdSpinBox(){
    _thresholdSlider->setValue( int(_thresholdSpinBox->value()/0.1) );
    _TOSParameters.threshold = _thresholdSpinBox->value();
}

void OptionWidget::setCompactnessSlider(){
    _compactnessSpinBox->setValue( _compactnessSlider->value()*0.1 );
    _TOSParameters.kappa = _compactnessSpinBox->value();
}

void OptionWidget::setCompactnessSpinBox(){
    _compactnessSlider->setValue( int(_compactnessSpinBox->value()/0.1) );
    _TOSParameters.kappa = _compactnessSpinBox->value();
}


void OptionWidget::setCompactnessDictSlider(){
    _compactnessDictSpinBox->setValue( _compactnessDictSlider->value()*0.1 );
    _dictionaryParameters.kappaDict = _compactnessDictSpinBox->value();
}

void OptionWidget::setCompactnessDictSpinBox(){
    _compactnessDictSlider->setValue( int(_compactnessDictSpinBox->value()/0.1) );
    _dictionaryParameters.kappaDict = _compactnessDictSpinBox->value();
}


void OptionWidget::setThetaSlider(){
    _thetaSpinBox->setValue( _thetaSlider->value()*0.1 );
    _TOSParameters.theta = _thetaSpinBox->value();
}

void OptionWidget::setThetaSpinBox(){
    _thetaSlider->setValue( int(_thetaSpinBox->value()/0.1) );
    _TOSParameters.theta = _thetaSpinBox->value();
}

void OptionWidget::setShiftSlider(){
    _shiftSpinBox->setValue( _shiftSlider->value()*0.1 );
    _TOSParameters.shift = _shiftSpinBox->value();
}

void OptionWidget::setShiftSpinBox(){
    _shiftSlider->setValue( int(_shiftSpinBox->value()/0.1) );
    _TOSParameters.shift = _shiftSpinBox->value();
}

void OptionWidget::setKerStdSlider(){
    _kerStdSpinBox->setValue( _kerStdSlider->value()*0.1 );
    _TOSParameters.kerStd = _kerStdSpinBox->value();
}

void OptionWidget::setKerStdSpinBox(){
    _kerStdSlider->setValue( int(_kerStdSpinBox->value()/0.1) );
    _TOSParameters.kerStd = _kerStdSpinBox->value();
}

void OptionWidget::setKerSizeSlider(){
    _kerSizeSpinBox->setValue( _kerSizeSlider->value() );
    _TOSParameters.kerSize = _kerSizeSpinBox->value();
}

void OptionWidget::setKerSizeSpinBox(){
    _kerSizeSlider->setValue( _kerSizeSpinBox->value() );
    _TOSParameters.kerSize = _kerSizeSpinBox->value();
}

void OptionWidget::setMedianSizeSlider(){
    _medianSizeSpinBox->setValue( _medianSizeSlider->value() );
    _TOSParameters.median = _medianSizeSpinBox->value();
}

void OptionWidget::setMedianSizeSpinBox(){
    _medianSizeSlider->setValue( _medianSizeSpinBox->value() );
    _TOSParameters.median = _medianSizeSpinBox->value();
}

void OptionWidget::setReliefOrientationSlider(){
    _reliefOrientationSpinBox->setValue( _reliefOrientationSlider->value() );
    _TOSParameters.reliefOrientation = _reliefOrientationSpinBox->value();
}

void OptionWidget::setReliefOrientationSpinBox(){
    _reliefOrientationSlider->setValue( _reliefOrientationSpinBox->value() );
    _TOSParameters.reliefOrientation = _reliefOrientationSpinBox->value();
}

void OptionWidget::setReliefHeightSlider(){
    _reliefHeightSpinBox->setValue( _reliefHeightSlider->value() );
    _TOSParameters.reliefHeight = _reliefHeightSpinBox->value();
}

void OptionWidget::setReliefHeightSpinBox(){
    _reliefHeightSlider->setValue( _reliefHeightSpinBox->value() );
    _TOSParameters.reliefHeight = _reliefHeightSpinBox->value();
}

void OptionWidget::setTransparencySlider(){
    _transparencySpinBox->setValue( _transparencySlider->value()*0.1 );
    _TOSParameters.alpha = _transparencySpinBox->value();
}

void OptionWidget::setTransparencySpinBox(){
    _transparencySlider->setValue( int(_transparencySpinBox->value()/0.1) );
    _TOSParameters.alpha = _transparencySpinBox->value();
}

void OptionWidget::setEpsSlider(){
    _epsSpinBox->setValue( _epsSlider->value() );
    _TOSParameters.eps = _epsSpinBox->value();
}

void OptionWidget::setEpsSpinBox(){
    _epsSlider->setValue( _epsSpinBox->value() );
    _TOSParameters.eps = (float)_epsSpinBox->value();
}

void OptionWidget::setMinAreaSlider(){
    _minAreaSpinBox->setValue( _minAreaSlider->value() );
    _TOSParameters.mpixel = _minAreaSpinBox->value();
}

void OptionWidget::setMinAreaSpinBox(){
    _minAreaSlider->setValue( _minAreaSpinBox->value() );
    _TOSParameters.mpixel = _minAreaSpinBox->value();
}

void OptionWidget::setLargeShapeSlider(){
    _largeShapeSpinBox->setValue( _largeShapeSlider->value() );
    _TOSParameters.maxarea = _largeShapeSpinBox->value();
}

void OptionWidget::setLargeShapeSpinBox(){
    _largeShapeSlider->setValue( _largeShapeSpinBox->value() );
    _TOSParameters.maxarea = _largeShapeSpinBox->value();
}

void OptionWidget::setSegSigmaSlider(){
    _segSigmaSpinBox->setValue( _segSigmaSlider->value()*0.1 );
    _segParameters.sigma = _segSigmaSpinBox->value();
}

void OptionWidget::setSegSigmaSpinBox(){
    _segSigmaSlider->setValue( int(_segSigmaSpinBox->value()/0.1) );
    _segParameters.sigma = _segSigmaSpinBox->value();
}

void OptionWidget::setSegMinSizeSlider(){
    _segMinSizeSpinBox->setValue( _segMinSizeSpinBox->value() );
    _segParameters.min_size = _segMinSizeSpinBox->value();
}

void OptionWidget::setSegMinSizeSpinBox(){
    _segMinSizeSlider->setValue( _segMinSizeSpinBox->value() );
    _segParameters.min_size = _segMinSizeSpinBox->value();
}


void OptionWidget::setSegConstantSlider(){
    _segConstantSpinBox->setValue( _segConstantSlider->value() );
    _segParameters.c = _segConstantSpinBox->value();
}

void OptionWidget::setSegConstantSpinBox(){
    _segConstantSlider->setValue( _segConstantSpinBox->value() );
    _segParameters.c = _segConstantSpinBox->value();
}

void OptionWidget::setMode( int mode ){
    _effectIntensitySlider->setEnabled(true);
    if( mode==0 ){
        updateFromParameters ( getAbstractionTOSParameters() );
        std::cout << "Abstraction " << std::endl;
    } else if( mode==1 ){
        std::cout << "Watercolor " << std::endl;
        updateFromParameters ( getWaterColorTOSParameters() );
    } else if( mode==2 ){
        std::cout << "Shaking " << std::endl;
        updateFromParameters ( getShapeShakingTOSParameters() );
    } else if( mode==3 ){
        std::cout << "Shape smoothing " << std::endl;
        updateFromParameters ( getShapeSmoothingTOSParameters() );
    } else if( mode==4 ){
        std::cout << "Style transfer " << std::endl;
        updateFromParameters ( getStyleTransferTOSParameters() );
        _effectIntensitySlider->setEnabled(false);
    }
    _effect_intensity = 5;
    _effectIntensitySlider->setValue(5);
}

void OptionWidget::setRelief( bool addRelief ){
    if( _reliefGroupBox->isChecked() )
        _TOSParameters.relief = 1;
    else
        _TOSParameters.relief = 0;
}

void OptionWidget::setShakingType( int smodel ){
    _TOSParameters.smodel = smodel;
}

void OptionWidget::setBlur( bool addBlur ){
    if( _blurGroupBox->isChecked() )
        _TOSParameters.blur = 1;
    else
        _TOSParameters.blur = 0;
}

void OptionWidget::setColorSketch( bool color_sketch ){
    if( _colorSketchRadioButton->isChecked() )
        _TOSParameters.color_sketch = 1;
    else
        _TOSParameters.color_sketch = 0;
}

void OptionWidget::setRenderModel( int renderModel ){
    _TOSParameters.model = renderModel;
    if( renderModel == 4 ){
        _dictionaryParameterGroupBox->setVisible(true);
        _dictionariesWidget->setCurrentRow(0);
        _dictionariesWidget->setVisible(true);
    } else {
        _dictionaryParameterGroupBox->setVisible(false);
        _dictionariesWidget->setVisible(false);
    }
}

void OptionWidget::setRenderOrder( int renderOrder ){
    _TOSParameters.order = renderOrder;  //rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2"
}

void OptionWidget::setSelectionModel( int selectionModel ){
    _dictionaryParameters.randS = selectionModel;  //     "selection model: randS=0, randomly select shapes;
    //              randS=1, select shapes according to elongation, compactness and scale;
    //             randS=2, select shapes according to elongation, compactness, scale and color",
}

void OptionWidget::setColorModel( int colorModel ){
    _dictionaryParameters.mcolor = colorModel;  //   "color model: mcolor=1, use the color of tranferred image, otherwise use the color of the orignal image",
}
void OptionWidget::setShapeScaling( int shapeScaling ){
    _dictionaryParameters.equal = shapeScaling;  //"scaling shape with equal aspect ratio or not, if equal=1, scaling shape with equal aspect ratio, otherwise scaling    shape with different aspect ratio on x-axis and y-axis",
}

void OptionWidget::setEffectIntensity( int effect_intensity ){
    int mode = _modeComboBox->currentIndex();

    if( mode==0 ){
        _TOSParameters.threshold = effect_intensity*1./(_effect_intensity_nb); //  "threshold for color filtering",
        _TOSParameters.mpixel = 5 + effect_intensity*95./(_effect_intensity_nb); //minimal area (in pixel) for FLST",
        _TOSParameters.kappa = 0.3*effect_intensity/(_effect_intensity_nb); //compactness parameter of the attribute filtering on the orignal image",
        _TOSParameters.eps= -10 + 10.*effect_intensity/_effect_intensity_nb; // "-log10(max number of false alarms)",
        _minAreaSpinBox->setValue(_TOSParameters.mpixel);
        _minAreaSlider->setValue(_TOSParameters.mpixel);
        _scaleRatioSpinBox->setValue(_TOSParameters.ns);
        _scaleRatioSlider->setValue( _TOSParameters.ns );
        _thresholdSpinBox->setValue(_TOSParameters.threshold);
        _thresholdSlider->setValue( int(_TOSParameters.threshold/0.1) );
        _compactnessSpinBox->setValue(_TOSParameters.kappa);
        _compactnessSlider->setValue( int(_TOSParameters.kappa/0.1) );
        if( _TOSParameters.color_sketch == 1)
            _colorSketchRadioButton->setChecked(true);
        else
            _colorSketchRadioButton->setChecked(false);
        _epsSpinBox->setValue(_TOSParameters.eps);
        _epsSlider->setValue( int(_TOSParameters.eps));
    } else if( mode==1 ){
        _TOSParameters.median = 6 + effect_intensity*14/(_effect_intensity_nb);
        //  _TOSParameters.kerSize = effect_intensity*10./(_effect_intensity_nb);
        _TOSParameters.kerStd = 0.5 + 0.5*effect_intensity/(_effect_intensity_nb);
        _medianSizeSpinBox->setValue(_TOSParameters.median);
        _medianSizeSlider->setValue(_TOSParameters.median);
        _kerSizeSpinBox->setValue(_TOSParameters.kerSize);
        _kerSizeSlider->setValue( _TOSParameters.kerSize );
        _kerStdSpinBox->setValue(_TOSParameters.kerStd);
        _kerStdSlider->setValue( int(_TOSParameters.kerStd/0.1) );
        _TOSParameters.shift = 5. + effect_intensity*10./(_effect_intensity_nb);
        _TOSParameters.theta = 0. + effect_intensity/(_effect_intensity_nb);
        _shiftSpinBox->setValue(_TOSParameters.shift);
        _shiftSlider->setValue(int(_TOSParameters.shift/0.1));
        _thetaSpinBox->setValue(_TOSParameters.theta);
        _thetaSlider->setValue( _TOSParameters.theta/0.1 );
    } else if( mode==2 ){
        _TOSParameters.shift = 5. + effect_intensity*10./(_effect_intensity_nb);
        _TOSParameters.theta = 0. + effect_intensity/(_effect_intensity_nb);
        _shiftSpinBox->setValue(_TOSParameters.shift);
        _shiftSlider->setValue(int(_TOSParameters.shift/0.1));
        _thetaSpinBox->setValue(_TOSParameters.theta);
        _thetaSlider->setValue( _TOSParameters.theta/0.1 );
    } else if( mode==3 ){
        _TOSParameters.mpixel = 30 + effect_intensity*340/(_effect_intensity_nb);
        _minAreaSpinBox->setValue(_TOSParameters.mpixel);
        _minAreaSlider->setValue(_TOSParameters.mpixel);
    }
    _effect_intensity = effect_intensity;
}

void OptionWidget::setNbShapes(int shapeNb){
    _largeShapeSpinBox->setMaximum(shapeNb);
    _largeShapeSlider->setMaximum(shapeNb);
    _largeShapeSlider->setValue(shapeNb);
    _largeShapeSpinBox->setValue(shapeNb);
}
