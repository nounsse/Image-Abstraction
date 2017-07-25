#include "ImageWidget.h"
#include <iostream>
#include <QFileInfo>

ImageWidget::ImageWidget(QWidget *parent, bool allow_crop) :
    QWidget(parent)
{

    _currentImageLabel = 0;
    _imageLabels.push_back(new ImageLabel);
    _imageLabels[_currentImageLabel]->setBackgroundRole(QPalette::Base);
    _imageLabels[_currentImageLabel]->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    _imageLabels[_currentImageLabel]->setScaledContents(true);

    _scrollAreas.push_back( new QScrollArea );
    _scrollAreas[_currentImageLabel]->setBackgroundRole(QPalette::Dark);
    _scrollAreas[_currentImageLabel]->setWidget(_imageLabels[_currentImageLabel]);
    _scrollAreas[_currentImageLabel]->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    _scaleFactors.push_back(1.);
    QVBoxLayout * vBoxLayout = new QVBoxLayout(this);
    _scrollAreasLayout = new QVBoxLayout;

    _fitToWindow = false;

    _scrollAreasLayout->addWidget( _scrollAreas[_currentImageLabel] );
    vBoxLayout->addLayout(_scrollAreasLayout);

    _allow_crop = allow_crop;
    if( _allow_crop ){
        _cropRadioButton = new QRadioButton("Crop");
        _cropRadioButton->setEnabled(false);
        vBoxLayout->addWidget(_cropRadioButton);

        connect(_cropRadioButton, SIGNAL(toggled(bool)), this, SLOT(activateCrop(bool) ) );
        connect(_cropRadioButton, SIGNAL(toggled(bool)), parent, SLOT(displayCropImage(bool))  );
    }

    _imageSlider = new QSlider( Qt::Horizontal);
    _imageSlider->setRange(0, 0);
    _imageSlider->setVisible( false );
    _imageSlider->setSingleStep(1);
    _imageSlider->setValue(0);
    vBoxLayout->addWidget(_imageSlider);

    connect(_imageSlider, SIGNAL(valueChanged(int)), this, SLOT(setCurrentImage(int)));

    connect(_imageLabels[_currentImageLabel], SIGNAL(selectionChanged()), this, SLOT(selectionChanged()));
    connect(_imageLabels[_currentImageLabel], SIGNAL(pointClicked()), this, SLOT(pointClicked()));
}

void ImageWidget::updateImage( const QImage & image ){

    if( _imageLabels.size() == 0 )
        return;

    _imageLabels[_currentImageLabel]->updateImage(image);

}

void ImageWidget::setImage( const QImage & image ){

    if( _imageLabels.size() > 1 )
        clear();

    _imageLabels[_currentImageLabel]->setImage(image);
    _scaleFactors[_currentImageLabel] = 1.0;

    _input_size = image.size();

    if( _allow_crop ){
        _cropRadioButton->setEnabled(true);
        _cropRadioButton->setChecked(false);
    }
}

void ImageWidget::setImages( const std::vector<QImage> & images ){

    if( images.size() > 0 ){

        setImage( images[0] );

        for( unsigned int i = 1 ; i < images.size(); i ++ ){
            addImage(images [i]);
        }
    }

}

void ImageWidget::addImage( const QImage & image ){

    _scrollAreas[_currentImageLabel]->setVisible(false);

    _currentImageLabel = _imageLabels.size();

    _scaleFactors.push_back(1.);

    _imageLabels.push_back(new ImageLabel);
    _imageLabels[_currentImageLabel]->setBackgroundRole(QPalette::Base);
    _imageLabels[_currentImageLabel]->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    _imageLabels[_currentImageLabel]->setScaledContents(true);
    _imageLabels[_currentImageLabel]->setImage(image);

    _scrollAreas.push_back( new QScrollArea );
    _scrollAreas[_currentImageLabel]->setBackgroundRole(QPalette::Dark);
    _scrollAreas[_currentImageLabel]->setWidget(_imageLabels[_currentImageLabel]);
    _scrollAreas[_currentImageLabel]->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);

    _scrollAreasLayout->addWidget( _scrollAreas[_currentImageLabel] );

    connect(_imageLabels[_currentImageLabel], SIGNAL(selectionChanged()), this, SLOT(selectionChanged()));

    _imageSlider->setRange(0, _currentImageLabel);
    _imageSlider->setValue(_currentImageLabel);
    _imageSlider->setVisible(true);

    _imageLabels[_currentImageLabel]->resize(_scrollAreas[_currentImageLabel]->size());
}

void ImageWidget::setCurrentImage(int currentImage){
    _scrollAreas[_currentImageLabel]->setVisible(false);
    _currentImageLabel = currentImage;
    _scrollAreas[_currentImageLabel]->setVisible(true);
}

void ImageWidget::adjustSize(){
    for( unsigned int i = 0 ; i < _imageLabels.size() ; i ++ ){
        _imageLabels[i]->adjustSize();
    }
}

void ImageWidget::scaleImages(double factor)
{
    for( unsigned int i = 0 ; i < _imageLabels.size() ; i ++ ){
        scaleImage(i, factor);
    }
}

void ImageWidget::scaleImage(int i, double factor)
{
    if( _imageLabels.size() > i ){
        Q_ASSERT(_imageLabels[i]->pixmap());
        _scaleFactors[i] *= factor;
        _imageLabels[i]->setScale(_scaleFactors[i]);

        adjustScrollBar(_scrollAreas[i]->horizontalScrollBar(), factor);
        adjustScrollBar(_scrollAreas[i]->verticalScrollBar(), factor);
    }

}


void ImageWidget::adjustScrollBar(QScrollBar *scrollBar, double factor)
{
    scrollBar->setValue(int(factor * scrollBar->value()
                            + ((factor - 1) * scrollBar->pageStep()/2)));
}

void ImageWidget::zoomIn()
{
    scaleImages(1.25);
}

void ImageWidget::zoomOut()
{
    scaleImages(0.75);
}

void ImageWidget::normalSize()
{
    for( unsigned int i = 0 ; i < _imageLabels.size() ; i ++ ){
        _imageLabels[i]->adjustSize();
        scaleImage(i, 1.0);
        _scaleFactors[i] = 1.;
    }
}

void ImageWidget::fitToWindow(){

    normalSize();

    for( unsigned int i = 0 ; i < _imageLabels.size() ; i ++ ){
        float scaleFactorToFit = float(this->width()-20)/_imageLabels[i]->pixmap()->width();

        if( _imageLabels[i]->pixmap()->height()*scaleFactorToFit >this->height() ){
            scaleFactorToFit = float(this->width() - 36)/_imageLabels[i]->pixmap()->width();
        }

        scaleImage(i, scaleFactorToFit);
    }
}

void ImageWidget::setImage( QString fileName ){

}

void ImageWidget::resizeEvent ( QResizeEvent * event ){

    if(_fitToWindow){
        fitToWindow();
    }
}

void ImageWidget::keyPressEvent(QKeyEvent *e)
{

    switch (e->key())
    {
    case Qt::Key_Plus : zoomIn();break;
    case Qt::Key_Minus : zoomOut();break;
    }

}

void ImageWidget::wheelEvent(QWheelEvent *event)
{

    if (event->modifiers()==Qt::ControlModifier){
        int delta = event->delta();
        if (delta > 0) {
            zoomIn();
        } else {
            zoomOut();
        }
    }
    event->accept();
}

QImage ImageWidget::crop(){
    return _imageLabels[_currentImageLabel]->crop();
}

void ImageWidget::activateCrop(bool activate){
    _imageLabels[_currentImageLabel]->setDrawSelection(activate);
}

void ImageWidget::clear(){

    _currentImageLabel = 0;
    for( unsigned int i = 1 ; i < _scrollAreas.size() ; i ++ ){
        _scrollAreas[i]->setVisible(false);
        _scrollAreasLayout->removeWidget(_scrollAreas[i]);
        delete _imageLabels[i];
        delete _scrollAreas[i];
    }
    QScrollArea * tmp = _scrollAreas[0];
    ImageLabel * tmp_img = _imageLabels[0];

    _imageLabels.clear();
    _scrollAreas.clear();
    _scaleFactors.clear();

    tmp->setVisible(true);
    _scrollAreas.push_back(tmp);
    _imageLabels.push_back(tmp_img);
    _scaleFactors.push_back(1.);
    _imageSlider->setVisible(false);
}

void ImageWidget::save( const QString &fileName ){

    _imageLabels[_currentImageLabel]->pixmap()->toImage().save(fileName);

}


void ImageWidget::saveAll( const QString &folderName ){

    QFileInfo fi(folderName);
    QString base_name = fi.baseName();

    for( unsigned int i=0 ; i < _imageLabels.size() ; i ++ ){
        QString fileName = folderName;
        fileName.append("/");
        fileName.append(base_name);
        fileName.append("_");
        fileName.append(QString::number(i));
        fileName.append(".png");
        _imageLabels[i]->pixmap()->toImage().save(fileName);
    }

}
