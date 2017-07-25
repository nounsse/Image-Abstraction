#ifndef IMAGEWIDGET_H
#define IMAGEWIDGET_H

#include <QWidget>
#include <QScrollArea>
#include <QScrollBar>
#include <QRadioButton>
#include <QLabel>
#include <QMouseEvent>
#include <QSlider>
#include <QVBoxLayout>

#include "ImageLabel.h"


class ImageWidget : public QWidget
{
    Q_OBJECT

protected:

    QVBoxLayout * _scrollAreasLayout;
    std::vector<QScrollArea *>_scrollAreas;
    std::vector<ImageLabel *> _imageLabels;
    int _currentImageLabel;
    QSlider * _imageSlider;

    std::vector<double> _scaleFactors;
    bool _fitToWindow;

    QSize _input_size;

    QRadioButton * _cropRadioButton;
    bool _allow_crop;

    void update();

    void scaleImages(double factor);
    void scaleImage(int i, double factor);
    void adjustScrollBar(QScrollBar *scrollBar, double factor);

    void resizeEvent ( QResizeEvent * event );
    void keyPressEvent(QKeyEvent *e);

    void wheelEvent(QWheelEvent *event);

    float scale;
public:
    explicit ImageWidget(QWidget *parent = 0, bool allow_crop=false);
    void zoomIn();
    void zoomOut();

    void normalSize();

    void setImage( const QImage & image );
    void updateImage( const QImage & image );
    void addImage( const QImage & image );
    void setImages( const std::vector<QImage> &images );

    void setImage( QString fileName );

    void adjustSize();

    QImage crop();
    QRect getSelection(){ return _imageLabels[_currentImageLabel]->getSelection(); }

    QPoint getClickedPixel(){ return _imageLabels[_currentImageLabel]->getClickedPixel(); }

    void clear();

    void save( const QString &fileName );
    void saveAll( const QString &folderName );

signals:
    void updateImage();
    void updatePointClicked();

public slots:
    void fitToWindow();
    void selectionChanged(){ emit updateImage();}
    void pointClicked(){ emit updatePointClicked();}

    void activateCrop(bool activate);
    void setCurrentImage( int currentImage );
};

#endif // IMAGEWIDGET_H
