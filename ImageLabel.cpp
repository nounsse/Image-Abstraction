#include "ImageLabel.h"
#include <QPainter>

#include <iostream>

ImageLabel::ImageLabel(QWidget *parent) :
    QLabel(parent)
{

    _draw_selection = false;
    _selection_active = false;
    _current_scale = 1.;
    _rubberBand = new QRubberBand(QRubberBand::Rectangle, this);
    _click_pos = QPoint(0,0);

}


void ImageLabel::mousePressEvent(QMouseEvent *event)
{
        _click_pos = event->pos();
    if( _draw_selection ){
         _move_selection = false;

        if( _selection_active
                && std::min( _origin.x(), _end.x() ) < _click_pos.x() && std::max( _origin.x(), _end.x() ) > _click_pos.x()
                && std::min( _origin.y(), _end.y() ) < _click_pos.y() && std::max( _origin.y(), _end.y() ) > _click_pos.y()){
            _move_selection = true;
            _previous_mouse_pos = _click_pos;
        } else {

            delete _rubberBand;
            _rubberBand = new QRubberBand(QRubberBand::Rectangle, this);

            _rubberBand->setGeometry(QRect(_click_pos, QSize()));
            _rubberBand->show();
            _selection_active = true;
        }
    }
    emit pointClicked();
}

void ImageLabel::mouseMoveEvent(QMouseEvent *event)
{
    if( _draw_selection ){
        if( _move_selection ){
            QPoint translation = event->pos() - _previous_mouse_pos;
            _origin = _origin+translation;
            _end = _end+translation;
            _previous_mouse_pos = event->pos();
        } else {
            _origin = _click_pos;
            _end = event->pos();

        }
        _rubberBand->setGeometry(QRect(_origin, _end).normalized());
        emit selectionChanged();
    }
}

void ImageLabel::mouseReleaseEvent(QMouseEvent *event)
{
    if( _draw_selection ){
        QPoint release_pos = event->pos();
        _selection_active = true;

        if(! _move_selection )
            if( abs(_click_pos.x() - release_pos.x()) > 1 && abs(_click_pos.y() - release_pos.y()) > 1 ){
                _origin = _click_pos;
                _end = release_pos;
                emit selectionChanged();
            } else {
                _selection_active = false;
            }

        _move_selection = false;
    }
    // if (rubberBand)
    //    rubberBand->hide();
    // determine selection, for example using QRect::intersects()
    // and QRect::contains().


}

void ImageLabel::updateImage( const QImage &image){
    this->setPixmap(QPixmap::fromImage(image));
}

void ImageLabel::setImage( const QImage &image){
    _current_scale = 1.;
    _input_size = image.size();
    this->setPixmap(QPixmap::fromImage(image));
}

void ImageLabel::setScale(float scale){

    this->resize(scale * this->pixmap()->size());

    _origin = (_origin/_current_scale)*scale;
    _end = (_end/_current_scale)*scale;

    if( _draw_selection && _selection_active )
        _rubberBand->setGeometry(QRect(_origin, _end).normalized());

    _current_scale = scale;
}


QImage ImageLabel::crop(){

    QImage image = this->pixmap()->toImage();

    if( _selection_active && _draw_selection ){
        QPoint A( std::min( _origin.x(), _end.x() ), std::min( _origin.y(), _end.y() ) );
        QPoint B( std::max( _origin.x(), _end.x() ), std::max( _origin.y(), _end.y() ) );

        A = A/_current_scale;
        B = B/_current_scale;

        A = QPoint ( std::min( std::max( 0, A.x() ), this->pixmap()->width()-1 ), std::min( std::max( 0, A.y() ), this->pixmap()->height()-1 ) );
        B = QPoint ( std::min( std::max( 0, B.x() ), this->pixmap()->width()-1 ), std::min( std::max( 0, B.y() ), this->pixmap()->height()-1 ) );

        image = image.copy(QRect( A, B ));
    }
    return image;
}


QRect ImageLabel::getSelection(){
    QPoint A( 0, 0 );
    QPoint B( this->pixmap()->width(), this->pixmap()->height() );

    if( _selection_active && _draw_selection ){
        A = QPoint( std::min( _origin.x(), _end.x() ), std::min( _origin.y(), _end.y() ) );
        B = QPoint( std::max( _origin.x(), _end.x() ), std::max( _origin.y(), _end.y() ) );

        A = A/_current_scale;
        B = B/_current_scale;

        A = QPoint ( std::min( std::max( 0, A.x() ), this->pixmap()->width() ), std::min( std::max( 0, A.y() ), this->pixmap()->height() ) );
        B = QPoint ( std::min( std::max( 0, B.x() ), this->pixmap()->width() ), std::min( std::max( 0, B.y() ), this->pixmap()->height() ) );
    }
    return QRect( A, B );

}
QPoint ImageLabel::getClickedPixel(){
    return _click_pos/_current_scale;
}

void ImageLabel::setDrawSelection(bool draw_selection) {
    _draw_selection = draw_selection;

    if(_draw_selection) {

        if( !_selection_active ){
            delete _rubberBand;
            _rubberBand = new QRubberBand(QRubberBand::Rectangle, this);

            QSize image_size = this->pixmap()->size();
            _origin = QPoint( image_size.width()*_current_scale*0.25, image_size.height()*_current_scale*0.25 );
            _end = QPoint( image_size.width()*0.75*_current_scale, image_size.height()*_current_scale*0.75 );
            _rubberBand->setGeometry(QRect(_origin, _end).normalized());

            _selection_active = true;
        }
        _rubberBand->show();
        emit selectionChanged();
    } else {
        _rubberBand->hide();
    }

}

void ImageLabel::paintEvent(QPaintEvent * e){
    QLabel::paintEvent(e);

    if(_selection_active)
    {

        //TODO draw center
    }
}
