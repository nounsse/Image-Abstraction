#ifndef IMAGELABEL_H
#define IMAGELABEL_H

#include <QLabel>
#include <QRubberBand>
#include <QMouseEvent>
#include <QImage>

class ImageLabel : public QLabel
{
    Q_OBJECT
protected:
    QRubberBand *_rubberBand;
    QPoint _click_pos;
    QPoint _previous_mouse_pos;
    QPoint _origin;
    QPoint _end;
    QRect _selection;

    float _current_scale;
    bool _selection_active;
    bool _move_selection;
    bool _draw_selection;

    QSize _input_size;
public:
    explicit ImageLabel(QWidget *parent = 0);
    
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);

    void paintEvent(QPaintEvent *e);

    void setScale( float scale );

    void setDrawSelection(bool draw_selection);
    QImage crop();

    QRect getSelection();

    QPoint getClickedPixel();

    void setImage( const QImage &image);
    void updateImage( const QImage &image);
signals:
    void selectionChanged();
    void pointClicked();
public slots:
    
};

#endif // IMAGELABEL_H
