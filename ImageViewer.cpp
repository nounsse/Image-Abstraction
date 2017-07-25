/****************************************************************************
**
** Copyright (C) 2012 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
**     the names of its contributors may be used to endorse or promote
**     products derived from this software without specific prior written
**     permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
** $QT_END_LICENSE$
**
****************************************************************************/

#include <QtGui>
#include <QDialog>
#include <QMessageBox>
#include <QFileDialog>
#include <QAction>
#include <QMenu>
#include <QMenuBar>
#include <QToolBar>
#include "OptionWidget.h"
#include "ImageViewer.h"

#include <fstream>
#include <iostream>

#define ZOOM_PMIN 5
#define ZOOM_COEFF 0.95

//! [0]
ImageViewer::ImageViewer()
{

    QWidget * widget = new QWidget;
    QHBoxLayout * imageLayout = new QHBoxLayout( widget );

    _imageWidget = new ImageWidget;
    _croppedImageWidget = new ImageWidget;
    _resultImageWidget = new ImageWidget;

    QImage color_image( _imageWidget->width(), _imageWidget->height(), QImage::Format_RGB32);
    color_image.fill(QColor(202,201,200));
    _imageWidget->setImage(color_image);
    _croppedImageWidget->setImage(color_image);
    _resultImageWidget->setImage(color_image);

    connect(_imageWidget, SIGNAL(updateImage()), this, SLOT(updateImage()));
    connect(_resultImageWidget, SIGNAL(updatePointClicked()), this, SLOT(updatePointClicked()));

    imageLayout->addWidget(_imageWidget);

    QVBoxLayout * resultLayout = new QVBoxLayout( );
    resultLayout->addWidget(_croppedImageWidget);
    resultLayout->addWidget(_resultImageWidget);

    imageLayout->addLayout(resultLayout);

    setCentralWidget(widget);

    createActions();
    createMenus();
    createToolBars();

    optionWidget = new OptionWidget(this);

    this->addDockWidget(Qt::RightDockWidgetArea, optionWidget);

    setWindowTitle(tr("Image abstraction"));
    resize(1200, 675);

    _croppedImageWidget->setVisible(false);

    _segmentation_displayed = false;
    _dictionary_loaded = false;
    _image_loaded = false;
    loadDictionaries();

}

void ImageViewer::loadDictionaries(){

#if 0

    addDictionary(QString("/home/nfaraj/Documents/Image_abstraction/ImageAbstractionInterface/dictionaries/delaunay_segmentation.png"),
                  QString("/home/nfaraj/Documents/Image_abstraction/ImageAbstractionInterface/dictionaries/delaunay.png"));
    _dictionary_loaded = true;
#else
    const QString folderPath ("./dictionaries");
    if(!folderPath.isEmpty())
    {
        QDir dir(folderPath);
        QStringList filter;
        filter << QLatin1String("*.png");
        filter << QLatin1String("*.jpeg");
        filter << QLatin1String("*.jpg");
        filter << QLatin1String("*.tiff");
        dir.setNameFilters(filter);
        QFileInfoList filelistinfo = dir.entryInfoList();
        foreach (const QFileInfo &fileinfo, filelistinfo) {
            QString dict_name = fileinfo.absoluteFilePath();

            QFileInfo dict_info (dict_name);
            QString basename = dict_info.baseName();

            QString seg ("_segmentation");
            if( ! basename.contains( seg )){
                basename.append(seg);
                //imageFile is the image path, just put your load image code here
                //look for segmentation
                QFileInfoList segmentationlistinfo = dir.entryInfoList();
                QString seg_name;
                bool segmentation_found = false;

                foreach (const QFileInfo &segmentationfileinfo, segmentationlistinfo) {
                    QFileInfo seg_info (segmentationfileinfo.absoluteFilePath());
                    QString seg_basename = seg_info.baseName();

                    //imageFile is the image path, just put your load image code here
                    if( seg_basename.contains( basename ) ){
                        segmentation_found = true;
                        seg_name = segmentationfileinfo.absoluteFilePath();
                    }
                }

                if( segmentation_found ){
                    addDictionary(seg_name, dict_name);
                    std::cout << "Found " << seg_name.toStdString() << " and " <<dict_name.toStdString() << std::endl;
                } else
                    addDictionary(dict_name);
                _dictionary_loaded = true;
            }
        }
    }
#endif
}

void ImageViewer::addDictionary(const QString &dict_name){
    QImage image_dict(dict_name);
    if (image_dict.isNull()) {
        QMessageBox::information(this, tr("Image Viewer"),
                                 tr("Cannot load %1.").arg(dict_name));
        return;
    }

    TreeOfShapes * dictionary = new TreeOfShapes(cfimages_from_qimage(image_dict));
    dictionary->compute_tree( getDefaultTOSParameters(), true);
    _dictionaries.push_back(dictionary);

    optionWidget->addDictionary( dict_name );
}

void ImageViewer::addDictionary(const QString &dict_name, const QString &texture_name){
    QImage image_dict(dict_name);
    QImage image_texture_dict(texture_name);
    if (image_dict.isNull() || image_texture_dict.isNull()) {
        QMessageBox::information(this, tr("Image Viewer"),
                                 tr("Cannot load %1.").arg(texture_name));
        return;
    }

    std::map<QRgb, int> color_map;
    std::vector<int> average_red;
    std::vector<int> average_green;
    std::vector<int> average_blue;
    std::vector<int> pixel_per_colors;
    int nb_of_segments = 0;
    for (int y = 0; y < image_dict.height(); y++) {
        for (int x = 0; x < image_dict.width(); x++) {
            QColor color = image_dict.pixel( x, y );
            QColor texture_color = image_texture_dict.pixel(x, y);
            QRgb color_id = color.rgb();
            std::map<QRgb, int>::iterator it = color_map.find(color_id);
            if(it == color_map.end()){
                //                            std::cout << color.red() << ", " <<
                //                                         color.green() << ", " <<
                //                                                 color.blue() << std::endl;
                //                                std::cout << color_id << " and " << average_red.size() << std::endl;
                color_map[color_id] = average_red.size();
                average_red.push_back(texture_color.red());
                average_green.push_back(texture_color.green());
                average_blue.push_back(texture_color.blue());
                pixel_per_colors.push_back(1);
            } else {
                // std::cout << "found "<< color_id << " and " << it->second << std::endl;
                average_red[it->second] += texture_color.red();
                average_green[it->second] += texture_color.green();
                average_blue[it->second] += texture_color.blue();
                pixel_per_colors[it->second]++;
            }
        }
    }

    QImage test (image_dict);

    for (int y = 0; y < image_dict.height(); y++) {
        for (int x = 0; x < image_dict.width(); x++) {
            QColor color = image_dict.pixel( x, y );
            int color_id = color_map[color.rgb()];
            QColor average_color (int(float(average_red[color_id])/pixel_per_colors[color_id]),
                                  int(float(average_green[color_id])/pixel_per_colors[color_id]),
                                  int(float(average_blue[color_id])/pixel_per_colors[color_id]));
            //            std::cout << color_id << std::endl;
            //            std::cout << float(average_red[color_id])/pixel_per_colors[color_id] << ", " <<
            //                         int(float(average_green[color_id])/pixel_per_colors[color_id])<< ", " <<
            //                                 int(float(average_blue[color_id])/pixel_per_colors[color_id]) << std::endl;

            //            std::cout << average_red[color_id]<< ", " <<
            //                         average_green[color_id]<< ", " <<
            //                                 average_blue[color_id] << std::endl;
            //            std::cout << pixel_per_colors[color_id] << std::endl;
            test.setPixel(x, y , qRgb(average_color.red(), average_color.green(), average_color.blue()));
        }
    }
  //  test.save(QString("./Test.png"));
    TreeOfShapes * dictionary = new TreeOfShapes(cfimages_from_qimage(test), cfimages_from_qimage(image_texture_dict));
    dictionary->compute_tree( getDictionaryTOSParameters(), true);
    _dictionaries.push_back(dictionary);

    optionWidget->addDictionary( texture_name );
}

Cfimage ImageViewer::cfimages_from_qimage( const QImage &input_image  ){

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


void ImageViewer::open()
//! [1] //! [2]
{

    QString selectedFilter;
    QString fileFilter = "Known Filetypes (*.tiff *.png *.jpg *.jpeg);;TIFF (*.tiff);;PNG (*.png);;JPEG (*.jpg *.jpeg)";

    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("Select an input image"),
                                                    QDir::currentPath(),
                                                    fileFilter,
                                                    &selectedFilter);

    if (!fileName.isEmpty()) {
        QImage image(fileName);
        if (image.isNull()) {
            QMessageBox::information(this, tr("Image Viewer"),
                                     tr("Cannot load %1.").arg(fileName));
            return;
        }

        //! [2] //! [3]
        _imageWidget->setImage(image);

        imageAbstractionProcess = new AbstractionProcess(image);

        segmentation = new Segmentation (image);
        //! [3] //! [4]

        printAct->setEnabled(true);
        _image_loaded = true;

        updateActions();

        fitToWindow();

        _croppedImageWidget->setImage(image);

    }
}
//! [4]

//! [5]
void ImageViewer::print()
//! [5] //! [6]
{
    /*
    Q_ASSERT(imageLabel->pixmap());
#ifndef QT_NO_PRINTER
//! [6] //! [7]
    QPrintDialog dialog(&printer, this);
//! [7] //! [8]
    if (dialog.exec()) {
        QPainter painter(&printer);
        QRect rect = painter.viewport();
        QSize size = imageLabel->pixmap()->size();
        size.scale(rect.size(), Qt::KeepAspectRatio);
        painter.setViewport(rect.x(), rect.y(), size.width(), size.height());
        painter.setWindow(imageLabel->pixmap()->rect());
        painter.drawPixmap(0, 0, *imageLabel->pixmap());
    }
#endif
*/
}
//! [8]

//! [9]
void ImageViewer::zoomIn()
//! [9] //! [10]
{
    _imageWidget->zoomIn();
}

void ImageViewer::zoomOut()
{
    _imageWidget->zoomOut();
}

//! [10] //! [11]
void ImageViewer::normalSize()
//! [11] //! [12]
{
    _imageWidget->normalSize();
}
//! [12]

//! [13]
void ImageViewer::fitToWindow()
//! [13] //! [14]
{

    _imageWidget->fitToWindow();

}
//! [14]


//! [15]
void ImageViewer::about()
//! [15] //! [16]
{
    QMessageBox::about(this, tr("About Image Viewer"),
                       tr("<p>The <b>Image Viewer</b> example shows how to combine QLabel "
                          "and QScrollArea to display an image. QLabel is typically used "
                          "for displaying a text, but it can also display an image. "
                          "QScrollArea provides a scrolling view around another widget. "
                          "If the child widget exceeds the size of the frame, QScrollArea "
                          "automatically provides scroll bars. </p><p>The example "
                          "demonstrates how QLabel's ability to scale its contents "
                          "(QLabel::scaledContents), and QScrollArea's ability to "
                          "automatically resize its contents "
                          "(QScrollArea::widgetResizable), can be used to implement "
                          "zooming and scaling features. </p><p>In addition the example "
                          "shows how to use QPainter to print an image.</p>"));
}
//! [16]

//! [17]
void ImageViewer::createActions()
//! [17] //! [18]
{
    openAct = new QAction(tr("&Open image..."), this);
    openAct->setShortcut(tr("Ctrl+O"));
    openAct->setIcon( QIcon("icons/open_icon.png") );
    connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

    saveAct = new QAction(tr("&Save current image..."), this);
    saveAct->setShortcut(tr("Ctrl+S"));
    saveAct->setIcon( QIcon("icons/save_icon.png") );
    saveAct->setEnabled(false);
    connect(saveAct, SIGNAL(triggered()), this, SLOT(save()));


    saveAllAct = new QAction(tr("&Save all images..."), this);
    saveAllAct->setShortcut(tr("Ctrl+Shift+S"));
    saveAllAct->setIcon( QIcon("icons/save_all_icon.png") );
    saveAllAct->setEnabled(false);
    connect(saveAllAct, SIGNAL(triggered()), this, SLOT(saveAll()));

    printAct = new QAction(tr("&Print..."), this);
    printAct->setShortcut(tr("Ctrl+P"));
    printAct->setEnabled(false);
    connect(printAct, SIGNAL(triggered()), this, SLOT(print()));

    exitAct = new QAction(tr("E&xit"), this);
    exitAct->setShortcut(tr("Ctrl+Q"));
    exitAct->setIcon( QIcon("icons/exit_icon.png") );

    connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));

    zoomInAct = new QAction(tr("Zoom &In (25%)"), this);
    //  zoomInAct->setShortcut(tr("Ctrl++"));
    zoomInAct->setEnabled(false);
    connect(zoomInAct, SIGNAL(triggered()), this, SLOT(zoomIn()));

    zoomOutAct = new QAction(tr("Zoom &Out (25%)"), this);
    //  zoomOutAct->setShortcut(tr("Ctrl+-"));
    zoomOutAct->setEnabled(false);
    connect(zoomOutAct, SIGNAL(triggered()), this, SLOT(zoomOut()));

    normalSizeAct = new QAction(tr("&Normal Size"), this);
    normalSizeAct->setShortcut(tr("Ctrl+N"));
    normalSizeAct->setEnabled(false);
    connect(normalSizeAct, SIGNAL(triggered()), this, SLOT(normalSize()));

    fitToWindowAct = new QAction(tr("&Reset size"), this);
    fitToWindowAct->setEnabled(false);
    fitToWindowAct->setShortcut(tr("Ctrl+F"));
    connect(fitToWindowAct, SIGNAL(triggered()), this, SLOT(fitToWindow()));

    aboutAct = new QAction(tr("&About"), this);
    connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));

    aboutQtAct = new QAction(tr("About &Qt"), this);
    connect(aboutQtAct, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}
//! [18]

//! [19]
void ImageViewer::createMenus()
//! [19] //! [20]
{
    fileMenu = new QMenu(tr("&File"), this);
    fileMenu->addAction(openAct);
    fileMenu->addAction(saveAct);
   // fileMenu->addAction(saveAllAct);
    fileMenu->addAction(printAct);
    fileMenu->addSeparator();
    fileMenu->addAction(exitAct);

    viewMenu = new QMenu(tr("&View"), this);
    viewMenu->addAction(zoomInAct);
    viewMenu->addAction(zoomOutAct);
    viewMenu->addAction(normalSizeAct);
    viewMenu->addSeparator();
    viewMenu->addAction(fitToWindowAct);

    helpMenu = new QMenu(tr("&Help"), this);
    helpMenu->addAction(aboutAct);
    helpMenu->addAction(aboutQtAct);

    menuBar()->addMenu(fileMenu);
    menuBar()->addMenu(viewMenu);
    menuBar()->addMenu(helpMenu);
}
//! [20]

void ImageViewer::createToolBars()
{
    fileToolBar = addToolBar(tr("File"));
    fileToolBar->addAction(openAct);
    fileToolBar->addAction(saveAct);
    //fileToolBar->addAction(saveAllAct);
    /*
    viewToolBar = addToolBar(tr("View"));
    viewToolBar->addAction(zoomOutAct);
    viewToolBar->addAction(zoomOutAct);
    viewToolBar->addAction(saveAllAct);
    */
}

//! [21]
void ImageViewer::updateActions()
//! [21] //! [22]
{
    zoomInAct->setEnabled(true);
    zoomOutAct->setEnabled(true);
    normalSizeAct->setEnabled(true);
    fitToWindowAct->setEnabled(true);
}
//! [22]

void ImageViewer::updateImage()
//! [23]
{

    QImage resulting_image = _imageWidget->crop();
    _croppedImageWidget->setImage(resulting_image);
    _croppedImageWidget->fitToWindow();

}
//! [22]

void ImageViewer::updatePointClicked()
//! [23]
{
    if( _segmentation_displayed ){
        QPoint clickedPixel = _resultImageWidget->getClickedPixel();
        QImage resulting_image = segmentation->removeRegionUnder(clickedPixel);
        _resultImageWidget->updateImage(resulting_image);
    }
}

void ImageViewer::displayCropImage(bool display){
    if( display ){
        _croppedImageWidget->resize( _resultImageWidget->size() );
        _croppedImageWidget->fitToWindow();
    }
    _croppedImageWidget->setVisible( display );
}

void ImageViewer::segment( ){

    if( _image_loaded ){
        QRect selection = _imageWidget->getSelection();

        QPoint A = selection.topLeft();
        QPoint B = selection.bottomRight();

        SegParameters segParameters = optionWidget->getSegParameters();

        QImage resulting_image = segmentation->segment( segParameters.sigma, segParameters.c, segParameters.min_size );
        _resultImageWidget->setImage(resulting_image);
        _resultImageWidget->fitToWindow();
        saveAct->setEnabled(true);
        saveAllAct->setEnabled(false);
        _segmentation_displayed = true;
    }
}


void ImageViewer::updateFromSegmentation(){
    if( _image_loaded & _segmentation_displayed ){
        QImage image = segmentation->getResult();
        _imageWidget->setImage(image);

        imageAbstractionProcess = new AbstractionProcess(image);

        segmentation = new Segmentation (image);
        //! [3] //! [4]

        _segmentation_displayed = false;
        printAct->setEnabled(true);

        updateActions();

        fitToWindow();

        QImage color_image( _imageWidget->width(), _imageWidget->height(), QImage::Format_RGB32);
        color_image.fill(QColor(202,201,200));
        _imageWidget->setImage(image);
        _resultImageWidget->setImage(color_image);
    }

}

void ImageViewer::run( ){



    if( _image_loaded ){
        QRect selection = _imageWidget->getSelection();

        QPoint A = selection.topLeft();
        QPoint B = selection.bottomRight();


        TOSParameters tosParameters = optionWidget->getTOSParameters();

        QImage resulting_image;

#if 0
        char * dictname = "../Data/birds_dict.tiff";
        resulting_image = imageAbstractionProcess->run(optionWidget->getProcessingMode(), dictname);
#else

        bool tree_recomputed;
        if( tosParameters.model == 4 ){

            resulting_image = imageAbstractionProcess->render(tosParameters, tree_recomputed, optionWidget->getDictionaryParameters(), _dictionaries[optionWidget->getSelectedDictionary()]);
        } else {
            resulting_image = imageAbstractionProcess->render(tosParameters, tree_recomputed);
        }

        if(tree_recomputed)
            optionWidget->setNbShapes(imageAbstractionProcess->getMaxArea());
#endif




        _resultImageWidget->setImage(resulting_image);
        _resultImageWidget->fitToWindow();
        saveAct->setEnabled(true);
        saveAllAct->setEnabled(true);

        _segmentation_displayed = false;
    }
}


int ImageViewer::computeNbLevels(const QImage & image){


    float coeff = ZOOM_COEFF;
    int pixel_min = ZOOM_PMIN;
    int size = std::min( image.width(), image.height() );

    return int( log10(float(pixel_min)/float(size))/log10(coeff) );

}

float ImageViewer::computeZoomCoefficient( float size, int nbLevels){

    float log_c = log10(float(ZOOM_PMIN)/float(size))/nbLevels;

    return pow(10, log_c);

}

void ImageViewer::saveTree(){
    QString fileName = QFileDialog::getSaveFileName(this, "Save current tree nodes file as ", QDir::currentPath(), "TREE (*.tree)");

    // In case of Cancel
    if ( fileName.isEmpty() ) {
        return;
    }

    if( ! fileName.endsWith( ".tree" ) )
        fileName.append(".tree");

    std::vector<QPoint> positions;
    std::vector<int> heights;
    std::vector<std::pair<int, int> > edges;
    std::vector< std::vector<std::pair<int, int> > > pixels;
    std::vector<QColor> colors;
    imageAbstractionProcess->getTreeInfo(positions, colors, pixels, heights, edges);

    std::ofstream out (fileName.toUtf8());
    if (!out)
        exit (EXIT_FAILURE);

    out << positions.size() << " " << edges.size() << std::endl;

    for(unsigned int i = 0 ; i < positions.size(); i++){
        const QPoint & p = positions[i];
        const QColor & color = colors[i];
        out << p.x() << " " << p.y()<< " " << heights[i] << " " << color.red() << " " << color.green() << " "<< color.blue()<< std::endl;
        out << pixels[i].size() << std::endl;
        for(unsigned int j = 0 ; j <  pixels[i].size(); j++){
            out << pixels[i][j].first << " " << pixels[i][j].second << std::endl;
        }
    }

    for(unsigned int i = 0 ; i < edges.size(); i++){
        out << edges[i].first << " " << edges[i].second << std::endl;
    }

    out << std::endl;
    out.close ();
}

void ImageViewer::save(){


    QString fileName = QFileDialog::getSaveFileName(this, "Save current image file as ", QDir::currentPath(), "TIFF (*.tiff);;PNG (*.png)");

    // In case of Cancel
    if ( fileName.isEmpty() ) {
        return;
    }

    if( ! fileName.endsWith( ".tiff" ) && ! fileName.endsWith(".png"))
        fileName.append(".png");

    _resultImageWidget->save( fileName );
    saveParameters(fileName);
}
void ImageViewer::saveParameters( const QString &filename){

    TOSParameters parameters = imageAbstractionProcess->getParameters();
    QFileInfo fileinfo (filename);
    QString name ( fileinfo.absolutePath() );
    name.append("/");
    name.append(fileinfo.baseName());
    name.append(".txt");
    std::ofstream myfile;
    myfile.open (name.toStdString());
    myfile << "order "<< parameters.order<< "\n";//rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2"
    myfile << "model "<< parameters.model<< "\n";//synthesis model: orignal shape: m=0; ellipse: m=1; rectangle: m=2; random: m=3;
    myfile << "alpha "<< parameters.alpha<< "\n";//alpha for transparent",
    myfile << "ns "<< parameters.ns<< "\n"; //  "scale ratio order for color filtering",
    myfile << "threshold "<< parameters.threshold<< "\n"; //  "threshold for color filtering",
    myfile << "smodel "<< parameters.smodel<< "\n"; //shaking type: uniform shaking: smode=0, dominant shaking: smode=1",
    myfile << "shift "<< parameters.shift<< "\n";//add a random shiftS to each shape, shiftS = shift*rand()",
    myfile << "theta "<< parameters.theta<< "\n";//add a random rotation to each shape, thetaS = theta*rand()",
    myfile << "mpixel "<< parameters.mpixel<< "\n"; //minimal area (in pixel) for FLST",
    myfile << "maxarea "<< parameters.maxarea<< "\n"; //minimal area (in pixel) for FLST",
    myfile << "kappa "<< parameters.kappa<< "\n";//compactness parameter of the attribute filtering on the orignal image",
    myfile << "relief "<< parameters.relief<< "\n";//     "add relief effects, if relief =1",
    myfile << "reliefOrientation "<< parameters.reliefOrientation<< "\n";//    "relief orentation, in degree",
    myfile << "reliefHeight "<< parameters.reliefHeight<< "\n"; //     "relief height",
    myfile << "blur "<< parameters.blur<< "\n";//     "add blur effects, if blur =1",
    myfile << "median "<< parameters.median<< "\n"; //kernel size for median filter",
    myfile << "kerSize "<< parameters.kerSize<< "\n";//kernel size for Gaussian blur",
    myfile << "kerStd "<< parameters.kerStd<< "\n";//std for the gaussian kernel",
    myfile << "color_sketch "<< parameters.color_sketch<< "\n";//     "compute the sketch: filter shapes based on contrast if=1",
    myfile << "eps "<< parameters.eps; // "-log10(max number of false alarms)",
    myfile.close();
    //imageAbstractionProcess->getParameters();
}
void ImageViewer::saveAll(){


#if 1
    QString folderName = QFileDialog::getExistingDirectory(this, tr("Save all images in selected directory"),
                                                           QDir::currentPath(),
                                                           QFileDialog::ShowDirsOnly
                                                           | QFileDialog::DontResolveSymlinks);

    // In case of Cancel
    if ( folderName.isEmpty() ) {
        return;
    }

    //imageAbstractionProcess->save_shapes( folderName, true );
    imageAbstractionProcess->save_shapes(folderName, optionWidget->getTOSParameters(), optionWidget->getDictionaryParameters(), _dictionaries[optionWidget->getSelectedDictionary()]);
    //_resultImageWidget->saveAll( folderName );
#else
        saveTree();
#endif
}
