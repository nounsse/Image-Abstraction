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

#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include <QMainWindow>
#include "ImageWidget.h"
#include "OptionWidget.h"
#include "Segmentation.h"

QT_BEGIN_NAMESPACE
class QAction;
class QMenu;
QT_END_NAMESPACE

#include "AbstractionProcess.h"
//! [0]
class ImageViewer : public QMainWindow
{
    Q_OBJECT

public:
    ImageViewer();

private slots:
    void open();
    void save();
    void saveAll();
    void print();
    void zoomIn();
    void zoomOut();
    void normalSize();
    void fitToWindow();
    void about();
    void updateImage();
    void updatePointClicked();
    void run();
    void segment( );
    void updateFromSegmentation();
    void displayCropImage(bool display);
    void saveTree();
private:
    void saveParameters(const QString &filename);
    int computeNbLevels(const QImage & image);
    float computeZoomCoefficient(float size, int nbLevels);
    void createActions();
    void createMenus();
    void createToolBars();
    void updateActions();
    Cfimage cfimages_from_qimage( const QImage &input_image  );
    void loadDictionaries();
    void addDictionary(const QString &dict_name);
    void addDictionary(const QString &dict_name, const QString &texture_name);

    ImageWidget * _imageWidget;
    ImageWidget * _croppedImageWidget;
    ImageWidget * _resultImageWidget;

    AbstractionProcess *imageAbstractionProcess;

    Segmentation *segmentation;

    OptionWidget * optionWidget;

    float _zoom_coeff;

    bool _segmentation_displayed;
    bool _image_loaded;
    bool _dictionary_loaded;

    std::vector<TreeOfShapes *>_dictionaries;

    QAction *openAct;
    QAction *saveAct;
    QAction *saveAllAct;
    QAction *printAct;
    QAction *exitAct;
    QAction *zoomInAct;
    QAction *zoomOutAct;
    QAction *normalSizeAct;
    QAction *fitToWindowAct;
    QAction *aboutAct;
    QAction *aboutQtAct;

    QToolBar *fileToolBar;

    QMenu *fileMenu;
    QMenu *viewMenu;
    QMenu *helpMenu;
};
//! [0]

#endif
