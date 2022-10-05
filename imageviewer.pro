MOC_DIR = ./moc
OBJECTS_DIR = ./obj

QT+= gui declarative

TARGET = image_abstraction

INCLUDEPATH += src/flst  \
    src/synth \
    src/mw3 \
    src/segment \
    src/kdtree \
    src/kdtree/ann/include
DEPENDPATH += src/flst \
    src/synth \
    src/segment \
    src/mw3 \
    src/kdtree \
    src/kdtree/ann/include

HEADERS       = ImageWidget.h \
    ImageLabel.h \
    ImageViewer.h \
    OptionWidget.h \
    src/flst/*.h \
    src/iio/*.h \
    src/mw3/*.h \
    src/synth/*.h \
    src/segment/*.h \
    src/kdtree/*.h \
    src/kdtree/ann/include/ANN/*.h \
    AbstractionProcess.h \
    Segmentation.h \
    TreeOfShapes.h
SOURCES       = main.cpp \
    ImageWidget.cpp \
    ImageLabel.cpp \
    ImageViewer.cpp \
    OptionWidget.cpp \
    src/flst/*.c \
    src/mw3/*.c \
    src/synth/*.c \
    src/kdtree/ann/src/*.cpp \
    src/segment/*.cpp \
    AbstractionProcess.cpp \
    Segmentation.cpp \
    TreeOfShapes.cpp
LIBS = -L/usr/lib/x86_64-linux-gnu/

QMAKE_CFLAGS += -std=c99
QMAKE_CXXFLAGS_RELEASE = -std=c++0x
QMAKE_CXXFLAGS_RELEASE += -Ofast

QMAKE_CXXFLAGS_DEBUG = -std=c++0x
QMAKE_CXXFLAGS_DEBUG += -Ofast
