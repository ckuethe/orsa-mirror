TEMPLATE = app
TARGET   = testGUI

CONFIG += qt
QT     += gui opengl

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib qwt_include qwt_lib spice_include spice_lib tbb_include tbb_lib

include(../../orsa.pri)

INCLUDEPATH += ../../src/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../bin/$${PLATFORM_NAME}

unix:LIBS += -L ../../lib/$${PLATFORM_NAME} -lorsa -lorsaOSG -lorsaQt -losgViewer -losgDB -losgdb_freetype -losgText -losgUtil -losgGA -losg -lOpenThreads -lgmp -lgmpxx -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads -lorsaUtil -losgDB -losg -lOpenThreads -lcurl

SOURCES = testGUI.cpp
HEADERS = 

