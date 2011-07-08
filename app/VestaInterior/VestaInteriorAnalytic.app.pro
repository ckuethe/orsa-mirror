TEMPLATE = app

CONFIG += qt
QT     -= gui opengl

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib qwt_include qwt_lib spice_include spice_lib

include(../../orsa.pri)

TARGET   = VestaInteriorAnalytic

INCLUDEPATH += ../../src/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../bin/$${PLATFORM_NAME}

unix:!macx {
	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaOSG -lorsaQt -lorsaSolarSystem -lorsaSPICE -lorsaPDS -lorsaUtil -losg -losgText -losgGA -losgViewer -lOpenThreads -losg -lGL -lGLU -losgDB -lOpenThreads -losg -lOpenThreads -losgUtil -losg -lOpenThreads
}


HEADERS += VestaInteriorAnalytic.h   vesta.h  CubicChebyshevMassDistribution.h
SOURCES += VestaInteriorAnalytic.cpp CubicChebyshevMassDistribution.cpp


