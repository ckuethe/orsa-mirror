TEMPLATE = app

CONFIG += qt
QT     -= gui opengl

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib spice_include spice_lib

include(../../orsa.pri)

TARGET   = dump3D

INCLUDEPATH += ../../src/ /home/tricaric/dislin/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = .

#unix:!macx {
#	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaEssentialOSG -lorsaSolarSystem -lorsaSPICE -lorsaPDS -lorsaUtil -lOpenThreads -lqd -L/home/tricaric/sqlite -lsqlite3 -L/home/tricaric/dislin/ -ldislin -lbsd
#}
LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaEssentialOSG -lorsaSolarSystem -lorsaSPICE -lorsaPDS -lorsaUtil -lOpenThreads -lqd -lsqlite3 -lcrypto -lssl -lgmp -lmpfr


HEADERS += CubicChebyshevMassDistribution.h shape.h simplex.h
SOURCES += dump3D.cpp CubicChebyshevMassDistribution.cpp


