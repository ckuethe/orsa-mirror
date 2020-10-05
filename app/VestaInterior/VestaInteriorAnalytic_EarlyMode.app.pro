TEMPLATE = app

CONFIG += qt
QT     -= gui opengl

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib spice_include spice_lib

include(../../orsa.pri)

TARGET   = VestaInteriorAnalytic_EarlyMode

INCLUDEPATH += ../../src/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = .

#unix:!macx {
#	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem -lorsaEssentialOSG -lorsaSPICE -lorsaPDS -lorsaUtil -lOpenThreads -lqd -lsqlite3 -lbsd
	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem -lorsaEssentialOSG -lorsaSPICE -lorsaPDS -lorsaUtil -lOpenThreads -lqd -lsqlite3 -lcrypto -lssl -lgmp -lmpfr
	#}

HEADERS += VestaInteriorAnalytic_EarlyMode.h   CCMD2SH.h shape.h  CubicChebyshevMassDistribution.h penalty.h
SOURCES += VestaInteriorAnalytic_EarlyMode.cpp CubicChebyshevMassDistribution.cpp 


