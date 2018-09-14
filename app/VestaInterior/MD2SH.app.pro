TEMPLATE = app

CONFIG += qt
QT     -= gui opengl

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib spice_include spice_lib

include(../../orsa.pri)

TARGET   = MD2SH

INCLUDEPATH += ../../src/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = .

#unix:!macx {
    INCLUDEPATH += /usr/local/opt/openssl/include/
	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaEssentialOSG -lorsaSolarSystem -lorsaSPICE -lorsaPDS -lorsaUtil -lOpenThreads -lqd -lsqlite3 -lcrypto -lssl -lgmp -lmpfr
	#}


HEADERS += CubicChebyshevMassDistribution.h CCMD2SH.h vesta.h gaskell.h simplex.h penalty.h
SOURCES += MD2SH.cpp CubicChebyshevMassDistribution.cpp


