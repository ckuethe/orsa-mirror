TEMPLATE = app

CONFIG += qt
QT     -= gui opengl

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib spice_include spice_lib

include(../../orsa.pri)

TARGET   = ThermalStressCraters

INCLUDEPATH += ../../src/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = .

unix:!macx {
	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaEssentialOSG -lorsaSolarSystem -lorsaSPICE -lOpenThreads
##-losg -losgText -losgGA -losgViewer
}

macx {
	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaQt -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads
}

HEADERS += ThermalStressCraters.h
SOURCES += ThermalStressCraters.cpp
