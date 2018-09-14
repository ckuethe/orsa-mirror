# version, for the libraries
# when using VERSION, libs on win32 are called liborsa1 instead of liborsa
# VERSION = 1.0.0

#CONFIG += thread debug warn_on
CONFIG += thread release warn_on

# enable this to build static libs, useful to create binaries for BOINC
#CONFIG += dll staticlib

#macx {
#	CONFIG += staticlib
#}

# less verbose compile
#CONFIG += silent

# not here... qt and gui are kept local on single directory's .pro files
#CONFIG += qt
#QT -= gui


#QMAKE_CXXFLAGS += -I/usr/local/opt/qt/include

# define ORSA_PRI as the absolute path to this file
unix:!macx {
	ORSA_PRI = /home/tricaric/orsa/orsa.pri
}
macx {
	ORSA_PRI = /Users/tricaric/orsa/orsa.pri
}
win32 {
	ORSA_PRI = C:\orsa\orsa.pri
}


# important tests, don't change
isEmpty(ORSA_PRI) {
	error(You need to set the ORSA_PRI variable correctly and to review all the other relevant settings in orsa.pri)
}
#
!exists($$ORSA_PRI) {
	error(ORSA_PRI is not set correctly [$$ORSA_PRI])
}

LIBS += -L/usr/local/opt/openssl/lib/ 

#QMAKE_MAC_SDK = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk

# platform name and other tweaking (unix,macx,win32)
unix:!macx {
	KERNEL_NAME   = $$system(uname -s)
	HARDWARE_NAME = $$system(uname -m)
	PLATFORM_NAME = $${KERNEL_NAME}_$${HARDWARE_NAME}

	DIR_SEP = "/"

#	QMAKE_CXXFLAGS += -std=c++0x

	QMAKE_CXXFLAGS_RELEASE +=
	QMAKE_LFLAGS_RELEASE   +=

	QMAKE_CXXFLAGS_DEBUG += -pg -ggdb
	QMAKE_LFLAGS_DEBUG   += -pg -ggdb
}
macx {
#	CONFIG += x86 # x86 ppc
	
	CONFIG -= app_bundle

#	CONFIG += staticlib

#	QMAKE_MAC_SDK = /Developer/SDKs/MacOSX10.5.sdk
#    QMAKE_MAC_SDK = macosx10.13

        KERNEL_NAME   = $$system(uname -s)
        HARDWARE_NAME = $$system(uname -m)
        PLATFORM_NAME = $${KERNEL_NAME}_$${HARDWARE_NAME}

	DIR_SEP = "/"

#        QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.5
#	QMAKE_LFLAGS += -undefined dynamic_lookup 

#	QMAKE_CXXFLAGS -= -mmacosx-version-min=10.5
#	QMAKE_LFLAGS   -= -mmacosx-version-min=10.5
		
#	QMAKE_CXXFLAGS += -mmacosx-version-min=10.6 -read_only_relocs suppress
#	QMAKE_LFLAGS   += -mmacosx-version-min=10.6 -read_only_relocs suppress

#	QMAKE_CFLAGS_X86_64 += -Xarch_x86_64 -mmacosx-version-min=10.13

	QMAKE_CXXFLAGS_DEBUG += -pg
	QMAKE_LFLAGS_DEBUG   += -pg
}
win32 {
	PLATFORM_NAME = "win32"
	DIR_SEP = "\\"

#	QMAKE_CXXFLAGS_RELEASE += -g
#	QMAKE_LFLAGS_RELEASE   += -g

#	QMAKE_CXXFLAGS_DEBUG += -pg -ggdb
#	QMAKE_LFLAGS_DEBUG   += -pg -ggdb

}
#message([$$PLATFORM_NAME])


ORSA_BASE = $$dirname(ORSA_PRI)
#message(ORSA_BASE = $$ORSA_BASE)

SUPPORT_DIR = $$DIR_SEP"support"$$DIR_SEP$$PLATFORM_NAME

ORSA_SUPPORT = $$ORSA_BASE$$SUPPORT_DIR
#message(ORSA_SUPPORT = $$ORSA_SUPPORT)

# 3rd party library configuration, platform by platform (unix,macx,win32)

unix:!macx {
	boinc_include {
		INCLUDEPATH += /home/tricaric/boinc/ /home/tricaric/boinc/api/ /home/tricaric/boinc/lib
	}
	boinc_lib {
		LIBS += /home/tricaric/boinc/api/libboinc_api.a /home/tricaric/boinc/lib/libboinc.a
	}
	gmp_include {
		INCLUDEPATH += 
##		INCLUDEPATH += /usr/local/include
	}
	gmp_lib {
		LIBS += -lgmp -lgmpxx
##		LIBS += /usr/local/lib/libgmp.a /usr/local/lib/libgmpxx.a
##		LIBS += -L/usr/local/lib -lgmp -lgmpxx
	}
	gsl_include {
		INCLUDEPATH += /home/tricaric/gsl-repository/install/include
	}
	gsl_lib {
		LIBS += -lgsl -lgslcblas
##		LIBS += /home/tricaric/gsl-1.14-installed/lib/libgsl.a /home/tricaric/gsl-1.14-installed/lib/libgslcblas.a
	}
	osg_include {
        	INCLUDEPATH += /home/tricaric/OpenSceneGraph/include		
##            INCLUDEPATH += /usr/local/opt/open-scene-graph/include
	}
	osg_lib {
		LIBS += -L/home/tricaric/OpenSceneGraph/lib -L/home/tricaric/OpenSceneGraph/lib/osgPlugins-3.0.1
##        LIBS += -L/usr/local/opt/open-scene-graph/lib -L/usr/local/opt/open-scene-graph/lib/osgPlugins-3.5.6
	}
	osg_src {
		OSG_SRC = /home/tricaric/OpenSceneGraph/src/
	}
	qwt_include {
		INCLUDEPATH += /home/tricaric/qwt-6.0/src
	}
	qwt_lib {
		LIBS += -L/home/tricaric/qwt-6.0/lib/ -lqwt
##		LIBS += -L/home/tricaric/qwt/lib/ -lqwt
	}
	spice_include {
		INCLUDEPATH += /home/tricaric/cspice/include
	}
	spice_lib {
		LIBS += /home/tricaric/cspice/lib/cspice.a /home/tricaric/cspice/lib/csupport.a
	}
	zlib_include {
        	INCLUDEPATH +=
	}
	zlib_lib {
		LIBS +=
	}
}

macx {
	boinc_include {
		INCLUDEPATH +=  /Users/tricaric/boinc/ /Users/tricaric/boinc/api/ /Users/tricaric/boinc/lib
	}
    boinc_lib {
#		LIBS += /Users/tricaric/boinc/mac_build/build/boinc.build/Deployment/libboinc.build/Objects-normal/i386/libboinc.a
#		LIBS += /Users/tricaric/boinc/mac_build/build/boinc.build/Deployment/api_libboinc.build/Objects-normal/i386/libboinc_api.a
	}
	gmp_include {
#	     	INCLUDEPATH += $$ORSA_SUPPORT/gmp/include
        INCLUDEPATH += /usr/local/include/
	}
	gmp_lib {
#		LIBS += $$ORSA_SUPPORT/gmp/lib/libgmpxx.a $$ORSA_SUPPORT/gmp/lib/libgmp.a
LIBS += -L/usr/local/lib/	-lgmp -lgmpxx
	}
	gsl_include {
#		INCLUDEPATH += $$ORSA_SUPPORT/gsl/include
	}
	gsl_lib {
#		LIBS += -L$$ORSA_SUPPORT/gsl/lib -lgsl -lgslcblas
LIBS += -lgsl -lgslcblas
	}
	osg_include {
      	INCLUDEPATH += /Users/tricaric/OpenSceneGraph/include
##        INCLUDEPATH += /usr/local/opt/open-scene-graph/include
	}
	osg_lib {
		LIBS += -L/Users/tricaric/OpenSceneGraph/lib/
#        LIBS += -L/usr/local/opt/open-scene-graph/lib -L/usr/local/opt/open-scene-graph/lib/osgPlugins-3.5.6
	}
	osg_src {
		OSG_SRC = /Users/tricaric/OpenSceneGraph/src/
	}
	qwt_include {
		INCLUDEPATH += /usr/local/qwt-6.1.3/lib/qwt.framework/Headers
	}
	qwt_lib {
#		LIBS += -L/usr/local/qwt-6.1.3/lib/ -lqwt
        LIBS += /usr/local/qwt-6.1.3/lib/qwt.framework/qwt
	}
	spice_include {
#		INCLUDEPATH += $$ORSA_SUPPORT/cspice/include
		INCLUDEPATH += /Users/tricaric/cspice/include
	}
	spice_lib {
#		LIBS += $$ORSA_SUPPORT/cspice/lib/cspice.a $$ORSA_SUPPORT/cspice/lib/csupport.a
		LIBS += /Users/tricaric/cspice/lib/cspice.a /Users/tricaric/cspice/lib/csupport.a
	}
#	zlib_include {
#        	INCLUDEPATH +=
#	}
#	zlib_lib {
#		LIBS +=
#	}
}

win32 {
	boinc_include {
		INCLUDEPATH += C:\boinc\boinc C:\boinc\boinc\api C:\boinc\boinc\lib	
	}
	boinc_lib {
		LIBS += C:\boinc\boinc\api\libboinc.a
	}
	gmp_include {
		INCLUDEPATH += $$ORSA_SUPPORT\gmp\include
	}
	gmp_lib {
		LIBS += $$ORSA_SUPPORT\gmp\lib\libgmpxx.a $$ORSA_SUPPORT\gmp\lib\libgmp.a 
	}
	gsl_include {
		INCLUDEPATH += $$ORSA_SUPPORT\gsl\include
	}
	gsl_lib {
		LIBS += -L$$ORSA_SUPPORT\gsl\lib -lgsl -lgslcblas
	}
	osg_include {
        INCLUDEPATH += C:\OpenSceneGraph\include
        QMAKE_CXXFLAGS += -DOSG_LIBRARY
	}
	osg_lib {
		LIBS += C:\OpenSceneGraph\lib\libOpenThreads.a C:\OpenSceneGraph\lib\libosg.a
	}
	osg_src {
		OSG_SRC = C:\OpenSceneGraph\src\
	}
	qwt_include {
		INCLUDEPATH += C:\qwt\src
	}
	qwt_lib {
		LIBS += -LC:\qwt\lib\ -lqwt5
	}
	spice_include {
		INCLUDEPATH += $$ORSA_SUPPORT\cspice\include
	}
	spice_lib {
		LIBS += $$ORSA_SUPPORT\cspice\lib\cspice.a $$ORSA_SUPPORT\cspice\lib\csupport.a
	}
	zlib_include {
        	INCLUDEPATH += C:\Qt\2010.02.1\qt\src\3rdparty\zlib
	}
	zlib_lib {
		LIBS +=
	}
}

