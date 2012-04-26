# ---------------------------------------------------------------------------- #
# PROJECT :         ORSA : Orbit Reconstruction, Simulation and Analysis
# AUTHOR :          Pasquale Tricarico <tricaric@psi.edu>
# CONTRIBUTORS :    Pasquale Tricarico <tricaric@psi.edu>
#                   Vincent Reverdy (ESA-SOCIS 2011) <vince.rev@gmail.com>
# ---------------------------------------------------------------------------- #
# ORSA - Orbit Reconstruction, Simulation and Analysis
# Copyright (C) 2002-2010 Pasquale Tricarico <tricaric@psi.edu>
# Portions copyright (C) by their respective authors.
#   
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#   
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#   
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, 51 Franklin Street, Suite 500, Boston, MA 02110-1335, USA.
# ---------------------------------------------------------------------------- #



# ----------------------------- CMAKE PARAMETERS ----------------------------- #
# VERSION
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.3)

# INCLUDE
INCLUDE(FindPackageHandleStandardArgs)

# USER SEARCH PATHS
SET(USER_INCLUDE_PATH /home/tricaric/OpenSceneGraph/src /home/tricaric/qwt/src /home/tricaric/cspice/lib)
SET(USER_LIBRARY_PATH /0/Qwt/qwt-5.2.2/lib/ /0/Cspice/cspice/lib/ /usr/local/lib/osgPlugins-3.0.1/)

# SEARCH DIRECTORIES
SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /usr/include/ ${USER_INCLUDE_PATH})
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /lib/i386-linux-gnu/ /usr/lib/i386-linux-gnu/ /usr/lib/ /usr/lib64/ /lib/ ${USER_LIBRARY_PATH})

# FLAGS
#SET(CMAKE_CXX_FLAGS "-std=c++0x -O3")
#SET(CMAKE_CXX_FLAGS "-std=c++0x -O2")
#SET(CMAKE_CXX_FLAGS "-std=c++0x")
#SET(CMAKE_CXX_FLAGS "-std=c++0x -O3 -W -Wall -fopenmp -g")
SET(CMAKE_CXX_FLAGS "-std=c++0x -O3 -W -Wall -g")
# ---------------------------------------------------------------------------- #



# ---------------------------- PROJECT PARAMETERS ---------------------------- #
# NAMES
PROJECT(ORSA)
SET(CMAKE_PROJECT_NAME ORSA)
SET(PLATFORM_NAME Linux_i686)
SET(LIBSTATIC_SUFFIX -static)
SET(LIBSHARED_SUFFIX -shared)

# ORSA ROOT DIRECTORY
SET(ORSA_DIR ${CMAKE_SOURCE_DIR})

# MAIN DIRECTORIES
SET(EXECUTABLE_OUTPUT_PATH ${ORSA_DIR}/bin/${PLATFORM_NAME})
SET(LIBRARY_OUTPUT_PATH ${ORSA_DIR}/lib/${PLATFORM_NAME})
SET(ORSA_DIR_SRC ${ORSA_DIR}/src)
SET(ORSA_DIR_APP ${ORSA_DIR}/app)
SET(ORSA_DIR_SUPPORT ${ORSA_DIR}/support)

# SEARCH DIRECTORIES
SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH})
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${LIBRARY_OUTPUT_PATH})

# APPS LINKING
SET(ORSA_LIBRARY_PREFIX lib)
SET(ORSA_LIBRARY_SUFFIX .so)
# ---------------------------------------------------------------------------- #


  
# -------------------------------- QT LIBRARY -------------------------------- #
FIND_PACKAGE(Qt4 REQUIRED)
SET(QT_USE_QTOPENGL TRUE)
SET(QT_USE_QTSVG TRUE)
INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})
# ---------------------------------------------------------------------------- # 



# ------------------------------- OSG LIBRARY -------------------------------- #
FIND_PACKAGE(OpenSceneGraph 3.0.1 REQUIRED osgViewer osgText osgUtil osgDB osgGA)
FIND_PATH(OSG_SOURCE_DIR NAMES /osg/ApplicationUsage.cpp /osg/DeleteHandler.cpp /osg/Matrixd.cpp /osg/Matrixf.cpp
          DOC "OpenSceneGraph src dir (where /osg/ApplicationUsage/cpp) is")
FIND_LIBRARY(OSG_LIBRARY_osgdbfreetype NAMES osgdb_freetype.a osgdb_freetype osgdb_freetype.so DOC "OSG plugin library osgdb_freetype (eg : /usr/local/lib/osgPlugins-3.0.1/osgdb_freetype.so)")
# ---------------------------------------------------------------------------- # 



# ------------------------------ BOINC LIBRARY ------------------------------- #
# INCLUDE DIRECTORIES
FIND_PATH(BOINC_INCLUDE_DIR NAMES boinc_api.h
          PATH_SUFFIXES boinc
          DOC "Boinc include directory (where boinc_api.h is)")
SET(BOINC_INCLUDE_DIR_api ${BOINC_INCLUDE_DIR}/api)
SET(BOINC_INCLUDE_DIR_lib ${BOINC_INCLUDE_DIR}/lib)

# LIBRARIES
FIND_LIBRARY(BOINC_LIBRARY NAMES libboinc.a boinc DOC "Boinc library (boinc)")
FIND_LIBRARY(BOINC_LIBRARY_api NAMES libboinc_api.a boinc_api DOC "Boinc api library (boinc_api)")

# COMPUTE RESULTS
FIND_PACKAGE_HANDLE_STANDARD_ARGS(BOINC 
                                  DEFAULT_MSG
                                  BOINC_INCLUDE_DIR
                                  BOINC_INCLUDE_DIR_api
                                  BOINC_INCLUDE_DIR_lib
                                  BOINC_LIBRARY
                                  BOINC_LIBRARY_api)

IF(BOINC_FOUND)
    SET(BOINC_INCLUDE_DIRS ${BOINC_INCLUDE_DIR} ${BOINC_INCLUDE_DIR_api} ${BOINC_INCLUDE_DIR_lib})
    SET(BOINC_LIBRARIES ${BOINC_LIBRARY} ${BOINC_LIBRARY_api})
ENDIF()

# MARK_AS_ADVANCED(BOINC_INCLUDE_DIR 
#                  BOINC_INCLUDE_DIR_api 
#                  BOINC_INCLUDE_DIR_lib 
#                  BOINC_LIBRARY 
#                  BOINC_LIBRARY_api)
# ---------------------------------------------------------------------------- # 



# ------------------------------- GMP LIBRARY -------------------------------- #
# INCLUDE DIRECTORIES
FIND_PATH(GMP_INCLUDE_DIR NAMES gmp.h gmpxx.h
          DOC "GMP include directory (optional)")

# LIBRARIES
FIND_LIBRARY(GMP_LIBRARY NAMES libgmp.a gmp DOC "GMP C library (gmp)")
FIND_LIBRARY(GMP_LIBRARY_gmpxx NAMES libgmpxx.a gmpxx DOC "GMP C++ library (gmpxx)")

# COMPUTE RESULTS
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMP
                                  DEFAULT_MSG
                                  GMP_INCLUDE_DIR
                                  GMP_LIBRARY
                                  GMP_LIBRARY_gmpxx)

IF(GMP_FOUND)
    SET(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
    SET(GMP_LIBRARIES ${GMP_LIBRARY} ${GMP_LIBRARY_gmpxx})
ENDIF()

# MARK_AS_ADVANCED(GMP_INCLUDE_DIR 
#                  GMP_LIBRARY 
#                  GMP_LIBRARY_gmpxx)
# ---------------------------------------------------------------------------- # 



# ------------------------------- GSL LIBRARY -------------------------------- #
# INCLUDE DIRECTORIES
FIND_PATH(GSL_INCLUDE_DIR NAMES gsl_version.h gsl_blas.h
          PATH_SUFFIXES gsl
          DOC "GSL include directory (optional)")

# LIBRARIES
FIND_LIBRARY(GSL_LIBRARY NAMES libgsl.a gsl DOC "GSL library (gsl)")
FIND_LIBRARY(GSL_LIBRARY_cblas NAMES libgslcblas.a gslcblas DOC "GSL CBLAS library (gslcblas)")

# COMPUTE RESULTS
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GSL
                                  DEFAULT_MSG
                                  GSL_INCLUDE_DIR
                                  GSL_LIBRARY
                                  GSL_LIBRARY_cblas)

IF(GSL_FOUND)
    SET(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR})
    SET(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_LIBRARY_cblas})
ENDIF()

# MARK_AS_ADVANCED(GSL_INCLUDE_DIR 
#                  GSL_LIBRARY 
#                  GSL_LIBRARY_cblas)
# ---------------------------------------------------------------------------- # 



# ------------------------------- QWT LIBRARY -------------------------------- #
# INCLUDE DIRECTORIES
FIND_PATH(QWT_INCLUDE_DIR NAMES qwt.h
          PATH_SUFFIXES qwt/src/
          DOC "QWT include/src directory (optional)")

# LIBRARIES
FIND_LIBRARY(QWT_LIBRARY NAMES libqwt.a qwt DOC "QWT library (qwt)")

# COMPUTE RESULTS
FIND_PACKAGE_HANDLE_STANDARD_ARGS(QWT
                                  DEFAULT_MSG
                                  QWT_INCLUDE_DIR
                                  QWT_LIBRARY)

IF(QWT_FOUND)
    SET(QWT_INCLUDE_DIRS ${QWT_INCLUDE_DIR})
    SET(QWT_LIBRARIES ${QWT_LIBRARY})
ENDIF()

# MARK_AS_ADVANCED(QWT_INCLUDE_DIR 
#                  QWT_LIBRARY)
# ---------------------------------------------------------------------------- # 



# ------------------------------ SPICE LIBRARY ------------------------------- #
# INCLUDE DIRECTORIES
FIND_PATH(SPICE_INCLUDE_DIR NAMES SpiceCel.h SpiceCK.h
          PATH_SUFFIXES cspice/include/
          DOC "SPICE include directory (where SpiceCel.h is)")

# LIBRARIES
FIND_LIBRARY(SPICE_LIBRARY NAMES cspice.a cspice DOC "SPICE library (cspice)")
FIND_LIBRARY(SPICE_LIBRARY_csupport NAMES csupport.a csupport DOC "SPICE support library (csupport)")

# COMPUTE RESULTS
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SPICE
                                  DEFAULT_MSG
                                  SPICE_INCLUDE_DIR
                                  SPICE_LIBRARY
                                  SPICE_LIBRARY_csupport)

IF(SPICE_FOUND)
    SET(SPICE_INCLUDE_DIRS ${SPICE_INCLUDE_DIR})
    SET(SPICE_LIBRARIES ${SPICE_LIBRARY} ${SPICE_LIBRARY_csupport})
ENDIF()

# MARK_AS_ADVANCED(SPICE_INCLUDE_DIR 
#                  SPICE_LIBRARY 
#                  SPICE_LIBRARY_csupport)
# ---------------------------------------------------------------------------- # 



# ------------------------------- ZLIB LIBRARY ------------------------------- #
# INCLUDE DIRECTORIES
FIND_PATH(ZLIB_INCLUDE_DIR NAMES zlib.h
          PATH_SUFFIXES src 3rdparty zlib
          DOC "ZLIB include directory (normally in [Qt_directory]/src/3rdparty/zlib)")

# LIBRARIES
FIND_LIBRARY(ZLIB_LIBRARY NAMES libz.a z libz libz.so DOC "ZLIB library (z,libz)")

# COMPUTE RESULTS
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ZLIB
                                  DEFAULT_MSG
                                  ZLIB_INCLUDE_DIR
                                  ZLIB_LIBRARY)

IF(ZLIB_FOUND)
    SET(ZLIB_INCLUDE_DIRS ${ZLIB_INCLUDE_DIR})
    SET(ZLIB_LIBRARIES ${ZLIB_LIBRARY})
ENDIF()

# MARK_AS_ADVANCED(ZLIB_INCLUDE_DIR 
#                  ZLIB_LIBRARY)
# ---------------------------------------------------------------------------- # 



# ----------------------------- SYSTEM LIBRARIES ----------------------------- #
FIND_LIBRARY(GLIB_LIBRARY NAMES glib-2.0 DOC "GLIB library (libglib-2.0.so)")
FIND_LIBRARY(GOBJECT_LIBRARY NAMES gobject-2.0 DOC "GOBJECT library (libgobject-2.0.so)")

SET(SYSTEM_LIBRARIES ${GLIB_LIBRARY} ${GOBJECT_LIBRARY})
# ---------------------------------------------------------------------------- # 



# ------------------------------- ALL LIBRARIES ------------------------------ #
# INCLUDES (QT not needed because of QT_USE_FILE)
SET(ALL_EXTERNAL_INCLUDES
                           ${OPENSCENEGRAPH_INCLUDE_DIRS}
                           ${BOINC_INCLUDE_DIRS}
		                   ${GMP_INCLUDE_DIRS}
                           ${GSL_INCLUDE_DIRS}
                           ${QWT_INCLUDE_DIRS}
                           ${SPICE_INCLUDE_DIRS}
                           ${ZLIB_INCLUDE_DIRS}
                           )

# LIBRARIES
SET(ALL_EXTERNAL_LIBRARIES 
                           ${QT_LIBRARIES} 
                           ${OPENSCENEGRAPH_LIBRARIES}
                           ${BOINC_LIBRARIES}
		           		   ${GMP_LIBRARIES}
                           ${GSL_LIBRARIES}
                           ${QWT_LIBRARIES}
                           ${SPICE_LIBRARIES}
                           ${ZLIB_LIBRARIES}
                           )
# ---------------------------------------------------------------------------- # 



# --------------------------------- INCLUDES --------------------------------- #
INCLUDE_DIRECTORIES(${ORSA_DIR_SRC})
INCLUDE_DIRECTORIES(${ALL_EXTERNAL_INCLUDES})
# ---------------------------------------------------------------------------- #



# ------------------------------- SOURCES NAMES ------------------------------ #
SET(ORSA_NAME_orsa orsa)
SET(ORSA_NAME_orsaEssentialOSG orsaEssentialOSG)
SET(ORSA_NAME_orsaInputOutput orsaInputOutput)
SET(ORSA_NAME_orsaOSG orsaOSG)
SET(ORSA_NAME_orsaPDS orsaPDS)
SET(ORSA_NAME_orsaQt orsaQt)
SET(ORSA_NAME_orsaSolarSystem orsaSolarSystem)
SET(ORSA_NAME_orsaSPICE orsaSPICE)
SET(ORSA_NAME_orsaUtil orsaUtil)
# ---------------------------------------------------------------------------- #



# ------------------------------- BUILD SOURCES ------------------------------ #
ADD_SUBDIRECTORY(${ORSA_DIR_SRC}) 
# ---------------------------------------------------------------------------- #



# ------------------------------- LIBRARY NAMES ------------------------------ #
SET(ORSA_LIBRARY_orsa ${LIBRARY_OUTPUT_PATH}/${ORSA_LIBRARY_PREFIX}${ORSA_NAME_orsa}${ORSA_LIBRARY_SUFFIX})
SET(ORSA_LIBRARY_orsaEssentialOSG ${LIBRARY_OUTPUT_PATH}/${ORSA_LIBRARY_PREFIX}${ORSA_NAME_orsaEssentialOSG}${ORSA_LIBRARY_SUFFIX})
SET(ORSA_LIBRARY_orsaInputOutput ${LIBRARY_OUTPUT_PATH}/${ORSA_LIBRARY_PREFIX}${ORSA_NAME_orsaInputOutput}${ORSA_LIBRARY_SUFFIX})
SET(ORSA_LIBRARY_orsaOSG ${LIBRARY_OUTPUT_PATH}/${ORSA_LIBRARY_PREFIX}${ORSA_NAME_orsaOSG}${ORSA_LIBRARY_SUFFIX})
SET(ORSA_LIBRARY_orsaPDS ${LIBRARY_OUTPUT_PATH}/${ORSA_LIBRARY_PREFIX}${ORSA_NAME_orsaPDS}${ORSA_LIBRARY_SUFFIX})
SET(ORSA_LIBRARY_orsaQt ${LIBRARY_OUTPUT_PATH}/${ORSA_LIBRARY_PREFIX}${ORSA_NAME_orsaQt}${ORSA_LIBRARY_SUFFIX})
SET(ORSA_LIBRARY_orsaSolarSystem ${LIBRARY_OUTPUT_PATH}/${ORSA_LIBRARY_PREFIX}${ORSA_NAME_orsaSolarSystem}${ORSA_LIBRARY_SUFFIX})
SET(ORSA_LIBRARY_orsaSPICE ${LIBRARY_OUTPUT_PATH}/${ORSA_LIBRARY_PREFIX}${ORSA_NAME_orsaSPICE}${ORSA_LIBRARY_SUFFIX})
SET(ORSA_LIBRARY_orsaUtil ${LIBRARY_OUTPUT_PATH}/${ORSA_LIBRARY_PREFIX}${ORSA_NAME_orsaUtil}${ORSA_LIBRARY_SUFFIX})
# ---------------------------------------------------------------------------- #



# -------------------------------- APPS NAMES -------------------------------- #
SET(ORSA_NAME_RendezvousWithVesta RendezvousWithVesta)
SET(ORSA_NAME_RendezvousWithVesta2 RendezvousWithVesta2)
SET(ORSA_NAME_GuiTest GuiTest)
SET(ORSA_NAME_CometGrainsDynamics CometGrainsDynamics)
SET(ORSA_NAME_dawnorbit dawnorbit)
# ---------------------------------------------------------------------------- #



# -------------------------------- BUILD APPS -------------------------------- #
ADD_SUBDIRECTORY(${ORSA_DIR_APP}) 
# ---------------------------------------------------------------------------- #
