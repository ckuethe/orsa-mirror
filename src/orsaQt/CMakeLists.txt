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



# ----------------------------------- NAMES ----------------------------------- #
GET_FILENAME_COMPONENT(ORSA_CURRENT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
SET(ORSA_CURRENT_LIBSTATIC ${ORSA_CURRENT_NAME}${LIBSTATIC_SUFFIX})
SET(ORSA_CURRENT_LIBSHARED ${ORSA_CURRENT_NAME}${LIBSHARED_SUFFIX})
# ---------------------------------------------------------------------------- #



# ---------------------------------- HEADERS --------------------------------- #
SET(ORSA_CURRENT_HEADERS
#about.h
#advancedpushbutton.h
#basicprojectwidget.h
#contentlistwidget.h
#datafilewidget.h
#displaymatrix.h
#displaymatrixelement.h
debug.h
#directories.h
event.h
#mainwindow.h
#objectwidget.h
#osg3dviewadapter.h
#osg3dviewer.h
#playbarwidget.h
#projectwidgetframe.h
#simulationparameterswidget.h
#simulationviewerwidget.h
)
# ---------------------------------------------------------------------------- #



# ---------------------------------- SOURCES --------------------------------- #
SET(ORSA_CURRENT_SOURCES
#about.cpp
#advancedpushbutton.cpp
#basicprojectwidget.cpp
#contentlistwidget.cpp
#datafilewidget.cpp
#displaymatrix.cpp
#displaymatrixelement.cpp
debug.cpp
#directories.cpp
#mainwindow.cpp
#objectwidget.cpp
#osg3dviewadapter.cpp
#osg3dviewer.cpp
#playbarwidget.cpp
#projectwidgetframe.cpp
#simulationparameterswidget.cpp
#simulationviewerwidget.cpp
)
# ---------------------------------------------------------------------------- #



# ------------------------------ LINK LIBRARIES ------------------------------ #
SET(ORSA_CURRENT_LINK_INCLUDE
${OPENSCENEGRAPH_INCLUDE_DIRS}
${GMP_INCLUDE_DIRS}
${GSL_INCLUDE_DIRS}
)

SET(ORSA_CURRENT_LINK_LIBRARIES
${QT_LIBRARIES} 
${OPENSCENEGRAPH_LIBRARIES}
${GMP_LIBRARIES}
${GSL_LIBRARIES}
)
# ---------------------------------------------------------------------------- #



# ----------------------------------- BUILD ---------------------------------- #
QT4_WRAP_CPP(ORSA_CURRENT_MOC ${ORSA_CURRENT_HEADERS})

ADD_LIBRARY(${ORSA_CURRENT_LIBSTATIC} STATIC ${ORSA_CURRENT_SOURCES} ${ORSA_CURRENT_MOC})
SET_TARGET_PROPERTIES(${ORSA_CURRENT_LIBSTATIC} PROPERTIES CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME ${ORSA_CURRENT_NAME})
ADD_LIBRARY(${ORSA_CURRENT_LIBSHARED} SHARED ${ORSA_CURRENT_SOURCES} ${ORSA_CURRENT_MOC})
SET_TARGET_PROPERTIES(${ORSA_CURRENT_LIBSHARED} PROPERTIES CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME ${ORSA_CURRENT_NAME})

TARGET_LINK_LIBRARIES(${ORSA_CURRENT_LIBSTATIC} ${ORSA_CURRENT_LINK_LIBRARIES})
TARGET_LINK_LIBRARIES(${ORSA_CURRENT_LIBSHARED} ${ORSA_CURRENT_LINK_LIBRARIES})
# ---------------------------------------------------------------------------- #
