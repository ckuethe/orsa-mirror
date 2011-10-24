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
SET(ORSA_CURRENT_EXECUTABLE ${ORSA_CURRENT_NAME})
# ---------------------------------------------------------------------------- #



# ---------------------------------- HEADERS --------------------------------- #
SET(ORSA_CURRENT_HEADERS
#BuildCounter.h
mainThread.h
multiminPhase.h
plot.h
RendezvousWithVesta.h
RendezvousWithVestaVersion.h
vesta.h
vestaViz.h
)
# ---------------------------------------------------------------------------- #



# ---------------------------------- SOURCES --------------------------------- #
SET(ORSA_CURRENT_SOURCES
#BuildCounter.cpp
main.cpp
mainThread.cpp
#mp.cpp
multiminPhase.cpp
plot.cpp
RendezvousWithVesta.cpp
RendezvousWithVestaVersion.cpp
)
# ---------------------------------------------------------------------------- #



# ------------------------------ LINK LIBRARIES ------------------------------ #
SET(ORSA_CURRENT_LINK_INCLUDE
${OPENSCENEGRAPH_INCLUDE_DIRS}
${GMP_INCLUDE_DIRS}
${GSL_INCLUDE_DIRS}
${QWT_INCLUDE_DIRS}
${SPICE_INCLUDE_DIRS}
)

SET(ORSA_CURRENT_LINK_LIBRARIES
${QT_LIBRARIES} 
${OPENSCENEGRAPH_LIBRARIES}
${OSG_LIBRARY_osgdbfreetype}
${GMP_LIBRARIES}
${GSL_LIBRARIES}
${QWT_LIBRARIES}
${SPICE_LIBRARIES}
${SYSTEM_LIBRARIES}
${ORSA_LIBRARY_orsa}
${ORSA_LIBRARY_orsaOSG}
${ORSA_LIBRARY_orsaSolarSystem}
${ORSA_LIBRARY_orsaSPICE}
)
# ---------------------------------------------------------------------------- #



# ----------------------------------- BUILD ---------------------------------- #
QT4_WRAP_CPP(ORSA_CURRENT_MOC ${ORSA_CURRENT_HEADERS})

ADD_EXECUTABLE(${ORSA_CURRENT_EXECUTABLE} ${ORSA_CURRENT_SOURCES} ${ORSA_CURRENT_MOC})

TARGET_LINK_LIBRARIES(${ORSA_CURRENT_EXECUTABLE} ${ORSA_CURRENT_LINK_LIBRARIES})
# ---------------------------------------------------------------------------- #
