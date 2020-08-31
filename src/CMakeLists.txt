#    UnBlocks-Gen: 3D rock mass generator and analyser
#    Copyright (C) 2020  Leandro Lima Rasmussen

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.6)
project(unblocks LANGUAGES CXX)

set(SOURCES
    Wrapper.cpp
    maths/Functions.cpp
    maths/Geometry.cpp
    
    dfn/DFN.cpp
    dfn/Fracture.cpp
    dfn/FractureSet.cpp
    dfn/Mapping.cpp
   
    generator/Block.cpp
    generator/Generator.cpp
    )

SET(CMAKE_CXX_FLAGS "-D DEBUG -O3 -std=c++11")
ADD_LIBRARY(${PROJECT_NAME} SHARED ${SOURCES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${PYTHON_LIBRARIES} Boost::python3 -lClpSolver -lClp -lCoinUtils -lboost_python)
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES SUFFIX .so)
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES PREFIX "")

FIND_PACKAGE(PythonLibs 3 REQUIRED)
FIND_PACKAGE(Boost REQUIRED COMPONENTS python3)
IF (Boost_FOUND AND PythonLibs_FOUND)
    INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}" ${PYTHON_INCLUDE_DIRS} ${PYTHON_LIBRARIES} ${CMAKE_SOURCE_DIR})	
ELSEIF(NOT Boost_FOUND OR NOT PythonLibs_FOUND)
    MESSAGE(FATAL_ERROR "Unable to find correct Boost library or Python library.")
ENDIF()
