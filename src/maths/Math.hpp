/*    
    UnBlocks-Gen: 3D rock mass generator and analyser
    Copyright (C) 2020  Leandro Lima Rasmussen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef MATH_H
#define MATH_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <map>
#include <boost/python.hpp>

#include "Assert.hpp"
#include "EigenTypes.hpp"
#include "TolControl.hpp"
#include "Geometry.hpp"
#include "Functions.hpp"

inline double randomize(double _value) {return _value*double((rand() % 1000000000))/1000000000;};

typedef boost::python::list PyList;
inline Vector3r pyListToVec3(boost::python::list _point){
    ASSERT(boost::python::len(_point) == 3);
    return {boost::python::extract<double>(_point[0]),boost::python::extract<double>(_point[1]), boost::python::extract<double>(_point[2])};
}
#endif //MATH_H
