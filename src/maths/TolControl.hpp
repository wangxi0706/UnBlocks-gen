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

#ifndef TOLCONTROL_H
#define TOLCONTROL_H

static double ND_tolerance = 100000*std::numeric_limits<double>::epsilon();
inline bool checkEquality(const double& _x, const double& _y) {return ((std::abs(_x - _y) <= ND_tolerance * std::max(1.0, std::max(std::abs(_x), std::abs(_y)))) ? true : false);}
inline bool checkEquality(const double& _x, const double& _y, double _tolerance) {return ((std::abs(_x - _y) <= _tolerance * std::max(1.0, std::max(std::abs(_x), std::abs(_y)))) ? true : false);}
inline double cZ(const double& _x) {return ((std::abs(_x) <= ND_tolerance * std::max(1.0, std::abs(_x))) ? 0 : _x);}

#endif //TOLCONTROL_H
