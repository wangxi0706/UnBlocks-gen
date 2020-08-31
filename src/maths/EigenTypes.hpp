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

#ifndef EIGENTYPES_H
#define EIGENTYPES_H

#include <Eigen/Core.c>
#include <Eigen/Geometry.c>
#include <Eigen/Eigenvalues.c>

using Quaternionr = Eigen::Quaternion<double>;
using AngleAxisr  = Eigen::AngleAxis<double>;

template<typename Scalar> using Vector3 = Eigen::Matrix<Scalar,3,1>;
using Vector3r = Vector3<double>;

template<typename Scalar> using Matrix3 = Eigen::Matrix<Scalar,3,3>;
using Matrix3r = Matrix3<double>;

using Quaternionr = Eigen::Quaternion<double>;
using AngleAxisr  = Eigen::AngleAxis<double>;

#endif //EIGENTYPES_H
