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

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "EigenTypes.hpp"
#include "Geometry.hpp"

class Functions{
public:
    static bool triangleTriangleIntersection(const Triangle& _triangle1,  const Triangle& _triangle2, Line& _resultLine);
    static bool pointVolumeIntersection(const Vector3r& _point, const std::vector<Triangle>& _modelBoundary);
    static bool lineTriangleIntersection(const Line& _line, const Triangle& _triangle, Vector3r& _resultPoint);
    static double vectorsAngle(Vector3r _centerPoint, Vector3r _unitVec, Vector3r _point1, Vector3r _point2);
    static bool calculate_threePlanesIntersection(const Plane& _p1, const Plane& _p2, const Plane& _p3, Vector3r& _intersectionPoint);
    static void organize_Vertices(Polygon& _polygon, std::vector<Vector3r>& _vertices);
    static void organize_Points(Vector3r& _center, Vector3r _unitVector, std::vector<Vector3r>& _points);
    static void eliminate_RedundantPlanes(std::vector<Plane>& _planes);
    static bool check_PointIntersection(const Vector3r& _point, const std::vector<Plane>& _planes);
};

#endif //FUNCTIONS_H
