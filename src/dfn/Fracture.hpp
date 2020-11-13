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

#ifndef FRACTURE_H
#define FRACTURE_H

#include "maths/Math.hpp"

class Fracture{
public:
    Fracture(int _id, Vector3r _center, Vector3r _unitVector, std::vector<Vector3r>& _borderPoints);
    ~Fracture(){};
    
    double get_Area();
    void export_FractureVtk(std::string _fileName);
    
    std::vector<Triangle> triangles;
    std::vector<Triangle> trimmedTriangles;
    std::vector<Vector3r> borderPoints;
    std::vector<Plane> borderPlanes;
    
    double boundingSphereRadius = 0;
    Vector3r unitVector;
    Vector3r center;
    int id;
};

#endif //FRACTURE_H
