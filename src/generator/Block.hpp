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

#ifndef BLOCK_H
#define BLOCK_H

#include "maths/Math.hpp"

class Block{
friend class Generator;
public:
    Block(){};
    Block(std::vector<Plane>& _planes, int _id):
        id(_id),
        planes(_planes)
        {
            calculate_BoundingSphere();
            calculate_InscribedSphere();
        };
    ~Block(){};
        
    std::vector<Plane> planes;
    
    std::vector<Polygon> polygons;
    std::vector<Edge> edges;
    std::vector<Vector3r> vertices;
       
    void export_BlockVtk(std::string _fileName);
    
    double get_Volume(){return volume;};
    double get_Alpha();
    double get_Beta();
    double get_InscribedSphereRadius(){return inscribedSphereRadius;};
    double get_AspectRatio(){return boundingSphereRadius/inscribedSphereRadius;};
    int get_Order(){return (int)polygons.size();};
    
    double volume;
    double alpha;
    double beta;
    int id;
        
private:
    void generate_Geometry(); 
    void calculate_Volume();
    void calculate_AlphaBeta();
    void calculate_BoundingSphere();
    void calculate_InscribedSphere();
    
    Vector3r inscribedSphereCenter;
    Vector3r boundingSphereCenter;
    double boundingSphereRadius = 0;
    double inscribedSphereRadius = 0;
};

#endif //BLOCK_H
