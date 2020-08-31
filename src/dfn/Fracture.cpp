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

#include "dfn/Fracture.hpp"
#include <boost/python.hpp>

using namespace boost::python;
void Fracture_Wrapper(){
    class_<Fracture, boost::noncopyable>("fractures", no_init)
        .def("get_Area", &Fracture::get_Area)
        .def("export_FractureVtk", &Fracture::export_FractureVtk)
    ;
};

Fracture::Fracture(int _id, Vector3r _center, Vector3r _unitVector, std::vector<Vector3r>& _borderPoints):
    id(_id),
    center(_center),
    unitVector(_unitVector),
    borderPoints(_borderPoints)
{
    ASSERT((int)borderPoints.size() >= 3);
    ASSERT(checkEquality(unitVector.norm(),1));
    for (int i = 0; i!=(int)borderPoints.size(); ++i){
        int j = i + 1;
        if (j == (int)borderPoints.size()) j = 0;
        Vector3r borderPlaneUnitVector = (borderPoints[j] - borderPoints[i]).cross(unitVector);
        borderPlaneUnitVector.normalize();
        ASSERT(checkEquality(borderPlaneUnitVector.norm(),1));
        ASSERT((borderPoints[i] - center).dot(borderPlaneUnitVector) > 0);
        Plane newBorderPlane(borderPoints[i], borderPlaneUnitVector);
        borderPlanes.push_back(newBorderPlane);
    }
    for (int i = 0; i!= (int)borderPoints.size(); ++i){
        int j = i+1;
        if (j == (int)borderPoints.size()) j = 0;
        Triangle newTri(center, borderPoints[i], borderPoints[j]);
        triangles.push_back(newTri);
    }
    
    for (int i = 0; i!= (int)borderPoints.size(); ++i){
        if ((borderPoints[i] - center).norm() > boundingSphereRadius) boundingSphereRadius = (borderPoints[i] - center).norm();
    }
}

double Fracture::get_Area(){
	double totalArea = 0;
    if ((int)trimmedTriangles.size() > 0){
        for (auto & Tri : trimmedTriangles){
            totalArea += Tri.area();
        }
    } else {
        for (auto & Tri : triangles){
            totalArea += Tri.area();
        }
    }
	return totalArea;
}

void Fracture::export_FractureVtk(std::string _fileName){
	std::ofstream out;
	out.open(_fileName + ".vtk");
	out << std::setprecision(15);
	
 	out << "# vtk DataFile Version 3.0" << std::endl;
	out << "Discrete fracture network" << std::endl;
	out << "ASCII" << std::endl;
 	int nFracs = 1;
 	int nBorderPoints = (int)borderPoints.size();
    int counter = 0;
		
	out << " " << std::endl;
	out << "DATASET POLYDATA" << std::endl;
	out << "POINTS " << nBorderPoints << " double" << std::endl;
    for (auto & p : borderPoints){
        out << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
		
	out << " " << std::endl;
	out << "POLYGONS " << nFracs << " " << (nFracs+nBorderPoints) << std::endl;
    out << (int)borderPoints.size();
    for (auto & p : borderPoints){
        out << " " << counter;
        counter++;
    }
    out << std::endl;

	out << " " << std::endl;
    out << "CELL_DATA " << nFracs << std::endl;
    out << "SCALARS FracId float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    out << id << std::endl;
    
    out << " " << std::endl;
    out << "SCALARS Area float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    out << get_Area() << std::endl;
        
    out.close();
	std::cout << "Fracture Vtk exported!" << std::endl;
}
