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

#include <fstream>
#include <iomanip>
#include "dfn/Mapping.hpp"
#include "dfn/DFN.hpp"
#include <boost/python.hpp>

using namespace boost::python;
void Mapping_Wrapper(){
    class_<LineMapping, boost::noncopyable>("LineMapping", no_init)
        .def("export_LineMappingVtk", &LineMapping::export_LineMappingVtk)
        .def("get_P10", &LineMapping::get_P10)
    ;
    
    class_<SurfaceMapping, boost::noncopyable>("SurfaceMapping", no_init)
        .def("export_SurfaceMappingVtk", &SurfaceMapping::export_SurfaceMappingVtk)
        .def("get_P20", &SurfaceMapping::get_P20)
        .def("get_P21", &SurfaceMapping::get_P21)
    ;
    
    class_<VolumeMapping, boost::noncopyable>("VolumeMapping", no_init)
        .def("export_VolumeMappingVtk", &VolumeMapping::export_VolumeMappingVtk)
        .def("get_P30", &VolumeMapping::get_P30)
        .def("get_P32", &VolumeMapping::get_P32)
    ;
};

double LineMapping::get_P10(int _fracSetId){
	Line line(pointA, pointB);
	int nFractures = 0;
	for (auto& frac : dfn->fractureSets[_fracSetId]->fractures){
		for (auto& tri : frac->triangles){
            Vector3r auxP = Vector3r::Zero();
			if (Functions::lineTriangleIntersection(line,tri,auxP)) {
                nFractures += 1;
                break;
            }
		}
	}
	return (double)nFractures/line.length();
}

void LineMapping::export_LineMappingVtk(std::string _fileName){
	std::ofstream out;
	out.open(_fileName + ".vtk");
	out << std::setprecision(15);
	
 	out << "# vtk DataFile Version 3.0" << std::endl;
	out << "Discrete fracture network" << std::endl;
	out << "ASCII" << std::endl;
	
	out << " " << std::endl;
	out << "DATASET POLYDATA" << std::endl;
	out << "POINTS " << 2 << " double" << std::endl;
	out << pointA[0] << " " << pointA[1] << " " << pointA[2] << std::endl;
    out << pointB[0] << " " << pointB[1] << " " << pointB[2] << std::endl;
	
	out << " " << std::endl;
	out << "LINES " << 1 << " " << 3 << std::endl;
	out << "2 0 1" << std::endl;
    out.close();
    
	std::cout << "Vtk Line Mapping Vtk exported!" << std::endl;    
}

void VolumeMapping::export_VolumeMappingVtk(std::string _fileName){
	dfn->export_RegionVtk(_fileName);
}

double SurfaceMapping::get_P21(int _fracSetId){
	double totalArea = 0; 
    double totalLineLength = 0;
	for (auto& surfTri : surfaceTriangles) totalArea += surfTri.area();
	for (auto& frac : dfn->fractureSets[_fracSetId]->fractures){
		for (auto & fracTri : frac->triangles){
			for (auto& surfTri : surfaceTriangles) {
				Line lineFromIntersect;
                if (Functions::triangleTriangleIntersection(fracTri,surfTri,lineFromIntersect)){
                    totalLineLength += lineFromIntersect.length();
                }
			}
		}
	}
	return totalLineLength/totalArea;
}

double SurfaceMapping::get_P20(int _fracSetId){
	double totalArea = 0; 
    int nFractures = 0;
	for (auto& surfTri : surfaceTriangles) totalArea += surfTri.area();
	for (auto& frac : dfn->fractureSets[_fracSetId]->fractures){
        bool intersecCheck = false;
		for (auto & fracTri : frac->triangles){
			for (auto& surfTri : surfaceTriangles) {
				Line lineFromIntersect;
				if (Functions::triangleTriangleIntersection(fracTri,surfTri,lineFromIntersect)){
                    nFractures++; 
                    intersecCheck = true;
                }
                if (intersecCheck) break;
			}
			if (intersecCheck) break;
		}	
	}
	return (double)nFractures/totalArea;
}

void SurfaceMapping::export_SurfaceMappingVtk(std::string _fileName){
	std::ofstream out;
	out.open(_fileName + ".vtk");
	out << std::setprecision(15);
	
 	out << "# vtk DataFile Version 3.0" << std::endl;
	out << "Discrete fracture network" << std::endl;
	out << "ASCII" << std::endl;
 	int count = 0;
 	int auxPolycount = 0;
	                 
    for (auto & Tri : surfaceTriangles){
        count+=3;
    }
    
	out << " " << std::endl;
	out << "DATASET POLYDATA" << std::endl;
	out << "POINTS " << count << " double" << std::endl;
    for (auto & Tri : surfaceTriangles){
        out << Tri.pointD[0] << " " << Tri.pointD[1] << " " << Tri.pointD[2] << std::endl;
        out << Tri.pointE[0] << " " << Tri.pointE[1] << " " << Tri.pointE[2] << std::endl;
        out << Tri.pointF[0] << " " << Tri.pointF[1] << " " << Tri.pointF[2] << std::endl; 
    }
	
	out << " " << std::endl;
	out << "POLYGONS " << count/3 << " " << 4*count/3 << std::endl;
    for (auto & Tri : surfaceTriangles){
        out << "3 " << auxPolycount << " " << auxPolycount+1 << " " << auxPolycount+2 << std::endl;
        auxPolycount += 3;
    }
    
	std::cout << "Surface Mapping Vtk exported!" << std::endl;
}

 double VolumeMapping::get_P32(int _fracSetId){
	double totalFracArea = 0;
	for (auto& frac : dfn->fractureSets[_fracSetId]->fractures){
		totalFracArea += frac->get_Area();
	}
	return totalFracArea/dfn->regionVolume; 
}

 double VolumeMapping::get_P30(int _fracSetId){
	int nFractures = (int)dfn->fractureSets[_fracSetId]->fractures.size();
	return (double)nFractures/dfn->regionVolume; 
 }
