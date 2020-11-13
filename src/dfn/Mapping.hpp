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

#ifndef MAPPING_H
#define MAPPING_H

#include "maths/Math.hpp"

class DFN;
  
class LineMapping{ 
	public:
        LineMapping(Vector3r _pointA, Vector3r _pointB, DFN* _dfn):
            pointA(_pointA),
            pointB(_pointB),
            dfn(_dfn)
            {};
		~LineMapping(){};
        
		Vector3r pointA = Vector3r::Zero();
		Vector3r pointB = Vector3r::Zero();
        
        void export_LineMappingVtk(std::string _fileName);
        double get_P10(int _fracSetId);
        
    private:
        DFN* dfn;
};

class SurfaceMapping{
	public:
        SurfaceMapping(DFN* _dfn):
            dfn(_dfn)
            {};
		~SurfaceMapping(){};
        
		std::vector<Triangle> surfaceTriangles;
        
        void export_SurfaceMappingVtk(std::string _fileName);
        double get_P21(int _fracSetId);
		double get_P20(int _fracSetId);
        
    private:
        DFN* dfn;
};

class VolumeMapping{
	public:
        VolumeMapping(DFN* _dfn):
            dfn(_dfn)
            {};
		~VolumeMapping(){};
        
        void export_VolumeMappingVtk(std::string _fileName);
        double get_P30(int _fracSetId);
		double get_P32(int _fracSetId);
        
    private:
        DFN* dfn;
};

#endif //MAPPING_H
