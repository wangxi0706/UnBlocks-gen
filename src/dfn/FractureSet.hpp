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

#ifndef FRACTURESET_H
#define FRACTURESET_H

#include "maths/Math.hpp"

class Fracture;
class DFN;

class FractureSet{ 
friend class DFN;
public:
    FractureSet(DFN* _dfn, int _id):
        dfn(_dfn),
        id(_id)
        {};
                
    ~FractureSet(){};
    
    void export_FractureSetVtk(std::string _fileName);
		
    void add_BaecherFracture(double _meanDipDirection, double _meanDipAngle, double _fisherConstant, std::string _sizeDistribution, double _meanFractureSize, double _sigmaFractureSize);
    void add_TriangularFracture(PyList _pointA, PyList _pointB, PyList _pointC);
    void add_CircularFracture(PyList _center, double _dipDirection, double _dipAngle, double _radius);
       
    std::vector<std::shared_ptr<Fracture>> fractures;
    int id;

private:
    void add_CircularFractureByUnitVec(Vector3r _center, Vector3r _unitVector, double _radius);
    void trim_FractureBorders(Fracture& _fracture);
    DFN* dfn;
};

#endif //FRACTURESET_H


