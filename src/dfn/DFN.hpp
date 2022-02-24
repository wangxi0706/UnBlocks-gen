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

#ifndef DFN_H
#define DFN_H

#include "maths/Math.hpp"
#include "dfn/Fracture.hpp"
#include "dfn/FractureSet.hpp"
#include "dfn/Mapping.hpp"

class DFN
{
    friend class FractureSet;
    friend class Generator;
    friend class VolumeMapping;

public:
    DFN();
    ~DFN()
    {
        srand(std::time(0));
        randomSeed = std::time(0);
    }

    void add_FractureSet();

    void add_LineMapping(PyList _pointA, PyList _pointB);
    void add_QuadrilateralMapping(PyList _pointA, PyList _pointB, PyList _pointC, PyList _pointD);
    void add_CircularMapping(PyList _center, double _dipDirection, double _dipAngle, double _radius);
    void add_VolumeMapping();

    void set_RandomSeed(int _seed);
    void set_RegionMaxCorner(PyList _point);
    // added
    void set_RegionMaxMinCorner(PyList _point1, PyList _point2);
    void set_FirstRecBlock(PyList _MinPoint, PyList _MaxPoint);

    void set_NumberOfBorderPoints(int _number) { nFracBorderPoints = _number; };

    void export_DFNVtk(std::string _fileName);
    void export_RegionVtk(std::string _fileName);

    std::vector<std::shared_ptr<FractureSet>> fractureSets;

    std::vector<std::shared_ptr<LineMapping>> linesMapping;
    std::vector<std::shared_ptr<SurfaceMapping>> surfacesMapping;
    std::vector<std::shared_ptr<VolumeMapping>> volumesMapping;

private:
    std::vector<Triangle> modelRegion;
    int randomSeed = 0;
    int nFracBorderPoints = 26;
    Vector3r offset = {0, 0, 0};
    Vector3r regionMinCorner = {0, 0, 0};
    Vector3r regionMaxCorner = {100, 100, 100};

    Vector3r firstBlkMin = {0, 0, 0};
    Vector3r firstBlkMax = {100, 100, 100};

    double regionVolume = 1000000;
};

#endif // DFN_H
