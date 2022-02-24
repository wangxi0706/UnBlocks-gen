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

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <fstream>
#include <iomanip>
#include "EigenTypes.hpp"
#include "Assert.hpp"
#include "TolControl.hpp"

// Added box region
struct Box
{
    Vector3r minCor;
    Vector3r maxCor;
};

class Plane
{
public:
    Plane(Vector3r _center, Vector3r _unitVector) : center(_center),
                                                    unitVector(_unitVector)
    {
        d = _center.dot(_unitVector);
    };
    Plane(Vector3r _unitVector, double _d) : unitVector(_unitVector),
                                             d(_d)
    {
        center = _d * _unitVector; //?? but the center may not on the unit vector??
    };

    Vector3r unitVector = Vector3r::Zero();
    Vector3r center = Vector3r::Zero();
    // can facilitate the point to plane distance, (a-c)*n=a*n-c*n=a*n-d
    double d;
};

class Polygon
{
public:
    Polygon(Vector3r _unitVector, std::vector<int> &_verticesId, std::vector<Vector3r> &_vertices);
    ~Polygon(){};

    std::vector<int> verticesId;
    Vector3r unitVector;
    Vector3r centroid;
    double area;
};

class Edge
{
public:
    Edge(int _verticeIdA, int _verticeIdB) : verticesIdA(_verticeIdA),
                                             verticesIdB(_verticeIdB){};

    ~Edge(){};

    int verticesIdA;
    int verticesIdB;
};

class Triangle
{
public:
    Triangle(Vector3r _pointD, Vector3r _pointE, Vector3r _pointF) : pointD(_pointD),
                                                                     pointE(_pointE),
                                                                     pointF(_pointF),
                                                                     unitVector(((_pointF - _pointD).cross(_pointF - _pointE)) / ((_pointF - _pointD).cross(_pointF - _pointE)).norm())
    {
        ASSERT(_pointD != _pointE && _pointD != _pointF && _pointE != _pointF);
        ASSERT(checkEquality(unitVector.norm(), 1));

        Vector3r minCorner = {INFINITY, INFINITY, INFINITY};
        Vector3r maxCorner = Vector3r::Zero();
        if (_pointD[0] < minCorner[0])
            minCorner[0] = _pointD[0];
        if (_pointE[0] < minCorner[0])
            minCorner[0] = _pointE[0];
        if (_pointF[0] < minCorner[0])
            minCorner[0] = _pointF[0];
        if (_pointD[1] < minCorner[1])
            minCorner[1] = _pointD[1];
        if (_pointE[1] < minCorner[1])
            minCorner[1] = _pointE[1];
        if (_pointF[1] < minCorner[1])
            minCorner[1] = _pointF[1];
        if (_pointD[2] < minCorner[2])
            minCorner[2] = _pointD[2];
        if (_pointE[2] < minCorner[2])
            minCorner[2] = _pointE[2];
        if (_pointF[2] < minCorner[2])
            minCorner[2] = _pointF[2];

        if (_pointD[0] > maxCorner[0])
            maxCorner[0] = _pointD[0];
        if (_pointE[0] > maxCorner[0])
            maxCorner[0] = _pointE[0];
        if (_pointF[0] > maxCorner[0])
            maxCorner[0] = _pointF[0];
        if (_pointD[1] > maxCorner[1])
            maxCorner[1] = _pointD[1];
        if (_pointE[1] > maxCorner[1])
            maxCorner[1] = _pointE[1];
        if (_pointF[1] > maxCorner[1])
            maxCorner[1] = _pointF[1];
        if (_pointD[2] > maxCorner[2])
            maxCorner[2] = _pointD[2];
        if (_pointE[2] > maxCorner[2])
            maxCorner[2] = _pointE[2];
        if (_pointF[2] > maxCorner[2])
            maxCorner[2] = _pointF[2];
        aabb.first = minCorner;
        aabb.second = maxCorner;
    };

    ~Triangle(){};

    double area() { return (((pointE - pointD).cross(pointF - pointD)).norm()) * 0.5; }

    Vector3r unitVector = Vector3r::Zero();
    Vector3r pointD = Vector3r::Zero();
    Vector3r pointE = Vector3r::Zero();
    Vector3r pointF = Vector3r::Zero();
    std::pair<Vector3r, Vector3r> aabb;

    std::pair<Vector3r, Vector3r> get_Aabb() { return aabb; };

    void export_TriangleVtk(std::string _fileName)
    {
        std::ofstream out;
        out.open(_fileName + ".vtk");
        out << std::setprecision(15);
        out << "# vtk DataFile Version 3.0" << std::endl;
        out << "Voronoi results" << std::endl;
        out << "ASCII" << std::endl;
        out << " " << std::endl;
        out << "DATASET POLYDATA" << std::endl;
        out << "POINTS "
            << "3"
            << " float" << std::endl;
        out << pointD[0] << " " << pointD[1] << " " << pointD[2] << std::endl;
        out << pointE[0] << " " << pointE[1] << " " << pointE[2] << std::endl;
        out << pointF[0] << " " << pointF[1] << " " << pointF[2] << std::endl;
        out << " " << std::endl;
        out << "POLYGONS " << 1 << " "
            << "4" << std::endl;
        out << "3 0 1 2" << std::endl;
        out.close();
    }
};

class Line
{
public:
    Line(){};
    Line(Vector3r _pointA, Vector3r _pointB) : pointA(_pointA),
                                               pointB(_pointB){};

    ~Line(){};

    double length() { return (pointA - pointB).norm(); }

    Vector3r pointA = Vector3r::Zero();
    Vector3r pointB = Vector3r::Zero();
};

#endif // GEOMETRY_H
