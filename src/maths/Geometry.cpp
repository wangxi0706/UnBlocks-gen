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

#include "Geometry.hpp"
#include "Functions.hpp"

Polygon::Polygon(Vector3r _unitVector, std::vector<int> &_verticesId, std::vector<Vector3r> &_vertices) : unitVector(_unitVector),
                                                                                                          verticesId(_verticesId)
{
    ASSERT((int)verticesId.size() >= 3);
    ASSERT(checkEquality(unitVector.norm(), 1));

    //calculate polygon center from average of vertices coordinates
    Vector3r center = Vector3r::Zero();
    for (auto &id : verticesId)
    {
        center += _vertices[id] / (double)verticesId.size();
    }

    //In this algorithm, the points with bigger angles are added first in a lowering angle sequence
    std::vector<double> angles;
    for (int i = 1; i != (int)verticesId.size(); i++)
    {
        double angle = Functions::vectorsAngle(center, unitVector, _vertices[verticesId[0]], _vertices[verticesId[i]]);
        angles.push_back(angle);
    }

    double previousBiggestAngle = M_PI * 2;
    std::vector<int> organizedVerticesIds;
    organizedVerticesIds.push_back(verticesId[0]);
    for (int i = 1; i != (int)verticesId.size(); i++)
    {
        int biggestAngleId = 0;
        double biggestAngle = 0;
        for (int j = 0; j != (int)angles.size(); j++)
        {
            if (angles[j] > biggestAngle)
            {
                biggestAngle = angles[j];
                biggestAngleId = j;
            }
        }
        ASSERT(!checkEquality(biggestAngle, previousBiggestAngle) && !checkEquality(biggestAngle, 0) && !checkEquality(biggestAngle, M_PI * 2));
        if (biggestAngle < previousBiggestAngle)
        {
            organizedVerticesIds.push_back(verticesId[biggestAngleId + 1]);
            previousBiggestAngle = biggestAngle;
        }
        angles[biggestAngleId] = 0;
    }

    ASSERT((int)verticesId.size() == (int)organizedVerticesIds.size());
    for (int i = 0; i != (int)verticesId.size(); i++)
        verticesId[i] = organizedVerticesIds[i];

    ASSERT((int)verticesId.size() >= 3);
    Vector3r unitVectorX = (center.cross(_unitVector)) / (center.cross(_unitVector)).norm();
    Vector3r unitVectorY = (_unitVector.cross(unitVectorX)) / (_unitVector.cross(unitVectorX)).norm();
    // ASSERT(checkEquality(unitVectorX.norm(), 1)); //if norm is zero?
    // ASSERT(checkEquality(unitVectorY.norm(), 1)); //
    Matrix3r Transf = Matrix3r::Zero();
    Transf(0, 0) = unitVectorX[0];
    Transf(0, 1) = unitVectorX[1];
    Transf(0, 2) = unitVectorX[2];
    Transf(1, 0) = unitVectorY[0];
    Transf(1, 1) = unitVectorY[1];
    Transf(1, 2) = unitVectorY[2];
    Transf(2, 0) = _unitVector[0];
    Transf(2, 1) = _unitVector[1];
    Transf(2, 2) = _unitVector[2];

    std::vector<Vector3r> verticesOf2DPolygon;
    for (auto &vId : verticesId)
    {
        Vector3r pointOn2DPlane = Transf * (_vertices[vId] - center);
        verticesOf2DPolygon.push_back(pointOn2DPlane);
    }

    area = 0;
    double Cx = 0;
    double Cy = 0;
    for (int i = 0; i != (int)verticesOf2DPolygon.size(); i++)
    {
        Vector3r Pi = verticesOf2DPolygon[i];
        Vector3r Pnext = Vector3r::Zero();
        if (i < (int)verticesOf2DPolygon.size() - 1)
            Pnext = verticesOf2DPolygon[i + 1];
        if (i == (int)verticesOf2DPolygon.size() - 1)
            Pnext = verticesOf2DPolygon[0];

        Cx += (Pi[0] + Pnext[0]) * (Pi[0] * Pnext[1] - Pnext[0] * Pi[1]);
        Cy += (Pi[1] + Pnext[1]) * (Pi[0] * Pnext[1] - Pnext[0] * Pi[1]);
        area += 0.5 * (Pi[0] * Pnext[1] - Pnext[0] * Pi[1]);
    }

    Cx *= double(1) / (6 * area);
    Cy *= double(1) / (6 * area);
    area = abs(area);
    // ASSERT(!std::isnan(area) && area > 0);

    Vector3r centroidIn2DPlaneCoordinates = {Cx, Cy, 0};
    centroid = Transf.transpose() * centroidIn2DPlaneCoordinates + center;
};
