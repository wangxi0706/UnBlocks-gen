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

#include "coin/ClpSimplex.hpp"
#include "Functions.hpp"
#include "Assert.hpp"
#include "TolControl.hpp"

bool Functions::triangleTriangleIntersection(const Triangle& _triangle1, const Triangle& _triangle2, Line& _resultLine){
	Line line1T1(_triangle1.pointD,_triangle1.pointE); Line line1T2(_triangle2.pointD,_triangle2.pointE); 
	Line line2T1(_triangle1.pointE,_triangle1.pointF); Line line2T2(_triangle2.pointE,_triangle2.pointF);
	Line line3T1(_triangle1.pointF,_triangle1.pointD); Line line3T2(_triangle2.pointF,_triangle2.pointD); 
    
    Vector3r auxVec;
	std::vector<Vector3r> resultingPoints;
	if (lineTriangleIntersection(line1T1,_triangle2, auxVec)) resultingPoints.push_back(auxVec);
	if (lineTriangleIntersection(line2T1,_triangle2, auxVec)) resultingPoints.push_back(auxVec);
	if (lineTriangleIntersection(line3T1,_triangle2, auxVec)) resultingPoints.push_back(auxVec);
	if (lineTriangleIntersection(line1T2,_triangle1, auxVec)) resultingPoints.push_back(auxVec);
	if (lineTriangleIntersection(line2T2,_triangle1, auxVec)) resultingPoints.push_back(auxVec);
	if (lineTriangleIntersection(line3T2,_triangle1, auxVec)) resultingPoints.push_back(auxVec);
	if ((int)resultingPoints.size() < 2) return false;
	
	for (int i = 0; i!= (int)resultingPoints.size()-1; i++){
		for (int j = i+1; j!= (int)resultingPoints.size(); j++){
			if (checkEquality((resultingPoints[i] - resultingPoints[j]).norm(), 0)) resultingPoints[j][0] = NAN;
		}
	}
    
    resultingPoints.erase(std::remove_if(std::begin(resultingPoints), std::end(resultingPoints), [](Vector3r& _p){return (std::isnan(_p[0]));}), std::end(resultingPoints));
	if ((int)resultingPoints.size() != 2 || checkEquality((resultingPoints[0] - resultingPoints[1]).norm(),0)) return false;
    
    ASSERT(!std::isnan(resultingPoints[0][0]) && !std::isnan(resultingPoints[1][0]));
    _resultLine.pointA = resultingPoints[0];
    _resultLine.pointB = resultingPoints[1];
    return true;
}

bool Functions::pointVolumeIntersection(const Vector3r& _point, const std::vector<Triangle>& _modelBoundary){
	ASSERT((int)_modelBoundary.size() >= 4);
	
	Line LineXDir(_point, {(rand()*9+1)*1e9,(rand()*9+1)*1e9,(rand()*9+1)*1e9});
    
	int nIntersections1 = 0;
	for (auto & Tri : _modelBoundary){
       		Vector3r auxVec;
		if (Functions::lineTriangleIntersection(LineXDir, Tri, auxVec)) nIntersections1++;
	}
	
	if (nIntersections1 % 2 == 0) return false;
	return true;
}

bool Functions::lineTriangleIntersection(const Line& _line, const Triangle& _triangle, Vector3r& _resultPoint){
	Vector3r vecAB =  _line.pointB - _line.pointA;
	Vector3r vecDE = _triangle.pointE - _triangle.pointD;
	Vector3r vecDF = _triangle.pointF - _triangle.pointD;
	Vector3r vecDA = _line.pointA - _triangle.pointD;
	Vector3r vecDB = _line.pointB - _triangle.pointD;
	Vector3r nD = vecDE.cross(vecDF)/(vecDE.cross(vecDF)).norm();
	double dA = vecDA.dot(nD);
	double dB = vecDB.dot(nD);
	Vector3r pointC = _line.pointA + vecAB * (abs(dA) / (abs(dA) + abs(dB)));
	Vector3r vecEF = _triangle.pointF - _triangle.pointE;
	Vector3r vecEC = pointC - _triangle.pointE;
	Vector3r vecDC = pointC - _triangle.pointD;
	Vector3r vecFD = _triangle.pointD - _triangle.pointF;
	Vector3r vecFC = pointC - _triangle.pointF;
	double SDEF = ((vecDE.cross(vecDF)).norm())/2;
	double SDEC = ((vecDE.cross(vecDC)).norm())/2;
	double SEFC = ((vecEF.cross(vecEC)).norm())/2;
	double SFDC = ((vecFD.cross(vecFC)).norm())/2;
	Vector3r pointM = _line.pointA + vecAB * 0.5;
	Vector3r vecDM = pointM - _triangle.pointD;
	double dM = vecDM.dot(nD);
	
	double areaCheck = SDEC+SEFC+SFDC - ND_tolerance*(SDEC+SEFC+SFDC);
	if (!(dA == 0 && dB == 0) && areaCheck <= SDEF && dA*dB == 0 && dM>0) {
		_resultPoint = pointC;
        return true;
	} else if (!(dA == 0 && dB == 0) && areaCheck <= SDEF && dA*dB < 0) {
		_resultPoint = pointC;
        return true;
	} 
    return false;
}

double Functions::vectorsAngle(Vector3r _centerPoint, Vector3r _unitVec, Vector3r _point1, Vector3r _point2){
 	ASSERT(checkEquality(_unitVec.norm(), 1));
    ASSERT(!checkEquality((_point1 - _centerPoint).norm(),0));
 	ASSERT(!checkEquality((_point2 - _centerPoint).norm(),0));
	ASSERT(!checkEquality((_point2 - _point1).norm(), 0));
	Vector3r P1 = (_point1 - _centerPoint); P1.normalize();
	Vector3r P2 = (_point2 - _centerPoint); P2.normalize();
    ASSERT(checkEquality(P1.norm(),1));
    ASSERT(checkEquality(P2.norm(),1));
    double aux = P1.dot(P2);
    if (aux > 1 && checkEquality(aux, 1)) aux = 1;
    if (aux < -1 && checkEquality(aux, -1)) aux = -1;
    ASSERT(aux >= -1 && aux <= 1);
	double angle = acos(aux);
    ASSERT(!std::isnan(angle));
	Vector3r crossP1P2 = P1.cross(P2);
 	if (crossP1P2.dot(_unitVec) > 0) angle = -angle + 2*M_PI;
    // if (crossP1P2.dot(_unitVec) < 0) angle = -angle + 2*M_PI;//modified, make sure outer normal 
    ASSERT(angle >= 0 && angle < 2*M_PI);
	return angle;
}


bool Functions::calculate_threePlanesIntersection(const Plane& _p1, const Plane& _p2, const Plane& _p3, Vector3r& _intersectionPoint){
    Vector3r u = _p2.unitVector.cross(_p3.unitVector);
    double denom = _p1.unitVector.dot(u);
    if (checkEquality(denom, 0)) return false;
    _intersectionPoint = (_p1.d * u + _p1.unitVector.cross(_p3.d * _p2.unitVector - _p2.d * _p3.unitVector)) / denom;
    return true;
}

bool Functions::check_PointIntersection(const Vector3r& _point, const std::vector<Plane>& _planes){
    bool isInside = true;
    for (int i = 0; i != (int)_planes.size(); ++i){
        double aux = _point.dot(_planes[i].unitVector) - _planes[i].d;
        if (aux > 0 && !checkEquality(aux, 0)) {
            isInside = false;
            break;
        }
    }
    return isInside;
}

void Functions::organize_Vertices(Polygon& _polygon, std::vector<Vector3r>& _vertices){
    ASSERT((int)_polygon.verticesId.size() >= 3);
    ASSERT(checkEquality(_polygon.unitVector.norm(), 1)); 
            
    //calculate polygon center from average of vertices coordinates
    Vector3r center = Vector3r::Zero();
    for (auto& id: _polygon.verticesId)
        center += _vertices[id] / (double)_polygon.verticesId.size();
    _polygon.centroid = center;
    
    //In this algorithm, the points with bigger angles are added first in a lowering angle sequence
    std::vector<double> angles; 
    for (int i = 1; i!= (int)_polygon.verticesId.size(); i++)
        angles.push_back(vectorsAngle(_polygon.centroid, _polygon.unitVector, _vertices[_polygon.verticesId[0]], _vertices[_polygon.verticesId[i]]));
            
    double previousBiggestAngle = M_PI*2;
    std::vector<int> organizedVerticesIds;
    organizedVerticesIds.push_back(_polygon.verticesId[0]);
    for (int i = 1; i!= (int)_polygon.verticesId.size(); i++){
        int biggestAngleId = 0; 
        double biggestAngle = 0;
        for (int j = 0; j!= (int)angles.size(); j++){
            if (angles[j] > biggestAngle) {
                biggestAngle = angles[j]; 
                biggestAngleId = j;
            }
        } 
        ASSERT(!checkEquality(biggestAngle, previousBiggestAngle) && !checkEquality(biggestAngle, 0) && !checkEquality(biggestAngle, M_PI*2));
        if (biggestAngle < previousBiggestAngle) {	
            organizedVerticesIds.push_back(_polygon.verticesId[biggestAngleId+1]);
            previousBiggestAngle = biggestAngle;
        }
        angles[biggestAngleId] = 0;	
    }
        
    ASSERT((int)_polygon.verticesId.size() == (int)organizedVerticesIds.size());
    for (int i = 0; i!= (int)_polygon.verticesId.size(); i++)
        _polygon.verticesId[i] = organizedVerticesIds[i];
}

void Functions::organize_Points(Vector3r& _center, Vector3r _unitVector, std::vector<Vector3r>& _points){
    ASSERT((int)_points.size() >= 3);
    ASSERT(checkEquality(_unitVector.norm(), 1)); 
    
    //In this algorithm, the points with bigger angles are added first in a lowering angle sequence
    std::vector<double> angles; 
    for (int i = 1; i!= (int)_points.size(); i++){
        double newAngle = vectorsAngle(_center, _unitVector, _points[0], _points[i]);
        ASSERT(!checkEquality(newAngle, M_PI*2));
        ASSERT(!checkEquality(newAngle, 0));
        angles.push_back(newAngle);
    }
    
    double previousBiggestAngle = M_PI*2;
    std::vector<Vector3r> organizedPoints;
    organizedPoints.push_back(_points[0]);
    for (int i = 1; i!= (int)_points.size(); i++){
        int biggestAngleId = 0; 
        double biggestAngle = 0;
        for (int j = 0; j!= (int)angles.size(); j++){
            if (angles[j] > biggestAngle) {
                biggestAngle = angles[j]; 
                biggestAngleId = j;
            }
        } 
        ASSERT(!checkEquality(biggestAngle, previousBiggestAngle));
        ASSERT(!checkEquality(biggestAngle, 0));
        ASSERT(!checkEquality(biggestAngle, M_PI*2));
        if (biggestAngle < previousBiggestAngle) {	
            organizedPoints.push_back(_points[biggestAngleId+1]);
            previousBiggestAngle = biggestAngle;
        }
        angles[biggestAngleId] = 0;	
    }
        
    ASSERT((int)_points.size() == (int)organizedPoints.size());
    for (int i = 0; i!= (int)_points.size(); i++)
        _points[i] = organizedPoints[i];
}

void Functions::eliminate_RedundantPlanes(std::vector<Plane>& _planes){
    for (int pId = 0; pId != (int)_planes.size(); ++pId){
        ClpSimplex model;
        model.setLogLevel(0);
        model.setPrimalTolerance(ND_tolerance*0.1);
        model.setOptimizationDirection(-1);
                
        int numberColumns = 3;
        int numberRows = (int)_planes.size();
        int numberElements = numberRows*numberColumns;
        model.resize(numberRows, numberColumns);
        
        double* objective = model.objective();
        double* columnUpper = model.columnUpper();
        double* columnLower = model.columnLower();
        double* rowLower = model.rowLower();
        double* rowUpper = model.rowUpper();
        for (int i = 0; i < numberColumns; i++) {
            objective[i] = cZ(_planes[pId].unitVector[i]);
            columnUpper[i] = COIN_DBL_MAX;
            columnLower[i] = -COIN_DBL_MAX;
        }
        for (int i = 0; i < numberRows; i++) {
            rowUpper[i] = cZ(_planes[i].d);
            rowLower[i] = -COIN_DBL_MAX;
        }
                
        double* elements = new double[numberElements];
        CoinBigIndex* starts = new CoinBigIndex[numberColumns+1];
        int* rows = new int[numberElements];;
        int* lengths = new int[numberColumns];
                
        // Now elements
        CoinBigIndex put = 0;
        for (int k = 0; k < numberColumns; k++) {
            starts[k] = put;
            lengths[k] = numberRows;
            for (int i = 0; i < numberRows; i++) {
                rows[put] = i;
                elements[put] = cZ(_planes[i].unitVector[k]);
                put++;
            }
        }
        starts[numberColumns] = put;    
                
        // assign to matrix
        CoinPackedMatrix * matrix = new CoinPackedMatrix(true, 0.0, 0.0);
        matrix->assignMatrix(true, numberRows, numberColumns, numberElements, elements, rows, starts, lengths);
        ClpPackedMatrix * clpMatrix = new ClpPackedMatrix(matrix);
        model.replaceMatrix(clpMatrix, true);
                    
        model.primal();
        ASSERT(model.status() == 0);
        double* columnPrimal = model.primalColumnSolution();
        Vector3r solPoint = {columnPrimal[0], columnPrimal[1], columnPrimal[2]};
                
        if (!(checkEquality(abs(_planes[pId].unitVector.dot(solPoint) - _planes[pId].d), 0))) _planes[pId].d = NAN;
    }
    _planes.erase(std::remove_if(std::begin(_planes), std::end(_planes), [](Plane& _p){return (std::isnan(_p.d));}), std::end(_planes));
}
