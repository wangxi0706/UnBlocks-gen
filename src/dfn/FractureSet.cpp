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

#include "dfn/FractureSet.hpp"
#include "dfn/DFN.hpp"
#include "dfn/Fracture.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/register_ptr_to_python.hpp>

using namespace boost::python;
void FractureSet_Wrapper()
{
    register_ptr_to_python<std::shared_ptr<Fracture>>();

    class_<std::vector<std::shared_ptr<Fracture>>>("fractures")
        .def(vector_indexing_suite<std::vector<std::shared_ptr<Fracture>>, true>());

    class_<FractureSet, boost::noncopyable>("fractureSets", no_init)
        .def("add_BaecherFracture", &FractureSet::add_BaecherFracture)
        .def("add_TriangularFracture", &FractureSet::add_TriangularFracture)
        .def("add_CircularFracture", &FractureSet::add_CircularFracture)
        .def("export_FractureSetVtk", &FractureSet::export_FractureSetVtk)

        .def_readwrite("fractures", &FractureSet::fractures);
};

void FractureSet::add_BaecherFracture(double _meanDipDirection, double _meanDipAngle, double _fisherConstant, std::string _sizeDistribution, double _meanFractureSize, double _sigmaFractureSize)
{
    ASSERT(_sizeDistribution == "log" || _sizeDistribution == "exp" || _sizeDistribution == "det");
    ASSERT(_fisherConstant >= 2 || _fisherConstant <= 0);

    _meanDipDirection *= M_PI / 180.0;
    _meanDipAngle *= M_PI / 180.0;

    boost::mt19937 gen;
    gen.seed(dfn->randomSeed);
    dfn->randomSeed++;
    boost::exponential_distribution<double> exponentialDist(double(1.0) / _meanFractureSize);
    boost::variate_generator<boost::mt19937, boost::exponential_distribution<double>> randomExponential(gen, exponentialDist);
    boost::lognormal_distribution<double> lognormalDist(_meanFractureSize, _sigmaFractureSize);
    boost::variate_generator<boost::mt19937, boost::lognormal_distribution<double>> randomLogNormal(gen, lognormalDist);

    Vector3r poleAllRotated = Vector3r::Zero();
    if (_fisherConstant >= 2)
    {
        double kapa = 0;
        double echis = (M_PI / 2) - _meanDipAngle;
        if (_meanDipDirection >= 0 && _meanDipDirection < M_PI)
            kapa = _meanDipDirection + M_PI;
        if (_meanDipDirection >= M_PI && _meanDipDirection <= 2 * M_PI)
            kapa = _meanDipDirection - M_PI;
        Vector3r poleMean = {cos(kapa) * cos(echis), -sin(kapa) * cos(echis), -sin(echis)};
        ASSERT(checkEquality(poleMean.norm(), 1));

        double fisherDipDevAngle = acos((_fisherConstant + log(1 - randomize(1))) / _fisherConstant);
        echis = (M_PI / 2) - (_meanDipAngle - fisherDipDevAngle);
        Vector3r poleDipRotated = {cos(kapa) * cos(echis), -sin(kapa) * cos(echis), -sin(echis)};
        double randomAngle = randomize(2 * M_PI);

        Quaternionr rotCorrectionUnitVec(AngleAxisr(randomAngle, poleMean));
        rotCorrectionUnitVec.normalize();
        poleAllRotated = rotCorrectionUnitVec * poleDipRotated;
        ASSERT(checkEquality(poleAllRotated.norm(), 1));

        //Add a fracture with random dip angle and random dip direction
    }
    else if (_fisherConstant <= 0)
    {
        double randomDipAngle = randomize(90);
        double randomDipDirection = randomize(2 * M_PI);
        double kapa = 0;
        double echis = (M_PI / 2) - randomDipAngle;
        if (randomDipDirection >= 0 && randomDipDirection < M_PI)
            kapa = randomDipDirection + M_PI;
        if (randomDipDirection >= M_PI && randomDipDirection <= 2 * M_PI)
            kapa = randomDipDirection - M_PI;
        poleAllRotated = {cos(kapa) * cos(echis), -sin(kapa) * cos(echis), -sin(echis)};
        ASSERT(checkEquality(poleAllRotated.norm(), 1));
    }

    //random fracture location (poisson process)
    double randCoordinateX = dfn->regionMinCorner[0] + randomize(dfn->regionMaxCorner[0] - dfn->regionMinCorner[0]);
    double randCoordinateY = dfn->regionMinCorner[1] + randomize(dfn->regionMaxCorner[1] - dfn->regionMinCorner[1]);
    double randCoordinateZ = dfn->regionMinCorner[2] + randomize(dfn->regionMaxCorner[2] - dfn->regionMinCorner[2]);
    ASSERT(randCoordinateX > 0 && randCoordinateY > 0 && randCoordinateZ > 0);
    ASSERT(!checkEquality(randCoordinateX, 0) && !checkEquality(randCoordinateY, 0) && !checkEquality(randCoordinateZ, 0));
    Vector3r newFracCenter = {randCoordinateX, randCoordinateY, randCoordinateZ};

    //add circular fracture
    if (_sizeDistribution == "log")
        add_CircularFractureByUnitVec({randCoordinateX, randCoordinateY, randCoordinateZ}, poleAllRotated, randomLogNormal());
    if (_sizeDistribution == "exp")
        add_CircularFractureByUnitVec({randCoordinateX, randCoordinateY, randCoordinateZ}, poleAllRotated, randomExponential());
    if (_sizeDistribution == "det")
        add_CircularFractureByUnitVec({randCoordinateX, randCoordinateY, randCoordinateZ}, poleAllRotated, _meanFractureSize);
}

void FractureSet::trim_FractureBorders(Fracture &_fracture)
{
    bool shouldTrim = false;
    std::cout << "region points "
              << dfn->regionMaxCorner[0] << " "
              << dfn->regionMaxCorner[1] << " "
              << dfn->regionMaxCorner[2] << "\n";
    for (auto &p : _fracture.borderPoints)
    {
        std::cout << "border points "
                  << p[0] << " " << p[1] << " " << p[2] << "\n";
        if (p[0] < dfn->regionMinCorner[0] || p[0] > dfn->regionMaxCorner[0])
        {
            shouldTrim = true;
            break;
        }
        if (p[1] < dfn->regionMinCorner[1] || p[1] > dfn->regionMaxCorner[1])
        {
            shouldTrim = true;
            break;
        }
        if (p[2] < dfn->regionMinCorner[2] || p[2] > dfn->regionMaxCorner[2])
        {
            shouldTrim = true;
            break;
        }
    }
    // ASSERT(shouldTrim == false);
    if (shouldTrim)
    {
        for (auto &fracTri : _fracture.triangles)
        {
            for (auto &borderTri : dfn->modelRegion)
            {
                Line intersecLine;
                if (Functions::triangleTriangleIntersection(fracTri, borderTri, intersecLine))
                {
                    _fracture.borderPoints.push_back(intersecLine.pointA);
                    _fracture.borderPoints.push_back(intersecLine.pointB);
                }
            }
        }
        for (int i = 0; i != (int)_fracture.borderPoints.size(); i++)
        {
            if (_fracture.borderPoints[i][0] < dfn->regionMinCorner[0] || _fracture.borderPoints[i][0] > dfn->regionMaxCorner[0])
                _fracture.borderPoints[i][0] = NAN;
            if (_fracture.borderPoints[i][1] < dfn->regionMinCorner[1] || _fracture.borderPoints[i][1] > dfn->regionMaxCorner[1])
                _fracture.borderPoints[i][0] = NAN;
            if (_fracture.borderPoints[i][2] < dfn->regionMinCorner[2] || _fracture.borderPoints[i][2] > dfn->regionMaxCorner[2])
                _fracture.borderPoints[i][0] = NAN;
        }
        _fracture.borderPoints.erase(std::remove_if(std::begin(_fracture.borderPoints), std::end(_fracture.borderPoints), [](Vector3r &_p) { return (std::isnan(_p[0])); }), std::end(_fracture.borderPoints));
        for (int i = 0; i < (int)_fracture.borderPoints.size() - 1; i++)
        {
            for (int j = i + 1; j < (int)_fracture.borderPoints.size(); j++)
            {
                if (checkEquality((_fracture.borderPoints[i] - _fracture.borderPoints[j]).norm(), 0))
                    _fracture.borderPoints[j][0] = NAN;
            }
        }
        _fracture.borderPoints.erase(std::remove_if(std::begin(_fracture.borderPoints), std::end(_fracture.borderPoints), [](Vector3r &_p) { return (std::isnan(_p[0])); }), std::end(_fracture.borderPoints));
        Functions::organize_Points(_fracture.center, _fracture.unitVector, _fracture.borderPoints);

        _fracture.trimmedTriangles.clear();
        for (int i = 0; i != (int)_fracture.borderPoints.size(); ++i)
        {
            int j = i + 1;
            if (j == (int)_fracture.borderPoints.size())
                j = 0;
            Triangle newTri(_fracture.center, _fracture.borderPoints[i], _fracture.borderPoints[j]);
            _fracture.trimmedTriangles.push_back(newTri);
        }
    }
}

void FractureSet::add_CircularFracture(PyList _center, double _dipDirection, double _dipAngle, double _radius)
{
    std::cout << "begin done\n";
    Vector3r center = pyListToVec3(_center);
    center[0] -= dfn->offset[0];
    center[1] -= dfn->offset[1];
    center[2] -= dfn->offset[2];
    std::cout << "offset done\n";
    double echis = (M_PI / 2) - _dipAngle * M_PI / 180.0;
    double kapa = 0;
    double dipDirection = _dipDirection * M_PI / 180.0;
    if (dipDirection >= 0 && dipDirection < M_PI)
        kapa = dipDirection + M_PI;
    if (dipDirection >= M_PI && dipDirection <= 2 * M_PI)
        kapa = dipDirection - M_PI;
    Vector3r unitVector = {cos(kapa) * cos(echis), -sin(kapa) * cos(echis), -sin(echis)};
    ASSERT(checkEquality(unitVector.norm(), 1));
    unitVector.normalize();
    // added to debug
    std::cout << "normalize done\n";
    std::vector<Vector3r> borderPoints;
    Vector3r vectorCenterToBorder = center.cross(unitVector);
    ASSERT(!checkEquality(vectorCenterToBorder.norm(), 0));
    vectorCenterToBorder.normalize();
    vectorCenterToBorder *= _radius;
    std::cout << "center " << center[0] << " " << center[1] << " " << center[2] << " "
              << "vectorCenterToBorder " << vectorCenterToBorder[0] << " " << vectorCenterToBorder[1]
              << " " << vectorCenterToBorder[2] << '\n';
    for (int i = 0; i != dfn->nFracBorderPoints; i++)
    {
        double rotAngle = (double)i * (2 * M_PI) / (double)dfn->nFracBorderPoints;
        Quaternionr q(AngleAxisr(rotAngle, unitVector));
        q.normalize();
        Vector3r pointToBorder = center + q * vectorCenterToBorder;
        borderPoints.push_back(pointToBorder);
    }
    std::cout << "begin make frac\n";
    std::shared_ptr<Fracture> fracture = std::make_shared<Fracture>((int)fractures.size(), center, unitVector, borderPoints);
    std::cout << "end make frac\n";
    trim_FractureBorders(*fracture);
    std::cout << "trim done\n";
    fractures.push_back(fracture);
    std::cout << "Fracture id " << fracture->id << " added!" << std::endl;
}

void FractureSet::add_CircularFractureByUnitVec(Vector3r _center, Vector3r _unitVector, double _radius)
{
    ASSERT(checkEquality(_unitVector.norm(), 1));

    std::vector<Vector3r> borderPoints;
    Vector3r vectorCenterToBorder = _center.cross(_unitVector);
    ASSERT(!checkEquality(vectorCenterToBorder.norm(), 0));
    vectorCenterToBorder.normalize();
    vectorCenterToBorder *= _radius;
    for (int i = 0; i != dfn->nFracBorderPoints; i++)
    {
        double rotAngle = i * (2 * M_PI) / (double)dfn->nFracBorderPoints;
        Quaternionr q(AngleAxisr(rotAngle, _unitVector));
        q.normalize();
        Vector3r pointToBorder = _center + q * vectorCenterToBorder;
        borderPoints.push_back(pointToBorder);
    }

    std::shared_ptr<Fracture> fracture = std::make_shared<Fracture>((int)fractures.size(), _center, _unitVector, borderPoints);
    trim_FractureBorders(*fracture);
    fractures.push_back(fracture);
    std::cout << "Fracture id " << fracture->id << " added!" << std::endl;
}

void FractureSet::add_TriangularFracture(PyList _pointA, PyList _pointB, PyList _pointC)
{
    Vector3r pointA = pyListToVec3(_pointA);
    Vector3r pointB = pyListToVec3(_pointB);
    Vector3r pointC = pyListToVec3(_pointC);
    pointA[0] -= dfn->offset[0];
    pointA[1] -= dfn->offset[1];
    pointA[2] -= dfn->offset[2];
    pointB[0] -= dfn->offset[0];
    pointB[1] -= dfn->offset[1];
    pointB[2] -= dfn->offset[2];
    pointC[0] -= dfn->offset[0];
    pointC[1] -= dfn->offset[1];
    pointC[2] -= dfn->offset[2];
    Vector3r center = (pointA + pointB + pointC) / 3.0;
    Vector3r unitVector = (pointB - pointA).cross(pointC - pointA);
    unitVector.normalize();
    ASSERT(checkEquality(unitVector.norm(), 1));

    std::vector<Vector3r> borderPoints;
    borderPoints.push_back(pointA);
    borderPoints.push_back(pointB);
    borderPoints.push_back(pointC);

    std::shared_ptr<Fracture> fracture = std::make_shared<Fracture>((int)fractures.size(), center, unitVector, borderPoints);
    trim_FractureBorders(*fracture); //trim if exceed the bounding region
    fractures.push_back(fracture);
    std::cout << "Triangular fracture id " << fracture->id << " added!" << std::endl;
}

void FractureSet::export_FractureSetVtk(std::string _fileName)
{
    std::ofstream out;
    out.open(_fileName + ".vtk");
    out << std::setprecision(15);

    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "Discrete fracture network" << std::endl;
    out << "ASCII" << std::endl;
    int nFracs = (int)fractures.size();
    int nBorderPoints = 0;
    int counter = 0;

    for (auto &Frac : fractures)
    {
        nBorderPoints += (int)Frac->borderPoints.size();
    }

    out << " " << std::endl;
    out << "DATASET POLYDATA" << std::endl;
    out << "POINTS " << nBorderPoints << " double" << std::endl;
    for (auto &Frac : fractures)
    {
        for (auto &p : Frac->borderPoints)
        {
            out << p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
    }

    out << " " << std::endl;
    out << "POLYGONS " << nFracs << " " << (nFracs + nBorderPoints) << std::endl;
    for (auto &Frac : fractures)
    {
        out << (int)Frac->borderPoints.size();
        for (auto &p : Frac->borderPoints)
        {
            out << " " << counter;
            counter++;
        }
        out << std::endl;
    }

    out << " " << std::endl;
    out << "CELL_DATA " << nFracs << std::endl;
    out << "SCALARS FracSetId float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (auto &Frac : fractures)
    {
        out << id << std::endl;
    }

    out << " " << std::endl;
    out << "SCALARS FracId float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (auto &Frac : fractures)
    {
        out << Frac->id << std::endl;
    }

    out << " " << std::endl;
    out << "SCALARS Area float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (auto &Frac : fractures)
    {
        out << Frac->get_Area() << std::endl;
    }

    out.close();
    std::cout << "DFN Vtk exported!" << std::endl;
}
