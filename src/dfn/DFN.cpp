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

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/register_ptr_to_python.hpp>
// #include "dfn/DFN.hpp"
#include "DFN.hpp"

using namespace boost::python;
void DFN_Wrapper()
{
	register_ptr_to_python<std::shared_ptr<FractureSet>>();
	register_ptr_to_python<std::shared_ptr<LineMapping>>();
	register_ptr_to_python<std::shared_ptr<SurfaceMapping>>();
	register_ptr_to_python<std::shared_ptr<VolumeMapping>>();

	class_<std::vector<std::shared_ptr<FractureSet>>>("fractureSets")
		.def(vector_indexing_suite<std::vector<std::shared_ptr<FractureSet>>, true>());
	class_<std::vector<std::shared_ptr<LineMapping>>>("linesMapping")
		.def(vector_indexing_suite<std::vector<std::shared_ptr<LineMapping>>, true>());
	class_<std::vector<std::shared_ptr<SurfaceMapping>>>("surfacesMapping")
		.def(vector_indexing_suite<std::vector<std::shared_ptr<SurfaceMapping>>, true>());
	class_<std::vector<std::shared_ptr<VolumeMapping>>>("volumesMapping")
		.def(vector_indexing_suite<std::vector<std::shared_ptr<VolumeMapping>>, true>());

	class_<DFN, boost::noncopyable>("DFN")
		.def("add_FractureSet", &DFN::add_FractureSet)
		.def("add_LineMapping", &DFN::add_LineMapping)
		.def("add_QuadrilateralMapping", &DFN::add_QuadrilateralMapping)
		.def("add_CircularMapping", &DFN::add_CircularMapping)
		.def("add_VolumeMapping", &DFN::add_VolumeMapping)
		.def("set_RegionMaxCorner", &DFN::set_RegionMaxCorner)

		// added
		.def("set_FirstRecBlock", &DFN::set_FirstRecBlock)

		.def("set_RandomSeed", &DFN::set_RandomSeed)
		.def("set_NumberOfBorderPoints", &DFN::set_NumberOfBorderPoints)
		.def("export_DFNVtk", &DFN::export_DFNVtk)
		.def("export_RegionVtk", &DFN::export_RegionVtk)

		.def_readwrite("fractureSets", &DFN::fractureSets)
		.def_readwrite("linesMapping", &DFN::linesMapping)
		.def_readwrite("surfacesMapping", &DFN::surfacesMapping)
		.def_readwrite("volumesMapping", &DFN::volumesMapping);
};

void DFN::set_RandomSeed(int _seed)
{
	if (_seed <= 0)
	{
		srand(std::time(0));
		randomSeed = std::time(0);
	}
	else
	{
		srand(_seed);
		randomSeed = _seed;
	}
};

void DFN::set_RegionMaxCorner(PyList _point)
{
	Vector3r point = pyListToVec3(_point);
	regionMaxCorner = point;
	regionVolume = (regionMaxCorner[0] - regionMinCorner[0]) * (regionMaxCorner[1] - regionMinCorner[1]) * (regionMaxCorner[2] - regionMinCorner[2]);

	Triangle triangle1({regionMinCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMaxCorner[0], regionMaxCorner[1], regionMinCorner[2]});
	Triangle triangle2({regionMinCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMaxCorner[0], regionMaxCorner[1], regionMinCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMinCorner[2]});
	Triangle triangle3({regionMinCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMinCorner[0], regionMinCorner[1], regionMaxCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMaxCorner[2]});
	Triangle triangle4({regionMinCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMaxCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMinCorner[2]});
	Triangle triangle5({regionMinCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMinCorner[0], regionMinCorner[1], regionMaxCorner[2]});
	Triangle triangle6({regionMinCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMinCorner[2]});
	Triangle triangle7({regionMaxCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMinCorner[0], regionMinCorner[1], regionMaxCorner[2]});
	Triangle triangle8({regionMaxCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMinCorner[0], regionMinCorner[1], regionMaxCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMaxCorner[2]});
	Triangle triangle9({regionMaxCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMaxCorner[0], regionMaxCorner[1], regionMinCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMinCorner[2]});
	Triangle triangle10({regionMaxCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMinCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMaxCorner[2]});
	Triangle triangle11({regionMaxCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMaxCorner[0], regionMaxCorner[1], regionMinCorner[2]});
	Triangle triangle12({regionMaxCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMaxCorner[2]});

	modelRegion.clear();
	modelRegion.push_back(triangle1);
	modelRegion.push_back(triangle2);
	modelRegion.push_back(triangle3);
	modelRegion.push_back(triangle4);
	modelRegion.push_back(triangle5);
	modelRegion.push_back(triangle6);
	modelRegion.push_back(triangle7);
	modelRegion.push_back(triangle8);
	modelRegion.push_back(triangle9);
	modelRegion.push_back(triangle10);
	modelRegion.push_back(triangle11);
	modelRegion.push_back(triangle12);
};

void DFN::set_FirstRecBlock(PyList _MinPoint, PyList _MaxPoint)
{
	firstBlkMax = pyListToVec3(_MaxPoint);
	firstBlkMin = pyListToVec3(_MinPoint);
}

DFN::DFN()
{
	regionVolume = (regionMaxCorner[0] - regionMinCorner[0]) * (regionMaxCorner[1] - regionMinCorner[1]) * (regionMaxCorner[2] - regionMinCorner[2]);

	Triangle triangle1({regionMinCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMaxCorner[0], regionMaxCorner[1], regionMinCorner[2]});
	Triangle triangle2({regionMinCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMaxCorner[0], regionMaxCorner[1], regionMinCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMinCorner[2]});
	Triangle triangle3({regionMinCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMinCorner[0], regionMinCorner[1], regionMaxCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMaxCorner[2]});
	Triangle triangle4({regionMinCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMaxCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMinCorner[2]});
	Triangle triangle5({regionMinCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMinCorner[0], regionMinCorner[1], regionMaxCorner[2]});
	Triangle triangle6({regionMinCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMinCorner[2]});
	Triangle triangle7({regionMaxCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMinCorner[0], regionMinCorner[1], regionMaxCorner[2]});
	Triangle triangle8({regionMaxCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMinCorner[0], regionMinCorner[1], regionMaxCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMaxCorner[2]});
	Triangle triangle9({regionMaxCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMaxCorner[0], regionMaxCorner[1], regionMinCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMinCorner[2]});
	Triangle triangle10({regionMaxCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMinCorner[2]}, {regionMinCorner[0], regionMaxCorner[1], regionMaxCorner[2]});
	Triangle triangle11({regionMaxCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMaxCorner[0], regionMaxCorner[1], regionMinCorner[2]});
	Triangle triangle12({regionMaxCorner[0], regionMaxCorner[1], regionMaxCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMinCorner[2]}, {regionMaxCorner[0], regionMinCorner[1], regionMaxCorner[2]});

	modelRegion.clear();
	modelRegion.push_back(triangle1);
	modelRegion.push_back(triangle2);
	modelRegion.push_back(triangle3);
	modelRegion.push_back(triangle4);
	modelRegion.push_back(triangle5);
	modelRegion.push_back(triangle6);
	modelRegion.push_back(triangle7);
	modelRegion.push_back(triangle8);
	modelRegion.push_back(triangle9);
	modelRegion.push_back(triangle10);
	modelRegion.push_back(triangle11);
	modelRegion.push_back(triangle12);
}

void DFN::add_FractureSet()
{
	std::shared_ptr<FractureSet> fractureSet = std::make_shared<FractureSet>(this, (int)fractureSets.size());
	fractureSets.push_back(fractureSet);
	std::cout << "Fracture Set number " << fractureSet->id << " added!" << std::endl;
}

void DFN::export_DFNVtk(std::string _fileName)
{
	std::ofstream out;
	out.open(_fileName + ".vtk");
	out << std::setprecision(15);

	out << "# vtk DataFile Version 3.0" << std::endl;
	out << "Discrete fracture network" << std::endl;
	out << "ASCII" << std::endl;
	int nFracs = 0;
	int nBorderPoints = 0;
	int counter = 0;

	for (auto &FracSet : fractureSets)
	{
		nFracs += (int)FracSet->fractures.size();
		for (auto &Frac : FracSet->fractures)
		{
			nBorderPoints += (int)Frac->borderPoints.size();
		}
	}
	out << " " << std::endl;
	out << "DATASET POLYDATA" << std::endl;
	out << "POINTS " << nBorderPoints << " double" << std::endl;
	for (auto &FracSet : fractureSets)
	{
		for (auto &Frac : FracSet->fractures)
		{
			for (auto &p : Frac->borderPoints)
			{
				out << p[0] << " " << p[1] << " " << p[2] << std::endl;
			}
		}
	}
	out << " " << std::endl;
	out << "POLYGONS " << nFracs << " " << (nFracs + nBorderPoints) << std::endl;
	for (auto &FracSet : fractureSets)
	{
		for (auto &Frac : FracSet->fractures)
		{
			out << (int)Frac->borderPoints.size();
			for (auto &p : Frac->borderPoints)
			{
				out << " " << counter;
				counter++;
			}
			out << std::endl;
		}
	}
	out << " " << std::endl;
	out << "CELL_DATA " << nFracs << std::endl;
	out << "SCALARS FracSetId float" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (auto &FracSet : fractureSets)
	{
		for (auto &Frac : FracSet->fractures)
		{
			out << FracSet->id << std::endl;
		}
	}
	out << " " << std::endl;
	out << "SCALARS FracId float" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (auto &FracSet : fractureSets)
	{
		for (auto &Frac : FracSet->fractures)
		{
			out << Frac->id << std::endl;
		}
	}
	out << " " << std::endl;
	out << "SCALARS Area float" << std::endl;
	out << "LOOKUP_TABLE default" << std::endl;
	for (auto &FracSet : fractureSets)
	{
		for (auto &Frac : FracSet->fractures)
		{
			out << Frac->get_Area() << std::endl;
		}
	}
	out.close();
	std::cout << "DFN Vtk exported!" << std::endl;
}

void DFN::export_RegionVtk(std::string _fileName)
{
	std::ofstream out;
	out.open(_fileName + ".vtk");
	out << std::setprecision(15);

	out << "# vtk DataFile Version 3.0" << std::endl;
	out << "Discrete fracture network" << std::endl;
	out << "ASCII" << std::endl;
	int count = 0;
	int auxPolycount = 0;

	for (auto &Tri : modelRegion)
	{
		count += 3;
	}

	out << " " << std::endl;
	out << "DATASET POLYDATA" << std::endl;
	out << "POINTS " << count << " double" << std::endl;
	for (auto &Tri : modelRegion)
	{
		out << Tri.pointD[0] << " " << Tri.pointD[1] << " " << Tri.pointD[2] << std::endl;
		out << Tri.pointE[0] << " " << Tri.pointE[1] << " " << Tri.pointE[2] << std::endl;
		out << Tri.pointF[0] << " " << Tri.pointF[1] << " " << Tri.pointF[2] << std::endl;
	}

	out << " " << std::endl;
	out << "POLYGONS " << count / 3 << " " << 4 * count / 3 << std::endl;
	for (auto &Tri : modelRegion)
	{
		out << "3 " << auxPolycount << " " << auxPolycount + 1 << " " << auxPolycount + 2 << std::endl;
		auxPolycount += 3;
	}

	std::cout << "DFN Region Vtk exported!" << std::endl;
}

void DFN::add_LineMapping(PyList _pointA, PyList _pointB)
{
	Vector3r pointA = pyListToVec3(_pointA);
	Vector3r pointB = pyListToVec3(_pointB);
	std::shared_ptr<LineMapping> lineMapping = std::make_shared<LineMapping>(pointA, pointB, this);
	linesMapping.push_back(lineMapping);
	std::cout << "Line mapping added!" << std::endl;
}

void DFN::add_QuadrilateralMapping(PyList _pointA, PyList _pointB, PyList _pointC, PyList _pointD)
{
	Vector3r pointA = pyListToVec3(_pointA);
	Vector3r pointB = pyListToVec3(_pointB);
	Vector3r pointC = pyListToVec3(_pointC);
	Vector3r pointD = pyListToVec3(_pointD);
	std::shared_ptr<SurfaceMapping> surfaceMapping = std::make_shared<SurfaceMapping>(this);
	Triangle tri1(pointA, pointB, pointC);
	Triangle tri2(pointC, pointD, pointA);
	if (tri1.unitVector.dot(tri2.unitVector) < 0)
		tri2.unitVector *= -1;
	surfaceMapping->surfaceTriangles.push_back(tri1);
	surfaceMapping->surfaceTriangles.push_back(tri2);
	surfacesMapping.push_back(surfaceMapping);
	std::cout << "Surface mapping added!" << std::endl;
}

void DFN::add_CircularMapping(PyList _center, double _dipDirection, double _dipAngle, double _radius)
{
	Vector3r center = pyListToVec3(_center);

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

	std::shared_ptr<SurfaceMapping> surfaceMapping = std::make_shared<SurfaceMapping>(this);
	std::vector<Vector3r> borderPoints;
	Vector3r vectorCenterToBorder = center.cross(unitVector);
	ASSERT(!checkEquality(vectorCenterToBorder.norm(), 0));
	vectorCenterToBorder.normalize();
	vectorCenterToBorder *= _radius;
	for (int i = 0; i != nFracBorderPoints; i++)
	{
		double rotAngle = i * (2 * M_PI) / (double)nFracBorderPoints;
		Quaternionr q(AngleAxisr(rotAngle, unitVector));
		q.normalize();
		Vector3r pointToBorder = center + q * vectorCenterToBorder;
		borderPoints.push_back(pointToBorder);
	}
	for (int i = 0; i != (int)borderPoints.size(); ++i)
	{
		int j = i + 1;
		if (j == (int)borderPoints.size())
			j = 0;
		Triangle newTri(center, borderPoints[i], borderPoints[j]);
		surfaceMapping->surfaceTriangles.push_back(newTri);
	}
	surfacesMapping.push_back(surfaceMapping);
	std::cout << "Surface mapping added!" << std::endl;
}

void DFN::add_VolumeMapping()
{
	std::shared_ptr<VolumeMapping> volumeMapping = std::make_shared<VolumeMapping>(this);
	volumesMapping.push_back(volumeMapping);
	std::cout << "Volume mapping added!" << std::endl;
}
