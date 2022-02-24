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
#include "Generator.hpp"
#include "dfn/DFN.hpp"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <ostream>

using namespace boost::python;
void Generator_Wrapper()
{
    register_ptr_to_python<std::shared_ptr<Block>>();

    class_<std::vector<std::shared_ptr<Block>>>("blocks")
        .def(vector_indexing_suite<std::vector<std::shared_ptr<Block>>, true>());

    class_<Generator, boost::noncopyable>("Generator")
        .def("set_MinInscribedSphereRadius", &Generator::set_MinInscribedSphereRadius)
        .def("set_MaxAspectRatio", &Generator::set_MaxAspectRatio)
        .def("export_BlocksVtk", &Generator::export_BlocksVtk)

        .def("export_BlocksDDA", &Generator::export_BlocksDDA)
        .def("export_BlocksDDAOpt", &Generator::export_BlocksDDAOpt)

        .def("export_ExcavationElementsVtk", &Generator::export_ExcavationElementsVtk)
        .def("import_ExcavationElementsObj", &Generator::import_ExcavationElementsObj)
        .def("generate_RockMass", &Generator::generate_RockMass)

        // added
        .def("generate_RockMass_Multi", &Generator::generate_RockMass_Multi)
        // added
        .def("add_Fixed_Region", &Generator::add_Fixed_Region)
        // added
        .def("input_BondedBlocks", &Generator::input_BondedBlocks)
        // added
        .def("input_ViscousBound", &Generator::input_ViscousBound)
        // added
        .def("input_RollerBound", &Generator::input_RollerBound)
        // added
        .def("input_InputBound", &Generator::input_InputBound)
        // added
        .def("input_FixFaceBound", &Generator::input_FixFaceBound)

        .def("excavate_RockMass", &Generator::excavate_RockMass)

        .def("get_Volumes", &Generator::get_Volumes)
        .def("get_AlphaValues", &Generator::get_AlphaValues)
        .def("get_BetaValues", &Generator::get_BetaValues)

        .def_readwrite("blocks", &Generator::blocks);
};

boost::python::list Generator::get_Volumes(bool _considerBorderBlocks)
{
    boost::python::list volumes;
    for (auto &b : blocks)
    {
        bool borderBlock = false;
        for (auto &v : b->vertices)
        {
            if (checkEquality(v[0], regionMinCorner[0]))
                borderBlock = true;
            if (checkEquality(v[1], regionMinCorner[1]))
                borderBlock = true;
            if (checkEquality(v[2], regionMinCorner[2]))
                borderBlock = true;
            if (checkEquality(v[0], regionMaxCorner[0]))
                borderBlock = true;
            if (checkEquality(v[1], regionMaxCorner[1]))
                borderBlock = true;
            if (checkEquality(v[2], regionMaxCorner[2]))
                borderBlock = true;
        }
        if (!borderBlock || _considerBorderBlocks)
            volumes.append(b->get_Volume());
    }
    return volumes;
}

boost::python::list Generator::get_AlphaValues(bool _considerBorderBlocks)
{
    boost::python::list alphaValues;
    for (auto &b : blocks)
    {
        bool borderBlock = false;
        for (auto &v : b->vertices)
        {
            if (checkEquality(v[0], regionMinCorner[0]))
                borderBlock = true;
            if (checkEquality(v[1], regionMinCorner[1]))
                borderBlock = true;
            if (checkEquality(v[2], regionMinCorner[2]))
                borderBlock = true;
            if (checkEquality(v[0], regionMaxCorner[0]))
                borderBlock = true;
            if (checkEquality(v[1], regionMaxCorner[1]))
                borderBlock = true;
            if (checkEquality(v[2], regionMaxCorner[2]))
                borderBlock = true;
        }
        if (!borderBlock || _considerBorderBlocks)
            alphaValues.append(b->get_Alpha());
    }
    return alphaValues;
}

boost::python::list Generator::get_BetaValues(bool _considerBorderBlocks)
{
    boost::python::list betaValues;
    for (auto &b : blocks)
    {
        bool borderBlock = false;
        for (auto &v : b->vertices)
        {
            if (checkEquality(v[0], regionMinCorner[0]))
                borderBlock = true;
            if (checkEquality(v[1], regionMinCorner[1]))
                borderBlock = true;
            if (checkEquality(v[2], regionMinCorner[2]))
                borderBlock = true;
            if (checkEquality(v[0], regionMaxCorner[0]))
                borderBlock = true;
            if (checkEquality(v[1], regionMaxCorner[1]))
                borderBlock = true;
            if (checkEquality(v[2], regionMaxCorner[2]))
                borderBlock = true;
        }
        if (!borderBlock || _considerBorderBlocks)
            betaValues.append(b->get_Beta());
    }
    return betaValues;
}

void Generator::generate_RockMass(DFN &_dfn)
{
    regionMaxCorner = _dfn.regionMaxCorner;
    regionMinCorner = _dfn.regionMinCorner;

    // added, to enable multi usage of generate_RockMass, but halfway, not tested
    blockNs.push_back((int)blocks.size());

    // First block generation
    Plane planeX0({-1, 0, 0}, _dfn.regionMinCorner[0]); // x=xmin的面，外法向朝x负方向
    Plane planeX1({1, 0, 0}, _dfn.regionMaxCorner[0]);  // x=xmax的面，外法向朝x正方向
    Plane planeY0({0, -1, 0}, _dfn.regionMinCorner[1]); // y=ymin的面，外法向朝y负方向
    Plane planeY1({0, 1, 0}, _dfn.regionMaxCorner[1]);  // y=ymax的面，外法向朝y正方向
    Plane planeZ0({0, 0, -1}, _dfn.regionMinCorner[2]); // z=zmin的面，外法向朝z负方向
    Plane planeZ1({0, 0, 1}, _dfn.regionMaxCorner[2]);  // z=zmax的面，外法向朝z正方向
    std::vector<Plane> firstBlockPlanes;
    firstBlockPlanes.push_back(planeX0);
    firstBlockPlanes.push_back(planeX1);
    firstBlockPlanes.push_back(planeY0);
    firstBlockPlanes.push_back(planeY1);
    firstBlockPlanes.push_back(planeZ0);
    firstBlockPlanes.push_back(planeZ1);
    std::shared_ptr<Block> firstBlock = std::make_shared<Block>(firstBlockPlanes, (int)blocks.size());
    blocks.push_back(firstBlock);

    for (auto &fracset : _dfn.fractureSets)
    {
        for (auto &frac : fracset->fractures)
        {
            std::cout << "Fracture id " << frac->id << " from Fracture Set " << fracset->id << " is being analysed!" << std::endl;
            int nBlocks = (int)blocks.size();
            for (int i = blockNs.back(); i != nBlocks; ++i)
            {
                // if bounding sphere intersects
                if ((blocks[i]->boundingSphereCenter - frac->center).norm() < blocks[i]->boundingSphereRadius + frac->boundingSphereRadius)
                {
                    if (check_BlockPlaneIntersection<Fracture>(*blocks[i], *frac))
                    {
                        Plane cuttingPlane(frac->center, frac->unitVector);
                        blocks[i]->planes.push_back(cuttingPlane);
                        blocks[i]->calculate_BoundingSphere();
                        blocks[i]->calculate_InscribedSphere();

                        std::vector<Plane> planesForNewBlock = blocks[i]->planes;
                        planesForNewBlock.back().unitVector *= -1; // invert its normal direction
                        planesForNewBlock.back().d *= -1;          // why -=1?
                        std::shared_ptr<Block> newblock = std::make_shared<Block>(planesForNewBlock, (int)blocks.size());

                        if (blocks[i]->boundingSphereRadius / blocks[i]->inscribedSphereRadius > maxAspectRatio || blocks[i]->inscribedSphereRadius < minInscribedSphereRadius)
                        {
                            blocks[i]->planes.pop_back();
                            blocks[i]->calculate_BoundingSphere();
                            blocks[i]->calculate_InscribedSphere();
                        }
                        else
                        {
                            if (newblock->boundingSphereRadius / newblock->inscribedSphereRadius > maxAspectRatio || newblock->inscribedSphereRadius < minInscribedSphereRadius)
                            {
                                blocks[i]->planes.pop_back();
                                blocks[i]->calculate_BoundingSphere();
                                blocks[i]->calculate_InscribedSphere();
                            }
                            else
                            {
                                blocks.push_back(newblock);
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = blockNs.back(); i != (int)blocks.size(); ++i)
    {
        blocks[i]->generate_Geometry();
        blocks[i]->offset = _dfn.offset;
        std::cout << blocks[i]->offset[0] << " "
                  << blocks[i]->offset[1] << " "
                  << blocks[i]->offset[2] << std::endl;
    }
}

void Generator::generate_RockMass_Multi(DFN &_dfn)
{
    int startIndex = (int)blocks.size();

    regionMaxCorner = _dfn.regionMaxCorner;
    regionMinCorner = _dfn.regionMinCorner;
    // First block generation
    //  Plane planeX0({-1, 0, 0}, _dfn.firstBlkMin[0]);
    //  Plane planeX1({1, 0, 0}, _dfn.firstBlkMax[0]);
    //  Plane planeY0({0, -1, 0}, _dfn.firstBlkMin[1]);
    //  Plane planeY1({0, 1, 0}, _dfn.firstBlkMax[1]);
    //  Plane planeZ0({0, 0, -1}, _dfn.firstBlkMin[2]);
    //  Plane planeZ1({0, 0, 1}, _dfn.firstBlkMax[2]);
    Plane planeX0(_dfn.firstBlkMin, {-1, 0, 0});
    Plane planeX1(_dfn.firstBlkMax, {1, 0, 0});
    Plane planeY0(_dfn.firstBlkMin, {0, -1, 0});
    Plane planeY1(_dfn.firstBlkMax, {0, 1, 0});
    Plane planeZ0(_dfn.firstBlkMin, {0, 0, -1});
    Plane planeZ1(_dfn.firstBlkMax, {0, 0, 1});
    std::vector<Plane> firstBlockPlanes;
    firstBlockPlanes.push_back(planeX0);
    firstBlockPlanes.push_back(planeX1);
    firstBlockPlanes.push_back(planeY0);
    firstBlockPlanes.push_back(planeY1);
    firstBlockPlanes.push_back(planeZ0);
    firstBlockPlanes.push_back(planeZ1);
    std::shared_ptr<Block> firstBlock = std::make_shared<Block>(firstBlockPlanes, (int)blocks.size());
    blocks.push_back(firstBlock);

    for (auto &fracset : _dfn.fractureSets)
    {
        for (auto &frac : fracset->fractures)
        {
            std::cout << "Fracture id " << frac->id << " from Fracture Set " << fracset->id << " is being analysed!" << std::endl;
            int nBlocks = (int)blocks.size();
            for (int i = startIndex; i != nBlocks; ++i)
            {
                if ((blocks[i]->boundingSphereCenter - frac->center).norm() < blocks[i]->boundingSphereRadius + frac->boundingSphereRadius)
                {
                    if (check_BlockPlaneIntersection<Fracture>(*blocks[i], *frac))
                    {
                        Plane cuttingPlane(frac->center, frac->unitVector);
                        blocks[i]->planes.push_back(cuttingPlane);
                        blocks[i]->calculate_BoundingSphere();
                        blocks[i]->calculate_InscribedSphere();

                        std::vector<Plane> planesForNewBlock = blocks[i]->planes;
                        planesForNewBlock.back().unitVector *= -1;
                        planesForNewBlock.back().d *= -1;
                        std::shared_ptr<Block> newblock = std::make_shared<Block>(planesForNewBlock, (int)blocks.size());

                        if (blocks[i]->boundingSphereRadius / blocks[i]->inscribedSphereRadius > maxAspectRatio || blocks[i]->inscribedSphereRadius < minInscribedSphereRadius)
                        {
                            blocks[i]->planes.pop_back();
                            blocks[i]->calculate_BoundingSphere();
                            blocks[i]->calculate_InscribedSphere();
                        }
                        else
                        {
                            if (newblock->boundingSphereRadius / newblock->inscribedSphereRadius > maxAspectRatio || newblock->inscribedSphereRadius < minInscribedSphereRadius)
                            {
                                blocks[i]->planes.pop_back();
                                blocks[i]->calculate_BoundingSphere();
                                blocks[i]->calculate_InscribedSphere();
                            }
                            else
                            {
                                blocks.push_back(newblock);
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = startIndex; i != (int)blocks.size(); ++i)
    {
        blocks[i]->generate_Geometry();
    }
}

void Generator::add_Fixed_Region(PyList _MinPoint, PyList _MaxPoint)
{
    Box box;
    box.minCor = pyListToVec3(_MinPoint);
    box.maxCor = pyListToVec3(_MaxPoint);
    fixRegion.push_back(box);
}

void Generator::input_BondedBlocks(PyList _bondedIndexes)
{
    // convert input list into vector of indexes
    vecInt Indexes = pyListToVecInt(_bondedIndexes);
    // insert input list into bondIndexes
    vBondIndex.insert(vBondIndex.end(), Indexes.begin(), Indexes.end());

    int l = (int)vBondIndex.size();
    ASSERT(l % 2 == 0);
}

void Generator::input_ViscousBound(int _visBlkId, PyList _visNor) // get viscous indexes and normals
{
    // store the block id
    vVisIndex.push_back(_visBlkId);

    // unify the normal vector
    Vector3r normalDir = pyListToVec3(_visNor);
    double norm = normalDir.norm();
    ASSERT(!checkEquality(norm, 0));

    // store the normal vector
    vVisNormal.push_back(normalDir / norm);
}

void Generator::input_InputBound(
    int _inputBlkId, PyList _inputNor)
{
    // store the input id
    vInputIndex.push_back(_inputBlkId);

    // unify normal vector
    Vector3r normalDir = pyListToVec3(_inputNor);
    double norm = normalDir.norm();
    ASSERT(!checkEquality(norm, 0));

    vInputNormal.push_back(normalDir);
}

void Generator::input_RollerBound(int _rollerBlkId, PyList _rollerNor)
{
    // store the block id
    vRollerIndex.push_back(_rollerBlkId);

    // unify the normal vector
    Vector3r normalDir = pyListToVec3(_rollerNor);
    double norm = normalDir.norm();
    ASSERT(!checkEquality(norm, 0));

    // store the normal vector
    vRollerNormal.push_back(normalDir / norm);
}

void Generator::input_FixFaceBound(int _fixFaceBlkId, PyList _fixFaceNor)
{
    // store the block id
    vFixFaceIndex.push_back(_fixFaceBlkId);

    // unify the normal vector
    Vector3r normalDir = pyListToVec3(_fixFaceNor);
    double norm = normalDir.norm();
    ASSERT(!checkEquality(norm, 0));

    // store the normal vector
    vFixFaceNormal.push_back(normalDir / norm);
}

void Generator::excavate_RockMass()
{
    for (auto &cElem : excavationElements)
    {
        std::cout << "Construction Element is being analysed!" << std::endl;
        Plane cuttingPlane(cElem.center, cElem.unitVector);
        int nBlocks = (int)blocks.size();
        for (int i = 0; i != nBlocks; ++i)
        {
            if ((blocks[i]->boundingSphereCenter - cElem.center).norm() < blocks[i]->boundingSphereRadius + cElem.boundingSphereRadius)
            {
                if (check_BlockPlaneIntersection<ExcavationElement>(*blocks[i], cElem))
                {
                    std::shared_ptr<Block> newblock = std::make_shared<Block>(*blocks[i]);
                    blocks[i]->planes.push_back(cuttingPlane);
                    newblock->id = blocks[i]->id;
                    newblock->planes.push_back(cuttingPlane);
                    newblock->planes.back().unitVector *= -1;
                    newblock->planes.back().d *= -1;
                    blocks.push_back(newblock);
                }
            }
        }
    }

    for (int i = 0; i != (int)blocks.size(); ++i)
        blocks[i]->generate_Geometry();

    std::vector<Triangle> excavationRegion;
    for (auto &cElem : excavationElements)
    {
        excavationRegion.push_back(cElem.triangle);
    }

    for (auto &b : blocks)
    {
        Vector3r center = Vector3r::Zero();
        for (auto &v : b->vertices)
        {
            center += v / (double)b->vertices.size();
        }
        if (Functions::pointVolumeIntersection(center, excavationRegion))
            b->id = -1;
    }
    blocks.erase(std::remove_if(std::begin(blocks), std::end(blocks), [](std::shared_ptr<Block> &_b)
                                { return (_b->id == -1); }),
                 std::end(blocks));
}

template <class X>
bool Generator::check_BlockPlaneIntersection(const Block &_block, const X &_plane)
{
    ClpSimplex model;
    model.setLogLevel(0);
    model.setPrimalTolerance(ND_tolerance * 0.1);
    model.setOptimizationDirection(1);

    int numberColumns = 4;
    int numberRows = (int)_block.planes.size() + (int)_plane.borderPlanes.size() + 1;
    int numberElements = numberRows * numberColumns;
    model.resize(numberRows, numberColumns);

    double *objective = model.objective();
    double *columnUpper = model.columnUpper();
    double *columnLower = model.columnLower();
    double *rowLower = model.rowLower();
    double *rowUpper = model.rowUpper();
    for (int i = 0; i < numberColumns; i++)
    {
        (i < numberColumns - 1) ? objective[i] = 0 : objective[i] = 1;
        columnUpper[i] = COIN_DBL_MAX;
        columnLower[i] = -COIN_DBL_MAX;
    }
    for (int i = 0; i < (int)_block.planes.size(); i++)
    {
        rowUpper[i] = cZ(_block.planes[i].d);
        rowLower[i] = -COIN_DBL_MAX;
    }
    for (int i = 0; i < (int)_plane.borderPlanes.size(); i++)
    {
        rowUpper[i + (int)_block.planes.size()] = cZ(_plane.borderPlanes[i].d);
        rowLower[i + (int)_block.planes.size()] = -COIN_DBL_MAX;
    }
    double d = cZ(_plane.center.dot(_plane.unitVector));
    rowUpper[numberRows - 1] = d;
    rowLower[numberRows - 1] = d;

    double *elements = new double[numberElements];
    CoinBigIndex *starts = new CoinBigIndex[numberColumns + 1];
    int *rows = new int[numberElements];
    ;
    int *lengths = new int[numberColumns];

    // Now elements
    CoinBigIndex put = 0;
    for (int k = 0; k < numberColumns; k++)
    {
        starts[k] = put;
        lengths[k] = numberRows;
        for (int i = 0; i < numberRows; i++)
        {
            rows[put] = i;
            if (k < numberColumns - 1 && i < (int)_block.planes.size())
                elements[put] = cZ(_block.planes[i].unitVector[k]);
            if (k == numberColumns - 1 && i < (int)_block.planes.size())
                elements[put] = -1;
            if (k < numberColumns - 1 && i >= (int)_block.planes.size() && i < numberRows - 1)
                elements[put] = cZ(_plane.borderPlanes[i - (int)_block.planes.size()].unitVector[k]);
            if (k == numberColumns - 1 && i >= (int)_block.planes.size() && i < numberRows - 1)
                elements[put] = -1;
            if (k < numberColumns - 1 && i == numberRows - 1)
                elements[put] = cZ(_plane.unitVector[k]);
            if (k == numberColumns - 1 && i == numberRows - 1)
                elements[put] = 0;
            put++;
        }
    }
    starts[numberColumns] = put;

    // assign to matrix
    CoinPackedMatrix *matrix = new CoinPackedMatrix(true, 0.0, 0.0);
    matrix->assignMatrix(true, numberRows, numberColumns, numberElements, elements, rows, starts, lengths);
    ClpPackedMatrix *clpMatrix = new ClpPackedMatrix(matrix);
    model.replaceMatrix(clpMatrix, true);

    model.primal();
    ASSERT(model.status() == 0);
    double *columnPrimal = model.primalColumnSolution();
    if (columnPrimal[3] < 0 && !checkEquality(columnPrimal[3], 0))
        return true;
    return false;
}

void Generator::export_BlocksVtk(std::string _fileName)
{
    int nTotalBlocks = (int)blocks.size();
    int nTotalFaces = 0;
    int nTotalVerts = 0;
    int nTotalPoints = 0;
    for (auto &b : blocks)
    {
        nTotalVerts += (int)b->vertices.size();
        nTotalFaces += (int)b->polygons.size();
        for (auto &p : b->polygons)
        {
            nTotalPoints += (int)p.verticesId.size();
        }
    }

    std::ofstream out;
    out.open(_fileName + ".vtk");
    out << std::setprecision(15);
    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "Voronoi results" << std::endl;
    out << "ASCII" << std::endl;
    out << " " << std::endl;
    out << "DATASET UNSTRUCTURED_GRID" << std::endl;
    out << "POINTS " << std::to_string(nTotalVerts) << " float" << std::endl;
    for (auto &b : blocks)
    {
        for (auto &v : b->vertices)
        {
            out
                << v[0] + b->offset[0] << " "
                << v[1] + b->offset[1] << " "
                << v[2] + b->offset[2]
                << std::endl;
        }
    }

    out << " " << std::endl;
    out << "CELLS " << nTotalBlocks << " " << std::to_string(2 * nTotalBlocks + nTotalFaces + nTotalPoints) << std::endl;
    int auxId = 0;
    for (auto &b : blocks)
    {
        int nFaces = (int)b->polygons.size();
        int nPoints = 0;
        for (auto &p : b->polygons)
            nPoints += (int)p.verticesId.size();
        out << 1 + nFaces + nPoints << " " << nFaces << " ";
        for (auto &p : b->polygons)
        {
            out << (int)p.verticesId.size() << " ";
            for (auto &vId : p.verticesId)
            {
                out << vId + auxId << " ";
            }
        }
        auxId += (int)b->vertices.size();
        out << std::endl;
    }
    out << std::endl;

    out << " " << std::endl;
    out << "CELL_TYPES " << nTotalBlocks << std::endl;
    ;
    for (auto &b : blocks)
    {
        out << 42 << std::endl;
    }

    out << " " << std::endl;
    out << "CELL_DATA " << nTotalBlocks << std::endl;
    out << "SCALARS id float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (auto &b : blocks)
    {
        out << b->id << std::endl;
    }

    out << " " << std::endl;
    out << "SCALARS volume float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (auto &b : blocks)
    {
        out << b->get_Volume() << std::endl;
    }

    out << " " << std::endl;
    out << "SCALARS alpha float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (auto &b : blocks)
    {
        out << b->get_Alpha() << std::endl;
    }

    out << " " << std::endl;
    out << "SCALARS beta float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (auto &b : blocks)
    {
        out << b->get_Beta() << std::endl;
    }

    out << " " << std::endl;
    out << "SCALARS order float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (auto &b : blocks)
    {
        out << b->get_Order() << std::endl;
    }

    out << " " << std::endl;
    out << "SCALARS aspectRatio float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (auto &b : blocks)
    {
        out << b->get_AspectRatio() << std::endl;
    }

    out << " " << std::endl;
    out << "SCALARS inscribedSphereRadius float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    for (auto &b : blocks)
    {
        out << b->get_InscribedSphereRadius() << std::endl;
    }

    out.close();
}

void Generator::export_BlocksDDA(std::string _fileName)
{
    int nTotalBlocks = (int)blocks.size();

    std::ofstream out;
    out.open(_fileName + ".dgi");
    out << std::setprecision(15);

    out << "dda geometrical input" << std::endl;
    out << "#Blocks: " << nTotalBlocks << std::endl;

    for (int i = 0; i < nTotalBlocks; i++)
    {
        out << "____________________________" << std::endl;
        out << "Block " << i << std::endl;
        auto &b = blocks[i];

        auto &vertices = b->vertices;
        int nv = (int)vertices.size();
        out << "#Vertices " << nv << std::endl;
        for (int j = 0; j < nv; j++)
        {
            out << vertices[j][0] << " " << vertices[j][1] << " " << vertices[j][2] << std::endl;
        }

        auto &faces = b->polygons;
        int nf = (int)faces.size();
        out << "#Faces " << nf << std::endl;
        for (int j = 0; j < nf; j++)
        {
            out << (int)faces[j].verticesId.size() << " ";
            for (auto &vId : faces[j].verticesId)
            {
                out << vId << " ";
            }
            out << std::endl;
        }

        auto &edges = b->edges;
        int ne = (int)edges.size();
        out << "#Edges " << ne << std::endl;
        for (int j = 0; j < ne; j++)
        {
            out << edges[j].verticesIdA << " " << edges[j].verticesIdB << std::endl;
        }
    }

    out.close();
}

// added, output in formats that can be easily processed
void Generator::export_BlocksDDAOpt(std::string _fileName)
{
    int nTotalBlocks = (int)blocks.size();
    int nTotalFaces = 0;
    int nTotalEdges = 0;
    int nTotalVerts = 0;
    int nTotalPoints = 0;
    for (auto &b : blocks)
    {
        nTotalVerts += (int)b->vertices.size();
        nTotalFaces += (int)b->polygons.size();
        nTotalEdges += (int)b->edges.size();
        for (auto &p : b->polygons)
        {
            nTotalPoints += (int)p.verticesId.size();
            // check if the vertices are arranged with counter-clock order
            if (p.unitVector.dot((b->vertices[p.verticesId[1]] - b->vertices[p.verticesId[0]]).cross(b->vertices[p.verticesId[2]] - b->vertices[p.verticesId[1]])) < 0)
                printf("Block %d, face , not outer normal\n", b->id);
        }
    }

    std::ofstream out;
    out.open(_fileName + ".dda");
    out << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
    // out << std::setprecision(15);
    out << "#_DDA_DataFile_Version_1.0" << std::endl;

    out << "POINTS_START_LENGTH " << std::to_string(nTotalBlocks) << std::endl;
    int auxId = 0;
    for (auto &b : blocks)
    {
        // output vertex start id and vertex length
        out << auxId << " " << (int)b->vertices.size() << std::endl;
        auxId += (int)b->vertices.size();
    }

    out << " " << std::endl;
    out << "POINTS " << std::to_string(nTotalVerts) << " double" << std::endl;
    for (auto &b : blocks)
    {
        for (auto &v : b->vertices)
        {
            out
                << v[0] + b->offset[0] << " "
                << v[1] + b->offset[1] << " "
                << v[2] + b->offset[2] << std::endl;
        }
    }
    out << " " << std::endl;

    struct index
    {
        int start;
        int len;
    };

    std::vector<index> faces;
    std::vector<index> Blockfaces;

    out << " " << std::endl;
    out << "FACES_NODELIST " << std::to_string(nTotalPoints) << std::endl;
    auxId = 0;
    int auxId_f = 0;
    int auxId_blkf = 0;
    for (auto &b : blocks)
    {
        index f = {auxId_blkf, (int)b->polygons.size()}; //每个块体face start, face size
        auxId_blkf += (int)b->polygons.size();

        Blockfaces.push_back(f);
        for (auto &p : b->polygons)
        {
            index ff = {auxId_f, (int)p.verticesId.size()}; //每个face vertex start, vertex size
            faces.push_back(ff);
            auxId_f += (int)p.verticesId.size();

            // out << (int)p.verticesId.size() << " ";
            // if (p.unitVector.dot((b->vertices[p.verticesId[1]] - b->vertices[p.verticesId[0]]).cross(b->vertices[p.verticesId[2]] - b->vertices[p.verticesId[1]])) > 0)
            // {
            for (auto &vId : p.verticesId)
            {
                out << vId << " "; // local id, global id=local id + start id
            }
            // }
            // else // if not counter-clock wise, revert it
            // {
            //     for (int i = p.verticesId.size() - 1; i >= 0; i--)
            //         out << p.verticesId[i] << " ";
            // }
        }
        // auxId += b->vertices.size();
        out << std::endl;
    }

    out << " " << std::endl;
    out << "FACES " << std::to_string(nTotalFaces) << std::endl;
    for (auto &f : faces)
    {
        out << f.start << " " << f.len << " " << std::endl;
    }
    out << " " << std::endl;

    out << " " << std::endl;
    out << "BLOCK_FACES " << std::to_string(nTotalBlocks) << std::endl;
    for (auto &bf : Blockfaces)
    {
        out << bf.start << " " << bf.len << " " << std::endl;
    }
    out << " " << std::endl;

    std::vector<index> Blockedges;
    out << " " << std::endl;
    out << "EDGES " << std::to_string(nTotalEdges) << std::endl;
    auxId = 0;
    for (auto &b : blocks)
    {
        int nEdges = (int)b->edges.size();
        index edge = {auxId, nEdges};
        Blockedges.push_back(edge);
        // out << auxId << " " << nEdges << " ";
        // out << 1 + nFaces + nPoints << " " << nFaces << " ";
        for (auto &e : b->edges)
        {
            out << e.verticesIdA << " " << e.verticesIdB << " "; // local id, global id=local id + start id
        }
        auxId += nEdges;
        out << std::endl;
    }
    out << std::endl;

    out << "EDGES " << std::to_string(nTotalBlocks) << std::endl;
    for (auto &e : Blockedges)
    {
        out << e.start << " " << e.len << " " << std::endl;
    }
    out << " " << std::endl;

    out << "SCALARS aspectRatio float" << std::endl;
    for (auto &b : blocks)
    {
        out << b->get_AspectRatio() << std::endl;
    }

    out << " " << std::endl;
    out << "SCALARS inscribedSphereRadius float" << std::endl;
    for (auto &b : blocks)
    {
        out << b->get_InscribedSphereRadius() << std::endl;
    }

    // out << " " << std::endl;
    // out << "SCALARS inscribedSphereRadius float" << std::endl;
    // for (auto &b : blocks)
    // {
    //     out << b->get_InscribedSphereRadius() << std::endl;
    // }

    out << " " << std::endl;
    out << "SCALARS fixedFlag int" << std::endl;
    for (auto &b : blocks)
    {
        int flag = 1;
        if (fixRegion.size() == 0)
        {
            flag = 0;
            goto to;
        }

        for (auto &box : fixRegion)
        {
            for (auto &v : b->vertices)
            {
                if (!inBox(box, v, b->offset))
                {
                    flag = 0; // not in box, so free
                    goto to;
                }
            }
        }
    to:
        out << flag << std::endl;
    }

    out << " " << std::endl;
    out << "SCALARS volume float" << std::endl;
    for (auto &b : blocks)
    {
        out << b->get_Volume() << std::endl;
    }

    ////1 BOND: should output the bonded blocks through bonded index pairs
    out << " " << std::endl;
    out << "/////////////Boundaries///////////" << std::endl;
    int nBond = (int)vBondIndex.size() / 2;
    ASSERT(vBondIndex.size() % 2 == 0);

    out << "Bond: " << std::to_string(nBond) << std::endl;
    for (size_t i = 0; i < nBond; i++)
    {
        // get block index and block
        int nB1 = vBondIndex[2 * i], nB2 = vBondIndex[2 * i + 1];
        std::vector<Polygon> &pos1 = blocks[nB1]->polygons;
        std::vector<Polygon> &pos2 = blocks[nB2]->polygons;
        // find the face id
        int ipoi(-1), jpoj(-1);
        for (size_t ipo = 0; ipo < pos1.size(); ipo++)
        {
            Vector3r &n1 = pos1[ipo].unitVector;
            for (size_t jpo = 0; jpo < pos2.size(); jpo++)
            {
                Vector3r &n2 = pos2[jpo].unitVector;
                if (abs(n1.dot(n2) + 1) < sin(TOLA / 180 * 3.1415926)) // found the opposite face
                {
                    // when opposite, calculate the distance
                    if (abs(pos2[jpo].unitVector.dot(
                            pos1[ipo].centroid - pos2[jpo].centroid)) < TOLD)
                    {
                        ipoi = ipo;
                        jpoj = jpo;
                    }
                }
            }
        }
        // ensure found the bond face, and the face has exact the same number of vertices
        ASSERT(ipoi >= 0 && jpoj >= 0 &&
               pos2[jpoj].verticesId.size() == pos1[ipoi].verticesId.size());
        out << nB1 << ' ' << nB2 << ' ' << ipoi << ' ' << jpoj
            << ' ' << pos2[jpoj].verticesId.size() << ' ' << std::endl;
    }

    ////2 VIS: should output viscous boundaries through block and orientation
    int nVis = (int)vVisIndex.size();
    out << std::endl
        << "Viscous: " << std::to_string(nVis) << std::endl;
    for (size_t i = 0; i < nVis; i++)
    {
        // get block index and block
        int nB = vVisIndex[i];
        // normal vector and polygons
        Vector3r &normal = vVisNormal[i];
        std::vector<Polygon> &pos = blocks[nB]->polygons;

        // to store the face index to be found
        int nF = -1;
        // find the face to be apply viscous boundary
        for (size_t j = 0; j < pos.size(); j++)
        {
            if (abs((pos[j].unitVector.dot(normal) - 1)) < sin(TOLA / 180 * 3.1415926))
            {
                nF = j;
            }
        }

        // debug
        // Vector3r &unit=pos[nF].unitVector;
        // printf("Vis: %e, %e, %e; %e, %e, %e\n",
        //     normal[0],normal[1],normal[2],  unit[0],unit[1],unit[2]);

        ASSERT(nF >= 0); // if fail, fail to find the face to be viscous boundary

        // Found the face id, output the block id and face id.
        out << nB << ' ' << nF << std::endl;
    }

    ////3 ROLL: should output roller boundaries through block and orientation
    int nRoller = (int)vRollerIndex.size();
    out << std::endl
        << "Roller: " << std::to_string(nRoller) << std::endl;
    for (size_t i = 0; i < nRoller; i++)
    {
        // get block index and block
        int nB = vRollerIndex[i];
        // normal vector and polygons
        Vector3r &normal = vRollerNormal[i];
        std::vector<Polygon> &pos = blocks[nB]->polygons;

        // to store the face index to be found
        int nF = -1;
        // find the face to be apply viscous boundary
        for (size_t j = 0; j < pos.size(); j++)
        {
            if (abs((pos[j].unitVector.dot(normal) - 1)) < sin(TOLA / 180 * 3.1415926))
            {
                nF = j;
            }
        }

        // debug
        // Vector3r &unit=pos[nF].unitVector;
        // printf("Vis: %e, %e, %e; %e, %e, %e\n",
        //     normal[0],normal[1],normal[2],  unit[0],unit[1],unit[2]);

        ASSERT(nF >= 0); // if fail, fail to find the face to be viscous boundary

        // Found the face id, output the block id and face id.
        out << nB << ' ' << nF << std::endl;
    }

    ////4 INPUT: should output input boundaries through block and orientation
    int nInput = (int)vInputIndex.size();
    out << std::endl
        << "Input: " << std::to_string(nInput) << std::endl;
    for (size_t i = 0; i < nInput; i++)
    {
        // get block index and block
        int nB = vInputIndex[i];
        // normal vector and polygons
        Vector3r &normal = vInputNormal[i];
        std::vector<Polygon> &pos = blocks[nB]->polygons;

        // to store the face index to be found
        int nF = -1;
        double dArea = 0;
        // find the face to be apply viscous boundary
        for (size_t j = 0; j < pos.size(); j++)
        {
            if (abs((pos[j].unitVector.dot(normal) - 1)) < sin(TOLA / 180 * 3.1415926))
            {
                nF = j;
                dArea = pos[j].area;
            }
        }
        // if fail, fail to find the face to be viscous boundary
        ASSERT(nF >= 0 && dArea > 0);

        // Found the face id, output the block id and face id.
        out << nB << ' ' << nF << ' ' << dArea << ' '
            << pos[nF].centroid[0] << ' '
            << pos[nF].centroid[1] << ' '
            << pos[nF].centroid[2] << std::endl;
    }

    ////5 FIXFACE: should output fixed faces
    int nFixFace = (int)vFixFaceIndex.size();
    out << std::endl
        << "FixFace: " << std::to_string(nFixFace) << std::endl;
    for (size_t i = 0; i < nFixFace; i++)
    {
        // get block index and block
        int nB = vFixFaceIndex[i];
        // normal vector and polygons
        Vector3r &normal = vFixFaceNormal[i];
        std::vector<Polygon> &pos = blocks[nB]->polygons;

        // to store the face index to be found
        int nF = -1;
        // find the face to be apply viscous boundary
        for (size_t j = 0; j < pos.size(); j++)
        { // find the face to be apply viscous boundary
            if (abs((pos[j].unitVector.dot(normal) - 1)) < sin(TOLA / 180 * 3.1415926))
            {
                nF = j;
            }
        }

        ASSERT(nF >= 0); // if fail, fail to find the face to be viscous boundary

        // Found the face id, output the block id and face id.
        out << nB << ' ' << nF << std::endl;
    }

    out.close();
}
// return false if do not find opposite faces
// bool get_BondInfo(Block& b1, Block& b2, Bond& bond){
// // get face id, point, direction, area, is mech needed?
// // first loop the faces
// std::vector<Vector3r> &vertices1=b1.vertices;
// std::vector<Vector3r> &vertices2=b2.vertices;
// for (size_t i = 0; i < b1.polygons; i++)
// {
//     for (size_t j = 0; j < b2.polygons; j++)
//     {
//         Polygon & po1=polygons[i], & po2=polygons[j];
//         // 0.5 °的限额，如果相差在0.5度以内，就算是对上了。
//         if(abs(po1.unitVector.dot(po2.unitVector)+1)
//                 <sin(TOLA/180*3.1415926)){
//             // 现在就算是找到了,然后存储点、normal、area等信息
//             bond.nFId1=i,bond.nFId2=j, bond.dArea=po1.area;
//             bond.nDir=(po2.unitVector-po1.unitVector)/2;

//             std::vector<int> &verticesId1=po1.verticesId;
//             std::vector<int> &verticesId2=po2.verticesId;
//             ASSERT(verticesId1.size()==verticesId2.size());//应该是两个相同的平面
//             ASSERT(verticesId1.size()==NBONDP);            //还是考虑长方体的粘结
//             bool flag=false;//denote whether found the same point or not
//             int nvi,nvj;
//             for (size_t vi = 0; i < verticesId1.size(); i++)
//             {
//                 for (size_t vj = 0; j < verticesId1.size(); j++)
//                 {
//                     Vector3r& pvi=vertices1[verticesId1[vi]],&pvj=vertices2[verticesId2[vj]];
//                     if((pvi-pvj).norm()<TOLD){//found the same point
//                         flag=true;
//                         nvi=vi,nvj=vj;
//                     }
//                 }
//             }
//             if(flag){//found the same point on the opposite plane
//                 for (size_t vi = 0; i < verticesId1.size(); i++)
//                 {
//                     int id1=(vi+nvi>=verticesId1.size())?
//                             (vi+nvi-verticesId1.size())
//                             :(vi+nvi);
//                     int id2=(vi+nvj>=verticesId2.size())?
//                             (vi+nvj-verticesId2.size())
//                             :(vi+nvj);
//                     Vector3r& pvi=vertices1[verticesId1[id1]],
//                             & pvj=vertices2[verticesId2[id2]];
//                     ASSERT((pvi-pvj).norm()<TOLD);  //should be the same point;
//                     bond.P1[vi]=pvi,bond.P2[vi]=pvj;//should record the adjacent points
//                 }
//                 return true;
//             }
//         }
//     }
// }
//     return false;
// }

void Generator::import_ExcavationElementsObj(std::string _fileName)
{
    std::ifstream in;
    in.open(_fileName + ".obj");
    ASSERT(!in.fail());

    // Read coordinates from mesh file
    std::cout.precision(16);
    std::vector<Vector3r> coordinates;
    std::string line;
    std::string buf;
    getline(in, line);
    while (line[0] == 'v')
    {
        std::stringstream ss(line);
        std::vector<std::string> tokens;
        while (ss >> buf)
            tokens.push_back(buf);
        Vector3r coord = Vector3r::Zero();
        coord[0] = stod(tokens[1]);
        coord[1] = stod(tokens[2]);
        coord[2] = stod(tokens[3]);
        coordinates.push_back(coord);
        getline(in, line);
    };

    // Read connectivity from mesh file and create triangles
    while (line[0] == 'f')
    {
        std::stringstream ss(line);
        std::vector<std::string> tokens;
        while (ss >> buf)
            tokens.push_back(buf);
        Triangle newTri(coordinates[stoi(tokens[1]) - 1], coordinates[stoi(tokens[2]) - 1], coordinates[stoi(tokens[3]) - 1]);
        ExcavationElement newcElem(newTri);
        excavationElements.push_back(newcElem);
        getline(in, line);
    }
    in.close();

    std::cout << "Construction Elements Imported!" << std::endl;
}

void Generator::export_ExcavationElementsVtk(std::string _fileName)
{
    std::ofstream out;
    out.open(_fileName + ".vtk");
    out << std::setprecision(15);
    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "Voronoi results" << std::endl;
    out << "ASCII" << std::endl;
    out << " " << std::endl;
    out << "DATASET POLYDATA" << std::endl;
    out << "POINTS " << 3 * (int)excavationElements.size() << " float" << std::endl;
    for (auto &cElem : excavationElements)
    {
        out << cElem.triangle.pointD[0] << " " << cElem.triangle.pointD[1] << " " << cElem.triangle.pointD[2] << std::endl;
        out << cElem.triangle.pointE[0] << " " << cElem.triangle.pointE[1] << " " << cElem.triangle.pointE[2] << std::endl;
        out << cElem.triangle.pointF[0] << " " << cElem.triangle.pointF[1] << " " << cElem.triangle.pointF[2] << std::endl;
    }
    int count = 0;
    out << " " << std::endl;
    out << "POLYGONS " << (int)excavationElements.size() << " " << 4 * (int)excavationElements.size() << std::endl;
    for (auto cElem : excavationElements)
    {
        out << "3";
        for (int i = 0; i != 3; ++i)
        {
            out << " " << count;
            count++;
        }
        out << std::endl;
    }
    out.close();
}

Generator::ExcavationElement::ExcavationElement(Triangle &_triangle) : triangle(_triangle)
{
    center = (triangle.pointD + triangle.pointE + triangle.pointF) / 3.0;
    unitVector = _triangle.unitVector;
    ASSERT(checkEquality(unitVector.norm(), 1));
    Vector3r borderPlaneUnitVector1 = (triangle.pointD - triangle.pointE).cross(unitVector);
    Vector3r borderPlaneUnitVector2 = (triangle.pointE - triangle.pointF).cross(unitVector);
    Vector3r borderPlaneUnitVector3 = (triangle.pointF - triangle.pointD).cross(unitVector);
    borderPlaneUnitVector1.normalize();
    borderPlaneUnitVector2.normalize();
    borderPlaneUnitVector3.normalize();
    ASSERT(checkEquality(borderPlaneUnitVector1.norm(), 1));
    ASSERT(checkEquality(borderPlaneUnitVector2.norm(), 1));
    ASSERT(checkEquality(borderPlaneUnitVector3.norm(), 1));
    if ((triangle.pointD - center).dot(borderPlaneUnitVector1) < 0)
        borderPlaneUnitVector1 *= -1;
    if ((triangle.pointE - center).dot(borderPlaneUnitVector2) < 0)
        borderPlaneUnitVector2 *= -1;
    if ((triangle.pointF - center).dot(borderPlaneUnitVector3) < 0)
        borderPlaneUnitVector3 *= -1;
    Plane newBorderPlane1(triangle.pointD, borderPlaneUnitVector1);
    Plane newBorderPlane2(triangle.pointE, borderPlaneUnitVector2);
    Plane newBorderPlane3(triangle.pointF, borderPlaneUnitVector3);
    borderPlanes.push_back(newBorderPlane1);
    borderPlanes.push_back(newBorderPlane2);
    borderPlanes.push_back(newBorderPlane3);
    if ((triangle.pointD - center).norm() > boundingSphereRadius)
        boundingSphereRadius = (triangle.pointD - center).norm();
    if ((triangle.pointE - center).norm() > boundingSphereRadius)
        boundingSphereRadius = (triangle.pointE - center).norm();
    if ((triangle.pointF - center).norm() > boundingSphereRadius)
        boundingSphereRadius = (triangle.pointF - center).norm();
}
