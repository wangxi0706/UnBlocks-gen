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

    //First block generation
    Plane planeX0({-1, 0, 0}, _dfn.regionMinCorner[0]);
    Plane planeX1({1, 0, 0}, _dfn.regionMaxCorner[0]);
    Plane planeY0({0, -1, 0}, _dfn.regionMinCorner[1]);
    Plane planeY1({0, 1, 0}, _dfn.regionMaxCorner[1]);
    Plane planeZ0({0, 0, -1}, _dfn.regionMinCorner[2]);
    Plane planeZ1({0, 0, 1}, _dfn.regionMaxCorner[2]);
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
            for (int i = 0; i != nBlocks; ++i)
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

    for (int i = 0; i != (int)blocks.size(); ++i)
        blocks[i]->generate_Geometry();
}

void Generator::generate_RockMass_Multi(DFN &_dfn)
{
    int startIndex = (int)blocks.size();

    regionMaxCorner = _dfn.regionMaxCorner;
    regionMinCorner = _dfn.regionMinCorner;
    //First block generation
    // Plane planeX0({-1, 0, 0}, _dfn.firstBlkMin[0]);
    // Plane planeX1({1, 0, 0}, _dfn.firstBlkMax[0]);
    // Plane planeY0({0, -1, 0}, _dfn.firstBlkMin[1]);
    // Plane planeY1({0, 1, 0}, _dfn.firstBlkMax[1]);
    // Plane planeZ0({0, 0, -1}, _dfn.firstBlkMin[2]);
    // Plane planeZ1({0, 0, 1}, _dfn.firstBlkMax[2]);
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
    blocks.erase(std::remove_if(std::begin(blocks), std::end(blocks), [](std::shared_ptr<Block> &_b) { return (_b->id == -1); }), std::end(blocks));
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
            out << v[0] << " " << v[1] << " " << v[2] << std::endl;
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
        }
    }

    std::ofstream out;
    out.open(_fileName + ".dda");
    out << std::setprecision(15);
    out << "#_DDA_DataFile_Version_1.0" << std::endl;

    out << "POINTS_START_LENGTH " << std::to_string(nTotalBlocks) << std::endl;
    int auxId = 0;
    for (auto &b : blocks)
    {
        //output vertex start id and vertex length
        out << auxId << " " << (int)b->vertices.size() << std::endl;
        auxId += (int)b->vertices.size();
    }

    out << " " << std::endl;
    out << "POINTS " << std::to_string(nTotalVerts) << " double" << std::endl;
    for (auto &b : blocks)
    {
        for (auto &v : b->vertices)
        {
            out << v[0] << " " << v[1] << " " << v[2] << std::endl;
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
            for (auto &vId : p.verticesId)
            {
                out << vId << " "; //local id, global id=local id + start id
            }
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
            out << e.verticesIdA << " " << e.verticesIdB << " "; //local id, global id=local id + start id
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
        for (auto &box : fixRegion)
        {
            for (auto &v : b->vertices)
            {
                if (!inBox(box, v))
                {
                    flag = 0; //not in box, so free
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

    out.close();
}

void Generator::import_ExcavationElementsObj(std::string _fileName)
{
    std::ifstream in;
    in.open(_fileName + ".obj");
    ASSERT(!in.fail());

    //Read coordinates from mesh file
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

    //Read connectivity from mesh file and create triangles
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
