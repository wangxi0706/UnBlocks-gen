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

#include "Block.hpp"
#include "coin/ClpSimplex.hpp"
#include <boost/python.hpp>

using namespace boost::python;
void Block_Wrapper()
{
    class_<Block, boost::noncopyable>("Block")
        .def("export_BlockVtk", &Block::export_BlockVtk)

        .def("get_Volume", &Block::get_Volume)
        .def("get_Alpha", &Block::get_Alpha)
        .def("get_Beta", &Block::get_Beta)
        .def("get_InscribedSphereRadius", &Block::get_InscribedSphereRadius)
        .def("get_AspectRatio", &Block::get_AspectRatio)
        .def("get_Order", &Block::get_Order);
};

double Block::get_Alpha()
{
    if (alpha > 10)
        return 10;
    if (alpha < 1)
        return 1;
    return alpha;
}

double Block::get_Beta()
{
    if (beta > 10)
        return 10;
    if (beta < 1)
        return 1;
    return beta;
}

void Block::generate_Geometry()
{
    Functions::eliminate_RedundantPlanes(planes);

    polygons.clear();
    edges.clear();
    vertices.clear();
    std::vector<std::vector<int>> planesVerticesId;

    int intCountCheck = 0;
    bool runAnew = true;
    while (runAnew)
    {
        ASSERT(intCountCheck < 100);
        intCountCheck++;
        vertices.clear();
        planesVerticesId.clear();
        for (int i = 0; i != (int)planes.size(); ++i)
        {
            std::vector<int> verticesId;
            planesVerticesId.push_back(verticesId);
            // printf("Plance center: %f, %f\n", planes[i].center.begin(), planes[i].center.end());
        }

        //Fill vertices vector
        if ((int)planes.size() < 4)
        {
            printf("plane size=%d\n", (int)planes.size());
        }
        ASSERT((int)planes.size() >= 4);
        for (int i = 0; i != (int)planes.size() - 2; ++i)
        {
            for (int j = i + 1; j != (int)planes.size() - 1; ++j)
            {
                for (int k = j + 1; k != (int)planes.size(); ++k)
                {
                    Vector3r intersectionPoint = Vector3r::Zero();
                    if (Functions::calculate_threePlanesIntersection(planes[i], planes[j], planes[k], intersectionPoint))
                    {
                        if (Functions::check_PointIntersection(intersectionPoint, planes))//if inside the halfspaces defined by block planes
                        {
                            bool verticeExist = false;
                            for (int l = 0; l != (int)vertices.size(); ++l)
                            {
                                if (checkEquality((intersectionPoint - vertices[l]).norm(), 0))
                                {
                                    bool addVertId = true;
                                    for (auto &id : planesVerticesId[i])
                                        if (id == l)
                                            addVertId = false;
                                    if (addVertId)
                                        planesVerticesId[i].push_back(l);
                                    addVertId = true;
                                    for (auto &id : planesVerticesId[j])
                                        if (id == l)
                                            addVertId = false;
                                    if (addVertId)
                                        planesVerticesId[j].push_back(l);
                                    addVertId = true;
                                    for (auto &id : planesVerticesId[k])
                                        if (id == l)
                                            addVertId = false;
                                    if (addVertId)
                                        planesVerticesId[k].push_back(l);
                                    verticeExist = true;
                                }
                            }
                            if (!verticeExist)
                            {
                                planesVerticesId[i].push_back((int)vertices.size());
                                planesVerticesId[j].push_back((int)vertices.size());
                                planesVerticesId[k].push_back((int)vertices.size());
                                vertices.push_back(intersectionPoint);
                            }
                        }
                    }
                }
            }
        }

        bool allPlanesOk = true;
        ASSERT((int)vertices.size() >= 4);
        for (int i = 0; i != (int)planes.size(); ++i)
        {
            bool isThisPlaneOk = false;
            if ((int)planesVerticesId[i].size() < 3)
                isThisPlaneOk = false;
            if ((int)planesVerticesId[i].size() >= 3)
            {
                for (int j = 1; j != (int)planesVerticesId[i].size() - 1; ++j)
                {
                    for (int k = j + 1; k != (int)planesVerticesId[i].size(); ++k)
                    {
                        double area = 0.5 * ((vertices[planesVerticesId[i][j]] - vertices[planesVerticesId[i][0]]).cross(vertices[planesVerticesId[i][k]] - vertices[planesVerticesId[i][0]])).norm();
                        if (!checkEquality(area, 0))
                            isThisPlaneOk = true;
                    }
                }
            }
            if (!isThisPlaneOk)
            {
                planes[i].d = NAN;
                allPlanesOk = false;
            }
        }
        if (!allPlanesOk)
        {
            planes.erase(std::remove_if(std::begin(planes), std::end(planes), [](Plane &_p) { return (std::isnan(_p.d)); }), std::end(planes));
            runAnew = true;
        }
        else
        {
            runAnew = false;
        }
    }

    // added, calculate the center, set all plane vectors as outer normal
    Vector3r center=Vector3r::Zero();
    for(auto v : vertices){
        center[0]+=v[0];center[1]+=v[1];center[2]+=v[2];
    }
    center[0]/=vertices.size();center[1]/=vertices.size();center[2]/=vertices.size();

    //Fill polygons vector
    for (int i = 0; i != (int)planes.size(); ++i)
    {
        // added, to make sure all plane vectors points the outer normal
        if(vertices[planesVerticesId[i][0]].dot(planes[i].unitVector)-
            center.dot(planes[i].unitVector)<0)
        {
            planes[i].unitVector[0]*=-1;
            planes[i].unitVector[1]*=-1;
            planes[i].unitVector[2]*=-1;
        }
        ASSERT((int)planesVerticesId[i].size() >= 3);
        Polygon newPolygon(planes[i].unitVector, planesVerticesId[i], vertices);
        polygons.push_back(newPolygon);
    }

    //Fill edges vector
    for (int i = 0; i != (int)polygons.size(); ++i)
    {
        for (int j = 0; j != (int)polygons[i].verticesId.size(); ++j)
        {
            int auxId1 = j;
            int auxId2 = j + 1;
            if (j == (int)polygons[i].verticesId.size() - 1)
                auxId2 = 0;
            bool isNewEdge = true;
            Vector3r VertBefore = vertices[polygons[i].verticesId[auxId1]];
            Vector3r VertAfter = vertices[polygons[i].verticesId[auxId2]];
            for (auto &edg : edges)
            {
                Vector3r edgVertA = vertices[edg.verticesIdA];
                Vector3r edgVertB = vertices[edg.verticesIdB];
                if ((VertBefore == edgVertA && VertAfter == edgVertB) || (VertBefore == edgVertB && VertAfter == edgVertA))
                {
                    isNewEdge = false;
                }
            }
            if (isNewEdge)
            {
                Edge newEdge(polygons[i].verticesId[auxId1], polygons[i].verticesId[auxId2]);
                edges.push_back(newEdge);
            }
        }
    }

    //Check Euler relation for polyhedron
    ASSERT((int)vertices.size() - (int)edges.size() + (int)polygons.size() == 2);

    //Calculate Volume
    calculate_Volume();

    //Calculate alpha and beta for shape characterization
    calculate_AlphaBeta();

    std::cout << "Generated Block " << id << " Geometry" << std::endl;
    // printf("Block v0: %f, %f, %f\n",
    //        vertices[0][0], vertices[0][1], vertices[0][2]);
}

void Block::calculate_Volume()
{
    volume = 0;
    for (auto &p : polygons)
    {
        ASSERT(checkEquality(p.unitVector.norm(), 1));
        volume += (p.centroid.dot(p.unitVector) * p.area) / 3.0;
    }
    volume = std::abs(volume);
}

void Block::calculate_AlphaBeta()
{
    std::vector<Vector3r> chords;
    std::vector<double> chordsLength;
    double averageChord = 0;
    for (int i = 0; i != (int)vertices.size() - 1; ++i)
    {
        for (int j = i + 1; j != (int)vertices.size(); ++j)
        {
            Vector3r newChord = vertices[i] - vertices[j];
            chords.push_back(newChord);
            chordsLength.push_back(newChord.norm());
            averageChord += newChord.norm();
        }
    }
    averageChord = averageChord / (double)chordsLength.size();

    //calculate median chord length
    double medianChord = 0;
    std::sort(chordsLength.begin(), chordsLength.end());
    if ((int)chordsLength.size() % 2 != 0)
        medianChord = chordsLength[((int)chordsLength.size() + 1) / 2 - 1];
    if ((int)chordsLength.size() % 2 == 0)
        medianChord = 0.5 * ((chordsLength[(int)chordsLength.size() / 2 - 1]) + (chordsLength[(int)chordsLength.size() / 2]));

    //erase chords with lengths smaller than the median value
    chords.erase(std::remove_if(std::begin(chords), std::end(chords), [&](Vector3r &_c) { return (_c.norm() < medianChord); }), std::end(chords));

    //calculate Beta
    double sumA = 0;
    double sumB = 0;
    for (int i = 0; i != (int)chords.size() - 1; ++i)
    {
        for (int j = i + 1; j != (int)chords.size(); ++j)
        {
            sumA += std::pow(chords[i].dot(chords[j]), 2);
            sumB += std::pow(chords[i].norm(), 2) * std::pow(chords[j].norm(), 2);
        }
    }
    beta = 10 * std::pow(sumA / sumB, 2);

    //calculate Alpha
    double surfaceArea = 0;
    for (auto &p : polygons)
        surfaceArea += p.area;
    alpha = surfaceArea * averageChord / (7.7 * volume);
}

void Block::export_BlockVtk(std::string _fileName)
{
    int nFaces = (int)polygons.size();
    int nVerts = (int)vertices.size();
    int nTotalPoints = 0;
    for (auto p : polygons)
        nTotalPoints += (int)p.verticesId.size();

    std::ofstream out;
    out.open(_fileName + ".vtk");
    out << std::setprecision(15);
    out << "# vtk DataFile Version 3.0" << std::endl;
    out << "Voronoi results" << std::endl;
    out << "ASCII" << std::endl;
    out << " " << std::endl;
    out << "DATASET UNSTRUCTURED_GRID" << std::endl;
    out << "POINTS " << std::to_string(nVerts) << " float" << std::endl;
    for (auto &v : vertices)
    {
        out << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }

    out << " " << std::endl;
    out << "CELLS " << 1 << " " << std::to_string((2 + nFaces + nTotalPoints)) << std::endl;
    out << 1 + nFaces + nTotalPoints << " " << nFaces << " ";
    for (auto &p : polygons)
    {
        out << (int)p.verticesId.size() << " ";
        for (auto &vId : p.verticesId)
        {
            out << vId << " ";
        }
    }
    out << std::endl;

    out << " " << std::endl;
    out << "CELL_TYPES " << 1 << std::endl;
    ;
    out << 42 << std::endl;

    out << " " << std::endl;
    out << "CELL_DATA " << 1 << std::endl;
    out << "SCALARS id float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    out << id << std::endl;

    out << " " << std::endl;
    out << "SCALARS volume float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    out << get_Volume() << std::endl;

    out << " " << std::endl;
    out << "SCALARS alpha float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    out << get_Alpha() << std::endl;

    out << " " << std::endl;
    out << "SCALARS beta float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    out << get_Beta() << std::endl;

    out << " " << std::endl;
    out << "SCALARS order float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    out << get_Order() << std::endl;

    out << " " << std::endl;
    out << "SCALARS aspectRatio float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    out << get_AspectRatio() << std::endl;

    out << " " << std::endl;
    out << "SCALARS inscribedSphereRadius float" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    out << get_InscribedSphereRadius() << std::endl;

    out.close();
};

void Block::calculate_BoundingSphere()
{
    double extremeCoordinates[6];
    for (int auxi = 0; auxi != 6; ++auxi)
    {
        ClpSimplex model;
        model.setLogLevel(0);
        model.setPrimalTolerance(ND_tolerance * 0.1);
        model.setOptimizationDirection(-1);

        int numberColumns = 3;
        int numberRows = (int)planes.size();
        int numberElements = numberRows * numberColumns;
        model.resize(numberRows, numberColumns);

        double *objective = model.objective();
        double *columnUpper = model.columnUpper();
        double *columnLower = model.columnLower();
        double *rowLower = model.rowLower();
        double *rowUpper = model.rowUpper();
        objective[0] = 0;
        objective[1] = 0;
        objective[2] = 0;
        if (auxi == 0)
            objective[0] = 1;
        if (auxi == 1)
            objective[0] = -1;
        if (auxi == 2)
            objective[1] = 1;
        if (auxi == 3)
            objective[1] = -1;
        if (auxi == 4)
            objective[2] = 1;
        if (auxi == 5)
            objective[2] = -1;
        for (int i = 0; i < numberColumns; i++)
        {
            columnUpper[i] = COIN_DBL_MAX;
            columnLower[i] = -COIN_DBL_MAX;
        }
        for (int i = 0; i < numberRows; i++)
        {
            rowUpper[i] = cZ(planes[i].d);
            rowLower[i] = -COIN_DBL_MAX;
        }

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
                elements[put] = cZ(planes[i].unitVector[k]);
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
        ASSERT(model.status() == 0); //???
        double *columnPrimal = model.primalColumnSolution();

        if (auxi == 0 || auxi == 1)
            extremeCoordinates[auxi] = columnPrimal[0];
        if (auxi == 2 || auxi == 3)
            extremeCoordinates[auxi] = columnPrimal[1];
        if (auxi == 4 || auxi == 5)
            extremeCoordinates[auxi] = columnPrimal[2];
    }

    Vector3r extremePosCoordinate = {extremeCoordinates[0], extremeCoordinates[2], extremeCoordinates[4]};
    Vector3r extremeNegCoordinate = {extremeCoordinates[1], extremeCoordinates[3], extremeCoordinates[5]};
    boundingSphereRadius = 0.5 * (extremePosCoordinate - extremeNegCoordinate).norm();
    boundingSphereCenter = 0.5 * (extremePosCoordinate + extremeNegCoordinate);
}

void Block::calculate_InscribedSphere()
{
    ClpSimplex model;
    model.setLogLevel(0);
    model.setPrimalTolerance(ND_tolerance * 0.1);
    model.setOptimizationDirection(1);

    int numberColumns = 4;
    int numberRows = (int)planes.size();
    int numberElements = numberRows * numberColumns;
    model.resize(numberRows, numberColumns);

    double *objective = model.objective();
    double *columnUpper = model.columnUpper();
    double *columnLower = model.columnLower();
    double *rowLower = model.rowLower();
    double *rowUpper = model.rowUpper();
    objective[0] = 0;
    objective[1] = 0;
    objective[2] = 0;
    objective[3] = 1;
    for (int i = 0; i < numberColumns; i++)
    {
        columnUpper[i] = COIN_DBL_MAX;
        columnLower[i] = -COIN_DBL_MAX;
    }
    for (int i = 0; i < numberRows; i++)
    {
        rowUpper[i] = cZ(planes[i].d);
        rowLower[i] = -COIN_DBL_MAX;
    }

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
            if (k < numberColumns - 1)
                elements[put] = cZ(planes[i].unitVector[k]);
            if (k == numberColumns - 1)
                elements[put] = -1;
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
    ASSERT(model.status() == 0); //optimization
    double *columnPrimal = model.primalColumnSolution();

    ASSERT(columnPrimal[3] < 0);
    inscribedSphereCenter[0] = columnPrimal[0];
    inscribedSphereCenter[1] = columnPrimal[1];
    inscribedSphereCenter[2] = columnPrimal[2];
    inscribedSphereRadius = -columnPrimal[3];
}
