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

#ifndef GENERATOR_H
#define GENERATOR_H

#include "maths/Math.hpp"
#include "Block.hpp"
#include <boost/python.hpp>

class DFN;
class Fracture;

class Generator
{
public:
    Generator(){};
    ~Generator(){};

    struct ExcavationElement
    {
        ExcavationElement(Triangle &_triangle);
        ~ExcavationElement(){};
        std::vector<Plane> borderPlanes;
        double boundingSphereRadius;
        Triangle triangle;
        Vector3r unitVector;
        Vector3r center;
    };

    void set_MinInscribedSphereRadius(double _minInscribedSphereRadius) { minInscribedSphereRadius = _minInscribedSphereRadius; };
    void set_MaxAspectRatio(double _maxAspectRatio) { maxAspectRatio = _maxAspectRatio; };
    void import_ExcavationElementsObj(std::string _fileName);
    void export_ExcavationElementsVtk(std::string _fileName);
    void export_BlocksVtk(std::string _fileName);

    /////////////////////////////////////////////////////////
    /// DDA: Added, extract DDA blocks
    typedef boost::python::list PyList;
#define TOLA 0.5  // tolerance for zero angle
#define TOLD 1e-6 // tolerance for zero dist
    void export_BlocksDDA(std::string _fileName);
    void export_BlocksDDAOpt(std::string _filename);

// Added, bond constraints
#define NBONDP 4 //# of bonded point
    // return false if do not find opposite faces
    // bool get_BondInfo(Block& b1, Block& b2,Bond& bond);
    vecInt vBondIndex;                              // to store bond indexes
    void input_BondedBlocks(PyList _bondedIndexes); // get bonded indexes

    // Added, viscous boundaries
    vecInt vVisIndex; // indexes of blocks that have viscous boudnaries
    vecP vVisNormal;  // normal direction of viscous boundaries
    // get viscous indexes and normals, do not need to be unit normals
    void input_ViscousBound(int _visBlkId, PyList _visNor);

    // Added, input boundaries
    vecInt vInputIndex;
    vecP vInputNormal;
    void input_InputBound(int _inputBlkId, PyList _inputNor);

    // Added, roller boundaries
    vecInt vRollerIndex;
    vecP vRollerNormal;
    // get roller indexes and normals, do not need to be unit normals
    void input_RollerBound(int _rollerBlkId, PyList _rollerNor);

    // Added, fixedFace with springs
    vecInt vFixFaceIndex;
    vecP vFixFaceNormal;
    void input_FixFaceBound(int _fixFaceBlkId, PyList _fixFaceNor);
    /////////////////////////////////////////////////////////

    boost::python::list get_Volumes(bool _considerBorderBlocks);
    boost::python::list get_AlphaValues(bool _considerBorderBlocks);
    boost::python::list get_BetaValues(bool _considerBorderBlocks);

    void generate_RockMass(DFN &_dfn);

    /////////////////////////////////////////////////////////
    /// DDA: Added, extract DDA blocks
    // add generate rock mass seperately
    void generate_RockMass_Multi(DFN &_dfn);
    void generate_Polyhedra(DFN &_dfn);

    // add fixed box region
    std::vector<Box> fixRegion;

    // add fixed box region
    void add_Fixed_Region(PyList _MinPoint, PyList _MaxPoint);
    inline bool inBox(Box &box, Vector3r &p, Vector3r &offset)
    {
        return box.minCor[0] < p[0] + offset[0] && box.minCor[1] < p[1] + offset[1] && box.minCor[2] < p[2] + offset[2] &&
               box.maxCor[0] > p[0] + offset[0] && box.maxCor[1] > p[1] + offset[1] && box.maxCor[2] > p[2] + offset[2];
    }
    /////////////////////////////////////////////////////////

    void excavate_RockMass();

    std::vector<std::shared_ptr<Block>> blocks;

    /////////////////////////////////////////////////////////
    /// DDA: added
    std::vector<int> blockNs;
    /////////////////////////////////////////////////////////

private:
    double calculate_BoundingSphereRadius(const Block &_block);
    double calculate_InscribedSphere(const Block &_block);
    template <class X>
    bool check_BlockPlaneIntersection(const Block &_block, const X &_plane);

    std::vector<ExcavationElement> excavationElements;
    Vector3r regionMinCorner = {0, 0, 0};
    Vector3r regionMaxCorner = {100, 100, 100};

    double minInscribedSphereRadius = 0;
    double maxAspectRatio = INFINITY;
};

#endif // GENERATOR_H
