
#ifndef POTTSTWOPOPULATIONSGENERATOR_HPP_
#define POTTSTWOPOPULATIONSGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "PottsMesh.hpp"

/**
 * Generator of regular Potts meshes, used as starting points for many simulations.
 *
 * This class takes in options such as width, height,
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
template<unsigned DIM>
class PottsTwoPopulationsGenerator
{
protected:

    /** A pointer to the mesh this class creates */
    PottsMesh<DIM>* mpMesh;

public:

    /**
     * Constructor.
     *
     * @param numNodesAcross  The number of columns of nodes in the mesh
     * @param numElementsAcross  The number of columns of elements in the mesh
     * @param elementWidth  Width of the elements
     * @param numNodesUp  The number of rows of nodes in the mesh (defaults to 1)
     * @param numElementsUp  The number of rows of elements in the mesh (defaults to 1)
     * @param elementHeight  Height of the elements (defaults to 1)
     * @param numNodesDeep  The number nodes deep for this mesh (defaults to 1)
     * @param numElementsDeep  The number of elements deep for this mesh (defaults to 1)
     * @param elementDepth  The number of rows of nodes in each element (defaults to 1)
     * @param startAtBottomLeft  If true then the mesh starts in the bottom left corner
     *     of the domain rather than the centre, used for simple tests (defaults to false)
     * @param isPeriodicInX  If true then the mesh is periodic in the x dimension (defaults to false)
     * @param isPeriodicInY  If true then the mesh is periodic in the y dimension (defaults to false)
     * @param isPeriodicInZ  If true then the mesh is periodic in the y dimension (defaults to false)
     */
    PottsTwoPopulationsGenerator(unsigned numNodesAcross,
                       unsigned numElementsAcross,
                       unsigned elementWidth,
                       unsigned numNodesUp=1u,
                       unsigned numElementsUp=1u,
                       unsigned elementHeight=1u,
                       unsigned numNodesDeep=1u,
                       unsigned numElementsDeep=1u,
                       unsigned elementDepth=1u,
                       bool startAtBottomLeft = false,
                       bool isPeriodicInX = false,
                       bool isPeriodicInY = false,
                       bool isPeriodicInZ = false);

    /**
     * Destructor - deletes the mesh object and pointer.
     */
    virtual ~PottsTwoPopulationsGenerator();

    /**
     * Helper method to calculate the Moore and Von Neumann Neighbourhoods of all nodes
     *
     * @param isPeriodicInX  If true then the mesh is periodic in the x dimension
     */
    void CaclulateNeighbouringNodeIndices(bool isPeriodicInX);

    /**
     * @return a Cuboid or rectangular Potts mesh.
     */
    virtual PottsMesh<DIM>* GetMesh();
};

#endif /*POTTSTWOPOPULATIONSGENERATOR_HPP_*/
