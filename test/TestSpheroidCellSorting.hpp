#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "HeterotypicBoundaryLengthWriter.hpp"

#include "OnLatticeSimulation.hpp"
#include "CellPopulationAdjacencyMatrixWriter.hpp"

#include "PottsBasedCellPopulation.hpp"
#include "PottsTwoPopulationsGenerator.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "GeneralisedChemotaxisPottsUpdateRule.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "RandomMotionForce.hpp"
#include "OffLatticeSimulation.hpp"

#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellLabelWriter.hpp"

#include "RadialCellKiller.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"

static const double M_TIME_TO_STEADY_STATE = 0.2; //10
static const double M_TIME_FOR_SIMULATION = 30.0; //100
static const double M_NUM_CELLS_ACROSS = 60; //20 // this ^2 cells
static const double M_CELL_FLUCTUATION = 1.0;

class TestSpheroidCellSorting : public AbstractCellBasedWithTimingsTestSuite
{
private:

    void RandomlyLabelCells(std::list<CellPtr>& rCells, boost::shared_ptr<AbstractCellProperty> pLabel, double labelledRatio)
    {
        for (std::list<CellPtr>::iterator cell_iter = rCells.begin();
             cell_iter != rCells.end();
             ++cell_iter)
        {
            if (RandomNumberGenerator::Instance()->ranf() < labelledRatio)
            {
                (*cell_iter)->AddCellProperty(pLabel);
            }
        }
    }

    void LabelCellsInConcentricCircles(AbstractCellPopulation<2>& rPopulation,
    		                           boost::shared_ptr<AbstractCellProperty> pLabel,
									   double R,
									   double x1,
									   double y1,
									   double x2,
									   double y2)
    {
        for (AbstractCellPopulation<2>::Iterator cell_iter = rPopulation.Begin();
             cell_iter != rPopulation.End();
             ++cell_iter)
        {
            c_vector<double,2> cell_centre = rPopulation.GetLocationOfCellCentre(*cell_iter);

            double dist_to_centre1_squared = pow(cell_centre[0]-x1,2) + pow(cell_centre[1]-y1,2);
            double dist_to_centre2_squared = pow(cell_centre[0]-x2,2) + pow(cell_centre[1]-y2,2);

            if (!(dist_to_centre1_squared < pow(R,2) || dist_to_centre2_squared < pow(R,2)))
            {
                (*cell_iter)->AddCellProperty(pLabel);
            }
        }
    }

public:

    void NoTestPottsMonolayerCellSorting() throw (Exception)
    {
        double R = 40;
        double x1 = 80;
        double y1 = 160;
        double x2 = 160;
        double y2 = 80;

        // Create a simple 2D PottsMesh
        unsigned element_size = 4;
        unsigned domain_size = M_NUM_CELLS_ACROSS * element_size;
        PottsTwoPopulationsGenerator<2> generator(domain_size, M_NUM_CELLS_ACROSS, element_size, domain_size, M_NUM_CELLS_ACROSS, element_size);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        for (std::vector<CellPtr>::iterator cell_iter = cells.begin();
             cell_iter != cells.end();
             ++cell_iter)
        {
            (*cell_iter)->SetApoptosisTime(0.02);
        }

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellLabelWriter>();

        // Set the temperature
        cell_population.SetTemperature(0.2); // Default is 0.1

        // Set up cell-based simulation and output directory
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellSorting/Potts");

        // Set time step and end time for simulation
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(25);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16); // i.e 4x4 cells
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.1);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_constraint_update_rule);
        p_surface_constraint_update_rule->SetMatureCellTargetSurfaceArea(16); // i.e 4x4 cells
        p_surface_constraint_update_rule->SetDeformationEnergyParameter(0.01);//0.01
        simulator.AddUpdateRule(p_surface_constraint_update_rule);

        MAKE_PTR(DifferentialAdhesionPottsUpdateRule<2>, p_differential_adhesion_update_rule);
        p_differential_adhesion_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.1);
        p_differential_adhesion_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.5); // 1.0
        p_differential_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.1); //0.1
        p_differential_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.2); // 1.0
        p_differential_adhesion_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(1.0); // 2.0
        simulator.AddUpdateRule(p_differential_adhesion_update_rule);

        c_vector<double,2> chemotactic_centre;
        chemotactic_centre[0] = domain_size/2.0;
        chemotactic_centre[1] = domain_size/2.0;
        MAKE_PTR_ARGS(GeneralisedChemotaxisPottsUpdateRule<2>, p_chemotaxis_update_rule, (chemotactic_centre));
        simulator.AddUpdateRule(p_chemotaxis_update_rule);

        MAKE_PTR_ARGS(RadialCellKiller<2>, p_killer, (&cell_population,0.5));
        simulator.AddCellKiller(p_killer);

        // Run simulation for one timestep just to remove cells outside the spheroid
        simulator.SetEndTime(0.01);
        simulator.Solve();
        simulator.RemoveAllCellKillers();

        // Run to get to steady state for a while
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE);
        simulator.Solve();

        // Now label some cells as epithelial
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        LabelCellsInConcentricCircles(cell_population, p_state, R-10, x1, y1, x2, y2);

        // Adjust parameters
        dynamic_cast <PottsBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation()))->SetTemperature(0.2*M_CELL_FLUCTUATION);

        // Run simulation
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), M_NUM_CELLS_ACROSS*M_NUM_CELLS_ACROSS);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }

    void TestNodeBasedMonolayerCellSorting() throw (Exception)
    {
        // Create a simple mesh
        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Extended to allow sorting for longer distances
        double cut_off_length = 2.5;

        // Delete any nodes outside two circles of specified centres and radii
        double R = 10.0;
        double x1 = 20.0;
        double y1 = 40.0;
        double x2 = 40.0;
        double y2 = 20.0;
        double dist_to_centre1_squared = 0.0;
        double dist_to_centre2_squared = 0.0;

        std::vector<Node<2>*> temp_nodes;
        for (AbstractMesh<2,2>::NodeIterator node_iter = p_generating_mesh->GetNodeIteratorBegin();
             node_iter != p_generating_mesh->GetNodeIteratorEnd();
             ++node_iter)
        {
            c_vector<double,2> location = node_iter->rGetLocation();

            dist_to_centre1_squared = pow(location[0]-x1,2) + pow(location[1]-y1,2);
            dist_to_centre2_squared = pow(location[0]-x2,2) + pow(location[1]-y2,2);
            if (dist_to_centre1_squared < pow(R,2) || dist_to_centre2_squared < pow(R,2))
            {
                temp_nodes.push_back(new Node<2>(node_iter->GetIndex(), location, node_iter->IsBoundaryNode()));
            }
        }
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(temp_nodes, cut_off_length);

        // Set up cells, one for each Node
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellLabelWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellSorting/Node");

        // Set time step and end time for simulation
//        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(50);

        // Create a force law and pass it to the simulation
        MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
        p_differential_adhesion_force->SetMeinekeSpringStiffness(50.0);
        p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1);
        p_differential_adhesion_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_differential_adhesion_force);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.05); // 0.1 causes dissasociation, 0.001 is not enough
        simulator.AddForce(p_random_force);

        // Run to get to steady state for a while
        simulator.SetEndTime(1.0);
        simulator.Solve();

        // Now label some cells as epithelial
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        LabelCellsInConcentricCircles(cell_population, p_state, R-2.5, x1, y1, x2, y2);

        // Adjust parameters
        p_random_force->SetMovementParameter(0.05*M_CELL_FLUCTUATION); // 0.1 causes dissociation

        // Run simulation
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);
        simulator.Solve();
    }
};
