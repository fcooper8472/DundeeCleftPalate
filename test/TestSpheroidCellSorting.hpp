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

#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellLabelWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

static const double M_TIME_TO_STEADY_STATE = 10; //10
static const double M_TIME_FOR_SIMULATION = 100; //100
static const double M_NUM_CELLS_ACROSS = 20; //20 // this ^2 cells
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

public:

    void TestPottsMonolayerCellSorting() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        unsigned element_size = 4;
        unsigned domain_size = M_NUM_CELLS_ACROSS * element_size * 3; // Three times the initial domain size
        PottsTwoPopulationsGenerator<2> generator(domain_size, M_NUM_CELLS_ACROSS, element_size, domain_size, M_NUM_CELLS_ACROSS, element_size);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellLabelWriter>();
        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();
        cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();

        // Set the Temperature
        cell_population.SetTemperature(0.2); //Default is 0.1

        // Set up cell-based simulation and output directory
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellSorting/Potts");

        // Set time step and end time for simulation
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE);

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

        // Run simulation
        simulator.Solve();

        // Now label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);

        // Adjust Parameters
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

};
