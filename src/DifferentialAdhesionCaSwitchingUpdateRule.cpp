
#include "DifferentialAdhesionCaSwitchingUpdateRule.hpp"

template<unsigned DIM>
DifferentialAdhesionCaSwitchingUpdateRule<DIM>::DifferentialAdhesionCaSwitchingUpdateRule()
    : AbstractCaSwitchingUpdateRule<DIM>(),
      mCellCellAdhesionEnergyParameter(0.2), // Educated guess
      mLabelledCellLabelledCellAdhesionEnergyParameter(0.2), // Educated guess
      mLabelledCellCellAdhesionEnergyParameter(0.1), // Educated guess
      mCellBoundaryAdhesionEnergyParameter(0.2), // Educated guess
      mLabelledCellBoundaryAdhesionEnergyParameter(0.2), // Educated guess
      mTemperature(0.1) // Educated guess
{
}

template<unsigned DIM>
DifferentialAdhesionCaSwitchingUpdateRule<DIM>::~DifferentialAdhesionCaSwitchingUpdateRule()
{
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::EvaluateHamiltonian(
    unsigned currentNodeIndex,
    unsigned neighbourNodeIndex,
    CaBasedCellPopulation<DIM>& rCellPopulation)
{
    // Energy before and after switch
    double H_0 = 0.0;
    double H_1 = 0.0;

    bool is_cell_on_node_1 = rCellPopulation.IsCellAttachedToLocationIndex(currentNodeIndex);
    bool is_cell_on_node_2 = rCellPopulation.IsCellAttachedToLocationIndex(neighbourNodeIndex);

    bool is_cell_1_labelled = false;
    if (is_cell_on_node_1)
    {
        CellPtr p_cell_1 = rCellPopulation.GetCellUsingLocationIndex(currentNodeIndex);
        is_cell_1_labelled = p_cell_1->HasCellProperty<CellLabel>();
    }
    bool is_cell_2_labelled = false;
    if (is_cell_on_node_2)
    {
        CellPtr p_cell_2 = rCellPopulation.GetCellUsingLocationIndex(neighbourNodeIndex);
        is_cell_2_labelled = p_cell_2->HasCellProperty<CellLabel>();
    }

    std::set<unsigned> node_1_neighbouring_node_indices = rCellPopulation.rGetMesh().GetVonNeumannNeighbouringNodeIndices(currentNodeIndex);

    // Remove node 2 from neighbouring indices
    node_1_neighbouring_node_indices.erase(neighbourNodeIndex);

    for (std::set<unsigned>::iterator iter = node_1_neighbouring_node_indices.begin();
         iter != node_1_neighbouring_node_indices.end();
         ++iter)
    {
        if (rCellPopulation.IsCellAttachedToLocationIndex(*iter))
        {
            // Cell is attached to neighbour node
            bool is_neighbour_labelled = rCellPopulation.GetCellUsingLocationIndex(*iter)->template HasCellProperty<CellLabel>();

            if (is_neighbour_labelled)
            {
                if (is_cell_on_node_1)
                {
                    if (is_cell_1_labelled)
                    {
                        H_0 += mLabelledCellLabelledCellAdhesionEnergyParameter;
                    }
                    else
                    {
                        H_0 += mLabelledCellCellAdhesionEnergyParameter;
                    }
                }
                else // no cell on node 1
                {
                    H_0 += mLabelledCellBoundaryAdhesionEnergyParameter;
                }

                if (is_cell_on_node_2)
                {
                    if (is_cell_2_labelled)
                    {
                        H_1 += mLabelledCellLabelledCellAdhesionEnergyParameter;
                    }
                    else
                    {
                        H_1 += mLabelledCellCellAdhesionEnergyParameter;
                    }
                }
                else // no cell on node 2
                {
                    H_1 += mLabelledCellBoundaryAdhesionEnergyParameter;
                }
            }
            else // Neighbour not labelled
            {
                if (is_cell_on_node_1)
                {
                    if (is_cell_1_labelled)
                    {
                        H_0 += mLabelledCellCellAdhesionEnergyParameter;
                    }
                    else
                    {
                        H_0 += mCellCellAdhesionEnergyParameter;
                    }
                }
                else // no cell on node 1
                {
                    H_0 += mCellBoundaryAdhesionEnergyParameter;
                }
                if (is_cell_on_node_2)
                {
                    if (is_cell_2_labelled)
                    {
                        H_1 += mLabelledCellCellAdhesionEnergyParameter;
                    }
                    else
                    {
                        H_1 += mCellCellAdhesionEnergyParameter;
                    }
                }
                else // no cell on node 2
                {
                    H_1 += mCellBoundaryAdhesionEnergyParameter;
                }
            }
        }
        else // No cell on neighbour node
        {
            if (is_cell_on_node_1)
            {
                if (is_cell_1_labelled)
                {
                    H_0 += mLabelledCellBoundaryAdhesionEnergyParameter;
                }
                else
                {
                    H_0 += mCellBoundaryAdhesionEnergyParameter;
                }
            }
            // else // no cell on node 1 or neighbour so no contribution
            if (is_cell_on_node_2)
            {
                if (is_cell_2_labelled)
                {
                    H_1 += mLabelledCellBoundaryAdhesionEnergyParameter;
                }
                else
                {
                    H_1 += mCellBoundaryAdhesionEnergyParameter;
                }
            }
            //else // no cell on node 2 or neighbour so no contribution
        }
    }

    std::set<unsigned> node_2_neighbouring_node_indices = rCellPopulation.rGetMesh().GetVonNeumannNeighbouringNodeIndices(neighbourNodeIndex);

    // Remove node 2 from neighbouring indices
    node_2_neighbouring_node_indices.erase(currentNodeIndex);

    for (std::set<unsigned>::iterator iter = node_2_neighbouring_node_indices.begin();
         iter != node_2_neighbouring_node_indices.end();
         ++iter)
    {
        if (rCellPopulation.IsCellAttachedToLocationIndex(*iter))
        {
            // Cell is attached to neighbour node
            bool is_neighbour_labelled = rCellPopulation.GetCellUsingLocationIndex(*iter)->template HasCellProperty<CellLabel>();

            if (is_neighbour_labelled)
            {
                if (is_cell_on_node_1)
                {
                    if (is_cell_1_labelled)
                    {
                        H_1 += mLabelledCellLabelledCellAdhesionEnergyParameter;
                    }
                    else
                    {
                        H_1 += mLabelledCellCellAdhesionEnergyParameter;
                    }
                }
                else // no cell on node 1
                {
                    H_1 += mLabelledCellBoundaryAdhesionEnergyParameter;
                }
                if (is_cell_on_node_2)
                {
                    if (is_cell_2_labelled)
                    {
                        H_0 += mLabelledCellLabelledCellAdhesionEnergyParameter;
                    }
                    else
                    {
                        H_0 += mLabelledCellCellAdhesionEnergyParameter;
                    }
                }
                else // no cell on node 2
                {
                    H_0 += mLabelledCellBoundaryAdhesionEnergyParameter;
                }
            }
            else // Neighbour not labelled
            {
                if (is_cell_on_node_1)
                {
                    if (is_cell_1_labelled)
                    {
                        H_1 += mLabelledCellCellAdhesionEnergyParameter;
                    }
                    else
                    {
                        H_1 += mCellCellAdhesionEnergyParameter;
                    }
                }
                else // no cell on node 1
                {
                   H_1 += mCellBoundaryAdhesionEnergyParameter;
                }
                if (is_cell_on_node_2)
                {
                    if (is_cell_2_labelled)
                    {
                        H_0 += mLabelledCellCellAdhesionEnergyParameter;
                    }
                    else
                    {
                        H_0 += mCellCellAdhesionEnergyParameter;
                    }
                }
                else // no cell on node 2
                {
                    H_0 += mCellBoundaryAdhesionEnergyParameter;
                }
            }
        }
        else // No cell on neighbour node
        {
            if (is_cell_on_node_1)
            {
                if (is_cell_1_labelled)
                {
                    H_1 += mLabelledCellBoundaryAdhesionEnergyParameter;
                }
                else
                {
                    H_1 += mCellBoundaryAdhesionEnergyParameter;
                }
            }
            // else // no cell on node 1 or neighbour so no contribution
            if (is_cell_on_node_2)
            {
                if (is_cell_2_labelled)
                {
                    H_0 += mLabelledCellBoundaryAdhesionEnergyParameter;
                }
                else
                {
                    H_0 += mCellBoundaryAdhesionEnergyParameter;
                }
            }
            // else // no cell on node 1 or neighbour so no contribution
        }
    }

    return H_1-H_0;
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::EvaluateSwitchingProbability(
    unsigned currentNodeIndex,
    unsigned neighbourNodeIndex,
    CaBasedCellPopulation<DIM>& rCellPopulation,
    double dt,
    double deltaX)
{
    // Check if cell will have a Moore neighbour after the move?
    bool cells_connected = false;

    bool is_cell_on_node_1 = rCellPopulation.IsCellAttachedToLocationIndex(currentNodeIndex);
    bool is_cell_on_node_2 = rCellPopulation.IsCellAttachedToLocationIndex(neighbourNodeIndex);

    // To get here need at least one cell
    assert(is_cell_on_node_1 || is_cell_on_node_2);

    if ( !is_cell_on_node_1 || !is_cell_on_node_2 )
    {
        if (is_cell_on_node_1)
        {
            assert(!is_cell_on_node_2);

            std::set<unsigned> node_2_neighbouring_node_indices = rCellPopulation.rGetMesh().GetMooreNeighbouringNodeIndices(neighbourNodeIndex);

            // Remove node 1 from neighbouring indices
            node_2_neighbouring_node_indices.erase(currentNodeIndex);

            for (std::set<unsigned>::iterator iter = node_2_neighbouring_node_indices.begin();
                 iter != node_2_neighbouring_node_indices.end();
                 ++iter)
            {
                if (rCellPopulation.IsCellAttachedToLocationIndex(*iter))
                {
                    cells_connected = true;
                }
            }
        }
        else
        {
            assert(is_cell_on_node_2);
            assert(!is_cell_on_node_1);

            std::set<unsigned> node_1_neighbouring_node_indices = rCellPopulation.rGetMesh().GetMooreNeighbouringNodeIndices(currentNodeIndex);

            // Remove node 1 from neighbouring indices
            node_1_neighbouring_node_indices.erase(neighbourNodeIndex);

            for (std::set<unsigned>::iterator iter = node_1_neighbouring_node_indices.begin();
                            iter != node_1_neighbouring_node_indices.end();
                            ++iter)
            {
                if (rCellPopulation.IsCellAttachedToLocationIndex(*iter))
                {
                    cells_connected = true;
                }
            }
        }
    }
    else
    {
        // Two cells so must be connected after move
        cells_connected = true;
    }

    double probability_of_switch = 0.0;
    if (cells_connected)
    {
        double hamiltonian_difference = EvaluateHamiltonian(currentNodeIndex,neighbourNodeIndex,rCellPopulation);

        if (hamiltonian_difference <= 0)
        {
            probability_of_switch =  dt;
        }
        else
        {
            probability_of_switch = dt*exp(-hamiltonian_difference/mTemperature);
        }
    }

   return probability_of_switch;
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::GetCellCellAdhesionEnergyParameter()
{
    return mCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::SetCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::GetLabelledCellLabelledCellAdhesionEnergyParameter()
{
    return mLabelledCellLabelledCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::SetLabelledCellLabelledCellAdhesionEnergyParameter(double labelledCellLabelledCellAdhesionEnergyParameter)
{
    mLabelledCellLabelledCellAdhesionEnergyParameter = labelledCellLabelledCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::GetLabelledCellCellAdhesionEnergyParameter()
{
    return mLabelledCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::SetLabelledCellCellAdhesionEnergyParameter(double labelledCellCellAdhesionEnergyParameter)
{
    mLabelledCellCellAdhesionEnergyParameter = labelledCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::GetCellBoundaryAdhesionEnergyParameter()
{
    return mCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::SetCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::GetLabelledCellBoundaryAdhesionEnergyParameter()
{
    return mLabelledCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::SetLabelledCellBoundaryAdhesionEnergyParameter(double labelledCellBoundaryAdhesionEnergyParameter)
{
    mLabelledCellBoundaryAdhesionEnergyParameter = labelledCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::GetTemperature()
{
    return mTemperature;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::SetTemperature(double temperature)
{
    mTemperature = temperature;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellCellAdhesionEnergyParameter>" << mCellCellAdhesionEnergyParameter << "</CellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<LabelledCellLabelledCellAdhesionEnergyParameter>" << mLabelledCellLabelledCellAdhesionEnergyParameter << "</LabelledCellLabelledCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<LabelledCellCellAdhesionEnergyParameter>" << mLabelledCellCellAdhesionEnergyParameter << "</LabelledCellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<CellBoundaryAdhesionEnergyParameter>" << mCellBoundaryAdhesionEnergyParameter << "</CellBoundaryAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<LabelledCellBoundaryAdhesionEnergyParameter>" << mLabelledCellBoundaryAdhesionEnergyParameter << "</LabelledCellBoundaryAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<Temperature>" << mTemperature << "</Temperature>\n";

    // Call method on direct parent class
    AbstractCaSwitchingUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

// Explicit instantiation
template class DifferentialAdhesionCaSwitchingUpdateRule<1u>;
template class DifferentialAdhesionCaSwitchingUpdateRule<2u>;
template class DifferentialAdhesionCaSwitchingUpdateRule<3u>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DifferentialAdhesionCaSwitchingUpdateRule)
