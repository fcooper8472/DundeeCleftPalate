
#include "GeneralisedChemotaxisPottsUpdateRule.hpp"
#include "UblasIncludes.hpp"

template<unsigned DIM>
GeneralisedChemotaxisPottsUpdateRule<DIM>::GeneralisedChemotaxisPottsUpdateRule(c_vector<double,DIM> chemotacticSourceLocation)
    : AbstractPottsUpdateRule<DIM>(),
	  mChemotacticSourceLocation(chemotacticSourceLocation)
{
}

template<unsigned DIM>
GeneralisedChemotaxisPottsUpdateRule<DIM>::~GeneralisedChemotaxisPottsUpdateRule()
{
}

template<unsigned DIM>
double GeneralisedChemotaxisPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(
    unsigned currentNodeIndex,
    unsigned targetNodeIndex,
    PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    // Note that we define these vectors before setting them as otherwise the profiling build will break (see #2367)
    c_vector<double, DIM> current_location;
    current_location = rCellPopulation.GetNode(currentNodeIndex)->rGetLocation();
    c_vector<double, DIM> target_location;
    target_location = rCellPopulation.GetNode(targetNodeIndex)->rGetLocation();

    c_vector<double, DIM> proposed_target_direction = target_location-current_location;

    c_vector<double, DIM> chemotactic_direction = current_location-mChemotacticSourceLocation;

    double scalar_product = inner_prod(proposed_target_direction/norm_2(proposed_target_direction),chemotactic_direction/norm_2(chemotactic_direction));

    double delta_H = scalar_product*0.1;
    return delta_H;
}

template<unsigned DIM>
void GeneralisedChemotaxisPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

// Explicit instantiation
template class GeneralisedChemotaxisPottsUpdateRule<1>;
template class GeneralisedChemotaxisPottsUpdateRule<2>;
template class GeneralisedChemotaxisPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GeneralisedChemotaxisPottsUpdateRule)
