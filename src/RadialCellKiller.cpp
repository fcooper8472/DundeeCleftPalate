
#include "RadialCellKiller.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsElement.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned DIM>
RadialCellKiller<DIM>::RadialCellKiller(AbstractCellPopulation<DIM>* pCellPopulation, double probabilityOfDeathInAnHour)
        : AbstractCellKiller<DIM>(pCellPopulation),
          mProbabilityOfDeathInAnHour(probabilityOfDeathInAnHour)
{
    if ((mProbabilityOfDeathInAnHour<0) || (mProbabilityOfDeathInAnHour>1))
    {
        EXCEPTION("Probability of death must be between zero and one");
    }
}

template<unsigned DIM>
double RadialCellKiller<DIM>::GetDeathProbabilityInAnHour() const
{
    return mProbabilityOfDeathInAnHour;
}

template<unsigned DIM>
void RadialCellKiller<DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
    double R = 40;

    double x_1 = 80;
    double y_1 = 160;

    double x_2 = 160;
    double y_2 = 80;

    double dist_to_centre1_squared = 0.0;
    double dist_to_centre2_squared = 0.0;

    if (dynamic_cast<PottsBasedCellPopulation<DIM>* >(this->mpCellPopulation)!=nullptr)
    {
		PottsBasedCellPopulation<DIM>* p_PottsPopulation=dynamic_cast<PottsBasedCellPopulation<DIM>* >(this->mpCellPopulation);
		PottsMesh<DIM>* p_PottsMesh=dynamic_cast<PottsMesh<DIM>* >(&(p_PottsPopulation->rGetMesh()));

		unsigned element_index = p_PottsPopulation->GetElementCorrespondingToCell(pCell)->GetIndex();
		c_vector<double,DIM> element_centre = p_PottsMesh->GetCentroidOfElement(element_index);
		 dist_to_centre1_squared = pow(element_centre[0]-x_1,2) + pow(element_centre[1]-y_1,2);
		 dist_to_centre2_squared = pow(element_centre[0]-x_2,2) + pow(element_centre[1]-y_2,2);
    }

	if (!pCell->HasApoptosisBegun() )
	{
		if (!(dist_to_centre1_squared < pow(R,2) || dist_to_centre2_squared < pow(R,2)))
		{
			pCell->StartApoptosis();
		}
	}
}

template<unsigned DIM>
void RadialCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        CheckAndLabelSingleCellForApoptosis(*cell_iter);
    }
}

template<unsigned DIM>
void RadialCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ProbabilityOfDeathInAnHour>" << mProbabilityOfDeathInAnHour << "</ProbabilityOfDeathInAnHour>\n";

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class RadialCellKiller<1>;
template class RadialCellKiller<2>;
template class RadialCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RadialCellKiller)
