
#ifndef GENERALISEDCHEMOTAXISUPDATERULE_HPP_
#define GENERALISEDCHEMOTAXISUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "UblasVectorInclude.hpp"
#include "AbstractPottsUpdateRule.hpp"
#include "PottsBasedCellPopulation.hpp"

/**
 * A simple update rule class to represent simple chemotaxis towards a chemotactic
 * centre (set in the constructor).
 *
 * Note this currently assumes cells don't grow, i.e the target volume is constant
 * for each cell over time.
 */
template<unsigned DIM>
class GeneralisedChemotaxisPottsUpdateRule : public AbstractPottsUpdateRule<DIM>
{
private:

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractPottsUpdateRule<DIM> >(*this);
        archive & mChemotacticSourceLocation;
    }

    /** \todo Document member */
    c_vector<double,DIM> mChemotacticSourceLocation;

    /**
     * Private default constructor for archiving.
     */
     GeneralisedChemotaxisPottsUpdateRule(){};

public:

    /**
     * Constructor.
     */
    GeneralisedChemotaxisPottsUpdateRule(c_vector<double,DIM> chemotacticSourceLocation);

    /**
     * Destructor.
     */
    ~GeneralisedChemotaxisPottsUpdateRule();

    /**
     * Overridden EvaluateHamiltonianContribution() method.
     *
     * Assigns a greater propensity for moving towards the chemotactic centre set in the constructor.
     *
     * @param currentNodeIndex The index of the current node/lattice site
     * @param targetNodeIndex The index of the target node/lattice site
     * @param rCellPopulation The cell population
     *
     * @return The difference in the Hamiltonian with the configuration of the target node
     * having the same spin as the current node with the current configuration. i.e H_1-H_0
     */
    double EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                           unsigned targetNodeIndex,
                                           PottsBasedCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GeneralisedChemotaxisPottsUpdateRule)

#endif /*GENERALISEDCHEMOTAXISUPDATERULE_HPP_*/
