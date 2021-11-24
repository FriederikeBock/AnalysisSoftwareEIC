
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/GhostedAreaSpec.hh>

bool _do_jetfinding = false;

// Uncomment if you have a sufficiently recent fj contrib.
//#include <fastjet/contrib/Centauro.hh>

/**
  * ANCHOR Perform jet finding on the given input four vectors.
  *
  * Note: We have to return the cluster sequence because the PseudoJet objects depend
  *       on the CS still being alive. Otherwise, you get:
  *       ```
  *       fastjet::Error:  you requested information about the internal structure of a jet, but its associated ClusterSequence has gone out of scope.
  *       ```
  *       Since the CS is handled with a shared_ptr, it should be cleaned up when it goes out of scope.
  *
  * @param[in] jetR Jet resolution parameter.
  * @param[in] jetAlgorithm Name of the jet algorithm to use. Only some are implemented here.
  * @param[in] px Px of the input four vectors.
  * @param[in] py Py of the input four vectors.
  * @param[in] pz Pz of the input four vectors.
  * @param[in] E Energy of the input four vectors.
  * @param[in] areaDefinition Jet area definition. Defaults should be fine for most purposes.
  *
  * @returns ClusterSequenceArea and jets found (as PseudoJet objects) in a tuple.
  */
template<typename T>
std::tuple<std::shared_ptr<fastjet::ClusterSequenceArea>, std::vector<fastjet::PseudoJet>> findJets(
        double jetR,
        std::string jetAlgorithm,
        // Particle 4 vectors
        std::vector<T> px,
        std::vector<T> py,
        std::vector<T> pz,
        std::vector<T> E,
        fastjet::AreaDefinition areaDefinition = fastjet::AreaDefinition(fastjet::AreaType::active_area, fastjet::GhostedAreaSpec(4, 1, 0.05))
)
{
    // Validation
    std::vector<std::string> implementedAlgorithms = {"anti-kt", "centauro"};
    if (std::find(implementedAlgorithms.begin(), implementedAlgorithms.end(), jetAlgorithm) == implementedAlgorithms.end())
    {
        throw std::out_of_range(std::string("Algorithm ") + jetAlgorithm + " is not implemented");
    }
    // Determine jet definition
    // NOTE: Implicitly uses the E scheme
    std::shared_ptr<fastjet::JetDefinition> jetDefinition = nullptr;
    if (jetAlgorithm == "anti-kt") {
        //jetDefinition = std::shared_ptr<fastjet::JetDefinition>(new fastjet::JetDefinition(fastjet::JetAlgorithm::antikt_algorithm, jetR));
        jetDefinition = std::make_shared<fastjet::JetDefinition>(fastjet::JetAlgorithm::antikt_algorithm, jetR);
    }
    // Uncomment if you have a sufficiently recent fj contrib.
    /*else if (jetAlgorithm == "centauro") {
        fastjet::contrib::CentauroPlugin centauroPlugin(jetR);
        jetDefinition = new fastjet::JetDefinition(centauroPlugin);
    }*/

    // Convert inputs to PseudoJet objects
    std::vector<fastjet::PseudoJet> inputPseudoJets;
    for (std::size_t i = 0; i < px.size(); i++) {
        inputPseudoJets.emplace_back(
            fastjet::PseudoJet(
                px.at(i),
                py.at(i),
                pz.at(i),
                E.at(i)
            )
        );
    }

    // Perform jet finding
    std::shared_ptr<fastjet::ClusterSequenceArea> cs = std::make_shared<fastjet::ClusterSequenceArea>(inputPseudoJets, *jetDefinition, areaDefinition);
    auto jets = fastjet::sorted_by_E(cs->inclusive_jets(0));

    return std::make_tuple(cs, jets);
}

/**
  * ANCHOR Trivial helper to print jets and consittuents.
  *
  * Based on fastjet example.
  *
  * @param[in] jets Jets found by the jet finder.
  */
void printJets(const std::vector<fastjet::PseudoJet> & jets) {
    for (std::size_t i = 0; i < jets.size(); i++) {
        std::cout << "jet " << i << ": "<< jets[i].pt() << " " << jets[i].rap() << " " << jets[i].phi() << std::endl;
        // vector<fastjet::PseudoJet> constituents = jets[i].constituents();
        // for (std::size_t j = 0; j < constituents.size(); j++) {
        //     std::cout << "    constituent " << j << "'s pt: " << constituents[j].pt() << std::endl;
        // }
    }
}
/**
  * ANCHOR Example of finding jets.
  *
  * Based on the fastjet example.
  *
  * Note: We have to return the cluster sequence because the PseudoJet objects depend
  *       on the CS still being alive. Otherwise, you get:
  *       ```
  *       fastjet::Error:  you requested information about the internal structure of a jet, but its associated ClusterSequence has gone out of scope.
  *       ```
  *       Since the CS is handled with a shared_ptr, it should be cleaned up when it goes out of scope.
  *
  * @returns ClusterSequenceArea and jets found (as PseudoJet objects) in a tuple.
  */
std::tuple<std::shared_ptr<fastjet::ClusterSequenceArea>, std::vector<fastjet::PseudoJet>> jet_finder()
{
    std::vector<float> px = {99.0, 4.0, -99.0};
    std::vector<float> py = {0.1, -0.1, 0};
    std::vector<float> pz = {0, 0, 0};
    std::vector<float> E = {100.0, 5.0, 99.0};

    // Perform the actual jet finding.
    auto res = findJets(0.5, "anti-kt", px, py, pz, E);
    // Reminder of how to get values out of a tuple.
    printJets(std::get<1>(res));

    return res;
}

