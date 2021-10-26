#include <cmath>
#include <iostream>

#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

/**
  * Event level utilities
  *
  * - Contains things like DIS kinematics
  *
  * \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
  */

 void printParticles()
 {
    std::cout << "Particles\n";
    for(unsigned int i=0; i<_nMCPart; ++i) {
        std::cout << i << "-> id=" << _mcpart_PDG[i] << ", barcode=" << _mcpart_BCID[i] << ", px=" << _mcpart_px[i]
                << ", py=" << _mcpart_py[i]
                << ", pz=" << _mcpart_pz[i]
                << ", E=" << _mcpart_E[i]
                << "\n";
    }
    std::cout << "HepMC\n";
    for (unsigned int i = 0; i < _nHepmcp; ++i) {
        std::cout << i << "-> id=" << _hepmcp_PDG[i] << ", barcode=" << _hepmcp_BCID[i] << ", px=" << _hepmcp_px[i]
                << ", py=" << _hepmcp_py[i]
                << ", pz=" << _hepmcp_pz[i]
                << ", E=" << _hepmcp_E[i]
                << ", status=" << _hepmcp_status[i]
                << "\n";
    }
    std::cout << "End particles" << std::endl;
 }

struct DISKinematics {
    DISKinematics() = default;

    double x{-1};
    double y{-1};
    double Q2{-1};
};

enum KinematicsErrors_t {
    kInvalid_x = 0,
    kInvalidOption,
    kIncomingProtonNotFound,
    kIncomingElectronNotFound,
    kStruckQuarkNotFound,
    kPartLevelElectronNotFound,
    kRecoElectronNotFound,
    kSuccess,
};

constexpr unsigned int STARTING_INDEX = 1e6;

/**
  * Find the index of the incoming proton in the HepMC truth particles.
  *
  * NOTE: It seems to be in position 0 in tests, but I don't know if that's absolutely required,
  *       especially for different generators.
  *
  * @returns Index of the proton.
  */
unsigned int findHepMCIndexOfIncomingProton()
{
    unsigned int index = STARTING_INDEX;
    for (unsigned int i = 0; i < _nHepmcp; ++i) {
        // Take the first proton, which should be the incoming proton.
        if (_hepmcp_PDG[i] == 2212) {
            index = i;
            break;
        }
    }
    // Validation
    if (index == STARTING_INDEX) {
        throw KinematicsErrors_t::kIncomingProtonNotFound;
    }
    return index;
}

/**
  * Find the index of the incoming electron in the HepMC truth particles.
  *
  * NOTE: It seems to be in index 1 in tests, but I don't know if that's absolutely required,
  *       especially for different generators.
  *
  * @returns Index of the electron.
  */
unsigned int findHepMCIndexOfIncomingElectron()
{
    unsigned int index = STARTING_INDEX;
    for (unsigned int i = 0; i < _nHepmcp; ++i) {
        // Take the first electron, which should be the incoming electron.
        if (_hepmcp_PDG[i] == 11) {
            index = i;
            break;
        }
    }
    // Validation
    if (index == STARTING_INDEX) {
        throw KinematicsErrors_t::kIncomingElectronNotFound;
    }
    return index;
}

/**
  * Find the quark struck by the electron based on the HepMC truth particles.
  *
  * @returns Index of the struck quark in the HepMC truth particles
  */
unsigned int findHepMCIndexOfStruckQuark()
{
    unsigned int index = STARTING_INDEX;
    for (unsigned int i = 0; i < _nHepmcp; ++i) {
        // We'll take any struck quark (ie. PDG code under 9)
        if (std::abs(_hepmcp_PDG[i]) < 9) {
            index = i;
            break;
        }
    }
    // Validation
    if (index == STARTING_INDEX) {
        throw KinematicsErrors_t::kStruckQuarkNotFound;
    }
    return index;
}

unsigned int findHepMCIndexOfOutgoingElectron()
{
    unsigned int index = STARTING_INDEX;
    // We'll take the first final state (==1) electron. Since it follows the pythia event record,
    // it seems like this is highly likely to be the final state electron.
    for (unsigned int i = 0; i < _nHepmcp; ++i) {
        if (_hepmcp_PDG[i] == 11 && _hepmcp_status[i] == 1) {
            index = i;
            break;
        }
    }
    // Validation
    if (index == STARTING_INDEX) {
        throw KinematicsErrors_t::kPartLevelElectronNotFound;
    }
    return index;
}

unsigned int findPartLevelIndexOfOutgoingElectron()
{
    unsigned int index = STARTING_INDEX;
    int barcode = 1e6;
    for (unsigned int i = 0; i < _nMCPart; ++i) {
        if (_mcpart_PDG[i] != 11) {
            continue;
        }

        // We have an electron. Take the lowest barcode index (meaning the earliest in the HepMC record),
        // which should correspond to the outgoing electron
        if (_mcpart_BCID[i] > -1 && (_mcpart_BCID[i] < barcode)) {
            index = i;
            barcode = _mcpart_BCID[i];
        }
    }
    // Validation
    if (index == STARTING_INDEX) {
        throw KinematicsErrors_t::kPartLevelElectronNotFound;
    }
    return index;
}

/**
 * @brief Finds the detector level outgoing electron using true (particle level) information.
 *
 * @param partLevelIndex Part level index.
 * @param primaryTrackSource Primary track source to be used for finding the electrons.
 * @return unsigned int Detector level index
 */
unsigned int findDetLevelIndexOfOutgoingElectron(unsigned int partLevelIndex, unsigned int primaryTrackSource)
{
    unsigned int index = STARTING_INDEX;
    for (unsigned int i = 0; i < static_cast<unsigned int>(_nTracks); ++i) {
        if (_track_source[i] != primaryTrackSource) {
            continue;
        }
        //if (_mcpart_PDG[static_cast<int>(_track_trueID[i])] == 11) {
        //    std::cout << i << ": PDG electron\n";
        //}
        if (_track_trueID[i] == partLevelIndex) {
            index = i;
        }
    }
    // Validation
    if (index == STARTING_INDEX) {
        throw KinematicsErrors_t::kRecoElectronNotFound;
    }
    return index;
}

DISKinematics CalculateKinematicsFromFourVectors(const ROOT::Math::PxPyPzEVector & incomingProton, const ROOT::Math::PxPyPzEVector & incomingElectron, const ROOT::Math::PxPyPzEVector & outgoingElectron)
{
    auto virtualPhoton = incomingElectron - outgoingElectron;
    double q2 = -virtualPhoton.mag2();

    double x = q2 / (2 * virtualPhoton.Dot(incomingProton));
    double y = (incomingProton.Dot(virtualPhoton)) / (incomingProton.Dot(incomingElectron));

    return DISKinematics{x, y, q2};
}

enum DirectKinematicsOptions_t {
    kHepMCStored = 0,
    kHepMCCalculated,
    kPartLevel,
    kDetLevel,
};

DISKinematics KinematicsUsingTrueInfo(DirectKinematicsOptions_t option, unsigned int primaryTrackSource = 0, bool isPythia6 = false, bool isPythia8 = false, unsigned int verbosity = 0)
{
    // Cross check. This seems the most relevant place, so we'll look here
    // NOTE: Scoped so it doesn't mess with anything else.
    {
        unsigned int hepmc_outgoingElectronIndex = findHepMCIndexOfOutgoingElectron();
        unsigned int part_outgoingElectronIndex = findPartLevelIndexOfOutgoingElectron();
        if (_mcpart_BCID[part_outgoingElectronIndex] != _hepmcp_BCID[hepmc_outgoingElectronIndex]) {
            std::cout << "WARNING: BCID mismatch for outgoing electron index, meaning different electrons were found at particle level and in HepMC. Check this out...\n";
            std::cout << "part level index: " << part_outgoingElectronIndex
                      << "\nhepmc index: " << hepmc_outgoingElectronIndex << "\n";
            std::cout << "mcpart BCID: " << _mcpart_BCID[part_outgoingElectronIndex]
                    << "\nhepmc BCID: " << _hepmcp_BCID[hepmc_outgoingElectronIndex] << "\n";
            // Print particles to help with debugging:
            printParticles();
        }
    }

    DISKinematics returnKinematics;
    if (option == DirectKinematicsOptions_t::kHepMCStored || option == DirectKinematicsOptions_t::kHepMCCalculated) {
        // Incoming proton
        unsigned int hepmc_incomingProtonIndex = findHepMCIndexOfIncomingProton();
        ROOT::Math::PxPyPzEVector hepmc_incomingProton(
            // NOTE: px and py are effectively 0 (usually 1e-13 or so).
            //       Here, I set them to be identically 0 to avoid numerical instabilities, especially
            //       since their values shouldn't matter.
            //_hepmcp_px[hepmc_incomingProtonIndex],
            //_hepmcp_py[hepmc_incomingProtonIndex],
            0, 0,
            // NOTE: The sign on pz is positive.
            _hepmcp_pz[hepmc_incomingProtonIndex],
            _hepmcp_E[hepmc_incomingProtonIndex]
        );
        // Incoming electron
        unsigned int hepmc_incomingElectronIndex = findHepMCIndexOfIncomingElectron();
        ROOT::Math::PxPyPzEVector hepmc_incomingElectron(
            // NOTE: px and py are effectively 0 (usually 1e-13 or so).
            //       Here, I set them to be identically 0 to avoid numerical instabilities, especially
            //       since their values shouldn't matter.
            //_hepmcp_px[hepmc_incomingElectronIndex],
            //_hepmcp_py[hepmc_incomingElectronIndex],
            0, 0,
            // NOTE: The sign on the pz is negative.
            _hepmcp_pz[hepmc_incomingElectronIndex],
            _hepmcp_E[hepmc_incomingElectronIndex]
        );
        //std::cout << "Proton"
        //    << ", px=" << _hepmcp_px[hepmc_incomingProtonIndex]
        //    << ", py=" << _hepmcp_py[hepmc_incomingProtonIndex]
        //    << ", pz=" << _hepmcp_pz[hepmc_incomingProtonIndex]
        //    << ", E=" << _hepmcp_E[hepmc_incomingProtonIndex] << "\n";
        //std::cout << "Electron"
        //    << ", px=" << _hepmcp_px[hepmc_incomingElectronIndex]
        //    << ", py=" << _hepmcp_py[hepmc_incomingElectronIndex]
        //    << ", pz=" << _hepmcp_pz[hepmc_incomingElectronIndex]
        //    << ", E=" << _hepmcp_E[hepmc_incomingElectronIndex] << "\n";
        // Outoging electron
        unsigned int hepmc_outgoingElectronIndex = findHepMCIndexOfOutgoingElectron();
        ROOT::Math::PxPyPzEVector hepmc_outgoingElectron(
            _hepmcp_px[hepmc_outgoingElectronIndex],
            _hepmcp_py[hepmc_outgoingElectronIndex],
            _hepmcp_pz[hepmc_outgoingElectronIndex],
            _hepmcp_E[hepmc_outgoingElectronIndex]
        );

        if (option == kHepMCStored) {
            auto virtualPhoton = hepmc_incomingElectron - hepmc_outgoingElectron;
            double q2 = -virtualPhoton.mag2();

            // Take the stored values, corrected as required for the different production quirks
            // These two values are taken by convention, and then are adpated as required.
            double xStored = _hepmcp_x1;
            double q2Stored = _hepmcp_Q2;
            // Pythia6:
            // - Stores the quark x in x2, so we swap it here.
            // - Properly stores q2 in hepMC.
            if (isPythia6) {
                xStored = _hepmcp_x2;
            }
            // Pythia8:
            // - Stores the quark x in x1.
            // - Stores Q, not Q2, in the hepMC, so we correct it here.
            //   Note that this Q2 still tends to be smaller than the other calculated Q2, perhaps because
            //   the stored Q2 is from the renormalization scale (?) rather than the factorization scale (?)
            if (isPythia8) {
                q2Stored = q2Stored * q2Stored;
            }
            //std::cout << "Calculated q2=" << q2 << ", stored=" << q2Stored << "\n";
            //std::cout << "stored x1=" << _hepmcp_x1 << ", stored x2=" << _hepmcp_x2 << ", xStored=" << xStored << "\n";

            // Here, we only have to calculate y, and then we return stored values for the rest.
            double y = (hepmc_incomingProton.Dot(virtualPhoton)) / (hepmc_incomingProton.Dot(hepmc_incomingElectron));
            returnKinematics = DISKinematics{xStored, y, q2Stored};
            // NOTE: Intentionally no newline so it blends with the line printout.
            if (verbosity > 1) { std::cout << "Stored HepMC "; }
        }
        else {
            returnKinematics = CalculateKinematicsFromFourVectors(
                hepmc_incomingProton, hepmc_incomingElectron, hepmc_outgoingElectron
            );
            // NOTE: Intentionally no newline so it blends with the line printout.
            if (verbosity > 1) { std::cout << "Calculated HepMC "; }
        }
    }
    else if (option == DirectKinematicsOptions_t::kPartLevel || option == DirectKinematicsOptions_t::kDetLevel) {
        // Handle kinematics in the lab frame.

        // Setup
        // Would be better to take this from an evaluator, but I'm not sure how to do that.
        // So, instead I'll just define it here.
        double crossingAngle = 0.025;

        // To start, we don't have the incoming beam particles in the lab frame (ie. we don't store them in mcpart).
        // So we'll do the best we can by rotating the proton incoming energy from HepMC (following John's example).
        unsigned int hepmc_incomingProtonIndex = findHepMCIndexOfIncomingProton();
        ROOT::Math::PxPyPzEVector lab_incomingProton(
            -_hepmcp_pz[hepmc_incomingProtonIndex] * std::sin(crossingAngle),
            0,
            _hepmcp_pz[hepmc_incomingProtonIndex] * std::cos(crossingAngle),
            _hepmcp_E[hepmc_incomingProtonIndex]
        );
        // For the electorn, we just take the energy
        unsigned int hepmc_incomingElectronIndex = findHepMCIndexOfIncomingElectron();
        ROOT::Math::PxPyPzEVector lab_incomingElectron(
            0, 0,
            _hepmcp_pz[hepmc_incomingElectronIndex],
            _hepmcp_E[hepmc_incomingElectronIndex]
        );

        unsigned int part_outgoingElectronIndex = findPartLevelIndexOfOutgoingElectron();

        // Here, we split paths depending on whether we want the particle level or detector level.
        if (option == DirectKinematicsOptions_t::kPartLevel) {
            // Outoging electron at particle level:
            ROOT::Math::PxPyPzEVector part_outgoingElectron(
                _mcpart_px[part_outgoingElectronIndex],
                _mcpart_py[part_outgoingElectronIndex],
                _mcpart_pz[part_outgoingElectronIndex],
                _mcpart_E[part_outgoingElectronIndex]
            );

            returnKinematics = CalculateKinematicsFromFourVectors(
                lab_incomingProton, lab_incomingElectron, part_outgoingElectron
            );
            // NOTE: Intentionally no newline so it blends with the line printout.
            if (verbosity > 1) { std::cout << "Part. level "; }
        }
        else {
            // At reconstructed level, using the reconstructed electron
            unsigned int det_outgoingElectronIndex = findDetLevelIndexOfOutgoingElectron(part_outgoingElectronIndex, primaryTrackSource);
            // NOTE: The lab proton and electron are the same because we're not mesasuring them
            if (det_outgoingElectronIndex < 1e6) {
                ROOT::Math::PxPyPzEVector reco_outgoingElectron(
                    _track_px[det_outgoingElectronIndex],
                    _track_py[det_outgoingElectronIndex],
                    _track_pz[det_outgoingElectronIndex],
                    std::sqrt(
                        std::pow(_track_px[det_outgoingElectronIndex], 2) +
                        std::pow(_track_py[det_outgoingElectronIndex], 2) +
                        std::pow(_track_pz[det_outgoingElectronIndex], 2)
                    )
                );
                returnKinematics = CalculateKinematicsFromFourVectors(
                    lab_incomingProton, lab_incomingElectron, reco_outgoingElectron
                );
                // NOTE: Intentionally no newline so it blends with the line printout.
                if (verbosity > 1) { std::cout << "Reco. level "; }
            }
            else {
                throw KinematicsErrors_t::kRecoElectronNotFound;
            }
        }
    }
    else {
        throw KinematicsErrors_t::kInvalidOption;
    }

    if (verbosity > 1) {
        std::cout << "calculated x: " << returnKinematics.x
                << ", y=" << returnKinematics.y
                << ", q2=" << returnKinematics.Q2 << "\n";
    }

    return returnKinematics;
}

/**
  * Determine DIS kinematics variables using the Jacquet-Blondel method
  *
  * Implemented using https://wiki.bnl.gov/eic/index.php/DIS_Kinematics
  */
DISKinematics JBKinematics(unsigned short primaryTrackSource = 0, double hadronMassHypothesis = 0.139)
{
    const unsigned int electronID = findDetLevelIndexOfOutgoingElectron(
        findPartLevelIndexOfOutgoingElectron(), primaryTrackSource
    );
    const double electronPt = std::sqrt(std::pow(_track_px[electronID], 2) + std::pow(_track_py[electronID], 2));

    // Loop over the final state hadrons to determine the kinematics.
    double px_h = 0;
    double py_h = 0;
    double e_minus_pz = 0;
    for(unsigned int i = 0; i < (unsigned int)_nTracks; ++i) {
        // Only consider tracks from the primary track source.
        if (_track_source[i] != primaryTrackSource) {
            continue;
        }
        // Skip the scattered electron.
        if (i == electronID) {
            continue;
        }
        px_h += _track_px[i];
        py_h += _track_py[i];
        double eTrack = std::sqrt(std::pow(_track_px[i],2) + std::pow(_track_py[i],2) + std::pow(_track_pz[i],2) + std::pow(hadronMassHypothesis, 2));
        e_minus_pz += (eTrack - _track_pz[i]);
    }

    // TODO: Add calorimeter towers.

    // Calculate kinematics
    double incomingElectronEnergy = _hepmcp_E[findHepMCIndexOfIncomingElectron()];
    double incomingHadronEnergy = _hepmcp_E[findHepMCIndexOfIncomingProton()];

    double sqrts = 4 * incomingElectronEnergy * incomingHadronEnergy;

    double y = e_minus_pz / (2 * incomingElectronEnergy);
    double Q2 = (std::pow(px_h, 2) + std::pow(py_h, 2)) / (1 - y);
    double x = Q2 / (sqrts * y);

    return DISKinematics{x, y, Q2};
}
