#include <cmath>

#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

/**
  * Event level utilities
  *
  * - Contains things like DIS kinematics
  *
  * \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
  */

struct DISKinematics {
    DISKinematics() = default;

    double x{-1};
    double y{-1};
    double Q2{-1};
};


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
    unsigned int index = 1e6;
    for (unsigned int i = 0; i < _nHepmcp; ++i) {
        // Take the first proton, which should be the incoming proton.
        if (_hepmcp_PDG[i] == 2212) {
            index = i;
            break;
        }
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
    unsigned int index = 1e6;
    for (unsigned int i = 0; i < _nHepmcp; ++i) {
        // Take the first electron, which should be the incoming electron.
        if (_hepmcp_PDG[i] == 11) {
            index = i;
            break;
        }
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
    unsigned int index = 1e6;
    for (unsigned int i = 0; i < _nHepmcp; ++i) {
        // We'll take any struck quark (ie. PDG code under 9)
        if (std::abs(_hepmcp_PDG[i]) < 9) {
            index = i;
            break;
        }
    }
    return index;
}

unsigned int findHepMCIndexOfOutgoingElectron()
{
    unsigned int index = 1e6;
    // We'll take the first final state (==1) electron
    for (unsigned int i = 0; i < _nHepmcp; ++i) {
        if (_hepmcp_PDG[i] == 11 && _hepmcp_status[i] == 1) {
            index = i;
            break;
        }
    }
    return index;
}

unsigned int findPartLevelIndexOfOutgoingElectron()
{
    unsigned int index = 1e6;
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
    unsigned int index = 1e6;
    for (unsigned int i = 0; i < static_cast<unsigned int>(_nTracks); ++i) {
        if (_track_source[i] != primaryTrackSource) {
            continue;
        }
        if (_mcpart_PDG[static_cast<int>(_track_trueID[i])] == 11) {
            std::cout << i << ": PDG electron\n";
        }
        if (_track_trueID[i] == partLevelIndex) {
            index = i;
        }
    }
    return index;
}

DISKinematics KinematicsUsingTrueInfo(unsigned int primaryTrackSource)
{
    // All based on HepMC because the mcpart info only contains final state particles
    unsigned int hepmc_incomingProtonIndex = findHepMCIndexOfIncomingElectron();
    unsigned int hepmc_incomingElectronIndex = findHepMCIndexOfIncomingElectron();
    unsigned int hepmc_struckQuarkIndex = findHepMCIndexOfStruckQuark();

    // Outgoing electron
    unsigned int part_outgoingElectronIndex = findPartLevelIndexOfOutgoingElectron();
    unsigned int hepmc_outgoingElectronIndex = findHepMCIndexOfOutgoingElectron();
    unsigned int det_outgoingElectronIndex = findDetLevelIndexOfOutgoingElectron(part_outgoingElectronIndex, primaryTrackSource);

    // Checks
    if (_mcpart_BCID[part_outgoingElectronIndex] != _hepmcp_BCID[hepmc_outgoingElectronIndex]) {
        std::cout << "BCID mismatch...\n";
        std::cout << "mcpart BCID: " << _mcpart_BCID[part_outgoingElectronIndex]
                  << "\nhepmc BCID: " << _hepmcp_BCID[hepmc_outgoingElectronIndex] << "\n";
    }
    /*
    std::cout << "Kinematics\n";
    std::cout << "mcpart"
              << ", px=" << _mcpart_px[part_outgoingElectronIndex]
              << ", py=" << _mcpart_py[part_outgoingElectronIndex]
              << ", pz=" << _mcpart_pz[part_outgoingElectronIndex]
              << ", E=" << _mcpart_E[part_outgoingElectronIndex]
              << "\n";
    std::cout << "hepmcp"
              << ", px=" << _hepmcp_px[hepmc_outgoingElectronIndex]
              << ", py=" << _hepmcp_py[hepmc_outgoingElectronIndex]
              << ", pz=" << _hepmcp_pz[hepmc_outgoingElectronIndex]
              << ", E=" << _hepmcp_E[hepmc_outgoingElectronIndex]
              << "\n";
    */

    // Electron kinematics
    std::cout << "hepmc_incomingElectronIndex=" << hepmc_incomingElectronIndex
              << ", det_outgoingElectronIndex=" << det_outgoingElectronIndex
              << "\n";
    ROOT::Math::PxPyPzEVector hepmc_incomingElectron(0, 0, _hepmcp_pz[hepmc_incomingElectronIndex], _hepmcp_E[hepmc_incomingElectronIndex]);
    //ROOT::Math::PxPyPzEVector hepmc_incomingElectron(0, 0, _hepmcp_pz[hepmc_incomingElectronIndex], _hepmcp_pz[hepmc_incomingElectronIndex]);
    /*ROOT::Math::PxPyPzEVector hepmc_outgoingElectron(
        _track_px[det_outgoingElectronIndex],
        _track_py[det_outgoingElectronIndex],
        _track_pz[det_outgoingElectronIndex],
        std::sqrt(
            std::pow(_track_px[det_outgoingElectronIndex], 2) +
            std::pow(_track_py[det_outgoingElectronIndex], 2) +
            std::pow(_track_pz[det_outgoingElectronIndex], 2)
        )
    );*/
    /*ROOT::Math::PxPyPzEVector hepmc_outgoingElectron(
        _mcpart_px[part_outgoingElectronIndex],
        _mcpart_py[part_outgoingElectronIndex],
        _mcpart_pz[part_outgoingElectronIndex],
        std::sqrt(
            std::pow(_mcpart_px[part_outgoingElectronIndex], 2) +
            std::pow(_mcpart_py[part_outgoingElectronIndex], 2) +
            std::pow(_mcpart_pz[part_outgoingElectronIndex], 2)
        )
    );*/
    ROOT::Math::PxPyPzEVector hepmc_outgoingElectron(
        _hepmcp_px[hepmc_outgoingElectronIndex],
        _hepmcp_py[hepmc_outgoingElectronIndex],
        _hepmcp_pz[hepmc_outgoingElectronIndex],
        std::sqrt(
            std::pow(_hepmcp_px[hepmc_outgoingElectronIndex], 2) +
            std::pow(_hepmcp_py[hepmc_outgoingElectronIndex], 2) +
            std::pow(_hepmcp_pz[hepmc_outgoingElectronIndex], 2)
        )
    );

    auto hepmc_virtualPhoton = hepmc_incomingElectron - hepmc_outgoingElectron;

    // Proton kinematics
    double crossingAngle = 0.025;
    ROOT::Math::PxPyPzEVector hepmc_incomingProton(
        -std::abs(_hepmcp_pz[hepmc_incomingProtonIndex]) * std::sin(crossingAngle),
        0,
        std::abs(_hepmcp_pz[hepmc_incomingProtonIndex]) * std::cos(crossingAngle),
        _hepmcp_E[hepmc_incomingProtonIndex]
    );

    double q2 = -hepmc_virtualPhoton.mag2();
    double x = q2 / (2 * hepmc_virtualPhoton.Dot(hepmc_incomingProton));
    double y = (hepmc_incomingProton.Dot(hepmc_incomingProton)) / (hepmc_incomingProton.Dot(hepmc_incomingElectron));

    std::cout << "Calculated: q2=" << q2 << ", x=" << x << ", y=" << y << "\n";
    std::cout << "From HepMC: q2=" << _hepmcp_Q2 << ", x=" << _hepmcp_x2 << "\n";

    return DISKinematics{x, y, q2};
}

/**
  * Determine DIS kinematics variables using the Jacquet-Blondel method
  *
  * Implemented using https://wiki.bnl.gov/eic/index.php/DIS_Kinematics
  */
DISKinematics JBKinematics(unsigned short primaryTrackSource, double incomingElectronEnergy, double incomingHadronEnergy, double hadronMassHypothesis = 0.139)
{
    // NOTE: We need electron ID here. Since we're don't yet have that available as of 16 June 2021, we use the true PID info.
    // NOTE: 1e6 is a sentinel value.
    unsigned int electronID = 1e6;
    double electronPt = -1e6;
    for (unsigned int i = 0; i < (unsigned int)_nTracks; ++i) {
        if (_track_source[i] != primaryTrackSource) {
            continue;
        }
        if (_mcpart_PDG[static_cast<int>(_track_trueID[i])] == 11) {
            // If it's an electron, take the highest pt one.
            // TODO: Update this when the electron PID is addressed.
            double pt = std::sqrt(std::pow(_track_px[i], 2) + std::pow(_track_py[i], 2));
            if (pt > electronPt) {
                electronID = i;
                electronPt = pt;
            }
        }
    }
    if (electronID == 1e6) {
        throw std::runtime_error("Unable to find scattered electron!");
    }

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

    double sqrts = 4 * incomingElectronEnergy * incomingHadronEnergy;

    double y = e_minus_pz / (2 * incomingElectronEnergy);
    double Q2 = (std::pow(px_h, 2) + std::pow(py_h, 2)) / (1 - y);
    double x = Q2 / (sqrts * y);

    return DISKinematics{x, y, Q2};
}
