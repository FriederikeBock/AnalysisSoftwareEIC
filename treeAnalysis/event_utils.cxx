/**
  * Event level utilities
  *
  * - Contains things like DIS kinematics
  *
  * \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
  */

struct DISKinematics {
    DISKinematics() = default;

    double x;
    double y;
    double Q2;
};

/**
  * Determine DIS kinematics variables using the Jacquet-Blondel method
  *
  * Implemented using https://wiki.bnl.gov/eic/index.php/DIS_Kinematics
  */
DISKinematics JBKinematics(unsigned int primaryTrackSource, double incomingElectronEnergy, double sqrts, double hadronMassHypothesis = 0.139)
{
    // NOTE: We need electron ID here. Since we're don't yet have that available as of 16 June 2021, we use the true PID info.
    // NOTE: 1e6 is a sentinel value.
    unsigned int electronID = 1e6;
    double electronPt = -1e6;
    for (unsigned int i = 0; i < _nTracks; ++i) {
        if (_track_source[i] != primaryTrackSource) {
            continue;
        }
        if (_mcpart_PDG[static_cast<int>(_track_trueID[i])] == 11) {
            // If it's an electron, take the highest pt one.
            double pt = std::sqrt(std::pow(_track_px[i], 2) + std::pow(_track_py[i], 2));
            if (pt > electronPt) {
                electronID = i;
                electronPt = pt;
            }
        }
    }
    if (electronID == 1e6) {
        throw std::runtime_error("Unable to find electron!");
    }

    // Loop over the final state hadrons to determine the kinematics.
    double px_h = 0;
    double py_h = 0;
    double e_minus_pz = 0;
    for(unsigned int i = 0; i < _nTracks; ++i) {
        if (_track_source[i] != primaryTrackSource) {
            continue;
        }
        // Don't include the electron
        // TODO: Update this when the electron ID is updated.
        if (_mcpart_PDG[static_cast<int>(_track_trueID[i])] == 11) {
            continue;
        }
        px_h += _track_px[i];
        py_h += _track_py[i];
        double eTrack = std::sqrt(std::pow(_track_px[i],2) + std::pow(_track_py[i],2) + std::pow(_track_pz[i],2) + std::pow(hadronMassHypothesis, 2));
        e_minus_pz += (eTrack - _track_pz[i]);
    }

    double y = e_minus_pz / (2 * incomingElectronEnergy);
    double Q2 = (std::pow(px_h, 2) + std::pow(py_h, 2)) / (1 - y);
    double x = Q2 / (sqrts * y);

    return DISKinematics{x, y, Q2};
}


