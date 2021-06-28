/**
  * Jet observables
  *
  * author: Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
  */

enum class EtaRegion_t {
    forward,
    midRapidity,
    backward,
};

struct RegionSpec {
    double etaMin;
    double etaMax;
    std::string name;
};

// Yeah, this is global. It's not great, but so it goes.
const std::map<EtaRegion_t, RegionSpec> regions = {
    {EtaRegion_t::forward, RegionSpec{1.5, 4.0, "forward"}},
    {EtaRegion_t::midRapidity, RegionSpec{-1.5, 1.5, "mid_rapidity"}},
    {EtaRegion_t::backward, RegionSpec{-4.0, -1.5, "backward"}},
};

RegionSpec findRegion(double eta) {
    for (auto && [region, r] : regions) {
        if (eta > r.etaMin && eta < r.etaMax) {
            return r;
        }
    }
    throw std::runtime_error("Could not find eta region");
}

enum class JetType_t {
    charged,
    calo,
    full,
};

const std::map<JetType_t, std::string> jetTypes = {
    {JetType_t::charged, "charged"},
    {JetType_t::calo, "calo"},
    {JetType_t::full, "full"},
};

struct EventObservables {
    bool initialized = false;
    TH1D Q2Measured = TH1D("q2Measured", "q2Measured", 200, 0, 200);

    void Init() {}

    void Write(std::string outputDir) {
        // define output file
        std::shared_ptr<TFile> fileOutput = std::make_shared<TFile>(TString::Format("%s/output_EventObservables.root", outputDir.c_str()), "RECREATE");

        // Write hists
        this->Q2Measured.Write();

        // Cleanup
        fileOutput->Close();
    }
};

std::string GetIdentifier(const double jetR, const JetType_t jetType, const RegionSpec & region, std::string name, std::string tag = "") {
    // Jet R
    std::string identifier = "jetR";
    identifier += jetR < 1.0 ? "0" : "";
    identifier += std::to_string(static_cast<int>(std::round(jetR * 100)));
    // tag
    if (tag != "") {
        identifier += "_" + tag;
    }
    // Type
    identifier += "_" + jetTypes.at(jetType);
    // Region
    identifier += "_" + region.name;
    // Name
    identifier += "_" + name;
    return identifier;
}

struct JetObservables {
    JetType_t jetType;
    std::string tag = "";
    std::map<std::string, TH1D> spectra;
    bool initialized = false;

    void Init(std::vector<double> jetRParameters)
    {
        std::string identifier = "";
        // Add types of jets, jet R. Make a string
        for (auto R : jetRParameters) {
            for (auto && [region, info] : regions) {
                // E spectra
                identifier = GetIdentifier(R, jetType, info, "spectra_E", tag);
                spectra[identifier] = TH1D(identifier.c_str(), identifier.c_str(), 150, 0, 150);
                spectra[identifier].Sumw2();
                // p spectra
                identifier = GetIdentifier(R, jetType, info, "spectra_p", tag);
                spectra[identifier] = TH1D(identifier.c_str(), identifier.c_str(), 150, 0, 150);
                spectra[identifier].Sumw2();

                // pt spectra
                identifier = GetIdentifier(R, jetType, info, "spectra_pt", tag);
                spectra[identifier] = TH1D(identifier.c_str(), identifier.c_str(), 150, 0, 150);
                spectra[identifier].Sumw2();
            }
        }
        initialized = true;
    }

    void Write(std::string outputDir)
    {
        // define output file
        std::shared_ptr<TFile> fileOutput = std::make_shared<TFile>(TString::Format("%s/output_JetObservables.root", outputDir.c_str()), "RECREATE");
        for (auto && [_, h] : spectra) {
            h.Write();
        }
    }
};

/**
  * Fill event level distributions
  */
void fillEventObservables(EventObservables & eventObservables, unsigned short primaryTrackSource)
{
    // NOTE: The cross section isn't available in the current test production, so set to 1 if not available.
    double cross_section = _cross_section ? _cross_section : 1;
    // Calculate the kinmeatics using J-B
    try {
        // TODO: Use args for e and p energies
        auto disKinematics = JBKinematics(primaryTrackSource, 10, 250);

        // Fill the hist.
        eventObservables.Q2Measured.Fill(disKinematics.Q2, cross_section);
        eventObservables.Q2Measured.Fill(disKinematics.Q2, cross_section);
    }
    catch (std::runtime_error & e) {
        //std::cout << e.what() << ". Continuing\n";
        return;
    }
}

/**
  * Fill jet spectra.
  */
void fillJetSpectra(JetObservables & observables, const std::vector<fastjet::PseudoJet> & jets, double jetR)
{
    // NOTE: The cross section isn't available in the current test production, so set to 1 if not available.
    // TODO: Grab last value for xsec to get best determination...
    double cross_section = _cross_section ? _cross_section : 1;
    for (auto & j : jets) {
        // Acceptance
        if (std::abs(j.eta()) > 4) {
            continue;
        }
        // TODO: Fiducial acceptance...
        auto region = findRegion(j.eta());
        observables.spectra[GetIdentifier(jetR, observables.jetType, region, "spectra_E", observables.tag)].Fill(j.e(), cross_section);
        observables.spectra[GetIdentifier(jetR, observables.jetType, region, "spectra_p", observables.tag)].Fill(j.modp(), cross_section);
        observables.spectra[GetIdentifier(jetR, observables.jetType, region, "spectra_pt", observables.tag)].Fill(j.perp(), cross_section);
    }
}

