/**
  * Jet observables
  *
  * author: Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
  */

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <TH1D.h>

#include <LHAPDF/LHAPDF.h>

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
    //for (auto && [region, r] : regions) {
    for (auto & v : regions) {
        if (eta > v.second.etaMin && eta < v.second.etaMax) {
            return v.second;
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

std::string GetIdentifier(const double jetR, const JetType_t jetType, const RegionSpec & region, std::string name, bool eA = false, std::string tag = "") {
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
    if (eA) {
        identifier += "_eA";
    }
    else {
        identifier += "_ep";
    }
    return identifier;
}

struct JetObservables {
    JetType_t jetType;
    std::string nPDFName = "";
    std::string tag = "";
    std::map<std::string, TH1D> spectra{};
    std::map<std::string, TH2D> angularity{};
    std::map<std::string, TH2D> jetHadronDPhi{};
    std::vector<TH1D> backwardHadrons{};
    std::unique_ptr<LHAPDF::PDF> pdf{};
    bool initialized{false};

    JetObservables(JetType_t _jetType, const std::string & _nPDFName = "", const std::string & _tag = ""):
        jetType(_jetType),
        nPDFName{_nPDFName},
        tag(_tag),
        spectra{},
        angularity{},
        jetHadronDPhi{},
        backwardHadrons{},
        pdf{},
        initialized{false}
    {}

    bool is_eA() { return nPDFName != "ep"; }

    void Init(std::vector<double> jetRParameters)
    {
        // Setup
        bool eA = is_eA();
        if (eA) {
            gSystem->Load("libLHAPDF");
            //this->pdf = std::make_unique<LHAPDF::PDF>(LHAPDF::mkPDF("EPPS16nlo_CT14nlo_Au197"));
            // We get the raw pointer for LHAPDF, so we want to encapsulate it in the unique_ptr
            this->pdf.reset(LHAPDF::mkPDF(nPDFName));
        }
        std::string identifier = "";
        // Log base 10 bins for angularity
        // From https://root-forum.cern.ch/t/how-to-define-a-log10-binning/11393
        // Double_t *xbins    = new Double_t[nbins+1];
        int nLogBins = 40;
        std::vector<double> logBins(nLogBins + 1);
        double xlogmin = std::log10(1e-4);
        double xlogmax = std::log10(1);
        double dlogx   = (xlogmax-xlogmin)/(static_cast<double>(nLogBins));
        for (int i=0; i<=nLogBins; i++) {
            double xlog = xlogmin + i*dlogx;
            logBins.at(i) = std::exp(std::log(10) * xlog);
            //std::cout << logBins[i] << " ";
        }
        //std::cout << "\n";
        //std::cout << "size: " << logBins.size() << "\n";
        // Add types of jets, jet R. Make a string
        for (auto R : jetRParameters) {
            //for (auto && [region, info] : regions) {
            for (auto & v : regions) {
                // E spectra
                identifier = GetIdentifier(R, jetType, v.second, "spectra_E", eA, tag);
                spectra[identifier] = TH1D(identifier.c_str(), identifier.c_str(), 150, 0, 150);
                spectra[identifier].Sumw2();
                // p spectra
                identifier = GetIdentifier(R, jetType, v.second, "spectra_p", eA, tag);
                spectra[identifier] = TH1D(identifier.c_str(), identifier.c_str(), 150, 0, 150);
                spectra[identifier].Sumw2();

                // pt spectra
                identifier = GetIdentifier(R, jetType, v.second, "spectra_pt", eA, tag);
                spectra[identifier] = TH1D(identifier.c_str(), identifier.c_str(), 150, 0, 150);
                spectra[identifier].Sumw2();

                // Angularity
                // a = 0 (mass)
                // As a function of E
                identifier = GetIdentifier(R, jetType, v.second, "angularity_a_0_E", eA, tag);
                angularity[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, nLogBins, logBins.data());
                angularity[identifier].Sumw2();
                // As a function of p
                identifier = GetIdentifier(R, jetType, v.second, "angularity_a_0_p", eA, tag);
                angularity[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, nLogBins, logBins.data());
                angularity[identifier].Sumw2();
                // a = 1 (girth)
                // As a function of E
                identifier = GetIdentifier(R, jetType, v.second, "angularity_a_1_E", eA, tag);
                angularity[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, nLogBins, logBins.data());
                angularity[identifier].Sumw2();
                // As a function of p
                identifier = GetIdentifier(R, jetType, v.second, "angularity_a_1_p", eA, tag);
                angularity[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, nLogBins, logBins.data());
                angularity[identifier].Sumw2();

                // Jet-hadron as function of E
                identifier = GetIdentifier(R, jetType, v.second, "jet_hadron_E", eA, tag);
                jetHadronDPhi[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, 72, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
                jetHadronDPhi[identifier].Sumw2();

                // Jet-hadron as function of p
                identifier = GetIdentifier(R, jetType, v.second, "jet_hadron_p", eA, tag);
                jetHadronDPhi[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, 72, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
                jetHadronDPhi[identifier].Sumw2();
            }
        }
        std::string name = "nBackwardHadrons";
        name += "_" + jetTypes.at(jetType);
        if (tag != "") {
            name += "_" + tag;
        }
        name += eA ? "_eA" : "_ep";
        backwardHadrons.emplace_back(TH1D(name.c_str(), name.c_str(), 100, 0, 100));
        backwardHadrons.back().Sumw2();
        initialized = true;
    }

    void Write(const std::string & outputDir, const std::string & openOption = "UPDATE")
    {
        // define output file
        std::string outputFilename = outputDir + "/output_JetObservables";
        if (is_eA()) {
            outputFilename += "_" + nPDFName;
        }
        outputFilename += ".root";
        std::shared_ptr<TFile> fileOutput = std::make_shared<TFile>(outputFilename.c_str(), openOption.c_str());
        //for (auto && [_, h] : spectra) {
        // Spectra
        for (auto & h : this->spectra) {
            h.second.Write();
        }
        // Angularity
        for (auto & h : this->angularity) {
            h.second.Write();
        }
        // Jet-hadron
        for (auto & h : this->jetHadronDPhi) {
            h.second.Write();
        }
        this->backwardHadrons[0].Write();
    }
};

/**
  * Fill event level distributions
  */
DISKinematics fillEventObservables(EventObservables & eventObservables, unsigned short primaryTrackSource)
{
    // NOTE: The cross section isn't available in the current test production, so set to 1 if not available.
    double cross_section = _cross_section ? _cross_section : 1;
    // TODO: Enable and replace the hard code below when available
    /*double incomingElectronEnergy = _hepmc_part_E[0];
    double incomingProtonEnergy = _hepmc_part_E[1];*/
    double incomingElectronEnergy = 10;
    double incomingProtonEnergy = 250;
    // Calculate the kinmeatics using J-B
    try {
        auto disKinematics = JBKinematics(primaryTrackSource, incomingElectronEnergy, incomingProtonEnergy);

        // Fill the hist.
        eventObservables.Q2Measured.Fill(disKinematics.Q2, cross_section);
        eventObservables.Q2Measured.Fill(disKinematics.Q2, cross_section);

        return disKinematics;
    }
    catch (std::runtime_error & e) {
        //std::cout << e.what() << ". Continuing\n";
        return DISKinematics{};
    }
}

/**
  * Fill jet spectra.
  */
//void fillJetObservables(JetObservables & observables, const std::vector<fastjet::PseudoJet> & jets, double jetR)
void fillJetObservables(JetObservables & observables, const std::vector<fastjet::PseudoJet> & jets, double jetR, DISKinematics & kinematics)
{
    // Setup
    bool eA = observables.is_eA();
    // NOTE: The cross section isn't available in the current test production, so set to 1 if not available.
    // TODO: Grab last value for xsec to get best determination...
    double weight = _cross_section ? _cross_section : 1;
    if (eA) {
        // TODO: Fully implement finding struck quark...
        int struckQuarkFlavor = 1;
        std::cout << "x=" << kinematics.x << ", y=" << kinematics.y << ", Q2=" << kinematics.Q2 << "\n";
        if (kinematics.x > -1 && kinematics.x < 2) {
            weight *= observables.pdf->xfxQ(struckQuarkFlavor, kinematics.x, std::sqrt(kinematics.Q2)) / kinematics.x;
        }
    }
    for (auto & j : jets) {
        // Acceptance
        //if (std::abs(j.eta()) > (4 - jetR)) {
        if (std::abs(j.eta()) > 4) {
            continue;
        }

        auto region = findRegion(j.eta());
        // Check fiducial acceptnace for the region. If outside, then continue
        // Skip for now, since I'm not sure it makes sense in the way that it does for HI collisions
        /*if (j.eta() > (region.etaMax - jetR) || j.eta() > (region.etaMin + jetR)) {
            continue;
        }*/
        observables.spectra[GetIdentifier(jetR, observables.jetType, region, "spectra_E", eA, observables.tag)].Fill(j.e(), weight);
        observables.spectra[GetIdentifier(jetR, observables.jetType, region, "spectra_p", eA, observables.tag)].Fill(j.modp(), weight);
        observables.spectra[GetIdentifier(jetR, observables.jetType, region, "spectra_pt", eA, observables.tag)].Fill(j.perp(), weight);

        double angularity_a_0 = 0;
        double angularity_a_1 = 0;
        for (auto constituent : j.constituents()) {
            // Angularity
            angularity_a_0 += (constituent.pt() * std::pow(j.delta_R(constituent), 2-0));
            angularity_a_1 += (constituent.pt() * std::pow(j.delta_R(constituent), 2-1));

            // Jet-hadron correlations
            // TODO: Use all hadrons, not just constituents. Duh...
            observables.jetHadronDPhi[GetIdentifier(jetR, observables.jetType, region, "jet_hadron_E", eA, observables.tag)].Fill(j.e(), j.delta_phi_to(constituent), weight);
            observables.jetHadronDPhi[GetIdentifier(jetR, observables.jetType, region, "jet_hadron_p", eA, observables.tag)].Fill(j.modp(), j.delta_phi_to(constituent), weight);
        }
        observables.angularity[GetIdentifier(jetR, observables.jetType, region, "angularity_a_0_E", eA, observables.tag)].Fill(j.e(), angularity_a_0 / j.pt(), weight);
        observables.angularity[GetIdentifier(jetR, observables.jetType, region, "angularity_a_0_p", eA, observables.tag)].Fill(j.modp(), angularity_a_0 / j.pt(), weight);
        observables.angularity[GetIdentifier(jetR, observables.jetType, region, "angularity_a_1_E", eA, observables.tag)].Fill(j.e(), angularity_a_1 / j.pt(), weight);
        observables.angularity[GetIdentifier(jetR, observables.jetType, region, "angularity_a_1_p", eA, observables.tag)].Fill(j.modp(), angularity_a_1 / j.pt(), weight);
    }
}

void fillHadronObservables(JetObservables & observables)
{
    // True hadrons going in the backwards direction
    unsigned int nHadrons = 0;
    for (std::size_t i = 0; i < (std::size_t)_nMCPart; ++i)
    {
        if (_mcpart_Eta[i] < -1.5 && _mcpart_E[i] > 0.2 && _mcpart_ID[i] != 11) {
            ++nHadrons;
        }
    }
    observables.backwardHadrons[0].Fill(nHadrons);
}
