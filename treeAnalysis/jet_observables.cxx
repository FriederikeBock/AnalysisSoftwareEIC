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
    {EtaRegion_t::forward, RegionSpec{1.5, 3.5, "forward"}},
    {EtaRegion_t::midRapidity, RegionSpec{-1.5, 1.5, "mid_rapidity"}},
    {EtaRegion_t::backward, RegionSpec{-3.5, -1.5, "backward"}},
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

struct NPDFInfo {
    std::string protonPDF;
    unsigned int nVariations;
};

struct PDFContainer {
    unsigned int nVariations{1};
    std::vector<std::unique_ptr<LHAPDF::PDF>> nPDF{};
    std::unique_ptr<LHAPDF::PDF> protonPDF{};
};

const std::map<std::string, NPDFInfo> nuclearPDFInfo = {
    {"EPPS16nlo_CT14nlo_Au197", NPDFInfo{"CT14nlo", 97}},
};

std::map<std::string, std::shared_ptr<PDFContainer>> createPDFContainers(const std::vector<std::string> & pdfNames, bool useNPDFVariations)
{
    gSystem->Load("libLHAPDF");
    //this->pdf = std::make_unique<LHAPDF::PDF>(LHAPDF::mkPDF("EPPS16nlo_CT14nlo_Au197"));
    std::map<std::string, std::shared_ptr<PDFContainer>> containers;
    for (auto name : pdfNames) {
        // Nothing to initialize for the "ep" case
        if (name == "ep") {
            continue;
        }

        // Otherwise, we need to retrieve the PDF info
        unsigned int nVariations = 1;
        auto npdfInfo = nuclearPDFInfo.at(name);
        if (useNPDFVariations) {
            nVariations = npdfInfo.nVariations;
        };
        //auto [tempName, container] = containers.emplace(name, PDFContainer{nVariations});
        auto container = std::make_shared<PDFContainer>();
        container->nVariations = nVariations;
        container->nPDF.resize(nVariations);
        for (unsigned int i = 0; i < nVariations; ++i) {
            // We get the raw pointer for LHAPDF, so we want to encapsulate it in the unique_ptr
            container->nPDF.at(i).reset(LHAPDF::mkPDF(name, i));
        }
        // We also need the proton PDF for reference
        container->protonPDF.reset(LHAPDF::mkPDF(npdfInfo.protonPDF));
        std::cout << "Done making PDFs for " << name << std::endl;
        containers.emplace(name, container);
        std::cout << "Done emplacing into map for " << name << std::endl;
    }
    return containers;
}


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

std::string GetIdentifier(const double jetR, const JetType_t jetType, const RegionSpec & region, std::string name, bool eA = false, int variation = 0, std::string tag = "") {
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
        identifier += "_eA_variation" + std::to_string(variation);
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
    std::map<std::string, TH2D> spectra2D{};
    std::map<std::string, TH2D> angularity{};
    std::map<std::string, TH2D> jetHadronDPhi{};
    std::vector<TH1D> backwardHadrons{};
    std::shared_ptr<PDFContainer> pdfContainer{};
    bool initialized{false};

    JetObservables(JetType_t _jetType, const std::string & _nPDFName = "", const std::string & _tag = ""):
        jetType(_jetType),
        nPDFName{_nPDFName},
        tag(_tag),
        spectra{},
        spectra2D{},
        angularity{},
        jetHadronDPhi{},
        backwardHadrons{},
        pdfContainer{},
        initialized{false}
    {}

    bool is_eA() { return nPDFName != "ep"; }

    void Init(std::vector<double> jetRParameters, std::map<std::string, std::shared_ptr<PDFContainer>> & pdfContainers)
    {
        // Setup
        bool eA = is_eA();
        // Setup pdfs if needed
        unsigned int nVariations = 1;
        if (eA) {
            this->pdfContainer = pdfContainers.at(this->nPDFName);
            nVariations = this->pdfContainer->nPDF.size();
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
        // Add types of jets, jet R, etc to make an identifier
        for (unsigned int i = 0; i < nVariations; ++i) {
            for (auto R : jetRParameters) {
                //for (auto && [region, info] : regions) {
                for (auto & v : regions) {
                    // E spectra
                    identifier = GetIdentifier(R, jetType, v.second, "spectra_E", eA, i, tag);
                    spectra[identifier] = TH1D(identifier.c_str(), identifier.c_str(), 150, 0, 150);
                    spectra[identifier].Sumw2();
                    // p spectra
                    identifier = GetIdentifier(R, jetType, v.second, "spectra_p", eA, i, tag);
                    spectra[identifier] = TH1D(identifier.c_str(), identifier.c_str(), 150, 0, 150);
                    spectra[identifier].Sumw2();
                    // pt spectra
                    identifier = GetIdentifier(R, jetType, v.second, "spectra_pt", eA, i, tag);
                    spectra[identifier] = TH1D(identifier.c_str(), identifier.c_str(), 150, 0, 150);
                    spectra[identifier].Sumw2();

                    // Save on huge memory costs by avoiding the variations for observables where we're
                    // not as interested at the moment
                    if (i == 0) {
                        // jet p vs Q2
                        // p spectra
                        identifier = GetIdentifier(R, jetType, v.second, "spectra_p_Q2", eA, i, tag);
                        spectra2D[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, 200, 0, 1000);
                        spectra2D[identifier].Sumw2();
                        // pt spectra
                        identifier = GetIdentifier(R, jetType, v.second, "spectra_pt_Q2", eA, i, tag);
                        spectra2D[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, 200, 0, 1000);
                        spectra2D[identifier].Sumw2();
                        // jet p vs x
                        // p spectra
                        identifier = GetIdentifier(R, jetType, v.second, "spectra_p_x", eA, i, tag);
                        spectra2D[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, nLogBins, logBins.data());
                        spectra2D[identifier].Sumw2();
                        // pt spectra
                        identifier = GetIdentifier(R, jetType, v.second, "spectra_pt_x", eA, i, tag);
                        spectra2D[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, nLogBins, logBins.data());
                        spectra2D[identifier].Sumw2();

                        // Angularity
                        // a = 0 (mass)
                        // As a function of E
                        identifier = GetIdentifier(R, jetType, v.second, "angularity_a_0_E", eA, i, tag);
                        angularity[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, nLogBins, logBins.data());
                        angularity[identifier].Sumw2();
                        // As a function of p
                        identifier = GetIdentifier(R, jetType, v.second, "angularity_a_0_p", eA, i, tag);
                        angularity[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, nLogBins, logBins.data());
                        angularity[identifier].Sumw2();
                        // a = 1 (girth)
                        // As a function of E
                        identifier = GetIdentifier(R, jetType, v.second, "angularity_a_1_E", eA, i, tag);
                        angularity[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, nLogBins, logBins.data());
                        angularity[identifier].Sumw2();
                        // As a function of p
                        identifier = GetIdentifier(R, jetType, v.second, "angularity_a_1_p", eA, i, tag);
                        angularity[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, nLogBins, logBins.data());
                        angularity[identifier].Sumw2();

                        // Jet-hadron as function of E
                        identifier = GetIdentifier(R, jetType, v.second, "jet_hadron_E", eA, i, tag);
                        jetHadronDPhi[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, 72, -1.0 * TMath::Pi(), TMath::Pi());
                        jetHadronDPhi[identifier].Sumw2();

                        // Jet-hadron as function of p
                        identifier = GetIdentifier(R, jetType, v.second, "jet_hadron_p", eA, i, tag);
                        jetHadronDPhi[identifier] = TH2D(identifier.c_str(), identifier.c_str(), 150, 0, 150, 72, -1.0 * TMath::Pi(), TMath::Pi());
                        jetHadronDPhi[identifier].Sumw2();
                    }
                }
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
        std::cout << "Done with init for " << jetTypes.at(jetType) << ", pdf name: " << nPDFName << std::endl;
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
        for (auto & h : this->spectra2D) {
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
void fillEventObservables(EventObservables & eventObservables, const DISKinematics & measuredKinematics)
{
    // NOTE: The cross section isn't available in the current test production, so set to 1 if not available.
    double cross_section = _cross_section ? _cross_section : 1;

    if (measuredKinematics.x < 0 || measuredKinematics.x > 1) {
        // Return early - the kinematics are invalid, so there's nothing to be done.
        return;
    }

    // Fill the hist.
    eventObservables.Q2Measured.Fill(measuredKinematics.Q2, cross_section);
    eventObservables.Q2Measured.Fill(measuredKinematics.Q2, cross_section);
}

/**
  * Fill jet spectra.
  */
void fillJetObservables(JetObservables & observables,
                        const std::vector<fastjet::PseudoJet> & jets, double jetR,
                        DISKinematics & kinematics,
                        const std::vector<float> particles_E,
                        const std::vector<float> particles_px,
                        const std::vector<float> particles_py,
                        const std::vector<float> particles_pz)
{
    // Setup
    bool eA = observables.is_eA();
    // Cross section for weighting
    // NOTE: The cross section somtimes isn't available in the productions, so set to 1 if not available.
    // NOTE: Even when available, we should be taking the last cross section to get the smallest error. But this isn't
    //       yet trivial to access. So keep it as 1 for now.
    //double baseWeight = _cross_section ? _cross_section : 1;
    double baseWeight = 1;

    // Setup information needed for evaluating the nPDF and protonPDF
    unsigned int struckQuarkIndexHepMC = -100;
    try {
        struckQuarkIndexHepMC = findHepMCIndexOfStruckQuark();
    }
    catch (KinematicsErrors_t & e) {
        std::cout << "Unable to find struck quark. Returning...\n";
        return;
    }
    int struckQuarkFlavor = _hepmcp_PDG[struckQuarkIndexHepMC];

    // Setup for variations
    unsigned int nVariations = 1;
    double weightPDF = 1.0;
    if (eA) {
        nVariations = observables.pdfContainer->nPDF.size();
        // Can only evaluator if the x is physically plausible
        if (kinematics.x > -1 && kinematics.x <= 1) {
            double weightPDF = observables.pdfContainer->protonPDF->xfxQ2(struckQuarkFlavor, kinematics.x, kinematics.Q2);
        }
    }

    for (auto & j : jets) {
        // Acceptance
        // Fiducial jet cut
        if (std::abs(j.eta()) > (3.5 - jetR)) {
            continue;
        }

        // Require at least two constituents
        // This is based on John's Centauro studies, which suggests that requiring two constituents will remove many beam jets
        if (j.constituents().size() < 2) {
            continue;
        }

        auto region = findRegion(j.eta());

        // Check fiducial acceptnace for the region. If outside, then continue
        // Skip for now, since I'm not sure it makes sense in the way that it does for say, ALICE
        /*if (j.eta() > (region.etaMax - jetR) || j.eta() > (region.etaMin + jetR)) {
            continue;
        }*/

        for (unsigned int i = 0; i < nVariations; ++i) {
            double weight = baseWeight;
            // Determine full weight (if needed)
            if (eA) {
                //std::cout << "x=" << kinematics.x << ", y=" << kinematics.y << ", Q2=" << kinematics.Q2 << "\n";

                // We want to reweight by the nPDF effects, so we reweight by nPDF / protonPDF
                // NOTE: In this case, we don't need to scale by 1/x because it will cancel in the weight ratio
                if (kinematics.x > -1 && kinematics.x <= 1) {
                    double weightNPDF = observables.pdfContainer->nPDF.at(i)->xfxQ2(struckQuarkFlavor, kinematics.x, kinematics.Q2);
                    double ratio = weightNPDF / weightPDF;
                    if (std::isfinite(ratio)) {
                        weight *= (weightNPDF / weightPDF);
                    }
                    else {
                        std::cout << "WARNING: Non finite value for weight, not modifying: weightNPDF=" << weightNPDF
                                << ", weightPDF=" << weightPDF
                                << ", ratio=" << ratio
                                << "\n";
                        std::cout << "struck quark flavor=" << struckQuarkFlavor
                                << ", weightNPDF (xf)=" << weightNPDF
                                << ", weightNDPF / x=" << weightNPDF / kinematics.x
                                << ", weightNDPF / weightPDF =" << weightNPDF / weightPDF
                                << ", weight=" << weight
                                << ", x=" << kinematics.x
                                << ", Q2=" << kinematics.Q2
                                << "\n";

                        // Print particles to help with debugging
                        // NOTE: Only do it once if there are many variations - otherewise it would get overwhelming quickly
                        if (i == 0) {
                            printParticles();
                        }
                    }
                }
            }
            observables.spectra[GetIdentifier(jetR, observables.jetType, region, "spectra_E", eA, i, observables.tag)].Fill(j.e(), weight);
            observables.spectra[GetIdentifier(jetR, observables.jetType, region, "spectra_p", eA, i, observables.tag)].Fill(j.modp(), weight);
            observables.spectra[GetIdentifier(jetR, observables.jetType, region, "spectra_pt", eA, i, observables.tag)].Fill(j.perp(), weight);

            // Save on huge memory costs by avoiding the variations for observables where we're
            // not as interested at the moment
            if (i == 0) {
                // 2D spectra
                // spectra vs Q2
                observables.spectra2D[GetIdentifier(jetR, observables.jetType, region, "spectra_p_Q2", eA, i, observables.tag)].Fill(j.modp(), kinematics.Q2, weight);
                observables.spectra2D[GetIdentifier(jetR, observables.jetType, region, "spectra_pt_Q2", eA, i, observables.tag)].Fill(j.perp(), kinematics.Q2, weight);
                // spectra vs x
                observables.spectra2D[GetIdentifier(jetR, observables.jetType, region, "spectra_p_x", eA, i, observables.tag)].Fill(j.modp(), kinematics.x, weight);
                observables.spectra2D[GetIdentifier(jetR, observables.jetType, region, "spectra_pt_x", eA, i, observables.tag)].Fill(j.perp(), kinematics.x, weight);

                double angularity_a_0 = 0;
                double angularity_a_1 = 0;
                for (auto constituent : j.constituents()) {
                    // Angularity
                    angularity_a_0 += (constituent.pt() * std::pow(j.delta_R(constituent), 2-0));
                    angularity_a_1 += (constituent.pt() * std::pow(j.delta_R(constituent), 2-1));
                }
                observables.angularity[GetIdentifier(jetR, observables.jetType, region, "angularity_a_0_E", eA, i, observables.tag)].Fill(j.e(), angularity_a_0 / j.pt(), weight);
                observables.angularity[GetIdentifier(jetR, observables.jetType, region, "angularity_a_0_p", eA, i, observables.tag)].Fill(j.modp(), angularity_a_0 / j.pt(), weight);
                observables.angularity[GetIdentifier(jetR, observables.jetType, region, "angularity_a_1_E", eA, i, observables.tag)].Fill(j.e(), angularity_a_1 / j.pt(), weight);
                observables.angularity[GetIdentifier(jetR, observables.jetType, region, "angularity_a_1_p", eA, i, observables.tag)].Fill(j.modp(), angularity_a_1 / j.pt(), weight);

                // Jet-hadron correlations
                // Separate loop because it uses all hadrons, not just constituents
                for (std::size_t k = 0; k < particles_E.size(); ++k) {
                    fastjet::PseudoJet part(particles_px[k], particles_py[k], particles_pz[k], particles_E[k]);
                    observables.jetHadronDPhi[GetIdentifier(jetR, observables.jetType, region, "jet_hadron_E", eA, i, observables.tag)].Fill(j.e(), j.delta_phi_to(part), weight);
                    observables.jetHadronDPhi[GetIdentifier(jetR, observables.jetType, region, "jet_hadron_p", eA, i, observables.tag)].Fill(j.modp(), j.delta_phi_to(part), weight);
                }
            }
        }
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
