#include <iostream>
#include <vector>
#include <cmath>
#include "TChain.h"
#include "TMath.h"

// Struct to hold calculated particle kinematics
struct Particle {
    double pt;
    double eta;
    double phi;
    int charge;
};

// Helper function to get charge from exact EPOS IDs (from idt.dt)
int GetChargeFromEposId(int epos_id) {
    int EPOS_PI_PLUS  = 120;   
    int EPOS_PI_MINUS = -120;  
    int EPOS_K_PLUS   = 130;   
    int EPOS_K_MINUS  = -130;  
    int EPOS_P        = 1120;  
    int EPOS_P_BAR    = -1120; 

    if (epos_id == EPOS_PI_PLUS || epos_id == EPOS_K_PLUS || epos_id == EPOS_P) return 1;
    if (epos_id == EPOS_PI_MINUS || epos_id == EPOS_K_MINUS || epos_id == EPOS_P_BAR) return -1;
    
    return 0; // Neutral or unlisted species
}

void run() {
    // ====================================================================
    // INPUT YOUR 30-40% CENTRALITY CUTS HERE (From DetermineCentrality.C)
    // ====================================================================
    double b_min = 7.71; // Replace with your actual 30% lower bound
    double b_max = 9.33; // Replace with your actual 40% upper bound
    // ====================================================================

    // 1. Load the EPOS ROOT files
    TChain chain("teposevent"); 
    chain.Add("z-auau_run_1.root");
    chain.Add("z-auau_run_2.root");
    chain.Add("z-auau_run_3.root");
    chain.Add("z-auau_run_4.root");
    chain.Add("z-auau_run_5.root");

    const int MAX_PARTICLES = 50000; 
    
    int np;
    float bim;
    float px[MAX_PARTICLES];
    float py[MAX_PARTICLES]; 
    float pz[MAX_PARTICLES];
    int id[MAX_PARTICLES];
    int ist[MAX_PARTICLES];

    chain.SetBranchAddress("np", &np);
    chain.SetBranchAddress("bim", &bim);
    chain.SetBranchAddress("px", px);
    chain.SetBranchAddress("py", py);
    chain.SetBranchAddress("pz", pz);
    chain.SetBranchAddress("id", id);
    chain.SetBranchAddress("ist", ist);

    // Correlator accumulators
    double sumGammaOS = 0.0;
    double sumGammaSS = 0.0;
    long long totalPairsOS = 0;
    long long totalPairsSS = 0;

    double sumSqGammaOS = 0.0;
    double sumSqGammaSS = 0.0;

    long nEvents = chain.GetEntries();
    std::cout << "Processing " << nEvents << " events..." << std::endl;
    std::cout << "Applying 30-40% Centrality cut: " << b_min << " < b < " << b_max << " fm" << std::endl;

    long eventsAnalyzed = 0;

    for (long iev = 0; iev < nEvents; ++iev) {
        chain.GetEntry(iev);

        // --- CENTRALITY SELECTION ---
        // Skip event if it is outside the 30-40% impact parameter range
        if (bim < b_min || bim >= b_max) continue; 
        
        eventsAnalyzed++;

        double Q2x = 0.0;
        double Q2y = 0.0;
        std::vector<Particle> posTracks;
        std::vector<Particle> negTracks;

        // Particle Loop
        for (int i = 0; i < np; ++i) {
            // ONLY keep hadrons of the last generation
            if (ist[i] != 0) continue; 

            int charge = GetChargeFromEposId(id[i]);
            if (charge == 0) continue; 

            // Calculate kinematics
            double pt_i = std::sqrt(px[i]*px[i] + py[i]*py[i]);
            double p_i = std::sqrt(pt_i*pt_i + pz[i]*pz[i]);
            double phi_i = std::atan2(py[i], px[i]);
            
            double eta_i = 0.0;
            if (p_i != pz[i] && p_i != -pz[i]) {
                eta_i = 0.5 * std::log((p_i + pz[i]) / (p_i - pz[i]));
            } else continue;

            // Event Plane selection (1.0 < |eta| < 3.0)
            if (std::abs(eta_i) > 1.0 && std::abs(eta_i) < 3.0) {
                Q2x += pt_i * std::cos(2.0 * phi_i);
                Q2y += pt_i * std::sin(2.0 * phi_i);
            }

            // Particle of Interest selection (|eta| < 1.0, pT > 0.2)
            if (std::abs(eta_i) < 1.0 && pt_i > 0.2) {
                Particle p_obj = {pt_i, eta_i, phi_i, charge};
                if (charge > 0) posTracks.push_back(p_obj);
                if (charge < 0) negTracks.push_back(p_obj);
            }
        }

        // Reconstruct Event Plane Angle
        double Psi2 = 0.5 * std::atan2(Q2y, Q2x);

        // --- Calculate OS Correlator (+-) ---
        for (size_t i = 0; i < posTracks.size(); ++i) {
            for (size_t j = 0; j < negTracks.size(); ++j) {
                double cosTerm = std::cos(posTracks[i].phi + negTracks[j].phi - 2.0 * Psi2);
                sumGammaOS += cosTerm;
                sumSqGammaOS += cosTerm * cosTerm;
                totalPairsOS++;
            }
        }

        // --- Calculate SS Correlator (++ and --) ---
        for (size_t i = 0; i < posTracks.size(); ++i) {
            for (size_t j = i + 1; j < posTracks.size(); ++j) {
                double cosTerm = std::cos(posTracks[i].phi + posTracks[j].phi - 2.0 * Psi2);
                sumGammaSS += cosTerm;
                sumSqGammaSS += cosTerm * cosTerm;
                totalPairsSS++;
            }
        }
        for (size_t i = 0; i < negTracks.size(); ++i) {
            for (size_t j = i + 1; j < negTracks.size(); ++j) {
                double cosTerm = std::cos(negTracks[i].phi + negTracks[j].phi - 2.0 * Psi2);
                sumGammaSS += cosTerm;
                sumSqGammaSS += cosTerm * cosTerm;
                totalPairsSS++;
            }
        }
    }

    // Final calculations
    double gammaOS = totalPairsOS > 0 ? sumGammaOS / totalPairsOS : 0.0;
    double gammaSS = totalPairsSS > 0 ? sumGammaSS / totalPairsSS : 0.0;

    double errOS = totalPairsOS > 0 ? std::sqrt((sumSqGammaOS / totalPairsOS - gammaOS * gammaOS) / totalPairsOS) : 0.0;
    double errSS = totalPairsSS > 0 ? std::sqrt((sumSqGammaSS / totalPairsSS - gammaSS * gammaSS) / totalPairsSS) : 0.0;
    
    double deltaGamma = gammaOS - gammaSS;
    double errDeltaGamma = std::sqrt(errOS * errOS + errSS * errSS);

    std::cout << "\n--- Final Results (30-40% Centrality) ---" << std::endl;
    std::cout << "Events actually analyzed: " << eventsAnalyzed << " out of " << nEvents << std::endl;
    std::cout << "Gamma OS: " << gammaOS << " +/- " << errOS << std::endl;
    std::cout << "Gamma SS: " << gammaSS << " +/- " << errSS << std::endl;
    std::cout << "Delta Gamma: " << deltaGamma << " +/- " << errDeltaGamma << std::endl;
}
