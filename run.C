#include <iostream>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TMath.h"

// Struct to store properties of Particles of Interest (POIs)
struct POI {
    double phi;
    int charge;
};

void run(const char* filename = "z-auau_run_1.root") {
    // 1. Open the EPOS ROOT file and extract the tree
    TFile *file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("teposevent"); 
    if (!tree) {
        std::cerr << "Error: Could not find tree 'tepos' in file." << std::endl;
        return;
    }

    // 2. Set up branch addresses
    const int MAXPTL = 150000; // Large enough for central heavy-ion collisions
    Int_t np;
    Int_t id[MAXPTL];
    Float_t px[MAXPTL], py[MAXPTL], pz[MAXPTL]; 
    // (We don't strictly need energy 'e' or mass 'm' for just eta/phi/pt)

    tree->SetBranchAddress("np", &np);
    tree->SetBranchAddress("id", id);
    tree->SetBranchAddress("px", px);
    tree->SetBranchAddress("py", py);
    tree->SetBranchAddress("pz", pz);

    // Initialize ROOT's particle database to look up charges
    TDatabasePDG *pdg = TDatabasePDG::Instance();

    // Variables for global statistics
    int valid_events = 0;
    double sum_delta_gamma = 0.0;
    std::vector<double> event_delta_gammas; // To compute variance later

    Long64_t n_events = tree->GetEntries();
    std::cout << "Processing " << n_events << " events..." << std::endl;

    // 3. Main Event Loop
    for (Long64_t ievt = 0; ievt < n_events; ++ievt) {
        tree->GetEntry(ievt);

        double Q2x = 0.0, Q2y = 0.0;
        std::vector<POI> pois;

        // --- Loop over particles in the current event ---
        for (int i = 0; i < np; ++i) {
            // Get particle charge using PDG database
            TParticlePDG *particle = pdg->GetParticle(id[i]);
            if (!particle) continue; // Skip unknown fragments/particles
            
            // TParticlePDG returns charge in units of |e|/3
            int charge = TMath::Nint(particle->Charge() / 3.0); 
            if (charge == 0) continue; // We only care about charged particles

            // Reconstruct kinematics
            TVector3 p3(px[i], py[i], pz[i]);
            double pt = p3.Pt();
            double eta = p3.Eta();
            double phi = p3.Phi();

            // Event Plane selection (Forward/Backward pseudo-rapidity)
            if (std::abs(eta) > 1.5 && std::abs(eta) < 3.0 && pt > 0.15) {
                Q2x += pt * std::cos(2.0 * phi); // pT-weighted flow vector
                Q2y += pt * std::sin(2.0 * phi);
            }

            // Particle of Interest selection (Mid-rapidity)
            if (std::abs(eta) < 1.0 && pt > 0.2 && pt < 2.0) {
                pois.push_back({phi, charge});
            }
        }

        // --- Calculate Event Plane ---
        if (Q2x == 0.0 && Q2y == 0.0) continue; // Skip event if no EP particles
        double Psi2 = 0.5 * std::atan2(Q2y, Q2x);

        // --- Calculate pairs ---
        double evt_gamma_OS = 0.0, evt_gamma_SS = 0.0;
        int n_OS = 0, n_SS = 0;

        for (size_t i = 0; i < pois.size(); ++i) {
            for (size_t j = i + 1; j < pois.size(); ++j) {
                double cos_val = std::cos(pois[i].phi + pois[j].phi - 2.0 * Psi2);
                
                if (pois[i].charge == pois[j].charge) {
                    evt_gamma_SS += cos_val;
                    n_SS++;
                } else {
                    evt_gamma_OS += cos_val;
                    n_OS++;
                }
            }
        }

        // --- Event-by-Event averages ---
        // Require at least one pair of each type to include the event
        if (n_OS > 0 && n_SS > 0) {
            double avg_OS = evt_gamma_OS / n_OS;
            double avg_SS = evt_gamma_SS / n_SS;
            double delta_gamma = avg_OS - avg_SS;

            sum_delta_gamma += delta_gamma;
            event_delta_gammas.push_back(delta_gamma);
            valid_events++;
        }
        
        if ((ievt + 1) % 100 == 0) {
            std::cout << "Processed " << (ievt + 1) << " events..." << "\r" << std::flush;
        }
    }
    std::cout << std::endl;

    // 4. Compute Final Statistics (Mean and Standard Error)
    if (valid_events == 0) {
        std::cerr << "No valid events found to compute correlator." << std::endl;
        file->Close();
        return;
    }

    double mean_delta_gamma = sum_delta_gamma / valid_events;
    
    // Compute variance of the event means
    double sum_sq_diff = 0.0;
    for (double dGamma : event_delta_gammas) {
        sum_sq_diff += (dGamma - mean_delta_gamma) * (dGamma - mean_delta_gamma);
    }
    
    // Standard error of the mean: sigma / sqrt(N)
    double std_dev = std::sqrt(sum_sq_diff / (valid_events - 1));
    double stat_error = std_dev / std::sqrt(valid_events);

    // 5. Print Results
    std::cout << "========================================" << std::endl;
    std::cout << "Analysis Complete." << std::endl;
    std::cout << "Valid Events Analyzed: " << valid_events << std::endl;
    std::cout << "Delta Gamma (OS - SS) : " << std::scientific << mean_delta_gamma 
              << " +/- " << stat_error << std::endl;
    std::cout << "========================================" << std::endl;

    file->Close();
}
