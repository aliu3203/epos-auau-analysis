#include <iostream>
#include <iomanip>
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TMath.h"

void determineCentrality() {
    // 1. Load the EPOS ROOT files
    // Replace "event_tree_name" with the actual name of your tree
    TChain chain("teposevent"); 
    chain.Add("z-auau_run_1.root");
    chain.Add("z-auau_run_2.root");
    chain.Add("z-auau_run_3.root");
    chain.Add("z-auau_run_4.root");
    chain.Add("z-auau_run_5.root");

    // Add the rest of your 5 files...

    // Variable to link to the impact parameter branch
    float bim;
    chain.SetBranchAddress("bim", &bim);

    // 2. Create a high-resolution histogram for the impact parameter (b)
    // Au-Au collisions rarely have b > 20 fm.
    TH1D* hBim = new TH1D("hBim", "Impact Parameter (b) Distribution;b (fm);Events", 15, 0.0, 20.0);

    long nEvents = chain.GetEntries();
    std::cout << "Processing " << nEvents << " events to determine centrality..." << std::endl;

    // 3. Fill the histogram
    for (long iev = 0; iev < nEvents; ++iev) {
        chain.GetEntry(iev);
        hBim->Fill(bim);
    }

    // 4. Define the Centrality percentile bins we want to find (e.g., 5%, 10%, 20%, ...)
    // In terms of quantiles, 0-5% most central means the bottom 5% of the b distribution.
    const int nQuants = 9;
    double percentiles[nQuants] = {0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80}; // 5%, 10%, 20%, etc.
    double b_cuts[nQuants]; // Array to store the calculated impact parameter cuts

    // 5. Use ROOT's GetQuantiles to find the exact b values for these percentiles
    hBim->GetQuantiles(nQuants, b_cuts, percentiles);

    // 6. Print out the Centrality Definitions
    std::cout << "\n==================================================" << std::endl;
    std::cout << "   Centrality Definition based on Impact Parameter" << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "Centrality Class |  Impact Parameter Range (fm)" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    
    std::cout << "      0 -  5%    |      0.00  < b <  " << std::fixed << std::setprecision(2) << b_cuts[0] << std::endl;
    std::cout << "      5 - 10%    |      " << b_cuts[0] << "  < b <  " << b_cuts[1] << std::endl;
    std::cout << "     10 - 20%    |      " << b_cuts[1] << "  < b <  " << b_cuts[2] << std::endl;
    std::cout << "     20 - 30%    |      " << b_cuts[2] << "  < b <  " << b_cuts[3] << std::endl;
    std::cout << "     30 - 40%    |      " << b_cuts[3] << "  < b <  " << b_cuts[4] << std::endl;
    std::cout << "     40 - 50%    |      " << b_cuts[4] << "  < b <  " << b_cuts[5] << std::endl;
    std::cout << "     50 - 60%    |      " << b_cuts[5] << "  < b <  " << b_cuts[6] << std::endl;
    std::cout << "     60 - 70%    |      " << b_cuts[6] << "  < b <  " << b_cuts[7] << std::endl;
    std::cout << "     70 - 80%    |      " << b_cuts[7] << "  < b <  " << b_cuts[8] << std::endl;
    std::cout << "     80 - 100%   |      " << b_cuts[8] << "  < b <  " << hBim->GetXaxis()->GetXmax() << " (or max b)" << std::endl;
    std::cout << "==================================================" << std::endl;

    // 7. Draw and SAVE the canvas
    TCanvas* c1 = new TCanvas("c1", "Impact Parameter Distribution", 800, 600);
    
    // Formatting the plot
    hBim->SetFillColor(kAzure+1);
    hBim->SetStats(0); // Turns off the stat box for a cleaner look
    
    // Draw only the raw histogram
    hBim->Draw();
    
    // SAVE THE PLOT TO DISK
    c1->SaveAs("Centrality_Distribution.pdf");
    c1->SaveAs("Centrality_Distribution.png");
    
    std::cout << "\nPlot saved as Centrality_Distribution.pdf and Centrality_Distribution.png" << std::endl;
}
