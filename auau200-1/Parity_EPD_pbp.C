using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

const int nHar = 2;
const int opt_boost = 0;
const float PI = TMath::Pi();
const float mp = 0.938272;
const float mk = 0.493667;
const float mpi= 0.13957;
const int ID_pip = 120, ID_pim = -120, ID_Kp = 130, ID_Km = -130, ID_pp = 1120, ID_pb = -1120;
const float pt_trig_up = 2;
const float pt_trig_lo = .2;  //0.2
const float EtaCut = 1.;
const float pt_asso_up = 2.0;//1;
const float pt_asso_lo = 0.2;//0.15;
const float EtaAssoCut = 5.0;
//const float cenDef[9] = {13.65, 12.76, 11.81, 10.78, 9.65, 8.36, 6.82, 4.83, 3.42};	// 200 GeV
//const float cenDef[9] = {13.14,12.29,11.38,10.39,9.29,8.05,6.57,4.65,3.29};   // 200 GeV Au+Au
const float cenDef[9] = {14.4,13.6, 12.8,11.6,10.8,9.2,8,6.8,5.2};   // 200 GeV Au+Au


void Parity_EPD_pbp(int cen = 5, int job = 0, const Char_t* inFileName = "test.list"){	//main_function

        

        float Eweight = 1;
        TChain* chain = new TChain("teposevent");
        chain->Add("z-auau_run_*.root");
        //std::cout << "Number of events in chain: " << chain->GetEntries() << endl;

        char fname_out[200];
        sprintf(fname_out,"cen%d.AVFD_Parity1%d%d_EPD_pT02_lab_job%d.root",cen,nHar-1,nHar,job);
        if(opt_boost) sprintf(fname_out,"cen%d.AVFD_Parity1%d%d_EPD_pT02_Boost_job%d.root",cen,nHar-1,nHar,job);
        TFile fout(fname_out,"RECREATE");
        //defining variables
        Int_t   Centrality, NPTracks, fParticles_, QB_index;	//
        Float_t psi, b, Charge, Eta, Theta, Phi, Pt, Mass, yy;		//track info	
        Float_t Px, Py, Pz, PID, energy;

        // Maps to store charge and names, keyed by Particle ID
        std::map<int, float> pCharge;
        std::map<int, std::string> pName;

        // Read the idt.dt file
        std::ifstream file("idt.dt");
        std::string line;
        while(std::getline(file, line)) {
                // Skip empty lines or comments starting with '!'
                if(line.empty() || line[0] == '!') continue;
                
                std::istringstream iss(line);
                int id_epos, id_pdg, id_qgsjet, id_gheisha, id_sibyll;
                std::string name;
                int ifl1, ifl2, ifl3, counter;
                float mass, charge;
                
                // Parse up to the 12th column (Charge)
                if(iss >> id_epos >> id_pdg >> id_qgsjet >> id_gheisha >> id_sibyll 
                        >> name >> ifl1 >> ifl2 >> ifl3 >> counter >> mass >> charge) {
                        
                        // NOTE: If your PIDAsso uses EPOS IDs rather than PDG IDs, 
                        // change `id_pdg` to `id_epos` inside the brackets below.
                        pCharge[id_epos] = charge; 
                        pName[id_epos] = name;
                }
        }

        // Create a histogram that automatically extends when given new string labels
        TH1F *hSpeciesAll = new TH1F("hSpeciesAll", "All Charged Particles", 1, 0, 1);
        hSpeciesAll->SetCanExtend(TH1::kAllAxes);

        //defining histograms
        TH1D *hCentrality = new TH1D("hCentrality","hCentrality",10,0,10);
        //TH1D *hMult = new T1D("hMult","hMult",1000,-0.5,999.5);
        TH1D *hMult = new TH1D("hMult","hMult",100,-0.5,99.5);

        TH1D *hMult_All = new TH1D("hMult_All", "Multiplicity (All Particles)", 1000, -0.5, 999.5); // Revert to 999.5 tmrw
        //TH1D *hPt_All = new TH1D("hPt_All", "Pt (All Particles)", 100, -0.5, 99.5);
        //TH1D *hPhi_All = new TH1D("hPt_All", )

        TH1D *hMult_ptCut = new TH1D("hMult_ptCut", "Multiplicity (pT > 0.05)", 100, -0.5, 99.5);
        TH1D *hMult_ptEtaCut = new TH1D("hMult_ptEtaCut", "Multiplicity (pT > 0.05 & |eta| < 0.5)", 100, -0.5, 99.5);

        TProfile *Hist_Parent_count = new TProfile("Hist_Parent_count","Hist_Parent_count",4,0.5,4.5,0,9999999);
        TH1D* Hist_Y = new TH1D("Hist_Y", "Hist_Y", 26, -6, 6);
        TH2D *hEtaPtDist = new TH2D("EtaPtDist","EtaPtDist",26, -6, 6,300,0,15);
        TH1D* Hist_Pt = new TH1D("Hist_Pt","Hist_Pt",60,0,3);
        TH1D* Hist_Pt_parent = new TH1D("Hist_Pt_parent","Hist_Pt_parent",300,0,15);
        TH1D* Hist_Eta_parent = new TH1D("Hist_Eta_parent","Hist_Eta_parent",60,-3,3);
        TH1D* Hist_Phi = new TH1D("Hist_Phi","Hist_Phi",72,-PI,PI);
        TH1D* Hist_Phi_parent = new TH1D("Hist_Phi_parent","Hist_Phi_parent",72,-PI,PI);
        TH1D* Hist_Mass_parent_SS = new TH1D("Hist_Mass_parent_SS","Hist_Mass_parent_SS",500,0,5);
        TH1D* Hist_Mass_parent_OS = new TH1D("Hist_Mass_parent_OS","Hist_Mass_parent_OS",500,0,5);
        TProfile *Hist_cos = new TProfile("Hist_cos","Hist_cos",6,0.5,6.5,-1,1,"");
        TH1D* hDpt   = new TH1D("hDpt","hDpt",200,0,2);

        TProfile *Hist_v2_pt_pip = new TProfile("Hist_v2_pt_pip","Hist_v2_pt_pip",100,0,5,-1,1,"");
        TProfile *Hist_v2_pt_pim = new TProfile("Hist_v2_pt_pim","Hist_v2_pt_pim",100,0,5,-1,1,"");
        TProfile *Hist_v2_pt_Kp = new TProfile("Hist_v2_pt_Kp","Hist_v2_pt_Kp",100,0,5,-1,1,"");
        TProfile *Hist_v2_pt_Km = new TProfile("Hist_v2_pt_Km","Hist_v2_pt_Km",100,0,5,-1,1,"");
        TProfile *Hist_v2_pt_pp = new TProfile("Hist_v2_pt_pp","Hist_v2_pt_pp",100,0,5,-1,1,"");
        TProfile *Hist_v2_pt_pm = new TProfile("Hist_v2_pt_pm","Hist_v2_pt_pm",100,0,5,-1,1,"");

        TProfile *Hist_a1_pt_pip = new TProfile("Hist_a1_pt_pip","Hist_a1_pt_pip",100,0,5,-1,1,"");
        TProfile *Hist_a1_pt_pim = new TProfile("Hist_a1_pt_pim","Hist_a1_pt_pim",100,0,5,-1,1,"");
        TProfile *Hist_a1_pt_Kp = new TProfile("Hist_a1_pt_Kp","Hist_a1_pt_Kp",100,0,5,-1,1,"");
        TProfile *Hist_a1_pt_Km = new TProfile("Hist_a1_pt_Km","Hist_a1_pt_Km",100,0,5,-1,1,"");
        TProfile *Hist_a1_pt_pp = new TProfile("Hist_a1_pt_pp","Hist_a1_pt_pp",100,0,5,-1,1,"");
        TProfile *Hist_a1_pt_pm = new TProfile("Hist_a1_pt_pm","Hist_a1_pt_pm",100,0,5,-1,1,"");

        TProfile *Hist_a3_pt_pip = new TProfile("Hist_a3_pt_pip","Hist_a3_pt_pip",100,0,5,-1,1,"");
        TProfile *Hist_a3_pt_pim = new TProfile("Hist_a3_pt_pim","Hist_a3_pt_pim",100,0,5,-1,1,"");
        TProfile *Hist_a3_pt_Kp = new TProfile("Hist_a3_pt_Kp","Hist_a3_pt_Kp",100,0,5,-1,1,"");
        TProfile *Hist_a3_pt_Km = new TProfile("Hist_a3_pt_Km","Hist_a3_pt_Km",100,0,5,-1,1,"");
        TProfile *Hist_a3_pt_pp = new TProfile("Hist_a3_pt_pp","Hist_a3_pt_pp",100,0,5,-1,1,"");
        TProfile *Hist_a3_pt_pm = new TProfile("Hist_a3_pt_pm","Hist_a3_pt_pm",100,0,5,-1,1,"");

        TProfile *Hist_v1_y_pip = new TProfile("Hist_v1_y_pip","Hist_v1_y_pip",20,-1,1,-1,1,"");
        TProfile *Hist_v1_y_pim = new TProfile("Hist_v1_y_pim","Hist_v1_y_pim",20,-1,1,-1,1,"");
        TProfile *Hist_v1_y_Kp = new TProfile("Hist_v1_y_Kp","Hist_v1_y_Kp",20,-1,1,-1,1,"");
        TProfile *Hist_v1_y_Km = new TProfile("Hist_v1_y_Km","Hist_v1_y_Km",20,-1,1,-1,1,"");
        TProfile *Hist_v1_y_pp = new TProfile("Hist_v1_y_pp","Hist_v1_y_pp",20,-1,1,-1,1,"");
        TProfile *Hist_v1_y_pm = new TProfile("Hist_v1_y_pm","Hist_v1_y_pm",20,-1,1,-1,1,"");

        TProfile *Hist_v2_y_pip = new TProfile("Hist_v2_y_pip","Hist_v2_y_pip",20,-1,1,-1,1,"");
        TProfile *Hist_v2_y_pim = new TProfile("Hist_v2_y_pim","Hist_v2_y_pim",20,-1,1,-1,1,"");
        TProfile *Hist_v2_y_Kp = new TProfile("Hist_v2_y_Kp","Hist_v2_y_Kp",20,-1,1,-1,1,"");
        TProfile *Hist_v2_y_Km = new TProfile("Hist_v2_y_Km","Hist_v2_y_Km",20,-1,1,-1,1,"");
        TProfile *Hist_v2_y_pp = new TProfile("Hist_v2_y_pp","Hist_v2_y_pp",20,-1,1,-1,1,"");
        TProfile *Hist_v2_y_pm = new TProfile("Hist_v2_y_pm","Hist_v2_y_pm",20,-1,1,-1,1,"");

        TProfile *Hist_a1_y_pip = new TProfile("Hist_a1_y_pip","Hist_a1_y_pip",20,-1,1,-1,1,"");
        TProfile *Hist_a1_y_pim = new TProfile("Hist_a1_y_pim","Hist_a1_y_pim",20,-1,1,-1,1,"");
        TProfile *Hist_a1_y_Kp = new TProfile("Hist_a1_y_Kp","Hist_a1_y_Kp",20,-1,1,-1,1,"");
        TProfile *Hist_a1_y_Km = new TProfile("Hist_a1_y_Km","Hist_a1_y_Km",20,-1,1,-1,1,"");
        TProfile *Hist_a1_y_pp = new TProfile("Hist_a1_y_pp","Hist_a1_y_pp",20,-1,1,-1,1,"");
        TProfile *Hist_a1_y_pm = new TProfile("Hist_a1_y_pm","Hist_a1_y_pm",20,-1,1,-1,1,"");

        TProfile *Hist_a2_y_pip = new TProfile("Hist_a2_y_pip","Hist_a2_y_pip",20,-1,1,-1,1,"");
        TProfile *Hist_a2_y_pim = new TProfile("Hist_a2_y_pim","Hist_a2_y_pim",20,-1,1,-1,1,"");
        TProfile *Hist_a2_y_Kp = new TProfile("Hist_a2_y_Kp","Hist_a2_y_Kp",20,-1,1,-1,1,"");
        TProfile *Hist_a2_y_Km = new TProfile("Hist_a2_y_Km","Hist_a2_y_Km",20,-1,1,-1,1,"");
        TProfile *Hist_a2_y_pp = new TProfile("Hist_a2_y_pp","Hist_a2_y_pp",20,-1,1,-1,1,"");
        TProfile *Hist_a2_y_pm = new TProfile("Hist_a2_y_pm","Hist_a2_y_pm",20,-1,1,-1,1,"");

        TProfile *Hist_a3_y_pip = new TProfile("Hist_a3_y_pip","Hist_a3_y_pip",20,-1,1,-1,1,"");
        TProfile *Hist_a3_y_pim = new TProfile("Hist_a3_y_pim","Hist_a3_y_pim",20,-1,1,-1,1,"");
        TProfile *Hist_a3_y_Kp = new TProfile("Hist_a3_y_Kp","Hist_a3_y_Kp",20,-1,1,-1,1,"");
        TProfile *Hist_a3_y_Km = new TProfile("Hist_a3_y_Km","Hist_a3_y_Km",20,-1,1,-1,1,"");
        TProfile *Hist_a3_y_pp = new TProfile("Hist_a3_y_pp","Hist_a3_y_pp",20,-1,1,-1,1,"");
        TProfile *Hist_a3_y_pm = new TProfile("Hist_a3_y_pm","Hist_a3_y_pm",20,-1,1,-1,1,"");

        TProfile *Hist_v2_pt = new TProfile("Hist_v2_pt","Hist_v2_pt",300,0,15,-100,100,"");
        TProfile *Hist_v2_pt_obs= new TProfile("Hist_v2_pt_obs","Hist_v2_pt_obs",300,0,15,-100,100,"");
        TProfile *Hist_v1_eta = new TProfile("Hist_v1_eta","Hist_v1_eta",300,-1.5,1.5,-100,100,"");
        TProfile *Hist_v2_eta = new TProfile("Hist_v2_eta","Hist_v2_eta",300,-1.5,1.5,-100,100,"");
        TProfile *Hist_v3_eta = new TProfile("Hist_v3_eta","Hist_v3_eta",300,-1.5,1.5,-100,100,"");
        TProfile *Hist_v2_eta_obs = new TProfile("Hist_v2_eta_obs1","Hist_v2_eta_obs1",300,-1.5,1.5,-100,100,"");

        TProfile *Hist_v2_pt_parent = new TProfile("Hist_v2_pt_parent","Hist_v2_pt_parent",300,0,15,-100,100,"");
        TProfile *Hist_v2_pt_parent_obs= new TProfile("Hist_v2_pt_parent_obs","Hist_v2_pt_parent_obs",300,0,15,-100,100,"");
        TProfile *Hist_v2_eta_parent = new TProfile("Hist_v2_eta_parent","Hist_v2_eta_parent",300,-1.5,1.5,-100,100,"");
        TProfile *Hist_v2_eta_parent_obs = new TProfile("Hist_v2_eta_parent_obs","Hist_v2_eta_parent_obs",300,-1.5,1.5,-100,100,"");

        TProfile *pParity_int_ss = new TProfile("Parity_int_ss","Parity_int_ss",8,0.5,8.5,-100,100,"");
        TProfile *pParity_int_ss_obs = new TProfile("Parity_int_ss_obs","Parity_int_ss_obs",8,0.5,8.5,-100,100,"");
        TProfile2D *pParity_eta_ss = new TProfile2D("Parity_eta_ss","Parity_eta_ss",12,0.5,12.5,20,-1,1,-100,100,"");
        TProfile2D *pParity_eta_ss_obs = new TProfile2D("Parity_eta_ss_obs","Parity_eta_ss_obs",12,0.5,12.5,20,-1,1,-100,100,"");
        TProfile2D *pParity_Deta_ss = new TProfile2D("Parity_Deta_ss","Parity_Deta_ss",12,0.5,12.5,20,0,2,-100,100,"");
        TProfile2D *pParity_Deta_ss_obs = new TProfile2D("Parity_Deta_ss_obs","Parity_Deta_ss_obs",12,0.5,12.5,20,0,2,-100,100,"");
        TProfile2D *pParity_pt_ss  = new TProfile2D("Parity_pt_ss","Parity_pt_ss",12,0.5,12.5,20,0,2.0,-100,100,"");
        TProfile2D *pParity_pt_ss_obs  = new TProfile2D("Parity_pt_ss_obs","Parity_pt_ss_obs",12,0.5,12.5,20,0,2.0,-100,100,"");
        TProfile2D *pParity_Dpt_ss = new TProfile2D("Parity_Dpt_ss","Parity_Dpt_ss",12,0.5,12.5,200,0,2.0,-100,100,"");
        TProfile2D *pParity_Dpt_ss_obs = new TProfile2D("Parity_Dpt_ss_obs","Parity_Dpt_ss_obs",12,0.5,12.5,200,0,2.0,-100,100,"");

        TProfile *pDelta_int_ss = new TProfile("Delta_int_ss","Delta_int_ss",4,0.5,4.5,-100,100,"");
        TProfile2D *pDelta_eta_ss = new TProfile2D("Delta_eta_ss","Delta_eta_ss",12,0.5,12.5,20,-1,1,-100,100,"");
        TProfile2D *pDelta_Deta_ss = new TProfile2D("Delta_Deta_ss","Delta_Deta_ss",12,0.5,12.5,20,0,2,-100,100,"");
        TProfile2D *pDelta_pt_ss  = new TProfile2D("Delta_pt_ss","Delta_pt_ss",12,0.5,12.5,20,0,2.0,-100,100,"");
        TProfile2D *pDelta_Dpt_ss  = new TProfile2D("Delta_Dpt_ss","Delta_Dpt_ss",12,0.5,12.5,20,0,2.0,-100,100,"");

        TProfile *pTemp_pT = new TProfile("pTemp_pT","pTemp_pT",2,0.5,2.5,0,10,"");
        TProfile *pTemp_a1 = new TProfile("pTemp_a1","pTemp_a1",2,0.5,2.5,-1,1,"");
        TProfile *pTemp_v2 = new TProfile("pTemp_v2","pTemp_v2",4,0.5,4.5,-100,100,"");
        TProfile *pTemp_v2_parent = new TProfile("pTemp_v2_parent","pTemp_v2_parent",2,0.5,2.5,-100,100,"");
        TProfile *pTemp_proton= new TProfile("pTemp_proton","pTemp_proton",8,0.5,8.5,-1,1,"");
        TProfile *pTemp_pion= new TProfile("pTemp_pion","pTemp_pion",8,0.5,8.5,-1,1,"");
        TProfile *pTemp_kaon= new TProfile("pTemp_kaon","pTemp_kaon",8,0.5,8.5,-1,1,"");
        TProfile *pTemp_parity = new TProfile("pTemp_parity","pTemp_parity",8,0.5,8.5,-100,100,"");
        TProfile *pTemp_parity2 = new TProfile("pTemp_parity2","pTemp_parity2",8,0.5,8.5,-100,100,"");
        TProfile *pTemp_delta = new TProfile("pTemp_delta","pTemp_delta",8,0.5,8.5,-100,100,"");

        TProfile *Hist_v2_v2parent = new TProfile("Hist_v2_v2parent","v2 parent, v2 correlation",200,-100,100,-100,100);
        TProfile *Hist_v2_v2parent_obs = new TProfile("Hist_v2_v2parent_obs","v2 parent, v2 correlation",200,-100,100,-100,100);
        TH1D* Hist_Q = new TH1D("Hist_Q","Hist_Q",250,0,25);
        TH1D* Hist_Q_test = new TH1D("Hist_Q_test","Hist_Q_test",250,0,5);
        TProfile *p_RefMult_Q = new TProfile("p_RefMult_Q","p_RefMult_Q",250,0,25,0,1000,"");
        TProfile *p_cos_Q = new TProfile("p_cos_Q","p_cos_Q",250,0,25,-1,1,"");
        TProfile *p_v2_Q = new TProfile("p_v2_Q","p_v2_Q",250,0,25,-100,100,"");
        TProfile *p_v2_Q_obs = new TProfile("p_v2_Q_obs","p_v2_Q_obs",250,0,25,-100,100,"");
        TProfile *p_v2_s_Q = new TProfile("p_v2_s_Q","p_v2_s_Q",250,0,25,-100,100,"");
        TProfile *p_v2_s_Q_obs = new TProfile("p_v2_s_Q_obs","p_v2_s_Q_obs",250,0,25,-100,100,"");
        TProfile2D *pParity_Q = new TProfile2D("Parity_Q","Parity_Q",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile2D *pParity_Q_obs = new TProfile2D("Parity_Q_obs","Parity_Q_obs",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile *pParity_SS_Q_obs = new TProfile("pParity_SS_Q_obs","pParity_SS_Q_obs",250, 0, 25, -100,100,"");
        TProfile *pParity_OS_Q_obs = new TProfile("pParity_OS_Q_obs","pParity_OS_Q_obs",250, 0, 25, -100,100,"");

        TProfile2D *pDelta_Q = new TProfile2D("Delta_Q","Delta_Q",4,0.5,4.5,250, 0, 25, -100,100,"");
        TH1D* Hist_Q2 = new TH1D("Hist_Q2","Hist_Q2",250,0,25);
        TH1D* Hist_Q4 = new TH1D("Hist_Q4","Hist_Q4",250,0,25);
        TProfile *p_RefMult_Q2 = new TProfile("p_RefMult_Q2","p_RefMult_Q2",250,0,25,0,1000,"");
        TProfile *p_cos_Q2 = new TProfile("p_cos_Q2","p_cos_Q2",250,0,25,-1,1,"");
        TProfile *p_v2_Q2 = new TProfile("p_v2_Q2","p_v2_Q2",250,0,25,-100,100,"");
        TProfile *p_v2_Q2_obs = new TProfile("p_v2_Q2_obs","p_v2_Q2_obs",250,0,25,-100,100,"");
        TProfile *p_v2_p_Q2 = new TProfile("p_v2_p_Q2","p_v2_p_Q2",250,0,25,-100,100,"");
        TProfile *p_v2_p_Q2_obs = new TProfile("p_v2_p_Q2_obs","p_v2_p_Q2_obs",250,0,25,-100,100,"");
        TProfile2D *pParity_Q2 = new TProfile2D("Parity_Q2","Parity_Q2",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile2D *pParity_Q2_obs = new TProfile2D("Parity_Q2_obs","Parity_Q2_obs",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile2D *pDelta_Q2 = new TProfile2D("Delta_Q2","Delta_Q2",4,0.5,4.5,250, 0, 25, -100,100,"");

        TH1D* Hist_QB	= new TH1D("Hist_QB","Hist_QB",250,0,10);
        TProfile* p_RefMult_QB = new TProfile("p_RefMult_QB","p_RefMult_QB",250,0,10,0,1000,"");
        TH2D* Hist_Q2_QB2 = new TH2D("Hist_Q2_QB2","Hist_Q2_QB2",250,0,10,250,0,10);
        TProfile* p_Q2_QB2= new TProfile("p_Q2_QB2","p_Q2_QB2",250,0,10,0,25);
        TH2D* Hist_Q_QB = new TH2D("Hist_Q_QB","Hist_Q_QB",250,0,10,250,0,10);
        TProfile* p_Q_QB= new TProfile("p_Q_QB","p_Q_QB",250,0,10,0,25);
        TH2D* Hist_PsiQ_PsiQB = new TH2D("Hist_PsiQ_PsiQB","Hist_PsiQ_PsiQB",72,-PI/2.,PI/2.,72,-PI/2.,PI/2.);
        TProfile2D* pMult_PsiQ_PsiQB = new TProfile2D("pMult_PsiQ_PsiQB","pMult_PsiQ_PsiQB",72,-PI/2.,PI/2.,72,-PI/2.,PI/2.,0,2000,"");
        TProfile2D* pv2_PsiQ_PsiQB = new TProfile2D("pv2_PsiQ_PsiQB","pv2_PsiQ_PsiQB",72,-PI/2.,PI/2.,72,-PI/2.,PI/2.,-100,100,"");

        TProfile* p_cos_QB = new TProfile("p_cos_QB","p_cos_QB",250,0,25,-1,1);
        TProfile *p_v2_QB = new TProfile("p_v2_QB","p_v2_QB",250,0,25,-100,100,"");
        TProfile *p_v2_QB_obs = new TProfile("p_v2_QB_obs","p_v2_QB_obs",250,0,25,-100,100,"");
        TProfile2D *pParity_QB = new TProfile2D("Parity_QB","Parity_QB",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile2D *pParity_QB_obs = new TProfile2D("Parity_QB_obs","Parity_QB_obs",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile2D *pDelta_QB = new TProfile2D("Delta_QB","Delta_QB",4,0.5,4.5,250, 0, 25, -100,100,"");

        TProfile* p_v2_QB_coarse = new TProfile("p_v2_QB_coarse","p_v2_QB_coarse",4,-0.5,3.5,-100,100);
        TProfile* p_v2_QB_coarse_obs = new TProfile("p_v2_QB_coarse_obs","p_v2_QB_coarse_obs",4,-0.5,3.5,-100,100);
        TProfile* p_cos_QB_coarse = new TProfile("p_cos_QB_coarse","p_cos_QB_coarse",4,-0.5,3.5,-1,1);
        TProfile2D *pParity_Q2_QB1 = new TProfile2D("Parity_Q2_QB1","Parity_Q2_QB1",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile2D *pParity_Q2_QB1_obs = new TProfile2D("Parity_Q2_QB1_obs","Parity_Q2_QB1_obs",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile2D *pParity_Q2_QB2 = new TProfile2D("Parity_Q2_QB2","Parity_Q2_QB2",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile2D *pParity_Q2_QB2_obs = new TProfile2D("Parity_Q2_QB2_obs","Parity_Q2_QB2_obs",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile2D *pParity_Q2_QB3 = new TProfile2D("Parity_Q2_QB3","Parity_Q2_QB3",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile2D *pParity_Q2_QB3_obs = new TProfile2D("Parity_Q2_QB3_obs","Parity_Q2_QB3_obs",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile2D *pParity_Q2_QB4 = new TProfile2D("Parity_Q2_QB4","Parity_Q2_QB4",8,0.5,8.5,250, 0, 25, -100,100,"");
        TProfile2D *pParity_Q2_QB4_obs = new TProfile2D("Parity_Q2_QB4_obs","Parity_Q2_QB4_obs",8,0.5,8.5,250, 0, 25, -100,100,"");

        TH1D* Hist_positive = new TH1D("Hist_positive","Hist_positive",1000,0.5,1000.5);
        TH1D* Hist_negative = new TH1D("Hist_negative","Hist_negative",1000,0.5,1000.5);

        TProfile* p_a1_v2_p = new TProfile("p_a1_v2_p","p_a1_v2_p",200,-1,1,-1,1);
        TProfile* p_a1_v2_n = new TProfile("p_a1_v2_n","p_a1_v2_n",200,-1,1,-1,1);
        TProfile* p_a1_Q2_p = new TProfile("p_a1_Q2_p","p_a1_Q2_p",250, 0, 25,-1,1);
        TProfile* p_a1_Q2_n = new TProfile("p_a1_Q2_n","p_a1_Q2_n",250, 0, 25,-1,1);
        TProfile* p_a1_QQ_p = new TProfile("p_a1_QQ_p","p_a1_QQ_p",250, 0, 25,-1,1);
        TProfile* p_a1_QQ_n = new TProfile("p_a1_QQ_n","p_a1_QQ_n",250, 0, 25,-1,1);
        TProfile* p_a1_QB_p = new TProfile("p_a1_QB_p","p_a1_QB_p",250, 0, 25,-1,1);
        TProfile* p_a1_QB_n = new TProfile("p_a1_QB_n","p_a1_QB_n",250, 0, 25,-1,1);
        TProfile* p_v2_Q2_p = new TProfile("p_v2_Q2_p","p_v2_Q2_p",250, 0, 25,-100,100);
        TProfile* p_v2_Q2_n = new TProfile("p_v2_Q2_n","p_v2_Q2_n",250, 0, 25,-100,100);
        TProfile* p_v4_Q2_p = new TProfile("p_v4_Q2_p","p_v4_Q2_p",250, 0, 25,-100,100);
        TProfile* p_v4_Q2_n = new TProfile("p_v4_Q2_n","p_v4_Q2_n",250, 0, 25,-100,100);
        TProfile* p_v2_Q4_p = new TProfile("p_v2_Q4_p","p_v2_Q4_p",250, 0, 25,-100,100);
        TProfile* p_v2_Q4_n = new TProfile("p_v2_Q4_n","p_v2_Q4_n",250, 0, 25,-100,100);
        TProfile* p_v4_Q4_p = new TProfile("p_v4_Q4_p","p_v4_Q4_p",250, 0, 25,-100,100);
        TProfile* p_v4_Q4_n = new TProfile("p_v4_Q4_n","p_v4_Q4_n",250, 0, 25,-100,100);

        TProfile* p_a1_M_p = new TProfile("p_a1_M_p","p_a1_M_p",1000,-0.5,999.5,-1,1);
        TProfile* p_a1_M_n = new TProfile("p_a1_M_n","p_a1_M_n",1000,-0.5,999.5,-1,1);
        TProfile* p_v2_M_p = new TProfile("p_v2_M_p","p_v2_M_p",1000,-0.5,999.5,-1,1);
        TProfile* p_v2_M_n = new TProfile("p_v2_M_n","p_v2_M_n",1000,-0.5,999.5,-1,1);
        TProfile* p_M_v2_p = new TProfile("p_M_v2_p","p_M_v2_p",2000,-100,100,0,1000);
        TProfile* p_M_v2_n = new TProfile("p_M_v2_n","p_M_v2_n",2000,-100,100,0,1000);

        TProfile* p_a1_v2_pp = new TProfile("p_a1_v2_pp","p_a1_v2_pp",200,-1,1,-1,1);
        TProfile* p_a1_v2_pm = new TProfile("p_a1_v2_pm","p_a1_v2_pm",200,-1,1,-1,1);
        TProfile* p_a1_v2_kp = new TProfile("p_a1_v2_kp","p_a1_v2_kp",200,-1,1,-1,1);
        TProfile* p_a1_v2_km = new TProfile("p_a1_v2_km","p_a1_v2_km",200,-1,1,-1,1);

        TProfile* p_v2_v4_pp = new TProfile("p_v2_v4_pp","p_v2_v4_pp",200,-1,1,-1,1);
        TProfile* p_v2_v4_pm = new TProfile("p_v2_v4_pm","p_v2_v4_pm",200,-1,1,-1,1);
        TProfile* p_v2_v4_kp = new TProfile("p_v2_v4_kp","p_v2_v4_kp",200,-1,1,-1,1);
        TProfile* p_v2_v4_km = new TProfile("p_v2_v4_km","p_v2_v4_km",200,-1,1,-1,1);
        TProfile* p_v2_v4_pip = new TProfile("p_v2_v4_pip","p_v2_v4_pip",200,-1,1,-1,1);
        TProfile* p_v2_v4_pim = new TProfile("p_v2_v4_pim","p_v2_v4_pim",200,-1,1,-1,1);

        TProfile* p_vn_pi = new TProfile("p_vn_pi","p_vn_pi",6,0.5,6.5,-1,1);
        TProfile* p_vn_k = new TProfile("p_vn_k","p_vn_k",6,0.5,6.5,-1,1);
        TProfile* p_vn_p = new TProfile("p_vn_p","p_vn_p",6,0.5,6.5,-1,1);

        TH1D* Hist_Ach = new TH1D("Hist_Ach","Hist_Ach",200,-1,1);
        TProfile* p_Mult_Ach = new TProfile("p_Mult_Ach","p_Mult_Ach",200,-1,1, 0,1000);
        TProfile* p_v2_Ach = new TProfile("p_v2_Ach","p_v2_Ach",200,-1,1, -1,1);
        TProfile* p_v2_pos_Ach = new TProfile("p_v2_pos_Ach","p_v2_pos_Ach",200,-1,1, -100,100);
        TProfile* p_v2_neg_Ach = new TProfile("p_v2_neg_Ach","p_v2_neg_Ach",200,-1,1, -100,100);
        TProfile* p_a1_pos_Ach = new TProfile("p_a1_pos_Ach","p_a1_pos_Ach",200,-1,1, -100,100);
        TProfile* p_a1_neg_Ach = new TProfile("p_a1_neg_Ach","p_a1_neg_Ach",200,-1,1, -100,100);
        TProfile* p_gamma_OS_Ach = new TProfile("p_gamma_OS_Ach","p_gamma_OS_Ach",200,-1,1, -100,100);
        TProfile* p_gamma_SS_Ach = new TProfile("p_gamma_SS_Ach","p_gamma_SS_Ach",200,-1,1, -100,100);
        TProfile* p_delta_OS_Ach = new TProfile("p_delta_OS_Ach","p_delta_OS_Ach",200,-1,1, -100,100);
        TProfile* p_delta_SS_Ach = new TProfile("p_delta_SS_Ach","p_delta_SS_Ach",200,-1,1, -100,100);

        TProfile *p_Nf_pos = new TProfile("p_Nf_pos","p_Nf_pos",500,0.5,500.5,0,500);
        TProfile *p_Nb_pos = new TProfile("p_Nb_pos","p_Nb_pos",500,0.5,500.5,0,500);
        TProfile *p_Nf2_pos= new TProfile("p_Nf2_pos","p_Nf2_pos",500,0.5,500.5,0,250000);
        TProfile *p_Nb2_pos= new TProfile("p_Nb2_pos","p_Nb2_pos",500,0.5,500.5,0,250000);
        TProfile *p_Nf_neg = new TProfile("p_Nf_neg","p_Nf_neg",500,0.5,500.5,0,500);
        TProfile *p_Nb_neg = new TProfile("p_Nb_neg","p_Nb_neg",500,0.5,500.5,0,500);
        TProfile *p_Nf2_neg= new TProfile("p_Nf2_neg","p_Nf2_neg",500,0.5,500.5,0,250000);
        TProfile *p_Nb2_neg= new TProfile("p_Nb2_neg","p_Nb2_neg",500,0.5,500.5,0,250000);
        TProfile *p_Nfb_pos_pos= new TProfile("p_Nfb_pos_pos","p_Nfb_pos_pos",500,0.5,500.5,0,250000);
        TProfile *p_Nfb_neg_neg= new TProfile("p_Nfb_neg_neg","p_Nfb_neg_neg",500,0.5,500.5,0,250000);
        TProfile *p_Nfb_pos_neg= new TProfile("p_Nfb_pos_neg","p_Nfb_pos_neg",500,0.5,500.5,0,250000);
        TProfile *p_Nfb_neg_pos= new TProfile("p_Nfb_neg_pos","p_Nfb_neg_pos",500,0.5,500.5,0,250000);

        TProfile *p_pT_Q = new TProfile("p_pT_Q","p_pT_Q",250,0,25,0,10,"");
        TProfile *p_pT_Q2 = new TProfile("p_pT_Q2","p_pT_Q2",250,0,25,0,10,"");
        TProfile *p_pT_parent_Q = new TProfile("p_pT_parent_Q","p_pT_parent_Q",250,0,25,0,10,"");
        TProfile *p_pT_parent_Q2 = new TProfile("p_pT_parent_Q2","p_pT_parent_Q2",250,0,25,0,10,"");

        TH2D* Hist_ParentPhi_Q2_InvM0309 = new TH2D("Hist_ParentPhi_Q2_InvM0309","Hist_ParentPhi_Q2_InvM0309",100,0,10,60,-PI,PI);
        TH2D* Hist_ParentPhi_QQ_InvM0309 = new TH2D("Hist_ParentPhi_QQ_InvM0309","Hist_ParentPhi_QQ_InvM0309",100,0,10,60,-PI,PI);
        TH2D* Hist_ParentPhi_QB_InvM0309 = new TH2D("Hist_ParentPhi_QB_InvM0309","Hist_ParentPhi_QB_InvM0309",100,0,10,60,-PI,PI);
        TH2D* Hist_ParentPhi_Q2_InvM_all = new TH2D("Hist_ParentPhi_Q2_InvM_all","Hist_ParentPhi_Q2_InvM_all",100,0,10,60,-PI,PI);
        TH2D* Hist_ParentPhi_QQ_InvM_all = new TH2D("Hist_ParentPhi_QQ_InvM_all","Hist_ParentPhi_QQ_InvM_all",100,0,10,60,-PI,PI);
        TH2D* Hist_ParentPhi_QB_InvM_all = new TH2D("Hist_ParentPhi_QB_InvM_all","Hist_ParentPhi_QB_InvM_all",100,0,10,60,-PI,PI);
        TH2D* Hist_CosTheta_Q2_InvM_all = new TH2D("Hist_CosTheta_Q2_InvM_all","Hist_CosTheta_Q2_InvM_all",100,0,10,50,0,1);
        TH2D* Hist_CosTheta_QQ_InvM_all = new TH2D("Hist_CosTheta_QQ_InvM_all","Hist_CosTheta_QQ_InvM_all",100,0,10,50,0,1);
        TH2D* Hist_CosTheta_QB_InvM_all = new TH2D("Hist_CosTheta_QB_InvM_all","Hist_CosTheta_QB_InvM_all",100,0,10,50,0,1);
        TH2D* Hist_CosTheta_Q2_InvM_rho = new TH2D("Hist_CosTheta_Q2_InvM_rho","Hist_CosTheta_Q2_InvM_rho",100,0,10,50,0,1);
        TH2D* Hist_CosTheta_QQ_InvM_rho = new TH2D("Hist_CosTheta_QQ_InvM_rho","Hist_CosTheta_QQ_InvM_rho",100,0,10,50,0,1);
        TH2D* Hist_CosTheta_QB_InvM_rho = new TH2D("Hist_CosTheta_QB_InvM_rho","Hist_CosTheta_QB_InvM_rho",100,0,10,50,0,1);
        TH2D* Hist_CosTheta_Q2_InvM0306 = new TH2D("Hist_CosTheta_Q2_InvM0306","Hist_CosTheta_Q2_InvM0306",100,0,10,50,0,1);
        TH2D* Hist_CosTheta_QQ_InvM0306 = new TH2D("Hist_CosTheta_QQ_InvM0306","Hist_CosTheta_QQ_InvM0306",100,0,10,50,0,1);
        TH2D* Hist_CosTheta_QB_InvM0306 = new TH2D("Hist_CosTheta_QB_InvM0306","Hist_CosTheta_QB_InvM0306",100,0,10,50,0,1);
        TH2D* Hist_CosTheta_Q2_InvM0914 = new TH2D("Hist_CosTheta_Q2_InvM0914","Hist_CosTheta_Q2_InvM0914",100,0,10,50,0,1);
        TH2D* Hist_CosTheta_QQ_InvM0914 = new TH2D("Hist_CosTheta_QQ_InvM0914","Hist_CosTheta_QQ_InvM0914",100,0,10,50,0,1);
        TH2D* Hist_CosTheta_QB_InvM0914 = new TH2D("Hist_CosTheta_QB_InvM0914","Hist_CosTheta_QB_InvM0914",100,0,10,50,0,1);

        Int_t nentries = chain->GetEntries();
        //loop through events
        for(int i = 0; i < nentries; i++){

                if((i+1)%1000==0) cout << "Processing entry == "<< i+1 <<" == out of "<<nentries<<".\n";
                chain->GetEntry(i);

                //		TLeaf* leaf_psi = chain->GetLeaf("Psi");
                TLeaf* leaf_b   = chain->GetLeaf("bim");
                TLeaf* leaf_NoTracks = chain->GetLeaf("np");
                TLeaf* leaf_Np_p= chain->GetLeaf("npartproj");
                TLeaf* leaf_Np_t= chain->GetLeaf("nparttarg");
                int Np = leaf_Np_p->GetValue(0) + leaf_Np_t->GetValue(0);
                //                if(Np<3) continue;

                psi = 0; //leaf_psi->GetValue(0);
                b = leaf_b->GetValue(0);
                NPTracks= (int)leaf_NoTracks->GetValue(0);
                //std::cout << NPTracks << endl;
                Centrality = 0;
                for(int j=0;j<9;j++) if(b<cenDef[j]) Centrality = j+1;
                hCentrality->Fill(Centrality,Eweight);
                //		if(cen && Centrality != cen) continue;   	

                TLeaf* leaf_PxV0       = chain->GetLeaf("px");
                TLeaf* leaf_PyV0       = chain->GetLeaf("py");
                TLeaf* leaf_PzV0       = chain->GetLeaf("pz");
                TLeaf* leaf_PIDV0    = chain->GetLeaf("id");
                TLeaf* leaf_e           = chain->GetLeaf("e");

                //TLeaf* leaf_MassV0	= chain->GetLeaf("Mass");

                //TPC EP reconstruction
                TVector2 mQ1, mQ2; 
                Double_t mQx=0., mQy=0., mQx1=0., mQy1=0., mQx2=0., mQy2=0., mQx_parent=0., mQy_parent=0.;
                double mQBx=0., mQBy=0.;
                int Fcount = 0, Ecount = 0, Wcount =0, Pcount=0, Parent_count = 0, Bcount = 0;

                int mult_all = 0;
                int mult_ptCut = 0;
                int mult_ptEtaCut = 0;

                int Npos = 0, Nneg = 0;
                for(int trk = 0; trk < NPTracks; trk++) {

                        float PxAsso    = leaf_PxV0->GetValue(trk);
                        float PyAsso    = leaf_PyV0->GetValue(trk); if((i+1)%2==0) PyAsso *= -1;
                        float PzAsso    = leaf_PzV0->GetValue(trk);
                        float PIDAsso   = leaf_PIDV0->GetValue(trk);
                        float MassAsso    = leaf_e->GetValue(trk);
                        //TLorentzVector part_vec;
                        //part_vec.SetPxPyPzE(PxAsso, PyAsso, PzAsso, energy);

                        // float MassAsso = part_vec.M();
                        //float MassAsso  = leaf_MassV0->GetValue(trk);
                        float ChargeAsso= 0;
                        if(PIDAsso== ID_pip || PIDAsso== ID_Kp || PIDAsso== ID_pp) ChargeAsso = 1;
                        if(PIDAsso== ID_pim || PIDAsso== ID_Km || PIDAsso== ID_pb) ChargeAsso =-1;
                        float PtAsso = sqrt(PxAsso*PxAsso+PyAsso*PyAsso);
                        float PhiAsso = atan2(PyAsso,PxAsso);
                        float ThetaAsso = atan2(PtAsso,PzAsso);
                        float EtaAsso = -log(tan(ThetaAsso/2.));
                        


                        if(ChargeAsso!=1 && ChargeAsso!=-1) continue;
                        if(PtAsso > pt_asso_up || PtAsso < pt_asso_lo) continue;
                        if(EtaAsso>EtaAssoCut || EtaAsso<-EtaAssoCut) continue;
                        if(MassAsso<0.1 || MassAsso>1) continue;

                        if(ChargeAsso>0) Npos++;
                        if(ChargeAsso<0) Nneg++;
                        Fcount++;

                        if(EtaAsso > EtaCut || EtaAsso < -EtaCut) continue;
                        Pcount++;

                }
                Hist_positive->Fill(Npos);
                Hist_negative->Fill(Nneg);

                //reshuffle charge
                int iCharge[Pcount];
                Pcount = 0;
                int match[NPTracks];

                int iTrack[Fcount], Scount = Fcount/2 -1;
                for(int q=0;q<Fcount;q++) iTrack[q] = q;
                random_shuffle(iTrack,iTrack+Fcount);
                Fcount = 0;
                float Psi_random = gRandom->Rndm()*2*PI;
                int Nf_pos = 0, Nf_neg = 0, Nb_pos = 0, Nb_neg = 0;
                float mQx4 = 0, mQy4 = 0;

                for(int trk = 0; trk < NPTracks; trk++) {

                        float PxAsso    = leaf_PxV0->GetValue(trk);
                        float PyAsso    = leaf_PyV0->GetValue(trk); if((i+1)%2==0) PyAsso *= -1;
                        float PzAsso    = leaf_PzV0->GetValue(trk);
                        float PIDAsso   = leaf_PIDV0->GetValue(trk);
                        float MassAsso    = leaf_e->GetValue(trk);
                        //TLorentzVector part_vec;
                        //part_vec.SetPxPyPzE(PxAsso, PyAsso, PzAsso, energy);
                        //std::cout << "First Particle ID: " << PIDAsso << std::endl;
                        //float MassAsso = part_vec.M();
                        //float MassAsso  = leaf_MassV0->GetValue(trk);
                        //float ChargeAsso = 0;
                        match[trk] = 0;
                        //if(PIDAsso== ID_pip || PIDAsso== ID_Kp || PIDAsso== ID_pp) ChargeAsso = 1;
                        //if(PIDAsso== ID_pim || PIDAsso== ID_Km || PIDAsso== ID_pb) ChargeAsso =-1;
                        float PtAsso = sqrt(PxAsso*PxAsso+PyAsso*PyAsso);
                        float PhiAsso = atan2(PyAsso,PxAsso);
                        float ThetaAsso = atan2(PtAsso,PzAsso);
                        float EtaAsso = -log(tan(ThetaAsso/2.));

                        int pid = (int)PIDAsso;

                        // Check if the particle exists in our table and if it has a non-zero charge
                        if (pCharge.find(pid) == pCharge.end() || pCharge[pid] == 0.0) {
                                continue; // Skip neutral particles or unknown IDs
                        }

                        // Dynamically assign the charge
                        float ChargeAsso = pCharge[pid]; 

                        // Fill the pie chart using the particle's string name!
                        hSpeciesAll->Fill(pName[pid].c_str(), 1.0);



                        if(false && ChargeAsso!=1 && ChargeAsso!=-1) {
                                continue;
                        }
                        mult_all++;
                        if(PtAsso > 0.05){
                                mult_ptCut++;
                                if(EtaAsso < 0.5 && EtaAsso > -0.5){
                                        mult_ptEtaCut++;
                                }
                        }
                        
                        if(PtAsso > pt_asso_up || PtAsso < pt_asso_lo) continue;
                        if(EtaAsso>EtaAssoCut || EtaAsso<-EtaAssoCut) continue;
                        if(MassAsso<0.1 || MassAsso>1) continue;

                        if(EtaAsso> 1.5) {mQx1 +=PtAsso*cos(PhiAsso*nHar); mQy1 +=PtAsso*sin(PhiAsso* nHar); Ecount++;}
                        if(EtaAsso<-1.5) {mQx2 +=PtAsso*cos(PhiAsso*nHar); mQy2 +=PtAsso*sin(PhiAsso*nHar); Wcount++;}
                        //		if(EtaAsso> 2) {mQx1 +=PtAsso*cos(PhiAsso*nHar); mQy1 +=PtAsso*sin(PhiAsso* nHar); Ecount++;}
                        //		if(EtaAsso<-2) {mQx2 +=PtAsso*cos(PhiAsso*nHar); mQy2 +=PtAsso*sin(PhiAsso*nHar); Wcount++;}
                        Fcount++;

                        //		if((EtaAsso> -1.5 && EtaAsso<-1) || (EtaAsso> 1 && EtaAsso<1.5)) {
                        if((EtaAsso> -2 && EtaAsso<-1) || (EtaAsso> 1 && EtaAsso<2)) {
                                mQBx += cos(PhiAsso*nHar); mQBy += sin(PhiAsso*nHar);
                                Bcount++;
                        }

                        if(EtaAsso > EtaCut || EtaAsso < -EtaCut) continue;
                        if(PIDAsso!= ID_pip && PIDAsso!= ID_pim) continue;
                        //		if(PIDAsso!= ID_Kp && PIDAsso!= ID_Km) continue;
                        //		mQx += PtAsso*cos(PhiAsso*nHar); mQy += PtAsso*sin(PhiAsso*nHar);
                        mQx += cos(PhiAsso*nHar); mQy += sin(PhiAsso*nHar);
                        mQx4 += cos(PhiAsso*4); mQy4 += sin(PhiAsso*4);
                        match[trk] = Pcount;
                        iCharge[Pcount] = ChargeAsso;
                        Pcount++;
                        float phi = PhiAsso + Psi_random;
                        if(phi > 2*PI) phi -= 2*PI;
                        if(phi < 0) phi += 2*PI;
                        if(phi>2*PI-PI/3. || phi<PI/3.) {if(ChargeAsso==1) Nf_pos++; if(ChargeAsso==-1) Nf_neg++;}
                        else if(phi<PI+PI/3. && phi>PI-PI/3.) {if(ChargeAsso==1) Nb_pos++; if(ChargeAsso==-1) Nb_neg++;}
                }
                hMult_All->Fill(mult_all);
                hMult_ptCut->Fill(mult_ptCut);
                hMult_ptEtaCut->Fill(mult_ptEtaCut);
                random_shuffle(iCharge,iCharge+Pcount);
                if(mQx1==0 || mQy1==0 || mQx2==0 || mQy2==0) continue;
                hMult->Fill(Pcount);

                //if(Pcount<100) continue; // ASK ABOUT THIS
                //	  if(Pcount<10) continue;
                p_Nf_pos->Fill(Pcount, Nf_pos);
                p_Nf_neg->Fill(Pcount, Nf_neg);
                p_Nb_pos->Fill(Pcount, Nb_pos);
                p_Nb_neg->Fill(Pcount, Nb_neg);
                p_Nf2_pos->Fill(Pcount, Nf_pos*Nf_pos);
                p_Nf2_neg->Fill(Pcount, Nf_neg*Nf_neg);
                p_Nb2_pos->Fill(Pcount, Nb_pos*Nb_pos);
                p_Nb2_neg->Fill(Pcount, Nb_neg*Nb_neg);
                p_Nfb_pos_pos->Fill(Pcount, Nf_pos*Nb_pos);
                p_Nfb_neg_neg->Fill(Pcount, Nf_neg*Nb_neg);
                p_Nfb_pos_neg->Fill(Pcount, Nf_pos*Nb_neg);
                p_Nfb_neg_pos->Fill(Pcount, Nf_neg*Nb_pos);

                //	  double Q2 = (mQx*mQx+mQy*mQy)/float(211.293+0.064*0.064*45719.4);
                double Q2 = (mQx*mQx+mQy*mQy)/float(Pcount+0.064*0.064*Pcount*Pcount);
                double Q4 = (mQx4*mQx4+mQy4*mQy4)/float(Pcount);
                //	  double Q2 = 2*(mQx*mQx+mQy*mQy)/float(Pcount+0.064*0.064*Pcount*Pcount);
                //	  double Q2 = (mQx*mQx+mQy*mQy)/float(Pcount+0.0535*0.0535*Pcount*Pcount);

                //	  double Q2 = (mQx*mQx+mQy*mQy)/float(Pcount+0.0505607*0.0505607*Pcount*Pcount);
                //	  double Q2 = (mQx*mQx+mQy*mQy)/float(Pcount+0.021212*0.021212*Pcount*Pcount);

                double QB2= (mQBx*mQBx+mQBy*mQBy)/float(Bcount+0.064*0.064*Bcount*Bcount);
                double QB = sqrt(QB2);
                //QB = QB2;
                if(QB<0.56) QB_index = 0;
                else if(QB<0.88) QB_index = 1;
                else if(QB<1.24) QB_index = 2;
                else if(QB<2)	QB_index = 3;
                else QB_index = 4;

                Hist_QB->Fill(QB);
                p_RefMult_QB->Fill(QB,Pcount);
                Hist_Q2_QB2->Fill(QB2,Q2);
                p_Q2_QB2->Fill(QB2,Q2);
                Hist_Q_QB->Fill(QB,sqrt(Q2));
                p_Q_QB->Fill(QB,sqrt(Q2));
                float Psi_B = atan2(mQBy,mQBx)/2.;
                float Psi_A = atan2(mQy,mQx)/2.;
                //	  Hist_PsiQ_PsiQB->Fill(Psi_B, Psi_A);
                pMult_PsiQ_PsiQB->Fill(Psi_B, Psi_A,Pcount);

                //	  if((fabs(Psi_A-Psi_B)>PI/4.)&&(fabs(Psi_A-Psi_B)<3*PI/4.)) continue;		
                //	  if(!((fabs(Psi_A-Psi_B)>PI/4.)&&(fabs(Psi_A-Psi_B)<3*PI/4.))) continue;
                //	  if(!((fabs(Psi_A-Psi_B)>PI/3.)&&(fabs(Psi_A-Psi_B)<2*PI/3.))) continue;
                Hist_PsiQ_PsiQB->Fill(Psi_B, Psi_A);

                mQ1.Set(mQx1, mQy1); mQ2.Set(mQx2, mQy2);
                float TPC_EP_east = mQ1.Phi()/nHar;
                float TPC_EP_west = mQ2.Phi()/nHar;
                float TPC_EP_east_new = TPC_EP_east;
                float TPC_EP_west_new = TPC_EP_west;
                float TPC_EP_full = gRandom->Gaus(0,0.3);
                Hist_cos->Fill(1,cos(nHar*TPC_EP_east_new-nHar*psi), Eweight);
                Hist_cos->Fill(2,cos(nHar*TPC_EP_east_new-nHar*TPC_EP_west_new), Eweight);
                p_cos_Q2->Fill(Q2,cos(nHar*TPC_EP_east_new-nHar*TPC_EP_west_new), Eweight);
                p_RefMult_Q2->Fill(Q2,Pcount);
                p_cos_QB->Fill(QB,cos(nHar*TPC_EP_east_new-nHar*TPC_EP_west_new));
                p_cos_QB_coarse->Fill(QB_index,cos(nHar*TPC_EP_east_new-nHar*TPC_EP_west_new));

                float Ach = float(Npos - Nneg)/(Npos + Nneg);
                //	  Ach = float(Npos - Nneg)/200.;
                Ach = float(Npos - Nneg)/sqrt(float(Npos + Nneg))/10.;
                Ach = float(Npos - Nneg)/(sqrt(float(Npos))+sqrt(float(Nneg)))/10.;
                Hist_Ach->Fill(Ach);

                //calculate parent Q first
                ///////////////////////////////////////////////////////////////////
                for(int trki = 0; trki < NPTracks; trki++){
                        Px    = leaf_PxV0->GetValue(trki);
                        Py    = leaf_PyV0->GetValue(trki); if((i+1)%2==0) Py *= -1;
                        Pz    = leaf_PzV0->GetValue(trki);
                        PID   = leaf_PIDV0->GetValue(trki);
                        //Mass  = leaf_MassV0->GetValue(trki);
                        float Mass    = leaf_e->GetValue(trki);
                        //TLorentzVector part_vec;
                        //part_vec.SetPxPyPzE(Px, Py, Pz, energy);

                        //float Mass = part_vec.M();
                        Charge= 0;
                        if(PID== ID_pip) Charge = 1;
                        if(PID== ID_pim) Charge =-1;
                        Pt = sqrt(Px*Px+Py*Py);
                        Phi = atan2(Py,Px);
                        Theta = atan2(Pt,Pz);
                        Eta = -log(tan(Theta/2.));

                        if(Mass<0.1 || Mass>1) continue;
                        if(Eta > EtaCut || Eta < -EtaCut) continue;
                        if(Charge!=1 && Charge!=-1) continue;
                        if(Pt < pt_trig_lo || Pt > pt_trig_up) continue;

                        for(int trkj = trki+1; trkj < NPTracks; trkj++) {
                                float Px2    = leaf_PxV0->GetValue(trkj);
                                float Py2    = leaf_PyV0->GetValue(trkj); if((i+1)%2==0) Py2 *= -1;
                                float Pz2    = leaf_PzV0->GetValue(trkj);
                                float PID2   = leaf_PIDV0->GetValue(trkj);
                                float Mass2    = leaf_e->GetValue(trkj);
                                //TLorentzVector part_vec2;
                                //part_vec2.SetPxPyPzE(Px2, Py2, Pz2, energy2);

                                //float Mass2 = part_vec.M();
                                float Charge2= 0;
                                if(PID2== ID_pip) Charge2 = 1;
                                if(PID2== ID_pim) Charge2 =-1;
                                float Pt2 = sqrt(Px2*Px2+Py2*Py2);
                                float Phi2 = atan2(Py2,Px2);
                                float Theta2 = atan2(Pt2,Pz2);
                                float Eta2 = -log(tan(Theta2/2.));

                                if(Charge2!=1 && Charge2!=-1) continue;
                                if(Mass2<0.1 || Mass2>1) continue;
                                if(Pt2 < pt_trig_lo || Pt2 > pt_trig_up) continue;
                                if(Eta2 > EtaCut || Eta2 < -EtaCut) continue;

                                // Boost back to the "parent's" rest frame
                                TLorentzVector lA; lA.SetPxPyPzE(Px,Py,Pz,sqrt(Mass*Mass+Pt*Pt+Pz*Pz));
                                TLorentzVector lB; lB.SetPxPyPzE(Px2,Py2,Pz2,sqrt(Mass2*Mass2+Pt2*Pt2+Pz2*Pz2));

                                float phi_parent = (lA+lB).Phi();
                                mQx_parent += cos(phi_parent*nHar);
                                mQy_parent += sin(phi_parent*nHar);
                                Parent_count++;
                        }
                }
                double QQ = (mQx_parent*mQx_parent + mQy_parent*mQy_parent)/float(Parent_count+0.05678*0.05678*Parent_count*Parent_count);
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                //loop through matched primary tracks
                for(int trki = 0; trki < NPTracks; trki++){
                        Px    = leaf_PxV0->GetValue(trki);
                        Py    = leaf_PyV0->GetValue(trki); if((i+1)%2==0) Py *= -1;
                        Pz    = leaf_PzV0->GetValue(trki);
                        PID   = leaf_PIDV0->GetValue(trki);
                        Mass    = leaf_e->GetValue(trki);
                        //TLorentzVector part_vec;
                        //part_vec.SetPxPyPzE(Px, Py, Pz, energy);

                        //Mass = part_vec.M();
                        //Mass  = leaf_MassV0->GetValue(trki);
                        Charge= 0;
                        if(PID== ID_pip) Charge = 1;
                        if(PID== ID_pim) Charge =-1;
			if(PID== ID_Kp) Charge = 1;
			if(PID== ID_Km) Charge =-1;
                        if(PID== ID_pp) Charge = 1;
			if(PID== ID_pb) Charge =-1;
                        Pt = sqrt(Px*Px+Py*Py);
                        Phi = atan2(Py,Px);
                        Theta = atan2(Pt,Pz);
                        Eta = -log(tan(Theta/2.));
                        yy = 0.5*log((sqrt(Mass*Mass+Px*Px+Py*Py+Pz*Pz)+Pz)/(sqrt(Mass*Mass+Px*Px+Py*Py+Pz*Pz)-Pz));
                        float eff = 1;

                        //if(Mass<0.1 || Mass>1) continue;
                        
                        //if(Eta > EtaCut || Eta < -EtaCut) continue;
                        

                        float Phi_new = Phi;
                        float Phi_new2 = Phi;
                        if((i+1)%2==0) Phi_new2 = atan2(-Py,Px);
                        if(PID== ID_pip) {
                                Hist_v2_pt_pip->Fill(Pt,cos(2*Phi_new-2*psi));
                                Hist_a1_pt_pip->Fill(Pt,sin(Phi_new2-psi));
                                Hist_a3_pt_pip->Fill(Pt,sin(3*Phi_new2-3*psi));
                                if(Pt>0.2 && Pt<2) {
                                        p_a1_v2_p->Fill(cos(2*Phi_new-2*psi),sin(Phi_new2-psi));
                                        p_v2_v4_pip->Fill(cos(4*Phi_new-4*psi),cos(2*Phi_new-2*psi));
                                        pTemp_pion->Fill(1,cos(2*Phi_new-2*psi));
                                        pTemp_pion->Fill(3,sin(Phi_new2-psi));
                                        pTemp_pion->Fill(5,cos(4*Phi_new-4*psi));
                                        pTemp_pion->Fill(7,cos(6*Phi_new-6*psi));

                                        p_a1_Q2_p->Fill(Q2,sin(Phi_new2-psi));
                                        p_a1_QQ_p->Fill(QQ,sin(Phi_new2-psi));
                                        p_a1_QB_p->Fill(QB,sin(Phi_new2-psi));
                                        p_v2_Q2_p->Fill(Q2,cos(2*Phi_new-2*psi));
                                        p_v4_Q2_p->Fill(Q2,cos(4*Phi_new-4*psi));
                                        p_v2_Q4_p->Fill(Q4,cos(2*Phi_new-2*psi));
                                        p_v4_Q4_p->Fill(Q4,cos(4*Phi_new-4*psi));

                                        Hist_v1_y_pip->Fill(yy,cos(Phi_new-psi));
                                        Hist_v2_y_pip->Fill(yy,cos(2*Phi_new-2*psi));
                                        Hist_a1_y_pip->Fill(yy,sin(Phi_new2-psi));
                                        Hist_a2_y_pip->Fill(yy,sin(2*Phi_new-2*psi));
                                        Hist_a3_y_pip->Fill(yy,sin(3*Phi_new2-3*psi));
                                }
                        }
                        if(PID== ID_pim) {
                                Hist_v2_pt_pim->Fill(Pt,cos(2*Phi_new-2*psi));
                                Hist_a1_pt_pim->Fill(Pt,sin(Phi_new2-psi));
                                Hist_a3_pt_pim->Fill(Pt,sin(3*Phi_new2-3*psi));
                                if(Pt>0.2 && Pt<2) {
                                        p_a1_v2_n->Fill(cos(2*Phi_new-2*psi),sin(Phi_new2-psi));
                                        p_v2_v4_pim->Fill(cos(4*Phi_new-4*psi),cos(2*Phi_new-2*psi));
                                        pTemp_pion->Fill(2,cos(2*Phi_new-2*psi));
                                        pTemp_pion->Fill(4,sin(Phi_new2-psi));
                                        pTemp_pion->Fill(6,cos(4*Phi_new-4*psi));
                                        pTemp_pion->Fill(8,cos(6*Phi_new-6*psi));

                                        p_a1_Q2_n->Fill(Q2,sin(Phi_new2-psi));
                                        p_a1_QQ_n->Fill(QQ,sin(Phi_new2-psi));
                                        p_a1_QB_n->Fill(QB,sin(Phi_new2-psi));
                                        p_v2_Q2_n->Fill(Q2,cos(2*Phi_new-2*psi));
                                        p_v4_Q2_n->Fill(Q2,cos(4*Phi_new-4*psi));
                                        p_v2_Q4_n->Fill(Q4,cos(2*Phi_new-2*psi));
                                        p_v4_Q4_n->Fill(Q4,cos(4*Phi_new-4*psi));

                                        Hist_v1_y_pim->Fill(yy,cos(Phi_new-psi));
                                        Hist_v2_y_pim->Fill(yy,cos(2*Phi_new-2*psi));
                                        Hist_a1_y_pim->Fill(yy,sin(Phi_new2-psi));
                                        Hist_a2_y_pim->Fill(yy,sin(2*Phi_new-2*psi));
                                        Hist_a3_y_pim->Fill(yy,sin(3*Phi_new2-3*psi));
                                }
                        }
                        if(PID== ID_Kp) {
                                Hist_v2_pt_Kp->Fill(Pt,cos(2*Phi_new-2*psi));
                                Hist_a1_pt_Kp->Fill(Pt,sin(Phi_new2-psi));
                                Hist_a3_pt_Kp->Fill(Pt,sin(3*Phi_new2-3*psi));
                                if(Pt>0.2 && Pt<2) {
                                        p_a1_v2_kp->Fill(cos(2*Phi_new-2*psi),sin(Phi_new2-psi));
                                        p_v2_v4_kp->Fill(cos(4*Phi_new-4*psi),cos(2*Phi_new-2*psi));
                                        pTemp_kaon->Fill(1,cos(2*Phi_new-2*psi));
                                        pTemp_kaon->Fill(3,sin(Phi_new2-psi));
                                        pTemp_kaon->Fill(5,cos(4*Phi_new-4*psi));
                                        pTemp_kaon->Fill(7,cos(6*Phi_new-6*psi));
                                        Hist_v1_y_Kp->Fill(yy,cos(Phi_new-psi));
                                        Hist_v2_y_Kp->Fill(yy,cos(2*Phi_new-2*psi));
                                        Hist_a1_y_Kp->Fill(yy,sin(Phi_new2-psi));
                                        Hist_a2_y_Kp->Fill(yy,sin(2*Phi_new-2*psi));
                                        Hist_a3_y_Kp->Fill(yy,sin(3*Phi_new2-3*psi));
                                }
                        }
                        if(PID== ID_Km) {
                                Hist_v2_pt_Km->Fill(Pt,cos(2*Phi_new-2*psi));
                                Hist_a1_pt_Km->Fill(Pt,sin(Phi_new2-psi));
                                Hist_a3_pt_Km->Fill(Pt,sin(3*Phi_new2-3*psi));
                                if(Pt>0.2 && Pt<2) {
                                        p_a1_v2_km->Fill(cos(2*Phi_new-2*psi),sin(Phi_new2-psi));
                                        p_v2_v4_km->Fill(cos(4*Phi_new-4*psi),cos(2*Phi_new-2*psi));
                                        pTemp_kaon->Fill(2,cos(2*Phi_new-2*psi));
                                        pTemp_kaon->Fill(4,sin(Phi_new2-psi));
                                        pTemp_kaon->Fill(6,cos(4*Phi_new-4*psi));
                                        pTemp_kaon->Fill(8,cos(6*Phi_new-6*psi));
                                        Hist_v1_y_Km->Fill(yy,cos(Phi_new-psi));
                                        Hist_v2_y_Km->Fill(yy,cos(2*Phi_new-2*psi));
                                        Hist_a1_y_Km->Fill(yy,sin(Phi_new2-psi));
                                        Hist_a2_y_Km->Fill(yy,sin(2*Phi_new-2*psi));
                                        Hist_a3_y_Km->Fill(yy,sin(3*Phi_new2-3*psi));
                                }
                        }
                        if(PID== ID_pp) {
                                Hist_v2_pt_pp->Fill(Pt,cos(2*Phi_new-2*psi));
                                Hist_a1_pt_pp->Fill(Pt,sin(Phi_new2-psi));
                                Hist_a3_pt_pp->Fill(Pt,sin(3*Phi_new2-3*psi));
                                if(Pt>0.2 && Pt<2) {
                                        p_a1_v2_pp->Fill(cos(2*Phi_new-2*psi),sin(Phi_new2-psi));
                                        p_v2_v4_pp->Fill(cos(4*Phi_new-4*psi),cos(2*Phi_new-2*psi));
                                        pTemp_proton->Fill(1,cos(2*Phi_new-2*psi));
                                        pTemp_proton->Fill(3,sin(Phi_new2-psi));
                                        pTemp_proton->Fill(5,cos(4*Phi_new-4*psi));
                                        pTemp_proton->Fill(7,cos(6*Phi_new-6*psi));
                                        Hist_v1_y_pp->Fill(yy,cos(Phi_new-psi));
                                        Hist_v2_y_pp->Fill(yy,cos(2*Phi_new-2*psi));
                                        Hist_a1_y_pp->Fill(yy,sin(Phi_new2-psi));
                                        Hist_a2_y_pp->Fill(yy,sin(2*Phi_new-2*psi));
                                        Hist_a3_y_pp->Fill(yy,sin(3*Phi_new2-3*psi));
                                }
                        }
                        if(PID== ID_pb) {
                                Hist_v2_pt_pm->Fill(Pt,cos(2*Phi_new-2*psi));
                                Hist_a1_pt_pm->Fill(Pt,sin(Phi_new2-psi));
                                Hist_a3_pt_pm->Fill(Pt,sin(3*Phi_new2-3*psi));
                                if(Pt>0.2 && Pt<2) {
                                        p_a1_v2_pm->Fill(cos(2*Phi_new-2*psi),sin(Phi_new2-psi));
                                        p_v2_v4_pm->Fill(cos(4*Phi_new-4*psi),cos(2*Phi_new-2*psi));
                                        pTemp_proton->Fill(2,cos(2*Phi_new-2*psi));
                                        pTemp_proton->Fill(4,sin(Phi_new2-psi));
                                        pTemp_proton->Fill(6,cos(4*Phi_new-4*psi));
                                        pTemp_proton->Fill(8,cos(6*Phi_new-6*psi));
                                        Hist_v1_y_pm->Fill(yy,cos(Phi_new-psi));
                                        Hist_v2_y_pm->Fill(yy,cos(2*Phi_new-2*psi));
                                        Hist_a1_y_pm->Fill(yy,sin(Phi_new2-psi));
                                        Hist_a2_y_pm->Fill(yy,sin(2*Phi_new-2*psi));
                                        Hist_a3_y_pm->Fill(yy,sin(3*Phi_new2-3*psi));
                                }
                        }

                        if(Charge!=1 && Charge!=-1) continue;
                        if(Pt > 0.05){
                                hEtaPtDist->Fill(Eta,Pt,Eweight);
                                Hist_Pt->Fill(Pt,Eweight);
                                Hist_Y->Fill(yy);
                                Hist_Phi->Fill(Phi);
                        }
                        if(Pt < pt_trig_lo || Pt > pt_trig_up) continue;

                        float v2a = cos(nHar*Phi_new-nHar*psi)*100;
                        float v2b1 = cos(nHar*Phi_new - nHar*TPC_EP_east_new)*100;
                        float v2b2 = cos(nHar*Phi_new - nHar*TPC_EP_west_new)*100;
                        Hist_v2_pt->Fill(Pt,v2a);
                        Hist_v2_pt_obs->Fill(Pt,v2b1);Hist_v2_pt_obs->Fill(Pt,v2b2);
                        Hist_v1_eta->Fill(Eta,cos(Phi_new-psi)*100);
                        Hist_v2_eta->Fill(Eta,v2a);
                        Hist_v3_eta->Fill(Eta,cos(3*Phi_new-3*psi)*100);
                        Hist_v2_eta_obs->Fill(Eta,v2b1);Hist_v2_eta_obs->Fill(Eta,v2b2);
                        pTemp_pT->Fill(1,Pt);
                        pTemp_v2->Fill(1,v2a);
                        pTemp_v2->Fill(2,v2b1);pTemp_v2->Fill(2,v2b2);

                        p_v2_Q2->Fill(Q2,v2a);
                        p_v2_s_Q->Fill(QQ,v2a);
                        p_v2_QB->Fill(QB,v2a);
                        p_v2_Q2_obs->Fill(Q2,v2b1); p_v2_Q2_obs->Fill(Q2,v2b2);
                        p_v2_s_Q_obs->Fill(QQ,v2b1); p_v2_s_Q_obs->Fill(QQ,v2b2);
                        p_v2_QB_obs->Fill(QB,v2b1); p_v2_QB_obs->Fill(QB,v2b2);

                        if((i+1)%2==0) Phi = atan2(-Py,Px);
                        float a1 = sin(Phi-psi);
                        float a3 = sin(3*Phi-3*psi);
                        if(Charge>0) {
                                Hist_cos->Fill(3,a1);
                                Hist_cos->Fill(5,a3);  
                                pTemp_a1->Fill(1,a1); 
                                pTemp_v2->Fill(3,v2a);
                        }
                        else {	
                                Hist_cos->Fill(4,a1); 
                                Hist_cos->Fill(6,a3);
                                pTemp_a1->Fill(2,a1); 
                                pTemp_v2->Fill(4,v2a);
                        }

                        for(int trkj = trki+1; trkj < NPTracks; trkj++) {
                                float Px2    = leaf_PxV0->GetValue(trkj);
                                float Py2    = leaf_PyV0->GetValue(trkj); if((i+1)%2==0) Py2 *= -1;
                                float Pz2    = leaf_PzV0->GetValue(trkj);
                                float PID2   = leaf_PIDV0->GetValue(trkj);
                                float Mass2    = leaf_e->GetValue(trkj);
                                //TLorentzVector part_vec;
                                //part_vec.SetPxPyPzE(Px2, Py2, Pz2, energy2);

                                //float Mass2 = part_vec.M();
                                //float Mass2  = leaf_MassV0->GetValue(trkj);
                                float Charge2= 0;
                                if(PID2== ID_pip) Charge2 = 1;
                                if(PID2== ID_pim) Charge2 =-1;
                                //				if(PID2== ID_Kp) Charge2 = 1;
                                //				if(PID2== ID_Km) Charge2 =-1;
                                float Pt2 = sqrt(Px2*Px2+Py2*Py2);
                                float Phi2 = atan2(Py2,Px2);
                                float Theta2 = atan2(Pt2,Pz2);
                                float Eta2 = -log(tan(Theta2/2.));

                                if(Charge2!=1 && Charge2!=-1) continue;
                                if(Mass2<0.1 || Mass2>1) continue;
                                if(Pt2 < pt_trig_lo || Pt2 > pt_trig_up) continue;
                                if(Eta2 > EtaCut || Eta2 < -EtaCut) continue;
                                float eff2 = 1;

                                //if(fabs(Eta-Eta2)<0.1) continue;

                                hDpt->Fill(fabs(Pt-Pt2),Eweight);

                                float Phi2_new = Phi2;
                                // Boost back to the "parent's" rest frame
                                TLorentzVector lA; lA.SetPxPyPzE(Px,Py,Pz,sqrt(Mass*Mass+Pt*Pt+Pz*Pz));
                                TLorentzVector lB; lB.SetPxPyPzE(Px2,Py2,Pz2,sqrt(Mass2*Mass2+Pt2*Pt2+Pz2*Pz2));
                                TVector3 boost = -(lA +lB).BoostVector();
                                if(opt_boost) {lA.Boost(boost); lB.Boost(boost); Phi_new = lA.Phi(); Phi2_new = lB.Phi();}
                                //if(gRandom->Rndm()>0.1) continue;
                                float phi_parent = (lA+lB).Phi();
                                float eta_parent = (lA+lB).Eta();
                                float pt_parent  = (lA+lB).Pt();
                                float mass_parent= (lA+lB).M();
                                //if(mass_parent>0.755 && mass_parent<0.775) continue;
                                //if(pt_parent<0.1) continue;

                                //bbost one pion to the parent's rest frame
                                lA.Boost(boost); 
                                TVector3 LL(0, -1, 0);
                                TVector3 pi_rest(lA.Px(),lA.Py(),lA.Pz());
                                double CosTheta = fabs(pi_rest.Unit().Dot(LL));

                                Hist_CosTheta_Q2_InvM_all->Fill(Q2,CosTheta,-Charge*Charge2);
                                Hist_CosTheta_QQ_InvM_all->Fill(QQ,CosTheta,-Charge*Charge2);
                                Hist_CosTheta_QB_InvM_all->Fill(QB,CosTheta,-Charge*Charge2);
                                if(mass_parent>0.76 && mass_parent<0.77) {
                                        Hist_CosTheta_Q2_InvM_rho->Fill(Q2,CosTheta,-Charge*Charge2);
                                        Hist_CosTheta_QQ_InvM_rho->Fill(QQ,CosTheta,-Charge*Charge2);
                                        Hist_CosTheta_QB_InvM_rho->Fill(QB,CosTheta,-Charge*Charge2);
                                }
                                if(mass_parent>0.3 && mass_parent<0.6) {
                                        Hist_CosTheta_Q2_InvM0306->Fill(Q2,CosTheta,-Charge*Charge2);
                                        Hist_CosTheta_QQ_InvM0306->Fill(QQ,CosTheta,-Charge*Charge2);
                                        Hist_CosTheta_QB_InvM0306->Fill(QB,CosTheta,-Charge*Charge2);
                                }
                                if(mass_parent>0.9 && mass_parent<1.4) {
                                        Hist_CosTheta_Q2_InvM0914->Fill(Q2,CosTheta,-Charge*Charge2);
                                        Hist_CosTheta_QQ_InvM0914->Fill(QQ,CosTheta,-Charge*Charge2);
                                        Hist_CosTheta_QB_InvM0914->Fill(QB,CosTheta,-Charge*Charge2);
                                }

                                if(Charge*Charge2>0) {
                                        Hist_ParentPhi_Q2_InvM_all->Fill(Q2,phi_parent,-1);
                                        Hist_ParentPhi_QQ_InvM_all->Fill(QQ,phi_parent,-1);
                                        Hist_ParentPhi_QB_InvM_all->Fill(QB,phi_parent,-1);
                                }
                                if(Charge*Charge2<0) {
                                        Hist_ParentPhi_Q2_InvM_all->Fill(Q2,phi_parent);
                                        Hist_ParentPhi_QQ_InvM_all->Fill(QQ,phi_parent);
                                        Hist_ParentPhi_QB_InvM_all->Fill(QB,phi_parent);
                                }
                                if(mass_parent>0.3 && mass_parent<0.9) {
                                        if(Charge*Charge2>0) {
                                                Hist_ParentPhi_Q2_InvM0309->Fill(Q2,phi_parent,-1);
                                                Hist_ParentPhi_QQ_InvM0309->Fill(QQ,phi_parent,-1);
                                                Hist_ParentPhi_QB_InvM0309->Fill(QB,phi_parent,-1);
                                        }
                                        if(Charge*Charge2<0) {
                                                Hist_ParentPhi_Q2_InvM0309->Fill(Q2,phi_parent);
                                                Hist_ParentPhi_QQ_InvM0309->Fill(QQ,phi_parent);
                                                Hist_ParentPhi_QB_InvM0309->Fill(QB,phi_parent);
                                        }

                                }
                                mQx_parent += cos(phi_parent*nHar);
                                mQy_parent += sin(phi_parent*nHar);
                                Parent_count++;

                                float v2_parent = cos(nHar*phi_parent-nHar*psi)*100;
                                float v2_parent_obs1 = cos(nHar*phi_parent - nHar*TPC_EP_east_new)*100;
                                float v2_parent_obs2 = cos(nHar*phi_parent - nHar*TPC_EP_west_new)*100;

                                Hist_Pt_parent->Fill(pt_parent);
                                Hist_Eta_parent->Fill(eta_parent);
                                Hist_Phi_parent->Fill(phi_parent);
                                pTemp_pT->Fill(2,pt_parent);
                                Hist_v2_pt_parent->Fill(pt_parent,v2_parent);
                                Hist_v2_pt_parent_obs->Fill(pt_parent,v2_parent_obs1);
                                Hist_v2_pt_parent_obs->Fill(pt_parent,v2_parent_obs2);
                                Hist_v2_eta_parent->Fill(eta_parent,v2_parent);
                                Hist_v2_eta_parent_obs->Fill(eta_parent,v2_parent_obs1);
                                Hist_v2_eta_parent_obs->Fill(eta_parent,v2_parent_obs2);

                                pTemp_v2_parent->Fill(1,v2_parent);
                                pTemp_v2_parent->Fill(2,v2_parent_obs1);
                                pTemp_v2_parent->Fill(2,v2_parent_obs2);

                                p_v2_Q->Fill(QQ,v2_parent); 
                                p_v2_p_Q2->Fill(Q2,v2_parent);
                                p_v2_Q_obs->Fill(QQ,v2_parent_obs1); p_v2_Q_obs->Fill(QQ,v2_parent_obs2);
                                p_v2_p_Q2_obs->Fill(Q2,v2_parent_obs1); p_v2_p_Q2_obs->Fill(Q2,v2_parent_obs2);

                                float correlator0 = 100*cos(Phi_new + (nHar-1)*Phi2_new - nHar*psi);
                                float correlator3 = 100*cos(Phi_new - Phi2_new);
                                float correlator4 = 50*(cos(Phi_new + (nHar-1)*Phi2_new - nHar*TPC_EP_east_new)+cos(Phi_new + (nHar-1)*Phi2_new - nHar*TPC_EP_west_new));
                                float correlator1 = 100*cos(Phi_new - (nHar+1)*Phi2_new + nHar*psi);
                                float correlator5 = 50*(cos(Phi_new - (nHar+1)*Phi2_new + nHar*TPC_EP_east_new)+cos(Phi_new - (nHar+1)*Phi2_new + nHar*TPC_EP_west_new));

                                if(Charge>0 && Charge2>0) {
                                        pParity_int_ss->Fill(1,correlator0);
                                        pParity_int_ss_obs->Fill(1,correlator4);
                                        pParity_int_ss->Fill(1+4,correlator1);
                                        pParity_int_ss_obs->Fill(1+4,correlator5);
                                        pParity_eta_ss->Fill(1,0.5*(Eta+Eta2),correlator0);
                                        pParity_eta_ss_obs->Fill(1,0.5*(Eta+Eta2),correlator4);
                                        pParity_Deta_ss->Fill(1,fabs(Eta-Eta2),correlator0);
                                        pParity_Deta_ss_obs->Fill(1,fabs(Eta-Eta2),correlator4);
                                        pParity_pt_ss->Fill(1,0.5*(Pt+Pt2),correlator0);
                                        pParity_pt_ss_obs->Fill(1,0.5*(Pt+Pt2),correlator4);
                                        pParity_Dpt_ss->Fill(1,fabs(Pt-Pt2),correlator0);
                                        pParity_Dpt_ss_obs->Fill(1,fabs(Pt-Pt2),correlator4);
                                        pTemp_parity->Fill(1,correlator0);
                                        pTemp_parity->Fill(1+4,correlator4);
                                        pTemp_parity2->Fill(1,correlator1);
                                        pTemp_parity2->Fill(1+4,correlator5);

                                        pParity_Q->Fill(1,QQ, correlator0);
                                        pParity_Q->Fill(1+4,QQ, correlator1);
                                        pParity_Q_obs->Fill(1,QQ, correlator4);
                                        pParity_Q_obs->Fill(1+4,QQ, correlator5);
                                        pParity_Q2->Fill(1,Q2, correlator0);
                                        pParity_Q2->Fill(1+4,Q2, correlator1);
                                        pParity_Q2_obs->Fill(1,Q2, correlator4);
                                        pParity_Q2_obs->Fill(1+4,Q2, correlator5);
                                        pParity_QB->Fill(1,QB, correlator0);
                                        pParity_QB->Fill(1+4,QB, correlator1);
                                        pParity_QB_obs->Fill(1,QB, correlator4);
                                        pParity_QB_obs->Fill(1+4,QB, correlator5);

                                        pDelta_int_ss->Fill(1,correlator3);
                                        pDelta_eta_ss->Fill(1,0.5*(Eta+Eta2),correlator3);
                                        pDelta_Deta_ss->Fill(1,fabs(Eta-Eta2),correlator3);
                                        pDelta_pt_ss->Fill(1,0.5*(Pt+Pt2),correlator3);
                                        pDelta_Dpt_ss->Fill(1,fabs(Pt-Pt2),correlator3);
                                        pTemp_delta->Fill(1,correlator3);
                                }
                                if(Charge<0 && Charge2<0) {
                                        pParity_int_ss->Fill(2,correlator0);
                                        pParity_int_ss_obs->Fill(2,correlator4);
                                        pParity_int_ss->Fill(2+4,correlator1);
                                        pParity_int_ss_obs->Fill(2+4,correlator5);
                                        pParity_eta_ss->Fill(2,0.5*(Eta+Eta2),correlator0);
                                        pParity_eta_ss_obs->Fill(2,0.5*(Eta+Eta2),correlator4);
                                        pParity_Deta_ss->Fill(2,fabs(Eta-Eta2),correlator0);
                                        pParity_Deta_ss_obs->Fill(2,fabs(Eta-Eta2),correlator4);
                                        pParity_pt_ss->Fill(2,0.5*(Pt+Pt2),correlator0);
                                        pParity_pt_ss_obs->Fill(2,0.5*(Pt+Pt2),correlator4);
                                        pParity_Dpt_ss->Fill(2,fabs(Pt-Pt2),correlator0);
                                        pParity_Dpt_ss_obs->Fill(2,fabs(Pt-Pt2),correlator4);
                                        pTemp_parity->Fill(2,correlator0);
                                        pTemp_parity->Fill(2+4,correlator4);
                                        pTemp_parity2->Fill(2,correlator1);
                                        pTemp_parity2->Fill(2+4,correlator5);

                                        pParity_Q->Fill(2,QQ, correlator0);
                                        pParity_Q->Fill(2+4,QQ, correlator1);
                                        pParity_Q_obs->Fill(2,QQ, correlator4);
                                        pParity_Q_obs->Fill(2+4,QQ, correlator5);
                                        pParity_Q2->Fill(2,Q2, correlator0);
                                        pParity_Q2->Fill(2+4,Q2, correlator1);
                                        pParity_Q2_obs->Fill(2,Q2, correlator4);
                                        pParity_Q2_obs->Fill(2+4,Q2, correlator5);
                                        pParity_QB->Fill(2,QB, correlator0);
                                        pParity_QB->Fill(2+4,QB, correlator1);
                                        pParity_QB_obs->Fill(2,QB, correlator4);
                                        pParity_QB_obs->Fill(2+4,QB, correlator5);

                                        pDelta_int_ss->Fill(2,correlator3);
                                        pDelta_eta_ss->Fill(2,0.5*(Eta+Eta2),correlator3);
                                        pDelta_Deta_ss->Fill(2,fabs(Eta-Eta2),correlator3);
                                        pDelta_pt_ss->Fill(2,0.5*(Pt+Pt2),correlator3);
                                        pDelta_Dpt_ss->Fill(2,fabs(Pt-Pt2),correlator3);
                                        pTemp_delta->Fill(2,correlator3);
                                        }
                                        if(Charge*Charge2>0) {
                                        Hist_Mass_parent_SS->Fill(mass_parent);
                                        pParity_int_ss->Fill(3,correlator0);
                                        pParity_int_ss_obs->Fill(3,correlator4);
                                        pParity_int_ss->Fill(3+4,correlator1);
                                        pParity_int_ss_obs->Fill(3+4,correlator5);
                                        pParity_eta_ss->Fill(3,0.5*(Eta+Eta2),correlator0);
                                        pParity_eta_ss_obs->Fill(3,0.5*(Eta+Eta2),correlator4);
                                        pParity_Deta_ss->Fill(3,fabs(Eta-Eta2),correlator0);
                                        pParity_Deta_ss_obs->Fill(3,fabs(Eta-Eta2),correlator4);
                                        pParity_pt_ss->Fill(3,0.5*(Pt+Pt2),correlator0);
                                        pParity_pt_ss_obs->Fill(3,0.5*(Pt+Pt2),correlator4);
                                        pParity_Dpt_ss->Fill(3,fabs(Pt-Pt2),correlator0);
                                        pParity_Dpt_ss_obs->Fill(3,fabs(Pt-Pt2),correlator4);
                                        pTemp_parity->Fill(3,correlator0);
                                        pTemp_parity->Fill(3+4,correlator4);
                                        pTemp_parity2->Fill(3,correlator1);
                                        pTemp_parity2->Fill(3+4,correlator5);

                                        pParity_Q->Fill(3,QQ, correlator0);
                                        pParity_Q->Fill(3+4,QQ, correlator1);
                                        pParity_Q_obs->Fill(3,QQ, correlator4);
                                        pParity_Q_obs->Fill(3+4,QQ, correlator5);
                                        pParity_Q2->Fill(3,Q2, correlator0);
                                        pParity_Q2->Fill(3+4,Q2, correlator1);
                                        pParity_Q2_obs->Fill(3,Q2, correlator4);
                                        pParity_Q2_obs->Fill(3+4,Q2, correlator5);
                                        pParity_QB->Fill(3,QB, correlator0);
                                        pParity_QB->Fill(3+4,QB, correlator1);
                                        pParity_QB_obs->Fill(3,QB, correlator4);
                                        pParity_QB_obs->Fill(3+4,QB, correlator5);

                                        pDelta_int_ss->Fill(3,correlator3);
                                        pDelta_eta_ss->Fill(3,0.5*(Eta+Eta2),correlator3);
                                        pDelta_Deta_ss->Fill(3,fabs(Eta-Eta2),correlator3);
                                        pDelta_pt_ss->Fill(3,0.5*(Pt+Pt2),correlator3);
                                        pDelta_Dpt_ss->Fill(3,fabs(Pt-Pt2),correlator3);
                                        pTemp_delta->Fill(3,correlator3);
                                        }
                                        if(Charge*Charge2<0) {
                                        Hist_Mass_parent_OS->Fill(mass_parent);
                                        pParity_int_ss->Fill(4,correlator0);
                                        pParity_int_ss_obs->Fill(4,correlator4);
                                        pParity_int_ss->Fill(4+4,correlator1);
                                        pParity_int_ss_obs->Fill(4+4,correlator5);
                                        pParity_eta_ss->Fill(4,0.5*(Eta+Eta2),correlator0);
                                        pParity_eta_ss_obs->Fill(4,0.5*(Eta+Eta2),correlator4);
                                        pParity_Deta_ss->Fill(4,fabs(Eta-Eta2),correlator0);
                                        pParity_Deta_ss_obs->Fill(4,fabs(Eta-Eta2),correlator4);
                                        pParity_pt_ss->Fill(4,0.5*(Pt+Pt2),correlator0);
                                        pParity_pt_ss_obs->Fill(4,0.5*(Pt+Pt2),correlator4);
                                        pParity_Dpt_ss->Fill(4,fabs(Pt-Pt2),correlator0);
                                        pParity_Dpt_ss_obs->Fill(4,fabs(Pt-Pt2),correlator4);
                                        pTemp_parity->Fill(4,correlator0);
                                        pTemp_parity->Fill(4+4,correlator4);
                                        pTemp_parity2->Fill(4,correlator1);
                                        pTemp_parity2->Fill(4+4,correlator5);

                                        pParity_Q->Fill(4,QQ, correlator0);
                                        pParity_Q->Fill(4+4,QQ, correlator1);
                                        pParity_Q_obs->Fill(4,QQ, correlator4);
                                        pParity_Q_obs->Fill(4+4,QQ, correlator5);
                                        pParity_Q2->Fill(4,Q2, correlator0);
                                        pParity_Q2->Fill(4+4,Q2, correlator1);
                                        pParity_Q2_obs->Fill(4,Q2, correlator4);
                                        pParity_Q2_obs->Fill(4+4,Q2, correlator5);
                                        pParity_QB->Fill(4,QB, correlator0);
                                        pParity_QB->Fill(4+4,QB, correlator1);
                                        pParity_QB_obs->Fill(4,QB, correlator4);
                                        pParity_QB_obs->Fill(4+4,QB, correlator5);

                                        pDelta_int_ss->Fill(4,correlator3);
                                        pDelta_eta_ss->Fill(4,0.5*(Eta+Eta2),correlator3);
                                        pDelta_Deta_ss->Fill(4,fabs(Eta-Eta2),correlator3);
                                        pDelta_pt_ss->Fill(4,0.5*(Pt+Pt2),correlator3);
                                        pDelta_Dpt_ss->Fill(4,fabs(Pt-Pt2),correlator3);
                                        pTemp_delta->Fill(4,correlator3);
                                }
                        } // 2nd track

                }  //1st Track

                Hist_Parent_count->Fill(1,Pcount);
                Hist_Parent_count->Fill(2,Pcount*Pcount/100.);
                Hist_Parent_count->Fill(3,Parent_count/100.);
                Hist_Parent_count->Fill(4,Parent_count*Parent_count/10000.);

                //		double QQ = (mQx_parent*mQx_parent + mQy_parent*mQy_parent)/float(Parent_count+0.05678*0.05678*Parent_count*Parent_count);

                //cout<<"mQx_parent = "<<mQx_parent<<" mQy_parent = "<<mQy_parent<<" Parent_count = "<<Parent_count<<endl;
                //cout<<" Q2 = "<<QQ<<endl; 

                float Temp_pT = pTemp_pT->GetBinContent(1);
                float Temp_pT_parent = pTemp_pT->GetBinContent(2);
                float Temp_v2 = pTemp_v2->GetBinContent(1);
                float Temp_v2b= pTemp_v2->GetBinContent(2);
                float Temp_v2_parent = pTemp_v2_parent->GetBinContent(1);
                float Temp_v2b_parent = pTemp_v2_parent->GetBinContent(2);
                Hist_v2_v2parent->Fill(Temp_v2_parent,Temp_v2);
                Hist_v2_v2parent_obs->Fill(Temp_v2b_parent,Temp_v2b);
                p_v2_QB_coarse->Fill(QB_index,Temp_v2);
                p_v2_QB_coarse_obs->Fill(QB_index,Temp_v2b);
                p_v2_Ach->Fill(Ach,Temp_v2);

                float Temp_v2p = pTemp_pion->GetBinContent(1);
                float Temp_v2n = pTemp_pion->GetBinContent(2);
                float Temp_a1p = pTemp_pion->GetBinContent(3);
                float Temp_a1n = pTemp_pion->GetBinContent(4);
                float Temp_v4p = pTemp_pion->GetBinContent(5);
                float Temp_v4n = pTemp_pion->GetBinContent(6);
                float Temp_v6p = pTemp_pion->GetBinContent(7);
                float Temp_v6n = pTemp_pion->GetBinContent(8);
                //		p_a1_v2_p->Fill(Temp_v2p,Temp_a1p);
                //		p_a1_v2_n->Fill(Temp_v2n,Temp_a1n);
                p_v2_v4_pip->Fill(Temp_v4p,Temp_v2p);
                p_v2_v4_pim->Fill(Temp_v4n,Temp_v2n);
                p_vn_pi->Fill(1,Temp_v2p); p_vn_pi->Fill(2,Temp_v4p);p_vn_pi->Fill(3,Temp_v6p);
                p_vn_pi->Fill(4,Temp_v2n); p_vn_pi->Fill(5,Temp_v4n);p_vn_pi->Fill(6,Temp_v6n);
                //		p_a1_Q2_p->Fill(Q2,Temp_a1p);
                //		p_a1_Q2_n->Fill(Q2,Temp_a1n);
                //                p_a1_QQ_p->Fill(QQ,Temp_a1p);
                //                p_a1_QQ_n->Fill(QQ,Temp_a1n);
                //		p_v2_Q2_p->Fill(Q2,Temp_v2p);
                //		p_v2_Q2_n->Fill(Q2,Temp_v2n);
                p_a1_M_p->Fill(Pcount,Temp_a1p);
                p_a1_M_n->Fill(Pcount,Temp_a1n);
                p_v2_M_p->Fill(Pcount,Temp_v2p);
                p_v2_M_n->Fill(Pcount,Temp_v2n);
                p_M_v2_p->Fill(Temp_v2p,Pcount);
                p_M_v2_n->Fill(Temp_v2n,Pcount);
                p_v2_pos_Ach->Fill(Ach,Temp_v2p*100);
                p_v2_neg_Ach->Fill(Ach,Temp_v2n*100);
                p_a1_pos_Ach->Fill(Ach,Temp_a1p*100);
                p_a1_neg_Ach->Fill(Ach,Temp_a1n*100);

                Temp_v2p = pTemp_proton->GetBinContent(1);
                Temp_v2n = pTemp_proton->GetBinContent(2);
                Temp_a1p = pTemp_proton->GetBinContent(3);
                Temp_a1n = pTemp_proton->GetBinContent(4);
                Temp_v4p = pTemp_proton->GetBinContent(5);
                Temp_v4n = pTemp_proton->GetBinContent(6);
                Temp_v6p = pTemp_proton->GetBinContent(7);
                Temp_v6n = pTemp_proton->GetBinContent(8);
                //		p_a1_v2_pp->Fill(Temp_v2p,Temp_a1p);
                //                p_a1_v2_pm->Fill(Temp_v2n,Temp_a1n);
                //		p_v2_v4_pp->Fill(Temp_v4p,Temp_v2p);
                //                p_v2_v4_pm->Fill(Temp_v4n,Temp_v2n);
                p_vn_p->Fill(1,Temp_v2p); p_vn_p->Fill(2,Temp_v4p);p_vn_p->Fill(3,Temp_v6p);
                p_vn_p->Fill(4,Temp_v2n); p_vn_p->Fill(5,Temp_v4n);p_vn_p->Fill(6,Temp_v6n);

                Temp_v2p = pTemp_kaon->GetBinContent(1);
                Temp_v2n = pTemp_kaon->GetBinContent(2);
                Temp_a1p = pTemp_kaon->GetBinContent(3);
                Temp_a1n = pTemp_kaon->GetBinContent(4);
                Temp_v4p = pTemp_kaon->GetBinContent(5);
                Temp_v4n = pTemp_kaon->GetBinContent(6);
                Temp_v6p = pTemp_kaon->GetBinContent(7);
                Temp_v6p = pTemp_kaon->GetBinContent(8);
                //                p_a1_v2_kp->Fill(Temp_v2p,Temp_a1p);
                //                p_a1_v2_km->Fill(Temp_v2n,Temp_a1n);
                //                p_v2_v4_kp->Fill(Temp_v4p,Temp_v2p);
                //                p_v2_v4_km->Fill(Temp_v4n,Temp_v2n);
                p_vn_k->Fill(1,Temp_v2p); p_vn_k->Fill(2,Temp_v4p);p_vn_k->Fill(3,Temp_v6p);
                p_vn_k->Fill(4,Temp_v2n); p_vn_k->Fill(5,Temp_v4n);p_vn_k->Fill(6,Temp_v6n);

                //g112{RP}
                float Temp_parity1 = pTemp_parity->GetBinContent(1);
                float Temp_parity2 = pTemp_parity->GetBinContent(2);
                float Temp_parity3 = pTemp_parity->GetBinContent(3);
                float Temp_parity4 = pTemp_parity->GetBinContent(4);
                //g112{PP}
                float Temp_parity1b = pTemp_parity->GetBinContent(1+4);
                float Temp_parity2b = pTemp_parity->GetBinContent(2+4);
                float Temp_parity3b = pTemp_parity->GetBinContent(3+4);
                float Temp_parity4b = pTemp_parity->GetBinContent(4+4);
                //g132{RP}
                float Temp2_parity1 = pTemp_parity2->GetBinContent(1);
                float Temp2_parity2 = pTemp_parity2->GetBinContent(2);
                float Temp2_parity3 = pTemp_parity2->GetBinContent(3);
                float Temp2_parity4 = pTemp_parity2->GetBinContent(4);
                //g132{PP}
                float Temp2_parity1b = pTemp_parity2->GetBinContent(1+4);
                float Temp2_parity2b = pTemp_parity2->GetBinContent(2+4);
                float Temp2_parity3b = pTemp_parity2->GetBinContent(3+4);
                float Temp2_parity4b = pTemp_parity2->GetBinContent(4+4);
                //delta
                float Temp_delta1 = pTemp_delta->GetBinContent(1);
                float Temp_delta2 = pTemp_delta->GetBinContent(2);
                float Temp_delta3 = pTemp_delta->GetBinContent(3);
                float Temp_delta4 = pTemp_delta->GetBinContent(4);

                p_pT_Q->Fill(QQ,Temp_pT);
                p_pT_parent_Q->Fill(QQ,Temp_pT_parent);
                p_pT_Q2->Fill(Q2,Temp_pT);
                p_pT_parent_Q2->Fill(Q2,Temp_pT_parent);

                Hist_Q->Fill(QQ);
                Hist_Q_test->Fill(sqrt(QQ));
                p_RefMult_Q->Fill(QQ,Pcount);
                p_cos_Q->Fill(QQ,cos(nHar*TPC_EP_east_new-nHar*TPC_EP_west_new), Eweight);
                if(Temp_parity3b!=0) pParity_SS_Q_obs->Fill(QQ, Temp_parity3b);
                if(Temp_parity4b!=0) pParity_OS_Q_obs->Fill(QQ, Temp_parity4b);

                if(Temp_delta1!=0) pDelta_Q->Fill(1,QQ, Temp_delta1);
                if(Temp_delta2!=0) pDelta_Q->Fill(2,QQ, Temp_delta2);
                if(Temp_delta3!=0) pDelta_Q->Fill(3,QQ, Temp_delta3);
                if(Temp_delta4!=0) pDelta_Q->Fill(4,QQ, Temp_delta4);

                Hist_Q2->Fill(Q2);
                Hist_Q4->Fill(Q4);
                if(Temp_v2!=0) {
                        pv2_PsiQ_PsiQB->Fill(Psi_B, Psi_A,Temp_v2);
                }
                if(Temp_parity1!=0) {
                        if(QB_index==0) pParity_Q2_QB1->Fill(1,Q2, Temp_parity1);
                        if(QB_index==1) pParity_Q2_QB2->Fill(1,Q2, Temp_parity1);
                        if(QB_index==2) pParity_Q2_QB3->Fill(1,Q2, Temp_parity1);
                        if(QB_index==3) pParity_Q2_QB4->Fill(1,Q2, Temp_parity1);
                }
                if(Temp_parity2!=0) {
                        if(QB_index==0) pParity_Q2_QB1->Fill(2,Q2, Temp_parity2);
                        if(QB_index==1) pParity_Q2_QB2->Fill(2,Q2, Temp_parity2);
                        if(QB_index==2) pParity_Q2_QB3->Fill(2,Q2, Temp_parity2);
                        if(QB_index==3) pParity_Q2_QB4->Fill(2,Q2, Temp_parity2);
                }
                if(Temp_parity3!=0) {
                        if(QB_index==0) pParity_Q2_QB1->Fill(3,Q2, Temp_parity3);
                        if(QB_index==1) pParity_Q2_QB2->Fill(3,Q2, Temp_parity3);
                        if(QB_index==2) pParity_Q2_QB3->Fill(3,Q2, Temp_parity3);
                        if(QB_index==3) pParity_Q2_QB4->Fill(3,Q2, Temp_parity3);
                        p_gamma_SS_Ach->Fill(Ach, Temp_parity3);
                }
                if(Temp_parity4!=0) {
                        if(QB_index==0) pParity_Q2_QB1->Fill(4,Q2, Temp_parity4);
                        if(QB_index==1) pParity_Q2_QB2->Fill(4,Q2, Temp_parity4);
                        if(QB_index==2) pParity_Q2_QB3->Fill(4,Q2, Temp_parity4);
                        if(QB_index==3) pParity_Q2_QB4->Fill(4,Q2, Temp_parity4);
                        p_gamma_OS_Ach->Fill(Ach, Temp_parity4);
                }
                if(Temp_parity1b!=0) {
                        if(QB_index==0) pParity_Q2_QB1_obs->Fill(1,Q2, Temp_parity1b);
                        if(QB_index==1) pParity_Q2_QB2_obs->Fill(1,Q2, Temp_parity1b);
                        if(QB_index==2) pParity_Q2_QB3_obs->Fill(1,Q2, Temp_parity1b);
                        if(QB_index==3) pParity_Q2_QB4_obs->Fill(1,Q2, Temp_parity1b);
                }
                if(Temp_parity2b!=0) {
                        if(QB_index==0) pParity_Q2_QB1_obs->Fill(2,Q2, Temp_parity2b);
                        if(QB_index==1) pParity_Q2_QB2_obs->Fill(2,Q2, Temp_parity2b);
                        if(QB_index==2) pParity_Q2_QB3_obs->Fill(2,Q2, Temp_parity2b);
                        if(QB_index==3) pParity_Q2_QB4_obs->Fill(2,Q2, Temp_parity2b);
                }
                if(Temp_parity3b!=0) {
                        if(QB_index==0) pParity_Q2_QB1_obs->Fill(3,Q2, Temp_parity3b);
                        if(QB_index==1) pParity_Q2_QB2_obs->Fill(3,Q2, Temp_parity3b);
                        if(QB_index==2) pParity_Q2_QB3_obs->Fill(3,Q2, Temp_parity3b);
                        if(QB_index==3) pParity_Q2_QB4_obs->Fill(3,Q2, Temp_parity3b);
                }
                if(Temp_parity4b!=0) {
                        if(QB_index==0) pParity_Q2_QB1_obs->Fill(4,Q2, Temp_parity4b);
                        if(QB_index==1) pParity_Q2_QB2_obs->Fill(4,Q2, Temp_parity4b);
                        if(QB_index==2) pParity_Q2_QB3_obs->Fill(4,Q2, Temp_parity4b);
                        if(QB_index==3) pParity_Q2_QB4_obs->Fill(4,Q2, Temp_parity4b);
                }
                if(Temp2_parity1!=0) {
                        if(QB_index==0) pParity_Q2_QB1->Fill(1+4,Q2, Temp2_parity1);
                        if(QB_index==1) pParity_Q2_QB2->Fill(1+4,Q2, Temp2_parity1);
                        if(QB_index==2) pParity_Q2_QB3->Fill(1+4,Q2, Temp2_parity1);
                        if(QB_index==3) pParity_Q2_QB4->Fill(1+4,Q2, Temp2_parity1);
                }
                if(Temp2_parity2!=0) {
                        if(QB_index==0) pParity_Q2_QB1->Fill(2+4,Q2, Temp2_parity2);
                        if(QB_index==1) pParity_Q2_QB2->Fill(2+4,Q2, Temp2_parity2);
                        if(QB_index==2) pParity_Q2_QB3->Fill(2+4,Q2, Temp2_parity2);
                        if(QB_index==3) pParity_Q2_QB4->Fill(2+4,Q2, Temp2_parity2);
                }
                if(Temp2_parity3!=0) {
                        if(QB_index==0) pParity_Q2_QB1->Fill(3+4,Q2, Temp2_parity3);
                        if(QB_index==1) pParity_Q2_QB2->Fill(3+4,Q2, Temp2_parity3);
                        if(QB_index==2) pParity_Q2_QB3->Fill(3+4,Q2, Temp2_parity3);
                        if(QB_index==3) pParity_Q2_QB4->Fill(3+4,Q2, Temp2_parity3);
                }
                if(Temp2_parity4!=0) {
                        if(QB_index==0) pParity_Q2_QB1->Fill(4+4,Q2, Temp2_parity4);
                        if(QB_index==1) pParity_Q2_QB2->Fill(4+4,Q2, Temp2_parity4);
                        if(QB_index==2) pParity_Q2_QB3->Fill(4+4,Q2, Temp2_parity4);
                        if(QB_index==3) pParity_Q2_QB4->Fill(4+4,Q2, Temp2_parity4);
                }
                if(Temp2_parity1b!=0) {
                        if(QB_index==0) pParity_Q2_QB1_obs->Fill(1+4,Q2, Temp2_parity1b);
                        if(QB_index==1) pParity_Q2_QB2_obs->Fill(1+4,Q2, Temp2_parity1b);
                        if(QB_index==2) pParity_Q2_QB3_obs->Fill(1+4,Q2, Temp2_parity1b);
                        if(QB_index==3) pParity_Q2_QB4_obs->Fill(1+4,Q2, Temp2_parity1b);
                }
                if(Temp2_parity2b!=0) {
                        if(QB_index==0) pParity_Q2_QB1_obs->Fill(2+4,Q2, Temp2_parity2b);
                        if(QB_index==1) pParity_Q2_QB2_obs->Fill(2+4,Q2, Temp2_parity2b);
                        if(QB_index==2) pParity_Q2_QB3_obs->Fill(2+4,Q2, Temp2_parity2b);
                        if(QB_index==3) pParity_Q2_QB4_obs->Fill(2+4,Q2, Temp2_parity2b);
                }
                if(Temp2_parity3b!=0) {
                        if(QB_index==0) pParity_Q2_QB1_obs->Fill(3+4,Q2, Temp2_parity3b);
                        if(QB_index==1) pParity_Q2_QB2_obs->Fill(3+4,Q2, Temp2_parity3b);
                        if(QB_index==2) pParity_Q2_QB3_obs->Fill(3+4,Q2, Temp2_parity3b);
                        if(QB_index==3) pParity_Q2_QB4_obs->Fill(3+4,Q2, Temp2_parity3b);
                }
                if(Temp2_parity4b!=0) {
                        if(QB_index==0) pParity_Q2_QB1_obs->Fill(4+4,Q2, Temp2_parity4b);
                        if(QB_index==1) pParity_Q2_QB2_obs->Fill(4+4,Q2, Temp2_parity4b);
                        if(QB_index==2) pParity_Q2_QB3_obs->Fill(4+4,Q2, Temp2_parity4b);
                        if(QB_index==3) pParity_Q2_QB4_obs->Fill(4+4,Q2, Temp2_parity4b);
                }
                if(Temp_delta1!=0) {
                        pDelta_Q2->Fill(1,Q2, Temp_delta1);
                        pDelta_QB->Fill(1,QB, Temp_delta1);
                }
                if(Temp_delta2!=0) {
                        pDelta_Q2->Fill(2,Q2, Temp_delta2);
                        pDelta_QB->Fill(2,QB, Temp_delta2);
                }
                if(Temp_delta3!=0) {
                        pDelta_Q2->Fill(3,Q2, Temp_delta3);
                        pDelta_QB->Fill(3,QB, Temp_delta3);
                        p_delta_SS_Ach->Fill(Ach, Temp_delta3);
                }
                if(Temp_delta4!=0) {
                        pDelta_Q2->Fill(4,Q2, Temp_delta4);
                        pDelta_QB->Fill(4,QB, Temp_delta4);
                        p_delta_OS_Ach->Fill(Ach, Temp_delta4);
                }

                p_Mult_Ach->Fill(Ach,Pcount);

                pTemp_pT->Reset();
                pTemp_a1->Reset();
                pTemp_v2->Reset();
                pTemp_v2_parent->Reset();
                pTemp_pion->Reset();
                pTemp_proton->Reset();
                pTemp_kaon->Reset();
                pTemp_parity->Reset();
                pTemp_parity2->Reset();
                pTemp_delta->Reset();

        } // Event



        fout.Write();
        return;
}
