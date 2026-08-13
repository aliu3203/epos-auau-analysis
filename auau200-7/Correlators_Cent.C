using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TMath.h"
#include "TProfile.h"
#include "EposPID.h"

// -----------------------------------------------------------------------------
// Centrality-differential CME/flow correlators vs the true reaction plane, for
// direct comparison with STAR (HEPData ins2928164, Au+Au 200 GeV):
//
//   v2        = <cos2(phi - Psi_RP)>                        (Fig. 5)
//   gamma112  = <cos(phi_a + phi_b - 2 Psi_RP)>             (Fig. 6)
//   gamma132  = <cos(phi_a - 3 phi_b + 2 Psi_RP)>           (Fig. 7)
//   delta     = <cos(phi_a - phi_b)>                        (Fig. 8)
//   kappa112  = Dgamma112 / (v2 * Ddelta)                   (Fig. 13)
//   kappa132  = Dgamma132 / (v2 * Ddelta)                   (Fig. 13)
// (the kappa definitions were verified against STAR's own published numbers:
//  their gamma/delta/v2 values reproduce their quoted kappa to 0.2%)
//
// All correlators are built from charge-separated Q-vectors -- O(N) per event,
// no pair loop. With Q1=Sum e^{i phi}, Q2=Sum e^{2i phi}, Q3=Sum e^{3i phi}:
//
//   gamma112_OS  = Re[Q1+ Q1- e^{-2iPsi}] / (N+ N-)
//   gamma112_SS  = Sum_c Re[(Q1c^2 - Q2c)/2 * e^{-2iPsi}] / Sum_c N_c(N_c-1)/2
//   gamma132_OS  = (Re[Q1+ conj(Q3-) e^{2iPsi}] + Re[Q1- conj(Q3+) e^{2iPsi}]) / (2 N+ N-)
//   gamma132_SS  = Sum_c Re[(Q1c conj(Q3c) - conj(Q2c)) e^{2iPsi}] / Sum_c N_c(N_c-1)
//   delta_OS     = Re[Q1+ conj(Q1-)] / (N+ N-)
//   delta_SS     = Sum_c (|Q1c|^2 - N_c)/2 / Sum_c N_c(N_c-1)/2
// gamma132 is asymmetric in (a,b), so ORDERED pairs are used for it (hence the
// factor 2 N+N- and N_c(N_c-1) denominators) while the symmetric 112/delta use
// unordered pairs.
//
// NOTE on the gamma132 pair convention (all of this verified by the brute-force pair
// loops in ValidateQvec.C).  Both correlators average over the SAME unordered pair set
// -- ~4751 OS and ~4729 SS pairs per event for charged mesons -- and both profiles are
// filled once per event with an equally-weighted per-event mean.  The difference is
// what each pair contributes.  Because cos(phi_a-3phi_b+2Psi) is not symmetric in
// (a,b), the symmetrised form above gives each pair the average of its two orderings,
// and that average obeys the exact identity
//
//     (1/2)[cos(pa-3pb+2Psi) + cos(pb-3pa+2Psi)] = cos(pa+pb-2Psi) * cos(2(pa-pb))
//
// i.e. the gamma132 pair term IS the gamma112 pair term times cos(2 dphi).  Since
// cos(2 dphi) has RMS 1/sqrt(2), the per-pair RMS drops from 0.7075 to 0.5003 and
// gamma132's error comes out ~sqrt(2) smaller -- an intrinsic property of the
// observable, not a filling mismatch.
//
// ORD1 = 1 (default) instead keeps ONE ordering per pair (alpha = the positive track
// for OS, the lower-index track for SS), so gamma132 and gamma112 sum exactly one
// cosine per pair.  That restores per-pair RMS 0.7071 ~ gamma112's 0.7075 and hence
// comparable errors.  The extra antisymmetric piece it carries is zero-mean, so the
// central value is unchanged -- only the error grows.  Both versions are computed and
// reported side by side.
//
// Centrality comes from refmult percentile cuts in cen_cuts.txt (see
// CentralityCalib.C). If that file is absent the macro runs in INCLUSIVE mode:
// one bin holding every event that passes, which is the correct behaviour for a
// b-restricted (single-centrality) sample like auau200-7.
//
// Statistical errors: gamma/delta/v2 from the event-to-event spread (single-bin
// TProfile GetBinError). kappa is a ratio of quantities measured on the SAME
// events, so its error uses sub-sampling (NSUB subgroups per centrality bin),
// which carries the Dgamma-Ddelta covariance.
// -----------------------------------------------------------------------------

const int ID_pip = 120, ID_pim = -120;
const float MPI_CHG = 0.13957;
const float yCut = 1.0, ptMin = 0.2, ptMax = 2.0;
// NSUB subgroups -> the sub-sample error estimate itself carries a relative
// uncertainty of 1/sqrt(2(NSUB-1)); 20 subgroups gives a sloppy +/-16% on the
// error bar, 100 tightens that to +/-7%.
const int NSUB = 100;
// gamma132 with one ordering per pair (identical filling to gamma112).  The SS part
// cannot be written with Q-vectors -- it depends on which track of the pair is called
// alpha -- so it costs an explicit ~4700-pair loop per event (~15 s over the full set).
const bool ORD1 = true;

// poiMode 0 = charged pions (the earlier auau200-7 results)
//         1 = all charged hadrons, including p/pbar (STAR's Fig.16 definition,
//             used only at 11.5-19.6 GeV -- NOT the 200 GeV figures)
//         2 = charged mesons, pi+/- and K+/- (STAR's POI for the 200 GeV
//             Figs. 5,6,7,8,13 -- "v2 of charged mesons"; p/pbar excluded).
//             THIS is the mode to use when comparing with the 200 GeV data.
void Correlators_Cent(int job = 0, int poiMode = 0, const char* filepat = "z-auau_run_*.root"){

    if(!EposPID::Load("idt.dt")) return;

    //------------------- centrality cuts (optional) -------------------
    int cut[NCENT]; bool useCent = false;
    {
        ifstream in("cen_cuts.txt");
        if(in){
            int n=0; string line;
            while(getline(in,line)){
                if(line.empty() || line[0]=='#') continue;
                if(n<NCENT) cut[n++] = atoi(line.c_str());
            }
            if(n==NCENT){ useCent=true;
                printf("centrality cuts from cen_cuts.txt: ");
                for(int i=0;i<NCENT;i++) printf("%d ",cut[i]);
                printf("\n");
            }
        }
    }
    const int NB = useCent ? NCENT : 1;
    if(!useCent)
        printf("\n*** cen_cuts.txt not found -> INCLUSIVE mode (single bin, all events).\n"
               "*** This is correct for a b-restricted sample; run CentralityCalib.C on a\n"
               "*** minimum-bias sample first to enable centrality-differential output.\n\n");

    //------------------- input -------------------
    TChain* chain = new TChain("teposevent");
    chain->Add(filepat);
    chain->SetBranchStatus("*",0);
    const char* act[] = {"np","bim","phi","px","py","pz","id","ist","npartproj","nparttarg"};
    for(int i=0;i<10;i++) chain->SetBranchStatus(act[i],1);

    char fname_out[200];
    sprintf(fname_out,"cen.correlators_cent_poi%d_job%d.root",poiMode,job);
    TFile fout(fname_out,"RECREATE");

    //------------------- output profiles, one set per centrality bin -------------------
    TProfile *pV2[NCENT], *pG112OS[NCENT], *pG112SS[NCENT], *pG132OS[NCENT], *pG132SS[NCENT];
    TProfile *pDOS[NCENT], *pDSS[NCENT], *pDg112[NCENT], *pDg132[NCENT], *pDd[NCENT];
    TProfile *pNpart[NCENT], *pRefmult[NCENT], *pBim[NCENT];
    // Cross-product profiles.  A TProfile carries Sum(x), Sum(x^2), N per bin, which
    // gives Var of each ingredient but NOT the covariance between two of them.  Filling
    // a profile with the per-event PRODUCT supplies the missing piece:
    //     Cov(X,Y) = <XY> - <X><Y>
    // With these, kappa's error follows analytically from the full covariance matrix --
    // no sub-sampling needed.  A1 = Dgamma112, A3 = Dgamma132, V = v2, D = Ddelta.
    // pA1A3 additionally carries cov(Dgamma112, Dgamma132), which is what decides
    // whether the observed kappa112 != kappa132 is a real effect: the two are measured
    // on the same events, so the difference cannot be errored by naive quadrature.
    TProfile *pA1V[NCENT], *pA1D[NCENT], *pA3V[NCENT], *pA3D[NCENT], *pVD[NCENT], *pA1A3[NCENT];
    // one-ordering gamma132 (identical filling to gamma112) + its own cross products
    TProfile *pG132OS1[NCENT], *pG132SS1[NCENT], *pDg132_1[NCENT];
    TProfile *pA31V[NCENT], *pA31D[NCENT], *pA1A31[NCENT];
    for(int ic=0; ic<NB; ic++){
        pV2[ic]    = new TProfile(Form("pV2_c%d",ic),    Form("v_{2} %s",CentLabel(ic)),1,0,1);
        pG112OS[ic]= new TProfile(Form("pG112OS_c%d",ic),Form("#gamma^{112}_{OS} %s",CentLabel(ic)),1,0,1);
        pG112SS[ic]= new TProfile(Form("pG112SS_c%d",ic),Form("#gamma^{112}_{SS} %s",CentLabel(ic)),1,0,1);
        pG132OS[ic]= new TProfile(Form("pG132OS_c%d",ic),Form("#gamma^{132}_{OS} %s",CentLabel(ic)),1,0,1);
        pG132SS[ic]= new TProfile(Form("pG132SS_c%d",ic),Form("#gamma^{132}_{SS} %s",CentLabel(ic)),1,0,1);
        pDOS[ic]   = new TProfile(Form("pDOS_c%d",ic),   Form("#delta_{OS} %s",CentLabel(ic)),1,0,1);
        pDSS[ic]   = new TProfile(Form("pDSS_c%d",ic),   Form("#delta_{SS} %s",CentLabel(ic)),1,0,1);
        pDg112[ic] = new TProfile(Form("pDg112_c%d",ic), Form("#Delta#gamma^{112} %s",CentLabel(ic)),1,0,1);
        pDg132[ic] = new TProfile(Form("pDg132_c%d",ic), Form("#Delta#gamma^{132} %s",CentLabel(ic)),1,0,1);
        pDd[ic]    = new TProfile(Form("pDd_c%d",ic),    Form("#Delta#delta %s",CentLabel(ic)),1,0,1);
        pNpart[ic] = new TProfile(Form("pNpart_c%d",ic), Form("N_{part} %s",CentLabel(ic)),1,0,1);
        pRefmult[ic]=new TProfile(Form("pRefmult_c%d",ic),Form("refmult %s",CentLabel(ic)),1,0,1);
        pBim[ic]   = new TProfile(Form("pBim_c%d",ic),   Form("b %s",CentLabel(ic)),1,0,1);
        pA1V[ic]   = new TProfile(Form("pA1V_c%d",ic),   Form("<#Delta#gamma^{112} v_{2}> %s",CentLabel(ic)),1,0,1);
        pA1D[ic]   = new TProfile(Form("pA1D_c%d",ic),   Form("<#Delta#gamma^{112} #Delta#delta> %s",CentLabel(ic)),1,0,1);
        pA3V[ic]   = new TProfile(Form("pA3V_c%d",ic),   Form("<#Delta#gamma^{132} v_{2}> %s",CentLabel(ic)),1,0,1);
        pA3D[ic]   = new TProfile(Form("pA3D_c%d",ic),   Form("<#Delta#gamma^{132} #Delta#delta> %s",CentLabel(ic)),1,0,1);
        pVD[ic]    = new TProfile(Form("pVD_c%d",ic),    Form("<v_{2} #Delta#delta> %s",CentLabel(ic)),1,0,1);
        pA1A3[ic]  = new TProfile(Form("pA1A3_c%d",ic),  Form("<#Delta#gamma^{112} #Delta#gamma^{132}> %s",CentLabel(ic)),1,0,1);
        pG132OS1[ic]=new TProfile(Form("pG132OS1_c%d",ic),Form("#gamma^{132,1ord}_{OS} %s",CentLabel(ic)),1,0,1);
        pG132SS1[ic]=new TProfile(Form("pG132SS1_c%d",ic),Form("#gamma^{132,1ord}_{SS} %s",CentLabel(ic)),1,0,1);
        pDg132_1[ic]=new TProfile(Form("pDg132_1_c%d",ic),Form("#Delta#gamma^{132,1ord} %s",CentLabel(ic)),1,0,1);
        pA31V[ic]  = new TProfile(Form("pA31V_c%d",ic),  Form("<#Delta#gamma^{132,1ord} v_{2}> %s",CentLabel(ic)),1,0,1);
        pA31D[ic]  = new TProfile(Form("pA31D_c%d",ic),  Form("<#Delta#gamma^{132,1ord} #Delta#delta> %s",CentLabel(ic)),1,0,1);
        pA1A31[ic] = new TProfile(Form("pA1A31_c%d",ic), Form("<#Delta#gamma^{112} #Delta#gamma^{132,1ord}> %s",CentLabel(ic)),1,0,1);
    }

    // sub-sampling accumulators
    // per-track trig, kept only when ORD1 needs the SS pair loop
    static vector<double> c1v[2], s1v[2], c3v[2], s3v[2];

    static double sV2[NCENT][NSUB], sG1[NCENT][NSUB], sG3[NCENT][NSUB], sDd[NCENT][NSUB];
    static double sG31[NCENT][NSUB];
    static long   cSub[NCENT][NSUB];
    memset(sV2,0,sizeof(sV2)); memset(sG1,0,sizeof(sG1));
    memset(sG3,0,sizeof(sG3)); memset(sDd,0,sizeof(sDd)); memset(cSub,0,sizeof(cSub));
    memset(sG31,0,sizeof(sG31));
    long nUsedBin[NCENT]; memset(nUsedBin,0,sizeof(nUsedBin));

    Long64_t nentries = chain->GetEntries();
    cout << nentries << " entries in chain\n";
    Long64_t nUsed = 0;

    for(Long64_t i=0;i<nentries;i++){
        if((i+1)%10000==0) cout<<"Processing entry == "<<i+1<<" == out of "<<nentries<<".\n";
        chain->GetEntry(i);

        int   NPTracks = (int)chain->GetLeaf("np")->GetValue(0);
        float PsiRP    = chain->GetLeaf("phi")->GetValue(0);
        float bim      = chain->GetLeaf("bim")->GetValue(0);
        int   Npart    = (int)chain->GetLeaf("npartproj")->GetValue(0)
                       + (int)chain->GetLeaf("nparttarg")->GetValue(0);
        double C2 = cos(2.*PsiRP), S2 = sin(2.*PsiRP);

        TLeaf* lpx=chain->GetLeaf("px"); TLeaf* lpy=chain->GetLeaf("py");
        TLeaf* lpz=chain->GetLeaf("pz"); TLeaf* lid=chain->GetLeaf("id");
        TLeaf* lst=chain->GetLeaf("ist");

        // charge-separated Q-vectors over the POI, plus refmult, in one track loop
        double pP=0,qP=0,p2P=0,q2P=0,p3P=0,q3P=0; int Np=0;
        double pM=0,qM=0,p2M=0,q2M=0,p3M=0,q3M=0; int Nm=0;
        int refmult = 0;
        if(ORD1) for(int k=0;k<2;k++){ c1v[k].clear(); s1v[k].clear(); c3v[k].clear(); s3v[k].clear(); }

        for(int trk=0; trk<NPTracks; trk++){
            if((int)lst->GetValue(trk) != 0) continue;              // final state
            int pid = (int)lid->GetValue(trk);
            int chg = EposPID::Charge(pid);
            if(chg == 0) continue;                                   // neutral: no use here

            float px=lpx->GetValue(trk), py=lpy->GetValue(trk), pz=lpz->GetValue(trk);
            float pt2 = px*px+py*py;

            // ---- refmult (centrality estimator) ----
            if(pt2 >= REFMULT_PTMIN*REFMULT_PTMIN && pz*pz < SINH_ETACUT_SQ*pt2) refmult++;

            // ---- POI selection ----
            int apid = abs(pid);
            bool isPOI;
            if      (poiMode==0) isPOI = (apid==120);                 // pi+/- only
            else if (poiMode==2) isPOI = (apid==120 || apid==130);    // charged mesons: pi+/-, K+/-
            else                 isPOI = true;                        // all charged hadrons
            if(!isPOI) continue;
            float pt = sqrt(pt2);
            if(pt < ptMin || pt > ptMax) continue;
            float m  = (poiMode==0) ? MPI_CHG : EposPID::Mass(pid);
            float E  = sqrt(m*m + pt2 + pz*pz);
            float y  = 0.5*log((E+pz)/(E-pz));
            if(fabs(y) > yCut) continue;

            float phi = atan2(py,px);
            double c1=cos(phi),   s1=sin(phi);
            double c2=cos(2.*phi),s2=sin(2.*phi);
            double c3=cos(3.*phi),s3=sin(3.*phi);
            if(chg > 0){ pP+=c1; qP+=s1; p2P+=c2; q2P+=s2; p3P+=c3; q3P+=s3; Np++; }
            else       { pM+=c1; qM+=s1; p2M+=c2; q2M+=s2; p3M+=c3; q3M+=s3; Nm++; }
            if(ORD1){ int k = (chg>0)?0:1;
                      c1v[k].push_back(c1); s1v[k].push_back(s1);
                      c3v[k].push_back(c3); s3v[k].push_back(s3); }
        }

        if(Np < 2 || Nm < 2) continue;

        int ic = useCent ? CentBinFromRefmult(refmult, cut) : 0;
        if(ic < 0) continue;                                        // more peripheral than 80%

        double NpOS = (double)Np*Nm;                                // unordered OS pairs
        double NpSS = 0.5*Np*(Np-1) + 0.5*Nm*(Nm-1);                // unordered SS pairs
        double NoOS = 2.0*Np*Nm;                                    // ordered OS pairs
        double NoSS = (double)Np*(Np-1) + (double)Nm*(Nm-1);        // ordered SS pairs

        //---- gamma112 ----
        double g112OS = ((pP*pM - qP*qM)*C2 + (pP*qM + qP*pM)*S2) / NpOS;
        double g112SS = ( 0.5*((pP*pP-qP*qP-p2P)*C2 + (2.*pP*qP-q2P)*S2)
                        + 0.5*((pM*pM-qM*qM-p2M)*C2 + (2.*pM*qM-q2M)*S2) ) / NpSS;

        //---- gamma132 : Re[Q1a conj(Q3b) e^{2iPsi}] = X*C2 - Y*S2 ----
        double Xpm = pP*p3M + qP*q3M,  Ypm = qP*p3M - pP*q3M;       // a in +, b in -
        double Xmp = pM*p3P + qM*q3P,  Ymp = qM*p3P - pM*q3P;       // a in -, b in +
        double g132OS = ((Xpm*C2 - Ypm*S2) + (Xmp*C2 - Ymp*S2)) / NoOS;

        double XppSelf = p2P, YppSelf = -q2P;                        // conj(Q2+)
        double XmmSelf = p2M, YmmSelf = -q2M;                        // conj(Q2-)
        double Xpp = pP*p3P + qP*q3P,  Ypp = qP*p3P - pP*q3P;
        double Xmm = pM*p3M + qM*q3M,  Ymm = qM*p3M - pM*q3M;
        double g132SS = ( ((Xpp-XppSelf)*C2 - (Ypp-YppSelf)*S2)
                        + ((Xmm-XmmSelf)*C2 - (Ymm-YmmSelf)*S2) ) / NoSS;

        //---- gamma132, ONE ordering per pair: exactly one cosine per unordered pair,
        //     the same pair set and the same denominators gamma112 uses.
        //     OS: alpha = the positive track.  SS: alpha = the lower-index track,
        //     which needs a real pair loop (no Q-vector form exists for it).
        double g132OS1 = (Xpm*C2 - Ypm*S2) / NpOS;
        double g132SS1 = 0.0;
        if(ORD1){
            double acc = 0.0;
            for(int k=0;k<2;k++){
                const vector<double>& C1=c1v[k]; const vector<double>& S1=s1v[k];
                const vector<double>& C3=c3v[k]; const vector<double>& S3=s3v[k];
                const int n = (int)C1.size();
                for(int a=0;a<n;a++){
                    const double c1a=C1[a], s1a=S1[a];
                    for(int b=a+1;b<n;b++){
                        // Re[ e^{i pa} conj(e^{3i pb}) e^{2i Psi} ]
                        double X = c1a*C3[b] + s1a*S3[b];
                        double Y = s1a*C3[b] - c1a*S3[b];
                        acc += X*C2 - Y*S2;
                    }
                }
            }
            g132SS1 = acc / NpSS;
        }

        //---- delta ----
        double dOS = (pP*pM + qP*qM) / NpOS;
        double dSS = ( 0.5*(pP*pP+qP*qP-Np) + 0.5*(pM*pM+qM*qM-Nm) ) / NpSS;

        //---- integrated v2 of the POI ----
        double v2ev = ((p2P+p2M)*C2 + (q2P+q2M)*S2) / (double)(Np+Nm);

        pV2[ic]    ->Fill(0.5, v2ev);
        pG112OS[ic]->Fill(0.5, g112OS);   pG112SS[ic]->Fill(0.5, g112SS);
        pG132OS[ic]->Fill(0.5, g132OS);   pG132SS[ic]->Fill(0.5, g132SS);
        pDOS[ic]   ->Fill(0.5, dOS);      pDSS[ic]   ->Fill(0.5, dSS);
        pDg112[ic] ->Fill(0.5, g112OS-g112SS);
        pDg132[ic] ->Fill(0.5, g132OS-g132SS);
        pDd[ic]    ->Fill(0.5, dOS-dSS);
        pNpart[ic] ->Fill(0.5, Npart);
        pRefmult[ic]->Fill(0.5, refmult);
        pBim[ic]   ->Fill(0.5, bim);

        // per-event cross products -> covariance matrix of the kappa ingredients
        double A1ev = g112OS-g112SS, A3ev = g132OS-g132SS, Ddev = dOS-dSS;
        pA1V[ic]->Fill(0.5, A1ev*v2ev);
        pA1D[ic]->Fill(0.5, A1ev*Ddev);
        pA3V[ic]->Fill(0.5, A3ev*v2ev);
        pA3D[ic]->Fill(0.5, A3ev*Ddev);
        pVD[ic] ->Fill(0.5, v2ev*Ddev);
        pA1A3[ic]->Fill(0.5, A1ev*A3ev);

        // one-ordering gamma132: filled in exactly the same place, same weight, same
        // pair set and denominators as gamma112
        double A31ev = g132OS1 - g132SS1;
        pG132OS1[ic]->Fill(0.5, g132OS1);
        pG132SS1[ic]->Fill(0.5, g132SS1);
        pDg132_1[ic]->Fill(0.5, A31ev);
        pA31V[ic] ->Fill(0.5, A31ev*v2ev);
        pA31D[ic] ->Fill(0.5, A31ev*Ddev);
        pA1A31[ic]->Fill(0.5, A1ev*A31ev);
        sG31[ic][nUsedBin[ic] % NSUB] += A31ev;

        int is = nUsedBin[ic] % NSUB;
        sV2[ic][is] += v2ev;
        sG1[ic][is] += (g112OS-g112SS);
        sG3[ic][is] += (g132OS-g132SS);
        sDd[ic][is] += (dOS-dSS);
        cSub[ic][is]++;
        nUsedBin[ic]++;
        nUsed++;
    }

    //------------------- kappa + sub-sampling errors -------------------
    TH1D* hK112 = new TH1D("hKappa112","#kappa^{112} vs centrality;centrality [%];#kappa^{112}",NB,0,NB);
    TH1D* hK132 = new TH1D("hKappa132","#kappa^{132} vs centrality;centrality [%];#kappa^{132}",NB,0,NB);
    TH1D* hNev  = new TH1D("hNev","events per centrality bin",NB,0,NB);

    printf("\n================ EPOS correlators%s ================\n",
           useCent? " vs centrality" : " (inclusive)");
    const char* poiName = (poiMode==0)? "charged pions (pi+/-)"
                        : (poiMode==2)? "charged mesons (pi+/-, K+/-)  [STAR 200 GeV POI]"
                                      : "all charged hadrons incl. p/pbar";
    printf("POI = %s,  |y|<1, %.1f<pT<%.1f,  true reaction plane\n", poiName, ptMin, ptMax);
    printf("events used = %lld\n\n", nUsed);

    for(int ic=0; ic<NB; ic++){
        double v2f = pV2[ic]->GetBinContent(1);
        double g1f = pDg112[ic]->GetBinContent(1);
        double g3f = pDg132[ic]->GetBinContent(1);
        double ddf = pDd[ic]->GetBinContent(1);
        double den = v2f*ddf;
        double k112 = den!=0 ? g1f/den : 0.0;
        double k132 = den!=0 ? g3f/den : 0.0;

        // sub-sampling errors on the two kappas.  We keep the per-subgroup value
        // of every ingredient too, so the sub-sample error on Dgamma/Ddelta/v2
        // can be compared against the independent TProfile (event-spread) error
        // -- if the two disagree, the sub-sampling machinery is broken.
        double kk1[NSUB], kk3[NSUB], vv[NSUB], gg1[NSUB], gg3[NSUB], ddv[NSUB], den_m[NSUB];
        double kk31[NSUB], gg31[NSUB];
        int M=0; double m1=0,m3=0,m31=0;
        for(int m=0;m<NSUB;m++){
            if(cSub[ic][m] <= 0) continue;
            double v=sV2[ic][m]/cSub[ic][m], g1=sG1[ic][m]/cSub[ic][m];
            double g3=sG3[ic][m]/cSub[ic][m], dd=sDd[ic][m]/cSub[ic][m];
            double g31=sG31[ic][m]/cSub[ic][m];
            double d = v*dd;
            vv[M]=v; gg1[M]=g1; gg3[M]=g3; ddv[M]=dd; den_m[M]=d; gg31[M]=g31;
            kk1[M] = d!=0? g1/d : 0.0;
            kk3[M] = d!=0? g3/d : 0.0;
            kk31[M]= d!=0? g31/d: 0.0;
            m1 += kk1[M]; m3 += kk3[M]; m31 += kk31[M]; M++;
        }
        double e112=0, e132=0, e132_1=0, sub_g31=0;
        // sub-sample errors on the ingredients + numerator/denominator correlations
        double sub_v2=0, sub_g1=0, sub_g3=0, sub_dd=0, sub_den=0, rho1=0, rho3=0;
        double naive1=0, naive3=0;
        if(M>=2){
            m1/=M; m3/=M; m31/=M;
            double s1=0,s3=0,s31=0,q31=0,ag31=0;
            for(int m=0;m<M;m++){ s1+=(kk1[m]-m1)*(kk1[m]-m1); s3+=(kk3[m]-m3)*(kk3[m]-m3);
                                  s31+=(kk31[m]-m31)*(kk31[m]-m31); ag31+=gg31[m]; }
            double dn=(double)M*(M-1);
            e112=sqrt(s1/dn); e132=sqrt(s3/dn); e132_1=sqrt(s31/dn);
            ag31/=M;
            for(int m=0;m<M;m++) q31+=(gg31[m]-ag31)*(gg31[m]-ag31);
            sub_g31=sqrt(q31/dn);

            // means of the ingredients across subgroups
            double av=0,ag1=0,ag3=0,add=0,aden=0;
            for(int m=0;m<M;m++){ av+=vv[m]; ag1+=gg1[m]; ag3+=gg3[m]; add+=ddv[m]; aden+=den_m[m]; }
            av/=M; ag1/=M; ag3/=M; add/=M; aden/=M;
            double qv=0,q1=0,q3=0,qd=0,qden=0,c1d=0,c3d=0;
            for(int m=0;m<M;m++){
                double dv=vv[m]-av, d1=gg1[m]-ag1, d3=gg3[m]-ag3, dd_=ddv[m]-add, dq=den_m[m]-aden;
                qv+=dv*dv; q1+=d1*d1; q3+=d3*d3; qd+=dd_*dd_; qden+=dq*dq;
                c1d+=d1*dq; c3d+=d3*dq;
            }
            sub_v2=sqrt(qv/dn); sub_g1=sqrt(q1/dn); sub_g3=sqrt(q3/dn);
            sub_dd=sqrt(qd/dn); sub_den=sqrt(qden/dn);
            rho1 = (q1*qden>0)? c1d/sqrt(q1*qden) : 0.0;
            rho3 = (q3*qden>0)? c3d/sqrt(q3*qden) : 0.0;
            // naive propagation IGNORING the numerator-denominator covariance
            if(g1f!=0 && den!=0) naive1 = fabs(k112)*sqrt(pow(sub_g1/g1f,2)+pow(sub_den/den,2));
            if(g3f!=0 && den!=0) naive3 = fabs(k132)*sqrt(pow(sub_g3/g3f,2)+pow(sub_den/den,2));
        }

        // ---------- analytic errors: TProfile variances + measured covariances ----------
        // kappa = A/(V D)  =>  (s_k/k)^2 = (sA/A)^2 + (sV/V)^2 + (sD/D)^2
        //                                 - 2cov(A,V)/(A V) - 2cov(A,D)/(A D) + 2cov(V,D)/(V D)
        // s_X is TProfile GetBinError (already the error ON THE MEAN); the covariances are
        // likewise covariances of the means, cov(X,Y) = (<XY> - <X><Y>)/(N-1).
        double Nev = pV2[ic]->GetBinEntries(1);
        double sA1 = pDg112[ic]->GetBinError(1), sA3 = pDg132[ic]->GetBinError(1);
        double sV  = pV2[ic]->GetBinError(1),    sD  = pDd[ic]->GetBinError(1);
        double cA1V=0,cA1D=0,cA3V=0,cA3D=0,cVD=0;
        if(Nev>1){
            double inv = 1.0/(Nev-1.0);
            cA1V = (pA1V[ic]->GetBinContent(1) - g1f*v2f)*inv;
            cA1D = (pA1D[ic]->GetBinContent(1) - g1f*ddf)*inv;
            cA3V = (pA3V[ic]->GetBinContent(1) - g3f*v2f)*inv;
            cA3D = (pA3D[ic]->GetBinContent(1) - g3f*ddf)*inv;
            cVD  = (pVD[ic] ->GetBinContent(1) - v2f*ddf)*inv;
        }
        double a112=0, a132=0;
        if(g1f!=0 && v2f!=0 && ddf!=0){
            double r2 = sA1*sA1/(g1f*g1f) + sV*sV/(v2f*v2f) + sD*sD/(ddf*ddf)
                      - 2.*cA1V/(g1f*v2f) - 2.*cA1D/(g1f*ddf) + 2.*cVD/(v2f*ddf);
            a112 = (r2>0)? fabs(k112)*sqrt(r2) : 0.0;
        }
        if(g3f!=0 && v2f!=0 && ddf!=0){
            double r2 = sA3*sA3/(g3f*g3f) + sV*sV/(v2f*v2f) + sD*sD/(ddf*ddf)
                      - 2.*cA3V/(g3f*v2f) - 2.*cA3D/(g3f*ddf) + 2.*cVD/(v2f*ddf);
            a132 = (r2>0)? fabs(k132)*sqrt(r2) : 0.0;
        }
        // correlation coefficients (scale factors cancel in the ratio)
        double rA1V = (sA1*sV>0)? cA1V/(sA1*sV):0.0, rA1D = (sA1*sD>0)? cA1D/(sA1*sD):0.0;
        double rA3V = (sA3*sV>0)? cA3V/(sA3*sV):0.0, rA3D = (sA3*sD>0)? cA3D/(sA3*sD):0.0;
        double rVD  = (sV*sD>0)?  cVD /(sV*sD) :0.0;

        hK112->SetBinContent(ic+1,k112); hK112->SetBinError(ic+1,e112);
        hK132->SetBinContent(ic+1,k132); hK132->SetBinError(ic+1,e132);
        hNev ->SetBinContent(ic+1,nUsedBin[ic]);

        printf("---- %s  (Nev=%ld, <refmult>=%.1f, <Npart>=%.1f, <b>=%.2f fm) ----\n",
               useCent? CentLabel(ic):"inclusive", nUsedBin[ic],
               pRefmult[ic]->GetBinContent(1), pNpart[ic]->GetBinContent(1), pBim[ic]->GetBinContent(1));
        printf("   %-12s %+12.5e +/- %10.3e\n","v2",       v2f, pV2[ic]->GetBinError(1));
        printf("   %-12s %+12.5e +/- %10.3e\n","g112_OS",  pG112OS[ic]->GetBinContent(1), pG112OS[ic]->GetBinError(1));
        printf("   %-12s %+12.5e +/- %10.3e\n","g112_SS",  pG112SS[ic]->GetBinContent(1), pG112SS[ic]->GetBinError(1));
        printf("   %-12s %+12.5e +/- %10.3e\n","Dgamma112",g1f, pDg112[ic]->GetBinError(1));
        printf("   %-12s %+12.5e +/- %10.3e\n","g132_OS",  pG132OS[ic]->GetBinContent(1), pG132OS[ic]->GetBinError(1));
        printf("   %-12s %+12.5e +/- %10.3e\n","g132_SS",  pG132SS[ic]->GetBinContent(1), pG132SS[ic]->GetBinError(1));
        printf("   %-12s %+12.5e +/- %10.3e\n","Dgamma132",g3f, pDg132[ic]->GetBinError(1));
        printf("   %-12s %+12.5e +/- %10.3e\n","delta_OS", pDOS[ic]->GetBinContent(1), pDOS[ic]->GetBinError(1));
        printf("   %-12s %+12.5e +/- %10.3e\n","delta_SS", pDSS[ic]->GetBinContent(1), pDSS[ic]->GetBinError(1));
        printf("   %-12s %+12.5e +/- %10.3e\n","Ddelta",   ddf, pDd[ic]->GetBinError(1));
        printf("   %-12s %+12.5f +/- %10.3f\n","kappa112", k112, e112);
        printf("   %-12s %+12.5f +/- %10.3f\n","kappa132", k132, e132);

        // ---- gamma132 with ONE ordering per pair: identical filling to gamma112 ----
        double g31f = pDg132_1[ic]->GetBinContent(1);
        double k132_1 = den!=0 ? g31f/den : 0.0;
        if(ORD1){
            printf("\n   [gamma132 with ONE ordering per pair -- exactly one cosine per\n"
                   "    unordered pair, same pair set/denominators as gamma112]\n");
            printf("   %-12s %+12.5e +/- %10.3e\n","g132_OS1", pG132OS1[ic]->GetBinContent(1), pG132OS1[ic]->GetBinError(1));
            printf("   %-12s %+12.5e +/- %10.3e\n","g132_SS1", pG132SS1[ic]->GetBinContent(1), pG132SS1[ic]->GetBinError(1));
            printf("   %-12s %+12.5e +/- %10.3e\n","Dgamma132_1",g31f, pDg132_1[ic]->GetBinError(1));
            printf("   %-12s %+12.5f +/- %10.3f\n","kappa132_1", k132_1, e132_1);
        }

        // ---- diagnostics: is the sub-sampling error self-consistent? ----
        printf("\n   [sub-sampling cross-check, M=%d subgroups]\n", M);
        printf("   %-12s %10s %10s %8s\n","quantity","TProfile","sub-samp","ratio");
        printf("   %-12s %10.3e %10.3e %8.3f\n","v2",       pV2[ic]->GetBinError(1),   sub_v2, sub_v2/pV2[ic]->GetBinError(1));
        printf("   %-12s %10.3e %10.3e %8.3f\n","Dgamma112",pDg112[ic]->GetBinError(1),sub_g1, sub_g1/pDg112[ic]->GetBinError(1));
        printf("   %-12s %10.3e %10.3e %8.3f\n","Dgamma132",pDg132[ic]->GetBinError(1),sub_g3, sub_g3/pDg132[ic]->GetBinError(1));
        printf("   %-12s %10.3e %10.3e %8.3f\n","Ddelta",   pDd[ic]->GetBinError(1),   sub_dd, sub_dd/pDd[ic]->GetBinError(1));
        printf("   denominator v2*Ddelta = %.5e +/- %.3e (rel %.3f%%)\n", den, sub_den, 100.*sub_den/den);
        printf("   rho(Dgamma112,den) = %+.3f   rho(Dgamma132,den) = %+.3f\n", rho1, rho3);
        printf("   kappa112 err: subsample %.4f  vs naive-uncorrelated %.4f\n", e112, naive1);
        printf("   kappa132 err: subsample %.4f  vs naive-uncorrelated %.4f\n", e132, naive3);
        printf("   rel err: Dg112 %.2f%%  Dg132 %.2f%%  k112 %.2f%%  k132 %.2f%%\n",
               100.*sub_g1/fabs(g1f), 100.*sub_g3/fabs(g3f), 100.*e112/fabs(k112), 100.*e132/fabs(k132));

        printf("\n   [analytic errors from TProfile variances + measured covariances]\n");
        printf("   kappa112  subsample %.4f   analytic %.4f   ratio %.3f\n", e112, a112, a112>0? e112/a112:0.0);
        printf("   kappa132  subsample %.4f   analytic %.4f   ratio %.3f\n", e132, a132, a132>0? e132/a132:0.0);
        printf("   rho(Dg112,v2)=%+.4f  rho(Dg112,Ddelta)=%+.4f\n", rA1V, rA1D);
        printf("   rho(Dg132,v2)=%+.4f  rho(Dg132,Ddelta)=%+.4f  rho(v2,Ddelta)=%+.4f\n", rA3V, rA3D, rVD);

        // ---------- is kappa112 != kappa132 significant? ----------
        // Both kappas are built from the SAME events and share the denominator v2*Ddelta,
        // so the difference must be errored with cov(Dg112,Dg132), not by quadrature.
        double cA1A3 = (Nev>1)? (pA1A3[ic]->GetBinContent(1) - g1f*g3f)/(Nev-1.0) : 0.0;
        double rA1A3 = (sA1*sA3>0)? cA1A3/(sA1*sA3) : 0.0;
        double dG    = g1f - g3f;
        double sdG2  = sA1*sA1 + sA3*sA3 - 2.*cA1A3;
        double sdG   = (sdG2>0)? sqrt(sdG2) : 0.0;
        double dK    = k112 - k132;
        // denominator uncertainty (1% here) adds in quadrature on the ratio; its
        // correlation with dG is negligible (see rho(Dgamma,den) above)
        double sdK   = (dG!=0 && den!=0)? fabs(dK)*sqrt(pow(sdG/dG,2)
                       + pow(sqrt(sV*sV/(v2f*v2f)+sD*sD/(ddf*ddf)),2)) : 0.0;
        printf("\n   [kappa112 vs kappa132 -- same events, shared denominator]\n");
        printf("   rho(Dg112,Dg132) = %+.4f   (event-by-event; ~0 means quadrature is fine)\n", rA1A3);
        printf("   Dg112-Dg132 = %+.4e +/- %.3e   (%.2f sigma)\n", dG, sdG, sdG>0? dG/sdG:0.0);
        printf("   k112 - k132 = %+.4f +/- %.4f   (%.2f sigma)\n", dK, sdK, sdK>0? dK/sdK:0.0);
        if(ORD1){
            double sA31  = pDg132_1[ic]->GetBinError(1);
            double cA1A31= (Nev>1)? (pA1A31[ic]->GetBinContent(1) - g1f*g31f)/(Nev-1.0) : 0.0;
            double dG1   = g1f - g31f;
            double sdG12 = sA1*sA1 + sA31*sA31 - 2.*cA1A31;
            double sdG1  = (sdG12>0)? sqrt(sdG12) : 0.0;
            double dK1   = k112 - k132_1;
            double sdK1  = (dG1!=0 && den!=0)? fabs(dK1)*sqrt(pow(sdG1/dG1,2)
                           + sV*sV/(v2f*v2f) + sD*sD/(ddf*ddf)) : 0.0;
            printf("\n   [identical filling: one cosine per pair for BOTH correlators]\n");
            printf("   err(Dg112) = %.3e   err(Dg132,1ord) = %.3e   ratio %.3f\n",
                   sA1, sA31, sA31>0? sA1/sA31 : 0.0);
            printf("   (for comparison, symmetrised Dg132 err = %.3e, ratio %.3f)\n",
                   sA3, sA3>0? sA1/sA3 : 0.0);
            printf("   Dg132,1ord - Dg132,sym = %+.4e  (same quantity, must agree)\n", g31f-g3f);
            printf("   k112 - k132,1ord = %+.4f +/- %.4f   (%.2f sigma)\n",
                   dK1, sdK1, sdK1>0? dK1/sdK1 : 0.0);
        }
        printf("\n");
    }

    fout.Write();
    printf("wrote %s\n", fname_out);
}
