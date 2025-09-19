#include "TFile.h"
#include "TTree.h"
#include <cmath>
#include <iostream>

float deltaR(float eta1, float phi1, float eta2, float phi2) {
    float dphi = std::fabs(phi1 - phi2);
    if (dphi > M_PI) dphi = 2 * M_PI - dphi;
    float deta = eta1 - eta2;
    return std::sqrt(deta * deta + dphi * dphi);
}

void select_signal(std::string path, const int genmass, const int genpdg) {
    
    // Open input file and tree
    TFile *infile = TFile::Open(path.c_str());
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error: cannot open input file!" << std::endl;
        return;
    }

    TTree *tree = (TTree*)infile->Get("Events");
    if (!tree) {
        std::cerr << "Error: cannot find TTree 'Events' in file!" << std::endl;
        return;
    }
    tree->SetMakeClass(0);

    // -----------------------------
    // Declare branches
    // -----------------------------

    // GenJet
    static const int MAXGENJET = 18;
    float genjet_pt[MAXGENJET], genjet_eta[MAXGENJET], genjet_phi[MAXGENJET], genjet_mass[MAXGENJET];
    int nGenJet = 0;

    // L1 jets (SC8)
    static const int MAXL1JET = 18;
    float l1jet_pt[MAXL1JET], l1jet_eta[MAXL1JET], l1jet_phi[MAXL1JET], l1jet_mass[MAXL1JET];
    int nL1Jet = 0;

    // GenPart
    static const int MAXGENPART = 10000;
    float GenPart_mass[MAXGENPART];
    int   GenPart_pdgId[MAXGENPART];
    float GenPart_pt[MAXGENPART];
    int   nGenPart = 0;

    // Set branch addresses
    tree->SetBranchAddress("GenJetAK8_pt",   genjet_pt);
    tree->SetBranchAddress("GenJetAK8_eta",  genjet_eta);
    tree->SetBranchAddress("GenJetAK8_phi",  genjet_phi);
    tree->SetBranchAddress("GenJetAK8_mass", genjet_mass);
    tree->SetBranchAddress("nGenJetAK8",     &nGenJet);

    tree->SetBranchAddress("L1puppiJetSC8_pt",   l1jet_pt);
    tree->SetBranchAddress("L1puppiJetSC8_eta",  l1jet_eta);
    tree->SetBranchAddress("L1puppiJetSC8_phi",  l1jet_phi);
    tree->SetBranchAddress("L1puppiJetSC8_mass", l1jet_mass);
    tree->SetBranchAddress("nL1puppiJetSC8",     &nL1Jet);

    tree->SetBranchAddress("GenPart_mass", GenPart_mass);
    tree->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
    tree->SetBranchAddress("GenPart_pt", GenPart_pt);
    tree->SetBranchAddress("nGenPart", &nGenPart);

    // -----------------------------
    // Output file + new tree
    // -----------------------------
    std::string filename = path.substr(path.find_last_of("/\\") + 1);
    filename = filename.substr(0, filename.find(".root"));
    TString outname = Form("./signal_samples/%s_signal.root", filename.c_str());
    TFile *outfile = new TFile(outname, "RECREATE");
    TTree *newtree = tree->CloneTree(0);

    // New branches for matched gen jet info
    float matchedGen_pt[MAXL1JET], matchedGen_mass[MAXL1JET], matchedGen_dr[MAXL1JET];
    newtree->Branch("L1puppiJetSC8_genpt",   matchedGen_pt,   "L1puppiJetSC8_genpt[nL1puppiJetSC8]/F");
    newtree->Branch("L1puppiJetSC8_gendr",   matchedGen_dr,   "L1puppiJetSC8_gendr[nL1puppiJetSC8]/F");
    newtree->Branch("L1puppiJetSC8_genmass", matchedGen_mass, "L1puppiJetSC8_genmass[nL1puppiJetSC8]/F");


    // -----------------------------
    // Event loop
    // -----------------------------
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        bool store = false;

        // loop over GenPart
        for (int j = 0; j < nGenPart; j++) {
            if (std::abs(GenPart_pdgId[j]) == genpdg && GenPart_pt[j] > 1.25 * GenPart_mass[j]) {
                store = true;
                break;
            }
        }

        // initialize matches
        for (int j = 0; j < nL1Jet; j++) {
            matchedGen_pt[j] = -1;
            // matchedGen_eta[j] = -99;
            // matchedGen_phi[j] = -99;
            matchedGen_dr[j] = -1;
            matchedGen_mass[j] = -1;

            float bestDR = 999.0;
            float bestPtDiff = 1e9;
            int bestIdx = -1;

            for (int k = 0; k < nGenJet; k++) {
                float dR = deltaR(l1jet_eta[j], l1jet_phi[j], genjet_eta[k], genjet_phi[k]);
                if (dR < 0.35) {
                    float dPt = std::fabs(l1jet_pt[j] - genjet_pt[k]);
                    if (dR < bestDR || (std::fabs(dR - bestDR) < 1e-6 && dPt < bestPtDiff)) {
                        bestDR = dR;
                        bestPtDiff = dPt;
                        bestIdx = k;
                    }
                }
            }

            if (bestIdx >= 0) {
                matchedGen_pt[j]   = genjet_pt[bestIdx];
                // matchedGen_eta[j]  = genjet_eta[bestIdx];
                // matchedGen_phi[j]  = genjet_phi[bestIdx];
                matchedGen_dr[j]   = bestDR;
                matchedGen_mass[j] = genjet_mass[bestIdx];
            }
        }

        // check jet + gen condition
        if (nGenJet > 0 && genjet_mass[0] > 0.9 * genmass && store) {
            newtree->Fill();
        }
    }

    // -----------------------------
    // Write output
    // -----------------------------
    newtree->Write();
    outfile->Close();
    infile->Close();

    std::cout << "Done! Filtered tree written to " << outname << std::endl;
}


// #include "TFile.h"
// #include "TTree.h"
// #include <cmath>
// #include <iostream>

// void select_signal(std::string path, const int genmass, const int genpdg) {
    
//     // Open input file and tree
//     TFile *infile = TFile::Open(path.c_str());
//     if (!infile || infile->IsZombie()) {
//         std::cerr << "Error: cannot open input file!" << std::endl;
//         return;
//     }

//     TTree *tree = (TTree*)infile->Get("Events");
//     if (!tree) {
//         std::cerr << "Error: cannot find TTree 'Events' in file!" << std::endl;
//         return;
//     }
//     tree->SetMakeClass(0);

//     // -----------------------------
//     // Declare branches
//     // -----------------------------

//     // GenJet
//     static const int MAXGENJET = 18;
//     float genjet_mass[MAXGENJET];
//     int nGenJet = 0;

//     // GenPart
//     static const int MAXGENPART = 10000;
//     float GenPart_mass[MAXGENPART];
//     int   GenPart_pdgId[MAXGENPART];
//     float GenPart_pt[MAXGENPART];
//     int   nGenPart = 0;

//     // Set branch addresses
//     tree->SetBranchAddress("GenJetAK8_mass", genjet_mass);
//     tree->SetBranchAddress("nGenJetAK8", &nGenJet);

//     tree->SetBranchAddress("GenPart_mass", GenPart_mass);
//     tree->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
//     tree->SetBranchAddress("GenPart_pt", GenPart_pt);
//     tree->SetBranchAddress("nGenPart", &nGenPart);

//     // -----------------------------
//     // Output file + new tree
//     // -----------------------------
//     std::string filename = path.substr(path.find_last_of("/\\") + 1);
//     filename = filename.substr(0, filename.find(".root"));
//     TString outname = Form("./signal_samples/%s_signal.root", filename.c_str());
//     TFile *outfile = new TFile(outname, "RECREATE");
//     TTree *newtree = tree->CloneTree(0);

//     // -----------------------------
//     // Event loop
//     // -----------------------------
//     Long64_t nentries = tree->GetEntries();
//     for (Long64_t i = 0; i < nentries; i++) {
//         tree->GetEntry(i);

//         bool store = false;

//         // loop over GenPart
//         for (int j = 0; j < nGenPart; j++) {
//             if (std::abs(GenPart_pdgId[j]) == genpdg && GenPart_pt[j] > 1.25 * GenPart_mass[j]) {
//                 store = true;
//                 break;
//             }
//         }

//         // check jet + gen condition
//         if (nGenJet > 0 && genjet_mass[0] > 0.9 * genmass && store) {
//             newtree->Fill();
//         }
//     }

//     // -----------------------------
//     // Write output
//     // -----------------------------
//     newtree->Write();
//     outfile->Close();
//     infile->Close();

//     std::cout << "Done! Filtered tree written to output.root" << std::endl;
// }

// // int main() {
// //     std::string filepath = "/eos/home-l/lroberts/Phase2-L1MenuTools/CMSSW_15_1_0_pre5/src/condor/jobs/XtoHH_MX-140to780_MH-20to90_PU200_NEW_1756385838/data/Phase2Spring24_XtoHH_MX-140to780_MH-20to90_PU200.root";
// //     const int genmass = 125; // Example mass
// //     const int genpdg = 25;   // Example PDG ID for Higgs
    
// //     select_signal(filepath, genmass, genpdg);
    
// //     return 0;
// // }