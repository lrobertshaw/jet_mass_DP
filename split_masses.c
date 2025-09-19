#include "TFile.h"
#include "TTree.h"
#include <cmath>
#include <iostream>
#include <vector>

void split_masses(const std::string path){
    // Input file
    TFile *infile = TFile::Open(path.c_str());
    TTree *tree = (TTree*)infile->Get("Events");

    // GenPart branches
    const int MAXGEN = 5000; // adjust if your nanoAOD has more
    int     nGenPart;
    float   GenPart_mass[MAXGEN];
    int     GenPart_pdgId[MAXGEN];
    Short_t GenPart_genPartIdxMother[MAXGEN];   // <-- fixed
    UShort_t GenPart_statusFlags[MAXGEN];       // <-- fixed


    tree->SetBranchAddress("nGenPart", &nGenPart);
    tree->SetBranchAddress("GenPart_mass", GenPart_mass);
    tree->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
    tree->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother);
    tree->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags);

    // Loop over Higgs mass windows
    for (int m = 5; m <= 5; m++) {
        float lower = (m*10) - 5;
        float mass_point = (m*10) + 0;
        float upper = (m*10) + 5;

        TString outname = Form("./light_higgs_samples/XtoHH_MX-350_MH-%d_qq_PU200.root", (int)mass_point);
        TFile *outfile = new TFile(outname, "RECREATE");
        TTree *newtree = tree->CloneTree(0); // copy structure, no entries

        Long64_t nentries = tree->GetEntries();
        for (Long64_t i = 0; i < nentries; i++) {
            tree->GetEntry(i);

            // 1) Find last-copy Higgs bosons
            std::vector<int> higgs_indices;
            for (int j = 0; j < nGenPart; j++) {
                if (std::abs(GenPart_pdgId[j]) == 25) {
                    // last copy check: statusFlags bit 13 (1<<13 == 8192)
                    if (GenPart_statusFlags[j] & (1 << 13)) {
                        higgs_indices.push_back(j);
                    }
                }
            }

            if (higgs_indices.size() < 2) continue; // need at least 2 Higgs

            // 2) Check decay products of each Higgs
            bool allHiggsToQuarksOrGluons = true;
            for (int higgsIdx : higgs_indices) {
                bool thisHiggsAllQuarksOrGluons = true;
                for (int k = 0; k < nGenPart; k++) {
                    if (GenPart_genPartIdxMother[k] == higgsIdx) {
                        int absId = std::abs(GenPart_pdgId[k]);
                        if (!((absId >= 1 && absId <= 8) || absId == 21)) {
                            thisHiggsAllQuarksOrGluons = false;
                            break;
                        }
                    }
                }
                if (!thisHiggsAllQuarksOrGluons) {
                    allHiggsToQuarksOrGluons = false;
                    break;
                }
            }

            if (!allHiggsToQuarksOrGluons) continue; // skip if any Higgs decays to something else


            // 3) Check Higgs mass window requirement
            bool pass_mass = false;
            for (int j = 0; j < nGenPart; j++) {
                if (std::abs(GenPart_pdgId[j]) == 25) {
                    float mH = GenPart_mass[j];
                    if (mH > lower && mH < upper) {
                        pass_mass = true;
                        break;
                    }
                }
            }

            if (pass_mass) newtree->Fill();
        }

        newtree->Write();
        outfile->Close();
        std::cout << "Wrote " << outname << std::endl;
    }

    infile->Close();
}
