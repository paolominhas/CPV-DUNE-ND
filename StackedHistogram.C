#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Vars/Vars.h" // Variables
#include "CAFAna/Cuts/TruthCuts.h" // Cuts
#include "StandardRecord/SRProxy.h" // A wrapper for the CAF format
#include "TCanvas.h" // Plots are drawn on a "canvas"
#include "TH1.h" // 1-dimensional histogram
#include "TPad.h" // Canvases are divided into pads. That could let you draw more than one plot on a canvas, if you wanted to, by using multiple pads. We will not bother with that today.
#include "TLegend.h" // Lets us draw a legend
#include "TMath.h" // I'll use some basic math functions
#include "TLatex.h"
#include <iostream>

const int MODE_QE = 1;
const int MODE_RES = 4;
const int MODE_DIS = 3;
const int MODE_MEC = 10;
const int PDG_MU=13;
const int PDG_E=11;
const int PDG_NUMU=14;
const int PDG_NUE=12;

using namespace ana;
void StackedHistogram()
{
    const std::string NDGAR_FHC = "/scratch/cpatrick/DuneSchool/CAFs/NDGAr/CAF_FHC_9**.root"; //ND-GAr FHC
    const std::string NDGAR_RHC = "/scratch/cpatrick/DuneSchool/CAFs/NDGAr/CAF_RHC_9**.root"; //ND-GAr RHC
    const std::string NDLAR_FHC = "/scratch/cpatrick/DuneSchool/CAFs/NDLAr/CAF_FHC_9**.root"; //ND-LAr FHC
    const std::string NDLAR_RHC = "/scratch/cpatrick/DuneSchool/CAFs/NDLAr/CAF_RHC_9**.root"; //ND-LAr RHC
    SpectrumLoader loader(NDLAR_FHC);
    const Binning binsTrue = Binning::Simple(80, 0, 10);

    const Var kTrueEnergy([](const caf::SRProxy* sr) { return sr->Ev; });

    // Define cuts for each interaction type
    const Cut kMEC([](const caf::SRProxy* sr) { return sr->mode == MODE_MEC; });
    const Cut kDIS([](const caf::SRProxy* sr) { return sr->mode == MODE_DIS; });
    const Cut kRES([](const caf::SRProxy* sr) { return sr->mode == MODE_RES; });
    const Cut kCCQE([](const caf::SRProxy* sr) { return sr->mode == MODE_QE; });

    // CC0pi condition: No pions, no other hadrons except protons and neutrons
    const Cut kCC0pi([](const caf::SRProxy* sr) {
        const int totalPions = sr->nipip + sr->nipim + sr->nipi0;
        const int otherHadrons = sr->nikp + sr->nikm + sr->nik0 + sr->niem + sr->nNucleus;
        return totalPions == 0 && otherHadrons == 0;
    });

    // Create histograms for MEC DIS RES CCQE (only events satisfying CC0pi)
    const HistAxis axTrue("True Neutrino Energy (GeV)", binsTrue, kTrueEnergy);

    Spectrum sMEC(loader, axTrue, kMEC && kCC0pi);
    Spectrum sDIS(loader, axTrue, kDIS && kCC0pi);
    Spectrum sRES(loader, axTrue, kRES && kCC0pi);
    Spectrum sCCQE(loader, axTrue, kCCQE && kCC0pi);

    loader.Go();  // Process the data

    // Convert Spectra to histograms
    TH1D* hMEC = sMEC.ToTH1(1e20);
    TH1D* hDIS = sDIS.ToTH1(1e20);
    TH1D* hRES = sRES.ToTH1(1e20);
    TH1D* hCCQE = sCCQE.ToTH1(1e20);

    // Assign colors
    hMEC->SetFillColor(kBlue);
    hDIS->SetFillColor(kRed);
    hRES->SetFillColor(kGreen);
    hCCQE->SetFillColor(kOrange);
    hMEC->SetLineColor(kBlack);
    hDIS->SetLineColor(kBlack);
    hRES->SetLineColor(kBlack);
    hCCQE->SetLineColor(kBlack);

    // Make the stacked historgram now
    THStack* hs = new THStack("hs", "Stacked CC0pi Interaction Types");

    hs->Add(hCCQE);
    hs->Add(hRES);
    hs->Add(hDIS);
    hs->Add(hMEC);

    // Details of the legend and canvas for plotting
    TCanvas* c1 = new TCanvas("c1", "Stacked Histogram", 1000, 750);
    hs->Draw("HIST");
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hMEC, "MEC", "f");
    legend->AddEntry(hDIS, "DIS", "f");
    legend->AddEntry(hRES, "RES", "f");
    legend->AddEntry(hCCQE, "CCQE", "f");
    legend->Draw();
    c1->SaveAs("StackedCC0pi.png");
}