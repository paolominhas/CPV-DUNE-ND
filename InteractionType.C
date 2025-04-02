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
void InteractionType()
{
  const std::string NDGAR_FHC = "/scratch/cpatrick/DuneSchool/CAFs/NDGAr/CAF_FHC_9**.root"; //ND-GAr FHC
  const std::string NDGAR_RHC = "/scratch/cpatrick/DuneSchool/CAFs/NDGAr/CAF_RHC_9**.root"; //ND-GAr RHC
  const std::string NDLAR_FHC = "/scratch/cpatrick/DuneSchool/CAFs/NDLAr/CAF_FHC_9**.root"; //ND-LAr FHC
  const std::string NDLAR_RHC = "/scratch/cpatrick/DuneSchool/CAFs/NDLAr/CAF_RHC_9**.root"; //ND-LAr RHC
  SpectrumLoader loader(NDLAR_FHC); // ***** Change this to use a different sample ***
  const Binning binsEnergy = Binning::Simple(100, 0, 10);
  const HistAxis axTrue("True neutrino energy (GeV)", binsEnergy, kTrueEnergy);

  const Cut kIsQE = SIMPLEVAR(mode) == MODE_QE; //The modes are defined at the top
  const Cut kIsMEC = SIMPLEVAR(mode) == MODE_MEC;
  const Cut kIsDIS = SIMPLEVAR(mode) == MODE_DIS;
  const Cut kIsRES = SIMPLEVAR(mode) == MODE_RES;
  const Cut kIsCCQE = kIsQE && kIsNumuCC && !kIsAntiNu;
  
  Spectrum sTrueEQE(loader, axTrue, kIsCCQE);
  Spectrum sTrueMEC(loader, axTrue, kIsMEC);
  Spectrum sTrueDIS(loader, axTrue, kIsDIS);
  Spectrum sTrueRES(loader, axTrue, kIsRES); 
  
  /*
   The CCQE final state is 1 proton and 1 muon.
   Pass: 1 proton and 1 muon, no other particles
   The input for this function is a CAF "Standard record"
   */
   const Cut kHasQEFinalState([](const caf::SRProxy* sr)
                             {
                               const int totOthers = sr->nN + sr->nipip + sr->nipim + sr->nipi0 + sr->nikp + sr->nikm + sr->nik0 + sr->niem + sr->nNucleus;     
                               return sr->LepPDG == PDG_MU && sr->nP == 1 && totOthers == 0;
                             });
  Spectrum sTrueEQEfs(loader, axTrue, kHasQEFinalState);

  const Cut kHasCC0PiFinalState([](const caf::SRProxy* sr)
                                {
                                  const int totPi = sr->nipip + sr->nipim + sr->nipi0;
                                  return sr->LepPDG == PDG_MU && sr->nP >= 1 && totPi == 0;
                                });
  Spectrum sTrueE0pifs(loader, axTrue, kHasCC0PiFinalState);

  loader.Go();
  const double pot = 1e20;
  TCanvas *canvas = new TCanvas; // Make a canvas
  
  TH1D *hTrueEQE = sTrueEQE.ToTH1(pot, 9);
  hTrueEQE->SetLineColor(kBlack);

  TH1D *hTrueEQEfs =sTrueEQEfs.ToTH1(pot, kRed-3);
  TH1D *hTrueE0pifs =sTrueE0pifs.ToTH1(pot, 3);
  
  // Set Axis Scale:
  double height= TMath::Max(hTrueEQE->GetMaximum(),TMath::Max(hTrueEQEfs->GetMaximum(),hTrueE0pifs->GetMaximum())); // height of the tallest histogram
  hTrueEQE->GetYaxis()->SetRangeUser(0,height * 1.1); // set the y axis range to 1.1 times the height
  // With error bars:
  hTrueEQE->Draw("HIST E");
  hTrueEQEfs->Draw("HIST E SAME");
  hTrueE0pifs->Draw("HIST E SAME");
  
  gPad->SetLogy(false);
  TLatex title;
  title.SetTextSize(0.04);
  title.DrawLatexNDC(0.2, 0.93, "Different Ways of Finding QE Interactions at NDLAr FHC"); 
  auto legend = new TLegend(0.65,0.65,0.9,0.9); // x and y coordinates of corners
  legend->SetHeader("Legend","C"); // option "C" to center the header
  legend->AddEntry(hTrueE0pifs,"CC0#pi","l");
  legend->AddEntry(hTrueEQE,"True CCQE","l");
  legend->AddEntry(hTrueEQEfs,"1 #mu^{-}, 1 p","l");
  legend->Draw();
  canvas->SaveAs("Task2LArFHCError.png");
 
  /*
  THStack * hs1 = new THStack("pot"," stacked");
  THStack * hs2 = new THStack("pot"," unstacked");

  THStack *hTrueEQE = sTrueEQE.ToTH1(pot, 9);


  
  TCanvas * c1 = new TCanvas("c1","stacked hists",700,900);
  c1->Divide(1, 2);
  c1->cd(1);
  hs1->Draw("");
  c1->cd(2);
  hs2->Draw("nostack");
 */
}
