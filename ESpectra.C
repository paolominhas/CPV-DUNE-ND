#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Vars/Vars.h"
#include "CAFAna/Cuts/TruthCuts.h" // Truth Cuts
#include "StandardRecord/SRProxy.h" // A wrapper for the CAF format
#include "TCanvas.h" 
#include "TH1.h" 
#include "TPad.h" 
#include "TLegend.h"
#include "TLatex.h"
#include <iostream>

const double M_P = .938; // Proton mass in GeV
const double M_N = .939; // Neutron mass in GeV
const double M_MU = .106; // Muon mass in GeV
const double E_B = .028; // Binding energy for nucleons in argon-40 in GeV
const int MODE_QE = 1;
const int MODE_RES = 4;
const int MODE_DIS = 3;
const int MODE_MEC = 10;
const int PDG_MU=13;
const int PDG_E=11;
const int PDG_NUMU=14;

using namespace ana;
using util::sqr; // Square

void ESpectra()
{
  const std::string NDGAR_FHC = "/scratch/cpatrick/DuneSchool/CAFs/NDGAr/CAF_FHC_9**.root"; //ND-GAr FHC
  const std::string NDGAR_RHC = "/scratch/cpatrick/DuneSchool/CAFs/NDGAr/CAF_RHC_9**.root"; //ND-GAr RHC
  const std::string NDLAR_FHC = "/scratch/cpatrick/DuneSchool/CAFs/NDLAr/CAF_FHC_9**.root"; //ND-LAr FHC
  const std::string NDLAR_RHC = "/scratch/cpatrick/DuneSchool/CAFs/NDLAr/CAF_RHC_9**.root"; //ND-LAr RHC
  SpectrumLoader loader(NDGAR_FHC);
  
  const Binning binsEnergy = Binning::Simple(100, 0, 10);
  const HistAxis axTrue("True neutrino energy (GeV)", binsEnergy, kTrueEnergy); // True Energy
  Spectrum sTrueENumu(loader, axTrue, kIsNumuCC && !kIsAntiNu); // Muon Neutrino CC
  Spectrum sTrueENumubar(loader, axTrue, kIsNumuCC && kIsAntiNu);  // Muon antineutrino CC interactions
  Spectrum sTrueENue(loader, axTrue, kIsBeamNue && !kIsAntiNu);    // Electron neutrino CC interactions
  Spectrum sTrueENuebar(loader, axTrue, kIsBeamNue && kIsAntiNu);  // Electron antineutrino CC interactions

  loader.Go();
  const double pot = 1e20; // Protons on target - this is a scaling factor to make the plot easier to read.
  TCanvas *canvas = new TCanvas; // Make a canvas
  TH1D *hTrueENumu = sTrueENumu.ToTH1(pot, kBlue);// Draw our spectrum in blue. ROOT colors are defined at https://root.cern.ch/doc/master/classTColor.html
  hTrueENumu->Draw("HIST"); // This time we turn our spectrum into a ROOT histogram, and draw that. It means we can use the histogram for other things - like a legend.
  TH1D *hTrueENumubar=sTrueENumubar.ToTH1(pot, kBlue, 7);//Antineutrinos are getting a dashed line.
  hTrueENumubar->Draw("HIST SAME"); // SAME canvas as the previous spectrum  
  TH1D *hTrueENue = sTrueENue.ToTH1(pot, kRed); // Red spectrum
  hTrueENue->Draw("HIST SAME"); // Drawn onto same canvas
  TH1D *hTrueENuebar = sTrueENuebar.ToTH1(pot, kRed, 7); // Red spectrum
  hTrueENuebar->Draw("HIST SAME"); // Drawn onto same canvas
  
  gPad->SetLogy();
  TLatex title;
  title.SetTextSize(0.04);
  title.DrawLatexNDC(0.2, 0.93, "Neutrino Beam Contents at the NDGAr FHC");
  auto legend = new TLegend(0.7,0.7,0.9,0.9); // x and y coordinates of corners
  legend->SetHeader("Legend","C"); // option "C" to center the header
  legend->AddEntry(hTrueENumu,"#nu_{#mu}","l");
  legend->AddEntry(hTrueENumubar,"#bar{#nu}_{#mu}","l");
  legend->AddEntry(hTrueENue,"#nu_{e}","l");
  legend->AddEntry(hTrueENuebar,"#bar{#nu}_{e}","l");
  legend->Draw();
  canvas->SaveAs("NDGArFHCSpectrum.png"); // Save the result
}
