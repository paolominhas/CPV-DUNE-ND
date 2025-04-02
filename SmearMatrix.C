#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/HistAxis.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Vars/Vars.h" // Variables
#include "CAFAna/Cuts/TruthCuts.h" // Cuts
#include "StandardRecord/SRProxy.h" // A wrapper for the CAF format
#include "TCanvas.h" // Plots are drawn on a "canvas"
#include "TH1.h" // 1-dimensional histogram
#include "TH2.h" // 2-dimensional histogram
#include "TPad.h" // Canvases are divided into pads.
#include "TLegend.h" // Lets us draw a legend
#include "TMath.h" // I'll use some basic math functions
#include "TLatex.h"
#include "TColor.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include <TH2D.h>

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
const int PDG_NUE=12;

using namespace ana;
using util::sqr; // Square

double QEFormula(double Emu, double cosmu) // Muon energy and cosine of muon angle
 {
   const double pmu = sqrt(sqr(Emu) - sqr(M_MU));
   const double num = sqr(M_P) - sqr(M_N - E_B) - sqr(M_MU) + 2 * (M_N - E_B) * Emu;
   const double denom = 2 * (M_N - E_B - Emu + pmu * cosmu);
  return num/denom;
}

void SmearMatrix()
{
  const std::string NDGAR_FHC = "/scratch/cpatrick/DuneSchool/CAFs/NDGAr/CAF_FHC_9**.root"; //ND-GAr FHC
  const std::string NDGAR_RHC = "/scratch/cpatrick/DuneSchool/CAFs/NDGAr/CAF_RHC_9**.root"; //ND-GAr RHC
  const std::string NDLAR_FHC = "/scratch/cpatrick/DuneSchool/CAFs/NDLAr/CAF_FHC_9**.root"; //ND-LAr FHC
  const std::string NDLAR_RHC = "/scratch/cpatrick/DuneSchool/CAFs/NDLAr/CAF_RHC_9**.root"; //ND-LAr RHC
  SpectrumLoader loader(NDGAR_FHC);
  const Binning binsEnergy = Binning::Simple(80, 0, 10);
  const Binning binsTrue = Binning::Simple(80, 0, 10);
  const Binning binsReco = Binning::Simple(80, 0, 10);

  const Var kConservedETrue([](const caf::SRProxy* sr)
                                { 
                            const double Emu = sr->LepE;
                            const double protonKE=sr->eP; 
                            return Emu + (protonKE + M_P) - (M_N - E_B);
  });

  const Var kConservedEReco([](const caf::SRProxy* sr)
                                {
                            const double Emu = sr->Elep_reco; // Final-state reconstructed lepton (muon) energy
                            const double protonKE=sr->eRecoP; // kinetic energy
                            return Emu + (protonKE + M_P) - (M_N - E_B);
  });

  const Var kRecoE([](const caf::SRProxy* sr)
                                {
                            if(std::isnan(sr->Ev_reco)) return 0.;
                            return sr->Ev_reco;
  });

  const Var kQEFormulaEnergy([](const caf::SRProxy* sr)
                                {
                            const double Emu = sr->Elep_reco;
                            const double cosmu = cos(sr->theta_reco);
                            if(Emu == 0) return 0.;
                            return QEFormula(Emu, cosmu);
  });

  const Cut kHasQEFinalState([](const caf::SRProxy* sr)
   {
      if(std::isnan(sr->Ev_reco)) return false;
      if(sr->Ev_reco<=0.) return false;
      if(sr->Elep_reco<=0.) return false;
      const int totOthers = sr->nN + sr->nipip + sr->nipim + sr->nipi0 + sr->nikp + sr->nikm + sr->nik0 + sr->niem + sr->nNucleus;
      return sr->LepPDG == PDG_MU && sr->nP == 1 && totOthers == 0 ;
   });

   const Var kTrueEnergy([](const caf::SRProxy* sr) { return sr->Ev; }); // True neutrino energy
   const Var kRecoEnergy([](const caf::SRProxy* sr) { return sr->Ev_reco; }); // Reco neutrino energy

   const HistAxis axTrue("True Energy (GeV)", binsTrue, kTrueEnergy);
   const HistAxis axReco("Reco Energy (GeV)", binsReco, kRecoEnergy);

   const Cut kValidReco([](const caf::SRProxy* sr) { return sr->Ev_reco > 0; });

   Spectrum sSmearing(loader, axTrue, axReco, kValidReco);

   loader.Go();  // Process the data

   TH2* hSmearing = sSmearing.ToTH2(1e20); // Normalize to POT
   TCanvas *c2 = new TCanvas("c2", "Smearing Matrix", 1000, 750);
   hSmearing->SetTitle("Smearing Matrix: True vs. Reconstructed Neutrino Energy - NDGAr FHC");
   gStyle->SetPalette(kWaterMelon); 
   hSmearing->Draw("COLZ"); // Plot as color map
   c2->SaveAs("SmearingMatrixGF.png"); 

   //TH2D *hTrueReco = sTrueReco.ToTH2(1e20); // Convert Spectrum2D to TH2D histogram
}
