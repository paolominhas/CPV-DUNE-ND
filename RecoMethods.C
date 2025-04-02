#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Vars/Vars.h" // Variables
#include "CAFAna/Cuts/TruthCuts.h" // Cuts
#include "StandardRecord/SRProxy.h" // A wrapper for the CAF format
#include "TCanvas.h" // Plots are drawn on a "canvas"
#include "TH1.h" // 1-dimensional histogram
#include "TPad.h" // Canvases are divided into pads.
#include "TLegend.h" // Lets us draw a legend
#include "TMath.h" // I'll use some basic math functions
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

void RecoMethods()
{
  const std::string NDGAR_FHC = "/scratch/cpatrick/DuneSchool/CAFs/NDGAr/CAF_FHC_90*.root"; //ND-GAr FHC
  const std::string NDGAR_RHC = "/scratch/cpatrick/DuneSchool/CAFs/NDGAr/CAF_RHC_90*.root"; //ND-GAr RHC
  const std::string NDLAR_FHC = "/scratch/cpatrick/DuneSchool/CAFs/NDLAr/CAF_FHC_90*.root"; //ND-LAr FHC
  const std::string NDLAR_RHC = "/scratch/cpatrick/DuneSchool/CAFs/NDLAr/CAF_RHC_90*.root"; //ND-LAr RHC
  SpectrumLoader loader(NDLAR_RHC);
  const Binning binsEnergy = Binning::Simple(40, 0, 10);

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


  const HistAxis axConservedETrue("E_#nu (conserve true energies) (GeV)", binsEnergy, kConservedETrue);
  const HistAxis axConservedEReco("E_#nu (conserve reco energies) (GeV)", binsEnergy, kConservedEReco);
  const HistAxis axEQE("E_#nu (QE formula) (GeV)", binsEnergy, kQEFormulaEnergy);
  const HistAxis axEReco("E_#nu reco (GeV)", binsEnergy, kRecoE);
  const HistAxis axETrue("True neutrino energy (GeV)", binsEnergy, kTrueEnergy);

  const Cut kHasQEFinalState([](const caf::SRProxy* sr)
   {
      if(std::isnan(sr->Ev_reco)) return false;
      if(sr->Ev_reco<=0.) return false;
      if(sr->Elep_reco<=0.) return false;
      const int totOthers = sr->nN + sr->nipip + sr->nipim + sr->nipi0 + sr->nikp + sr->nikm + sr->nik0 + sr->niem + sr->nNucleus;
      return sr->LepPDG == PDG_MU && sr->nP == 1 && totOthers == 0 ;
   });

   Spectrum sConservedETrue (loader, axConservedETrue, kHasQEFinalState);
   Spectrum sConservedEReco (loader, axConservedEReco, kHasQEFinalState);
   Spectrum sEQE (loader, axEQE, kHasQEFinalState);
   Spectrum sEReco (loader, axEReco, kHasQEFinalState);
   Spectrum sETrue (loader, axETrue, kHasQEFinalState);

   loader.Go();
   const double pot = 1e20;
   TCanvas *canvas = new TCanvas;

   TH1D *hConservedETrue = sConservedETrue.ToTH1(pot, kRed-10);
   hConservedETrue->SetFillColorAlpha(kRed-10, 0.4);  // 40% opacity fill
   hConservedETrue->SetLineColor(kRed-10);
   hConservedETrue->SetLineWidth(3);

   TH1D *hETrue = sETrue.ToTH1(pot, kBlack);
   //hETrue->SetFillColorAlpha(kRed+2, 0.4);
   hETrue->SetLineColor(kBlack);
   hETrue->SetLineWidth(3);

  TH1D *hConservedEReco = sConservedEReco.ToTH1(pot, kRed+2);
   hConservedEReco->SetLineColor(kRed+2);
   hConservedEReco->SetLineStyle(3);  // Dashed
   hConservedEReco->SetLineWidth(3);

   TH1D *hEQE = sEQE.ToTH1(pot, kBlue);
   hEQE->SetLineColor(kBlue);
   hEQE->SetLineStyle(3);  // Dashed
   hEQE->SetLineWidth(3);

   TH1D *hEReco = sEReco.ToTH1(pot, kGreen-4);
   hEReco->SetLineColor(kGreen-4);
   hEReco->SetLineStyle(3);  // Dashed
   hEReco->SetLineWidth(3);

   double height= TMath::Max(hConservedETrue->GetMaximum(),hConservedEReco->GetMaximum());
   height=max(height,hEQE->GetMaximum());
   height=max(height,hEReco->GetMaximum());
   height=max(height,hETrue->GetMaximum());

   hConservedEReco->SetLineStyle(2);
   hEQE->SetLineStyle(2);
   hEReco->SetLineStyle(2);

   hConservedETrue->GetYaxis()->SetRangeUser(0,height * 1.1); // set the y axis range to 1.1 times the height
   hConservedETrue->GetXaxis()->SetTitle("Energy calculated various ways (GeV)");

   hConservedETrue->Draw("HIST");
   hETrue->Draw("HIST SAME");
   hConservedEReco->Draw("HIST SAME");
   hEQE->Draw("HIST SAME");
   hEReco->Draw("HIST SAME");

   gPad->SetLogy(false);

   hConservedETrue->GetXaxis()->SetTitleSize(0.04);
   hConservedETrue->GetXaxis()->SetLabelSize(0.03);
   hConservedETrue->GetYaxis()->SetTitleSize(0.04);
   hConservedETrue->GetYaxis()->SetLabelSize(0.03);

   auto legend = new TLegend(0.65,0.65,0.9,0.9); // x and y coordinates of corners
   TLatex title;
   title.SetTextSize(0.04);
   title.DrawLatexNDC(0.2, 0.93, "Neutrino Energy Reconstruction Methods NDLAr RHC");
   legend->SetFillStyle(0);  // Transparent background
   legend->SetBorderSize(0);
   legend->SetTextSize(0.035);
   legend->AddEntry(hConservedETrue, "True Energy Conservation", "f");  // Filled
   legend->AddEntry(hETrue, "True Energy", "f");
   legend->AddEntry(hConservedEReco, "Reco Energy Conservation", "l");
   legend->AddEntry(hEQE, "QE Formula", "l");
   legend->AddEntry(hEReco, "Reconstructed from CAF", "l");
   legend->Draw();
   canvas->SaveAs("RecoMethodLR.png");

}
