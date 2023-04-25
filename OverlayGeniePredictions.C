#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TSpline.h>
#include <TProfile.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

void PrettyPlot(TH1D* h,int LineWidth = 2, int FontStyle = 132, int Ndivisions = 6, double TextSize = 0.06) {

	// ----------------------------------------------------------------------------------------------------------------

	h->SetLineWidth(LineWidth);

	// ----------------------------------------------------------------------------------------------------------------

	// X-axis

	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelSize(TextSize);
	h->GetXaxis()->SetTitleSize(TextSize);
	h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetNdivisions(Ndivisions);

	// ----------------------------------------------------------------------------------------------------------------

	// Y-axis

	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleSize(TextSize); 
	h->GetYaxis()->SetTickSize(0.02);
	h->GetYaxis()->SetLabelSize(TextSize);
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleOffset(1.05);
	h->GetYaxis()->SetNdivisions(Ndivisions);
	h->GetYaxis()->SetTitle("Cross Section [10^{-38}/cm^{2}]");

	return;	

}

void OverlayGeniePredictions() {

	// ----------------------------------------------------------------------------------------------------------------------------------------

	int LineWidth = 3;
	int FontStyle = 132;
	int Ndivisions = 6; 
	double TextSize = 0.06; 

	gStyle->SetOptStat(0);	

	// ----------------------------------------------------------------------------------------------------------------------------------------

	const std::vector<int> Colors{kAzure-4,kOrange+7,kGreen+1,kRed+1,kGreen+3,kBlue};

	// -----------------------------------------------------------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	vector<TString> FileNames; FileNames.clear();
	vector<TString> Label; Label.clear();
	vector<TH1D*> Plots; Plots.clear();

	// -----------------------------------------------------------------------------------------------------------------------------

	vector<TString> PlotName;
	PlotName.push_back("TrueDeltaPTPlot");
	PlotName.push_back("TrueDeltaAlphaTPlot");
	PlotName.push_back("TrueMuonCosThetaPlot");
	PlotName.push_back("TrueProtonCosThetaPlot");

	int nplots = PlotName.size();

	// -----------------------------------------------------------------------------------------------------------------------------

	FileNames.push_back("v3_2_0_G21_11b_00_000"); Label.push_back("v3.2.0 G21_11b_00_000");
	FileNames.push_back("v3_4_0_G21_11b_00_000"); Label.push_back("v3.4.0 G21_11b_00_000");
	FileNames.push_back("v3_2_0_G21_11b_00_000_1GeV"); Label.push_back("v3.2.0 G21_11b_00_000 1GeV");		 
	FileNames.push_back("SuSav2"); Label.push_back("feature branch (in paper)");
//	FileNames.push_back("Genie_v3_0_6_Nominal"); Label.push_back("G18");
//	FileNames.push_back("v3_4_0_G18_02a_02_11a"); Label.push_back("G18_02a_02_11a");		

	const int NFiles = FileNames.size();

	// -----------------------------------------------------------------------------------------------------------------------------

	// Loop over the plots

	for (int iplot = 0; iplot < nplots; iplot++) {

		Plots.clear();

		TString CanvasName = "OverlayCanvas"+TString(std::to_string(iplot));
		TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		can->SetBottomMargin(0.17);
		can->SetLeftMargin(0.15);

		// -----------------------------------------------------------------------------------------------------------------------------

		TLegend* leg = new TLegend(0.2,0.9,0.7,0.98);
		leg->SetNColumns(1);

		// -----------------------------------------------------------------------------------------------------------------------------

		for (int WhichFile = 0; WhichFile < NFiles; WhichFile++) {

			TFile* f = TFile::Open("OutputFiles/STVAnalysis_"+FileNames[WhichFile]+".root");
			TH1D* h = (TH1D*)(f->Get(PlotName.at(iplot)));
			Plots.push_back(h);
				
			Plots[WhichFile]->SetLineColor(Colors[WhichFile]);
			PrettyPlot(Plots[WhichFile]);
			Plots[WhichFile]->Draw("hist same");

			leg->AddEntry(Plots[WhichFile],Label[WhichFile],"l");

		}

		// -----------------------------------------------------------------------------------------------------------------------------

		leg->SetBorderSize(0);
		leg->SetTextSize(TextSize);
		leg->SetTextFont(FontStyle);
		leg->Draw();

	} // End of the loop over the plots

	// -----------------------------------------------------------------------------------------------------------------------------


} // End of the program
