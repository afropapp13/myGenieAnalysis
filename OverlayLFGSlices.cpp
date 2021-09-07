#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPaletteAxis.h>

#include <iostream>
#include <vector>

#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

// -------------------------------------------------------------------------------------------------------------------------------------------------

void PrettyPlot(TH1D* h, int Color) {

	int font = 132;
	double size = 0.05;

	h->GetXaxis()->SetRangeUser(0.01,0.79);

	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetTitleFont(font);
	h->GetXaxis()->SetTitleSize(size);	
	h->GetXaxis()->SetLabelFont(font);
	h->GetXaxis()->SetLabelSize(size);
	h->GetXaxis()->SetTitleOffset(1.1);	
	
	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleFont(font);
	h->GetYaxis()->SetTitleSize(size);	
	h->GetYaxis()->SetLabelFont(font);
	h->GetYaxis()->SetLabelSize(size);

	h->SetLineColor(Color);
	h->SetMarkerColor(Color);
	h->SetMarkerSize(1.5);
	h->SetMarkerStyle(20);

	h->Scale(1./h->GetMaximum());
	h->GetYaxis()->SetRangeUser(0.01,1.1);	

}

// -------------------------------------------------------------------------------------------------------------------------------------------------

void OverlayLFGSlices() {

	// -----------------------------------------------------------------------------

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();	

	// -------------------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetTitleSize(0.0,"t");
	gStyle->SetTitleFont(132,"t");

	gStyle->SetTitleOffset(0.9,"xyz");
	gStyle->SetTitleSize(0.05,"xyz");
	gStyle->SetTitleFont(132,"xyz");
	gStyle->SetLabelSize(0.05,"xyz");
	gStyle->SetLabelFont(132,"xyz");
	gStyle->SetNdivisions(8,"xyz");

	gStyle->SetOptStat(0);
	gStyle->SetPadRightMargin(0.13);

	gROOT->ForceStyle();

	// -----------------------------------------------------------------------------------------------------------------------------------------

	TFile* OverlayFile = TFile::Open("OutputFiles/STVAnalysis_Genie_v3_0_6_Out_Of_The_Box.root","readonly");

	// -----------------------------------------------------------------------------------------------------------------------------------------

	TString CanvasNameAll = "LFG_Pn_AllSlices";
	TCanvas* canAll = new TCanvas(CanvasNameAll,CanvasNameAll,205,34,1024,768);

	for (int i = 0; i < NRanges+1; i++) {

		TString CanvasName = "LFG_Pn_Slice"+TString(std::to_string(i));
		TCanvas* can = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		can->SetBottomMargin(0.12);

		TH1D* LFG = (TH1D*)(OverlayFile->Get("TrueLFGPnSlicePlot_"+TString(std::to_string(i))));
		PrettyPlot(LFG,kBlue+2);
		//LFG->SetTitle(RangeLabel[i]);
		LFG->Draw("p hist same");

		TLatex *text = new TLatex();
		text->SetTextFont(FontStyle);
		text->SetTextSize(0.05);
		if (i == 4) { text->DrawLatexNDC(0.65,0.85,RangeLabel[i]); }
		else { text->DrawLatexNDC(0.5,0.5,RangeLabel[i]); }			

		TH1D* kMiss = (TH1D*)(OverlayFile->Get("TruekMissSlicePlot_"+TString(std::to_string(i))));
		PrettyPlot(kMiss,kOrange+7);
		kMiss->Draw("p0 hist same");

		TH1D* kMissNoFSI = (TH1D*)(OverlayFile->Get("TruekMissSliceNoFSIPlot_"+TString(std::to_string(i))));
		PrettyPlot(kMissNoFSI,kRed+2);
		kMissNoFSI->Draw("p0 hist same");			

		TH1D* Pn = (TH1D*)(OverlayFile->Get("TruePnProxySlicePlot_"+TString(std::to_string(i))));
		PrettyPlot(Pn,kGreen+2);
		Pn->Draw("p0 hist same");	

		TH1D* PnNoFSI = (TH1D*)(OverlayFile->Get("TruePnProxySliceNoFSIPlot_"+TString(std::to_string(i))));
		PrettyPlot(PnNoFSI,kMagenta+1);
		PnNoFSI->Draw("p0 hist same");						

		TLegend* leg = new TLegend(0.7,0.6,0.87,0.9);
		leg->AddEntry(LFG,"LFG","p");
		leg->AddEntry(kMiss,"k_{Miss}","p");
		leg->AddEntry(kMissNoFSI,"k_{Miss}^{No FSI}","p");		
		leg->AddEntry(Pn,"P_{n,proxy}","p");
		leg->AddEntry(PnNoFSI,"P_{n,proxy}^{No FSI}","p");		
		leg->SetTextFont(FontStyle);
		leg->SetTextSize(0.05);			
		if (i ==0 || i == 1) { leg->Draw(); }

		TLegend* legPanel = new TLegend(0.8,0.7,1.,0.9);
		legPanel->AddEntry(LFG,"LFG","p");
		legPanel->AddEntry(kMiss,"k_{Miss}","p");
		legPanel->AddEntry(Pn,"P_{n,proxy}","p");
		legPanel->SetTextFont(FontStyle);
		legPanel->SetTextSize(0.05);	

		can->SaveAs("myPlots/"+CanvasName+".pdf");
		delete can;

		// -------------------------------------------------------------------------

		if (i > 0) {

			canAll->cd();
			TString PadName = "Pad";
			TPad* pad = new TPad(PadName,PadName,0.,0.5,0.5,1.,21); pad->SetBottomMargin(0.0); pad->SetRightMargin(0.0); 
			if (i == 2) { pad = new TPad(PadName,PadName,0.5,0.5,1.,1.,21); pad->SetBottomMargin(0.0); pad->SetLeftMargin(0.0); }
			if (i == 3) { pad = new TPad(PadName,PadName,0.,0.,0.5,0.5,21); pad->SetTopMargin(0.0); pad->SetRightMargin(0.0); pad->SetBottomMargin(0.12); }
			if (i == 4) { pad = new TPad(PadName,PadName,0.5,0.,1.,0.5,21); pad->SetTopMargin(0.0); pad->SetLeftMargin(0.0); pad->SetBottomMargin(0.12); }

			//pad->SetTopMargin(0.01);
			//pad->SetRightMargin(0.01);				

			pad->SetFillColor(kWhite);			
			pad->Draw();
			pad->cd();
			LFG->Draw("p hist same");
			kMiss->Draw("p0 hist same");	
			Pn->Draw("p0 hist same");

			if (i == 4) { text->DrawLatexNDC(0.65,0.85,RangeLabel[i]); }
			else { text->DrawLatexNDC(0.5,0.5,RangeLabel[i]); }					

			if (i == 1) { legPanel->Draw(); }

			gPad->RedrawAxis();											

		}

	}

	canAll->SaveAs("myPlots/"+CanvasNameAll+".pdf");
	delete canAll;		

	// -----------------------------------------------------------------------------------------------------------------------------------------

} // End of the program 
