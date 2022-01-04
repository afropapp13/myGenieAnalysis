#define GenieAnalysis_cxx
#include "GenieAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <iomanip>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>

#include "../myClasses/Tools.h"
#include "../myClasses/STV_Tools.h"
#include "../myClasses/Util.h"

using namespace std;

//----------------------------------------//

void GenieAnalysis::Loop() {

	if (fChain == 0) return;
	TH1D::SetDefaultSumw2();
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	// ---------------------------------------------------------------------------------------------------------------------------------	

	TString FileNameAndPath = "OutputFiles/STVAnalysis_"+OutFileName+".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");
	std::cout << std::endl << std::endl;

	// ---------------------------------------------------------------------------------------------------------------------------------

	// Plots for reco LFG pn, kMiss & pn,proxy in slices
	// Index 0 corresponds to all the events

	TH1D* TrueLFGPnSlicePlot[NRanges+1]; 
	TH1D* TruekMissSlicePlot[NRanges+1];
	TH1D* TruePnProxySlicePlot[NRanges+1];
	TH1D* TruekMissSliceNoFSIPlot[NRanges+1];
	TH1D* TruePnProxySliceNoFSIPlot[NRanges+1];		

	for (int i = 0; i < NRanges + 1; i++) {

		TrueLFGPnSlicePlot[i] = new TH1D("TrueLFGPnSlicePlot_"+TString(std::to_string(i)),";p_{n} [GeV/c]",NBins,LowEdge,HighEdge);
		TruekMissSlicePlot[i] = new TH1D("TruekMissSlicePlot_"+TString(std::to_string(i)),";p_{n} [GeV/c]",NBins,LowEdge,HighEdge);
		TruePnProxySlicePlot[i] = new TH1D("TruePnProxySlicePlot_"+TString(std::to_string(i)),";p_{n} [GeV/c]",NBins,LowEdge,HighEdge);
		TruekMissSliceNoFSIPlot[i] = new TH1D("TruekMissSliceNoFSIPlot_"+TString(std::to_string(i)),";p_{n} [GeV/c]",NBins,LowEdge,HighEdge);
		TruePnProxySliceNoFSIPlot[i] = new TH1D("TruePnProxySliceNoFSIPlot_"+TString(std::to_string(i)),";p_{n} [GeV/c]",NBins,LowEdge,HighEdge);		

	}

	// ---------------------------------------------------------------------------------------------------------------------------------	

	TH1D* TrueDeltaPTPlot = new TH1D("TrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* TrueDeltaAlphaTPlot = new TH1D("TrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* TrueDeltaPhiTPlot = new TH1D("TrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
	TH1D* TrueDeltaPLPlot = new TH1D("TrueDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
	TH1D* TrueDeltaPnPlot = new TH1D("TrueDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
	TH1D* TrueDeltaPtxPlot = new TH1D("TrueDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
	TH1D* TrueDeltaPtyPlot = new TH1D("TrueDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
	TH1D* TrueAPlot = new TH1D("TrueAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);
	TH1D* TrueMuonMomentumPlot = new TH1D("TrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* TrueMuonPhiPlot = new TH1D("TrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* TrueMuonCosThetaPlot = new TH1D("TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
	TH1D* TrueMuonCosThetaSingleBinPlot = new TH1D("TrueMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,-1.,1.);
	TH1D* TrueProtonMomentumPlot = new TH1D("TrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* TrueProtonPhiPlot = new TH1D("TrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* TrueProtonCosThetaPlot = new TH1D("TrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
	TH1D* TrueECalPlot = new TH1D("TrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueECalLowPTPlot = new TH1D("TrueECalLowPTPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueECalMidPTPlot = new TH1D("TrueECalMidPTPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueECalHighPTPlot = new TH1D("TrueECalHighPTPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueEQEPlot = new TH1D("TrueEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
	TH1D* TrueQ2Plot = new TH1D("TrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);
	TH1D* TrueCCQEMuonMomentumPlot = new TH1D("TrueCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
	TH1D* TrueCCQEMuonPhiPlot = new TH1D("TrueCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);
	TH1D* TrueCCQEMuonCosThetaPlot = new TH1D("TrueCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);	
	TH1D* TrueCCQEProtonMomentumPlot = new TH1D("TrueCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
	TH1D* TrueCCQEProtonPhiPlot = new TH1D("TrueCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);
	TH1D* TrueCCQEProtonCosThetaPlot = new TH1D("TrueCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
	TH1D* TrueCCQEECalPlot = new TH1D("TrueCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
	TH1D* TrueCCQEQ2Plot = new TH1D("TrueCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);	
	TH1D* TruekMissPlot = new TH1D("TruekMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
	TH1D* TruePMissMinusPlot = new TH1D("TruePMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
	TH1D* TruePMissPlot = new TH1D("TruePMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);
	TH1D* TruePMissNoFSIPlot = new TH1D("TruePMissNoFSIPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);
	TH1D* TruekMissNoFSIPlot = new TH1D("TruekMissNoFSIPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
	TH1D* TrueDeltaPnNoFSIPlot = new TH1D("TrueDeltaPnNoFSIPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);	 	 	 

	// For now and until box opening
	int NBins2DAnalysis = 4;
		
	TH2D* TrueCosThetaMuPmuPlot = new TH2D("TrueCosThetaMuPmuPlot",LabelXAxisMuonCosTheta+LabelXAxisMuonMomentum
			,NBins2DAnalysis,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			,NBins2DAnalysis,ArrayNBinsMuonMomentum[0],ArrayNBinsMuonMomentum[NBinsMuonMomentum]);
			
	TH2D* TrueCosThetaPPpPlot = new TH2D("TrueCosThetaPPpPlot",LabelXAxisProtonCosTheta+LabelXAxisProtonMomentum
			,NBins2DAnalysis,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			,NBins2DAnalysis,ArrayNBinsProtonMomentum[0],ArrayNBinsProtonMomentum[NBinsProtonMomentum]);	

	// ------------------------------------------------------------------------------------------------------------------------------------------------

	TH2D* CVWeightsVsEvPlot = new TH2D("CVWeightsVsEvPlot",";Ev [GeV];cv weights",200,0.,2.,300,0.,3.);
	TH2D* CVWeightsVsInteractionModePlot = new TH2D("CVWeightsVsInteractionModePlot",";CVWeightsVsInteractionModePlot;cv weights",4,-0.5,3.5,300,0.,3.);

	// ------------------------------------------------------------------------------------------------------------------------------------------------

	Tools tools;

	//--------------------------------------------------//

	// Now let's get serious with the 2D analysis

	// Ecal in DeltaPT and DeltaAlphaT bins		

	TH1D* TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];	

	//--------------------------------------------------//

	// DeltaAlphaT in DeltaPT slices

	TH1D* TrueDeltaAlphaT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];	

	//--------------------------------------------------//

	// DeltaPhiT in DeltaPT slices

	TH1D* TrueDeltaPhiT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];	

	//--------------------------------------------------//

	// DeltaPn in DeltaPT slices

	TH1D* TrueDeltaPn_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];	

	//--------------------------------------------------//	

	for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

		//------------------------------//

		// DeltaAlphaT in DeltaPT slices

		TString DeltaAlphaTTwoDInDeltaPTLabel = "DeltaAlphaT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			

		TrueDeltaAlphaT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("True"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);

		//------------------------------//

		// DeltaPhiT in DeltaPT slices

		TString DeltaPhiTTwoDInDeltaPTLabel = "DeltaPhiT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
		TrueDeltaPhiT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("True"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);	

		//------------------------------//

		// DeltaPn in DeltaPT slices

		TString DeltaPnTwoDInDeltaPTLabel = "DeltaPn_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
		TrueDeltaPn_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("True"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);	

		//------------------------------//

		for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {	

			TString ECalTwoDInDeltaPTDeltaAlphaTLabel = "ECal_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";
			TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[WhichDeltaPT][WhichDeltaAlphaT] = new TH1D("True"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);

		}

	}

	//--------------------------------------------------//

	// DeltaAlphaT in muon cos theta slices

	TH1D* TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];

	// DeltaPT in muon cos theta slices

	TH1D* TrueDeltaPT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];	

	// Pmu in muon cos theta slices

	TH1D* TrueMuonMomentum_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];

	for (int WhichMuonCosTheta = 0; WhichMuonCosTheta < TwoDNBinsMuonCosTheta; WhichMuonCosTheta++) {

		TString DeltaAlphaTTwoDInMuonCosThetaLabel = "DeltaAlphaT_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
		TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("True"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);

		TString DeltaPTTwoDInMuonCosThetaLabel = "DeltaPT_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
		TrueDeltaPT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("True"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);		

		TString MuonMomentumTwoDInMuonCosThetaLabel = "MuonMomentum_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
		TrueMuonMomentum_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("True"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);

	}	

	//--------------------------------------------------//

	// Pp in proton cos theta slices

	TH1D* TrueProtonMomentum_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];

	for (int WhichProtonCosTheta = 0; WhichProtonCosTheta < TwoDNBinsProtonCosTheta; WhichProtonCosTheta++) {

		TString ProtonMomentumTwoDInProtonCosThetaLabel = "ProtonMomentum_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
		TrueProtonMomentum_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("True"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);	

	}

	//--------------------------------------------------//

	// DeltaPty in DeltaPtx slices

	TH1D* TrueDeltaPty_InDeltaPtxTwoDPlot[TwoDNBinsDeltaPtx];

	for (int WhichDeltaPtx = 0; WhichDeltaPtx < TwoDNBinsDeltaPtx; WhichDeltaPtx++) {

		TString DeltaPtyTwoDInDeltaPtxLabel = "DeltaPty_DeltaPtx_"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx+1])+"Plot";			
		TrueDeltaPty_InDeltaPtxTwoDPlot[WhichDeltaPtx] = new TH1D("True"+DeltaPtyTwoDInDeltaPtxLabel,LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);	

	}	

	//--------------------------------------------------//		

	// ---------------------------------------------------------------------------------------------------------------------------------

	// Counters STV analysis

	int CounterSTVEventsPassedSelection = 0;
	int CounterSTVQEEventsPassedSelection = 0;
	int CounterSTVMECEventsPassedSelection = 0;
	int CounterSTVRESEventsPassedSelection = 0;
	int CounterSTVDISEventsPassedSelection = 0;
	int CounterSTVCOHEventsPassedSelection = 0;		

	// ---------------------------------------------------------------------------------------------------------------------------------	

	// T2K Weights to make GENIE consistent with the MicroBooNE tune v2 (or v1)
	
	TFile* fweights = nullptr;
	TTree* tweights = nullptr;
	float cv_weight = -99.;

	if (OutFileName == "Genie_v3_0_6_uB_Tune_1") {

		fweights = TFile::Open("mySamples/myWeights_uB_Tune_v1.root");
		tweights = (TTree*)fweights->Get("ub_tune_cv");
		tweights->SetBranchAddress("cv_weight", &cv_weight);

	}

	if (OutFileName == "Genie_v3_0_6_Nominal" || OutFileName == "Genie_v3_0_6_NoFSI" || OutFileName == "Genie_v3_0_6_NoRPA" || OutFileName == "Genie_v3_0_6_NoCoulomb" 
	|| OutFileName == "Genie_v3_0_6_hN2018" || OutFileName == "Genie_v3_0_6_RFG" || OutFileName == "Genie_v3_0_6_EffSF") {

		if (OutFileName == "Genie_v3_0_6_Nominal" || OutFileName == "Genie_v3_0_6_NoFSI" || OutFileName == "Genie_v3_0_6_hN2018") 
			{ fweights = TFile::Open("mySamples/myWeights_uB_Tune_Nominal.root"); }
		if (OutFileName == "Genie_v3_0_6_NoRPA") { fweights = TFile::Open("mySamples/myWeights_uB_Tune_NoRPA.root"); }
		if (OutFileName == "Genie_v3_0_6_NoCoulomb") { fweights = TFile::Open("mySamples/myWeights_uB_Tune_NoCoulomb.root"); }
		if (OutFileName == "Genie_v3_0_6_RFG") { fweights = TFile::Open("mySamples/myWeights_uB_Tune_RFG.root"); }
		if (OutFileName == "Genie_v3_0_6_EffSF") { fweights = TFile::Open("mySamples/myWeights_uB_Tune_EffSF.root"); }

		tweights = (TTree*)fweights->Get("GenericVectors__VARS");
		tweights->SetBranchAddress("Weight", &cv_weight);

	}

	// ---------------------------------------------------------------------------------------------------------------------------------	

	int SumCVWeights = 0;
	int SumSelectedCVWeights = 0;	

	// ---------------------------------------------------------------------------------------------------------------------------------	

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry); if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

		// ---------------------------------------------------------------------------------------------------------------------------------

		double weight = 1.;

		if (OutFileName == "Genie_v3_0_6_uB_Tune_1" || OutFileName == "Genie_v3_0_6_Nominal" || OutFileName == "Genie_v3_0_6_NoFSI" || OutFileName == "Genie_v3_0_6_hN2018" 
		|| OutFileName == "Genie_v3_0_6_NoRPA" || OutFileName == "Genie_v3_0_6_NoCoulomb" || OutFileName == "Genie_v3_0_6_RFG" || OutFileName == "Genie_v3_0_6_EffSF") 
			{ tweights->GetEntry(jentry); weight = cv_weight; }

		if (weight <= 0 || weight >= 30) { continue; }
		SumCVWeights += weight;
		
		// ---------------------------------------------------------------------------------------------------------------------------------

		// CV weight studies

		CVWeightsVsEvPlot->Fill(Ev,cv_weight);

		double genie_mode = -1.;

		if (qel==1) { genie_mode = 0; }
		if (mec==1) { genie_mode = 1; }
		if (res==1) { genie_mode = 2; }
		if (dis==1) { genie_mode = 3; }
		if (coh==1) { genie_mode = 4; }		
		
		CVWeightsVsInteractionModePlot->Fill(genie_mode,cv_weight);

		// ---------------------------------------------------------------------------------------------------------------------------------

		if (!cc) { continue; }

		// -----------------------------------------------------------------------------------------------------------------------------------

		// Muon

		TVector3 Muon3Vector(pxl,pyl,pzl); // GeV
		TLorentzVector Muon4Vector(pxl,pyl,pzl,El); // GeV
		double MuonMomentum = pl; // GeV / c
		double MuonTheta = TMath::ACos(cthl) * 180. / TMath::Pi(); // deg
		double MuonCosTheta = cthl;
		double MuonPhi = Muon4Vector.Phi() * 180. / TMath::Pi(); // deg		

		// ----------------------------------------------------------------------------------------------------------------------------------

		// True struck nucleon momentum (LFG / RFG et al)

		double LFG_pn = TMath::Sqrt(pxn*pxn+pyn*pyn+pzn*pzn);	

		int Index = -1;

		for (int i = 0; i < NRanges; i++) {

			// Don't forget the offste by 1, 0 corresponds to all the events
			if (LFG_pn > range[i] && LFG_pn < range[i+1]) { Index = i+1; }

		}			
		
		// ---------------------------------------------------------------------------------------------------------------------------------	

		int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0, TrueHeavierMesonCounter = 0;
		vector <int> ProtonID; ProtonID.clear();

		// --------------------------------------------------------------------

		// No FSI case, use initial state particles

		int ProtonTaggingNoFSI = 0, ChargedPionTaggingNoFSI = 0, NeutralPionTaggingNoFSI = 0, TrueHeavierMesonCounterNoFSI = 0;
		vector <int> ProtonIDNoFSI; ProtonIDNoFSI.clear();

		for (int i = 0; i < ni; i++) {

			double pi = TMath::Sqrt( TMath::Power(pxi[i],2.) + TMath::Power(pyi[i],2.) + TMath::Power(pzi[i],2.));

			if (pdgi[i] == ProtonPdg && pi > ArrayNBinsProtonMomentum[0]) {

				ProtonTaggingNoFSI ++;
				ProtonIDNoFSI.push_back(i);

			}

			if (fabs(pdgi[i]) == AbsChargedPionPdg && pi > ChargedPionMomentumThres)  {

				ChargedPionTaggingNoFSI ++;

			}

			if (pdgi[i] == NeutralPionPdg)  {

				NeutralPionTaggingNoFSI ++;

			}

			if ( pdgi[i] != NeutralPionPdg && fabs(pdgi[i]) != AbsChargedPionPdg && tools.is_meson_or_antimeson(pdgi[i]) ) { TrueHeavierMesonCounterNoFSI++; }

		} // End of the loop over the initial state particles, because there are no final state particles in a No FSI sample


		if ( ProtonTaggingNoFSI == 1 && ChargedPionTaggingNoFSI == 0 && NeutralPionTaggingNoFSI == 0 && TrueHeavierMesonCounterNoFSI == 0) { 

			// Proton

			int ProtonIndexNoFSI = ProtonIDNoFSI.at(0);
			TVector3 Proton3Vector(pxi[ProtonIndexNoFSI],pyi[ProtonIndexNoFSI],pzi[ProtonIndexNoFSI]); // GeV
			STV_Tools stv_tool(Muon3Vector,Proton3Vector,El,Ei[ProtonIndexNoFSI]);	

			TruekMissSliceNoFSIPlot[0]->Fill(stv_tool.ReturnkMiss(),weight);
			TruePnProxySliceNoFSIPlot[0]->Fill(stv_tool.ReturnPn(),weight);

			TruekMissSliceNoFSIPlot[Index]->Fill(stv_tool.ReturnkMiss(),weight);
			TruePnProxySliceNoFSIPlot[Index]->Fill(stv_tool.ReturnPn(),weight);							

		}

		// --------------------------------------------------------------------

		for (int i = 0; i < nf; i++) {

			if (pdgf[i] == ProtonPdg && pf[i] > ArrayNBinsProtonMomentum[0]) {

				ProtonTagging ++;
				ProtonID.push_back(i);

			}

			if (fabs(pdgf[i]) == AbsChargedPionPdg && pf[i] > ChargedPionMomentumThres)  {

				ChargedPionTagging ++;

			}

			if (pdgf[i] == NeutralPionPdg)  {

				NeutralPionTagging ++;

			}

			if ( pdgf[i] != NeutralPionPdg && fabs(pdgf[i]) != AbsChargedPionPdg && tools.is_meson_or_antimeson(pdgf[i]) ) { TrueHeavierMesonCounter++; }

		} // End of the loop over the final state particles


		// -----------------------------------------------------------------------------------------------------------------------------------

		// Impose signal
		// 1 muon > 0.1 GeV/c
		// 1 proton > 0.3 GeV/c
		// no charged pions > 0.07 GeV/c
		// no neutral pions of any momenta
		// any number of neutrons & photons

		if ( ProtonTagging != 1 || ChargedPionTagging != 0 || NeutralPionTagging != 0 ) { continue; }
		if ( TrueHeavierMesonCounter != 0 ) { continue; }

		// ------------------------------------------------------------------------------------------------------------------------------------

		// Proton

		int ProtonIndex = ProtonID.at(0);
		TVector3 Proton3Vector(pxf[ProtonIndex],pyf[ProtonIndex],pzf[ProtonIndex]); // GeV
		TLorentzVector Proton4Vector(pxf[ProtonIndex],pyf[ProtonIndex],pzf[ProtonIndex],Ef[ProtonIndex]); // GeV
		double ProtonMomentum = pf[ProtonIndex]; // GeV / c
//		double ProtonTheta = TMath::ACos(cthf[ProtonIndex]) * 180. / TMath::Pi(); // deg
		double ProtonCosTheta = cthf[ProtonIndex];
		double ProtonPhi = Proton4Vector.Phi() * 180. / TMath::Pi(); // deg
		double ProtonE = Proton4Vector.E(); // GeV

		// ---------------------------------------------------------------------------------------------------------------------------------

		// Relative Proton - Muon Angles

		double DeltaThetaProtonMuon = Proton3Vector.Angle(Muon3Vector) * 180. / TMath::Pi();
		if (DeltaThetaProtonMuon < 0.) { DeltaThetaProtonMuon += 180.; }
		if (DeltaThetaProtonMuon > 180.) { DeltaThetaProtonMuon -= 180.; }

		double DeltaPhiProtonMuon = Muon3Vector.DeltaPhi(Proton3Vector)* 180. / TMath::Pi();
		if (DeltaPhiProtonMuon < 0.) { DeltaPhiProtonMuon += 360.; }
		if (DeltaPhiProtonMuon > 360.) { DeltaPhiProtonMuon -= 360.; }

		// ---------------------------------------------------------------------------------------------------------------------------------

		STV_Tools stv_tool(Muon3Vector,Proton3Vector,El,ProtonE);

		double PTmissMomentum = stv_tool.ReturnPt();
		double TrueDeltaAlphaT = stv_tool.ReturnDeltaAlphaT();
		double TrueDeltaPhiT = stv_tool.ReturnDeltaPhiT();
		double ECal = stv_tool.ReturnECal();
		double EQE = stv_tool.ReturnEQE();
		double TrueQ2 = stv_tool.ReturnQ2();	

		double TruekMiss = stv_tool.ReturnkMiss();
		double TruePMissMinus = stv_tool.ReturnPMissMinus();
		double TrueMissMomentum = stv_tool.ReturnPMiss();

		double TruePL = stv_tool.ReturnPL();
		double TruePn = stv_tool.ReturnPn();
		double TruePtx = stv_tool.ReturnPtx();
		double TruePty = stv_tool.ReturnPty();
		double TrueA = stv_tool.ReturnA();

		// -----------------------------------------------------------------------------------------------------------------
		// -----------------------------------------------------------------------------------------------------------------	

		// Overflow bins
		// Affects

		// DeltaPT
		// DeltaPtx
		// DeltaPty
		// DeltaPL
		// DeltaPn
		// Q2
		// ECal
		// EQE
		// alpha
		// kMiss
		// PMiss
		// PMissMinus


		if (PTmissMomentum > ArrayNBinsDeltaPT[NBinsDeltaPT]) { PTmissMomentum = 0.5 * (ArrayNBinsDeltaPT[NBinsDeltaPT] + ArrayNBinsDeltaPT[NBinsDeltaPT-1]); }
		if (TruePtx > ArrayNBinsDeltaPtx[NBinsDeltaPtx]) { TruePtx = 0.5 * (ArrayNBinsDeltaPtx[NBinsDeltaPtx] + ArrayNBinsDeltaPtx[NBinsDeltaPtx-1]); }
		if (TruePty > ArrayNBinsDeltaPty[NBinsDeltaPty]) { TruePty = 0.5 * (ArrayNBinsDeltaPty[NBinsDeltaPty] + ArrayNBinsDeltaPty[NBinsDeltaPty-1]); }
		if (TruePL > ArrayNBinsDeltaPL[NBinsDeltaPL]) { TruePL = 0.5 * (ArrayNBinsDeltaPL[NBinsDeltaPL] + ArrayNBinsDeltaPL[NBinsDeltaPL-1]); }						
		if (TruePn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { TruePn = 0.5 * (ArrayNBinsDeltaPn[NBinsDeltaPn] + ArrayNBinsDeltaPn[NBinsDeltaPn-1]); }

		if (ECal > ArrayNBinsECal[NBinsECal]) { ECal = 0.5 * (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1]); }
		if (EQE > ArrayNBinsEQE[NBinsEQE]) { EQE = 0.5 * (ArrayNBinsEQE[NBinsEQE] + ArrayNBinsEQE[NBinsEQE-1]); }
		if (TrueQ2 > ArrayNBinsQ2[NBinsQ2]) { TrueQ2 = 0.5 * (ArrayNBinsQ2[NBinsQ2] + ArrayNBinsQ2[NBinsQ2-1]); }

		if (TrueA > ArrayNBinsA[NBinsA]) { TrueA = 0.5 * (ArrayNBinsA[NBinsA] + ArrayNBinsA[NBinsA-1]); }	
		if (TruekMiss > ArrayNBinskMiss[NBinskMiss]) { TruekMiss = 0.5 * (ArrayNBinskMiss[NBinskMiss] + ArrayNBinskMiss[NBinskMiss-1]); }														
		if (TrueMissMomentum > ArrayNBinsPMiss[NBinsPMiss]) { TrueMissMomentum = 0.5 * (ArrayNBinsPMiss[NBinsPMiss] + ArrayNBinsPMiss[NBinsPMiss-1]); }
		if (TruePMissMinus > ArrayNBinsPMissMinus[NBinsPMissMinus]) { TruePMissMinus = 0.5 * (ArrayNBinsPMissMinus[NBinsPMissMinus] + ArrayNBinsPMissMinus[NBinsPMissMinus-1]); }

		// ---------------------------------------------------------------------------------------------------------------------------

		// Underflow bins
		// Affects

		// ECal
		// EQE
		// DeltaPtx
		// DeltaPty
		// DeltaPL
		// alpha
		// PMissMinus
			
		if (ECal < ArrayNBinsECal[0]) { ECal = 0.5 * (ArrayNBinsECal[0] + ArrayNBinsECal[1]); }			
		if (EQE < ArrayNBinsEQE[0]) { EQE = 0.5 * (ArrayNBinsEQE[0] + ArrayNBinsEQE[1]); }			
		if (TruePtx < ArrayNBinsDeltaPtx[0]) { TruePtx = 0.5 * (ArrayNBinsDeltaPtx[0] + ArrayNBinsDeltaPtx[1]); }
		if (TruePty < ArrayNBinsDeltaPty[0]) { TruePty = 0.5 * (ArrayNBinsDeltaPty[0] + ArrayNBinsDeltaPty[1]); }
		if (TruePL < ArrayNBinsDeltaPL[0]) { TruePL = 0.5 * (ArrayNBinsDeltaPL[0] + ArrayNBinsDeltaPL[1]); }						
		if (TrueA < ArrayNBinsA[0]) { TrueA = 0.5 * (ArrayNBinsA[0] + ArrayNBinsA[1]); }
		if (TruePMissMinus < ArrayNBinsPMissMinus[0]) { TruePMissMinus = 0.5 * (ArrayNBinsPMissMinus[0] + ArrayNBinsPMissMinus[1]); }		

		// ----------------------------------------------------------------------------------------------------------------------------
		// ---------------------------------------------------------------------------------------------------------------------------									

		if (
			MuonMomentum > ArrayNBinsMuonMomentum[0]
		     && ProtonMomentum > ArrayNBinsProtonMomentum[0]
		) {

			if (		
			    // Same events fill all the plots
			        
			    PTmissMomentum > ArrayNBinsDeltaPT[0] //&& PTmissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
			    && TrueDeltaAlphaT > ArrayNBinsDeltaAlphaT[0] && TrueDeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
		 	    && TrueDeltaPhiT > ArrayNBinsDeltaPhiT[0] && TrueDeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]
		 	    
			    && MuonMomentum < ArrayNBinsMuonMomentum[NBinsMuonMomentum]  
			    && ProtonMomentum < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
			    
			    && MuonCosTheta > ArrayNBinsMuonCosTheta[0]
			    && MuonCosTheta < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			    && ProtonCosTheta > ArrayNBinsProtonCosTheta[0]
			    && ProtonCosTheta < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			    
//			    && ECal > ArrayNBinsECal[0] && ECal < ArrayNBinsECal[NBinsECal]
//			    && EQE > ArrayNBinsEQE[0] && EQE < ArrayNBinsEQE[NBinsEQE]
//			    && TrueQ2 > ArrayNBinsQ2[0] && TrueQ2 < ArrayNBinsQ2[NBinsQ2]
			    
			) {

				SumSelectedCVWeights += weight;

				TrueDeltaPTPlot->Fill(PTmissMomentum,weight);
				TrueDeltaAlphaTPlot->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaPhiTPlot->Fill(TrueDeltaPhiT,weight);

				TrueECalPlot->Fill(ECal,weight);
				if (PTmissMomentum > LowPT[0] && PTmissMomentum < HighPT[0]) { TrueECalLowPTPlot->Fill(ECal,weight); }
				if (PTmissMomentum > LowPT[1] && PTmissMomentum < HighPT[1]) { TrueECalMidPTPlot->Fill(ECal,weight); }
				if (PTmissMomentum > LowPT[2] && PTmissMomentum < HighPT[2]) { TrueECalHighPTPlot->Fill(ECal,weight); }
				TrueEQEPlot->Fill(EQE,weight);
				TrueQ2Plot->Fill(TrueQ2,weight);

				CounterSTVEventsPassedSelection++;
				if (qel==1) { CounterSTVQEEventsPassedSelection++; }
				if (mec==1) { CounterSTVMECEventsPassedSelection++; }
				if (res==1) { CounterSTVRESEventsPassedSelection++; }
				if (dis==1) { CounterSTVDISEventsPassedSelection++; }
				if (coh==1) { CounterSTVCOHEventsPassedSelection++; }				

				TrueMuonMomentumPlot->Fill(MuonMomentum,weight);
				TrueMuonPhiPlot->Fill(MuonPhi,weight);
				TrueMuonCosThetaPlot->Fill(MuonCosTheta,weight);
				TrueMuonCosThetaSingleBinPlot->Fill(MuonCosTheta,weight);

				TrueProtonMomentumPlot->Fill(ProtonMomentum,weight);
				TrueProtonPhiPlot->Fill(ProtonPhi,weight);
				TrueProtonCosThetaPlot->Fill(ProtonCosTheta,weight);			

				TruekMissPlot->Fill(TruekMiss,weight);
				TruePMissMinusPlot->Fill(TruePMissMinus,weight);
				TruePMissPlot->Fill(TrueMissMomentum,weight);

				TrueDeltaPLPlot->Fill(TruePL,weight);
				TrueDeltaPnPlot->Fill(TruePn,weight);
				TrueDeltaPtxPlot->Fill(TruePtx,weight);
				TrueDeltaPtyPlot->Fill(TruePty,weight);
				TrueAPlot->Fill(TrueA,weight);

				//-----------------------------------------------------------------------------------------

				TrueLFGPnSlicePlot[0]->Fill(LFG_pn,weight);				
				TruekMissSlicePlot[0]->Fill(TruekMiss,weight);
				TruePnProxySlicePlot[0]->Fill(TruePn,weight);

				TrueLFGPnSlicePlot[Index]->Fill(LFG_pn,weight);				
				TruekMissSlicePlot[Index]->Fill(TruekMiss,weight);
				TruePnProxySlicePlot[Index]->Fill(TruePn,weight);	

				//----------------------------------------//

				// Indices for 2D analysis

				int DeltaPTTwoDIndex = tools.ReturnIndex(PTmissMomentum, TwoDArrayNBinsDeltaPT);
				int DeltaAlphaTTwoDIndex = tools.ReturnIndex(TrueDeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
				int MuonCosThetaTwoDIndex = tools.ReturnIndex(MuonCosTheta, TwoDArrayNBinsMuonCosTheta);
				int ProtonCosThetaTwoDIndex = tools.ReturnIndex(ProtonCosTheta, TwoDArrayNBinsProtonCosTheta);
				int DeltaPtxTwoDIndex = tools.ReturnIndex(TruePtx, TwoDArrayNBinsDeltaPtx);								

				//----------------------------------------//

				// 2D analysis	

				TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
				TrueDeltaAlphaT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaPhiT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(TrueDeltaPhiT,weight);	
				TrueDeltaPn_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(TruePn,weight);
				TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaPT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(PTmissMomentum,weight);
				TrueMuonMomentum_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(MuonMomentum,weight);
				TrueProtonMomentum_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(ProtonMomentum,weight);
				TrueDeltaPty_InDeltaPtxTwoDPlot[DeltaPtxTwoDIndex]->Fill(TruePty,weight);

				//----------------------------------------//

				if (		
					// CCQElike measurement
						
					PTmissMomentum < 0.35
					&& TMath::Abs(DeltaPhiProtonMuon - 180) < 35
					&& TMath::Abs(DeltaThetaProtonMuon - 90) < 55

					&& MuonMomentum > CCQEArrayNBinsMuonMomentum[0]  
					&& ProtonMomentum > CCQEArrayNBinsProtonMomentum[0]					
					&& MuonMomentum < CCQEArrayNBinsMuonMomentum[CCQENBinsMuonMomentum]  
					&& ProtonMomentum < CCQEArrayNBinsProtonMomentum[CCQENBinsProtonMomentum]
					
					&& MuonCosTheta > CCQEArrayNBinsMuonCosTheta[0]
					&& MuonCosTheta < CCQEArrayNBinsMuonCosTheta[CCQENBinsMuonCosTheta]
					&& ProtonCosTheta > CCQEArrayNBinsProtonCosTheta[0]
					&& ProtonCosTheta < CCQEArrayNBinsProtonCosTheta[CCQENBinsProtonCosTheta]

				) {

					TrueCCQEMuonMomentumPlot->Fill(MuonMomentum,weight);
					TrueCCQEMuonPhiPlot->Fill(MuonPhi,weight);
					TrueCCQEMuonCosThetaPlot->Fill(MuonCosTheta,weight);

					TrueCCQEProtonMomentumPlot->Fill(ProtonMomentum,weight);
					TrueCCQEProtonPhiPlot->Fill(ProtonPhi,weight);
					TrueCCQEProtonCosThetaPlot->Fill(ProtonCosTheta,weight);	

					TrueCCQEECalPlot->Fill(ECal,weight);
					TrueCCQEQ2Plot->Fill(TrueQ2,weight);
													
				}	

			}

		} // End of the selection cuts && the demand that we fill the plots with the same events

		// ---------------------------------------------------------------------------------------------------------------------------------

	} // End of the loop over the events

	// -----------------------------------------------------------------------------------------------------------------------------------------

	// STV analysis
	
	std::cout << std::endl << "----------------------- STV analysis -------------------------" << std::endl << std::endl;
	std::cout << "Percetage of events passing the selection cuts = " << 
	double(CounterSTVEventsPassedSelection)/ double(nentries)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting QE events = " << 
	double(CounterSTVQEEventsPassedSelection)/ double(CounterSTVEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting MEC events = " << 
	double(CounterSTVMECEventsPassedSelection)/ double(CounterSTVEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting RES events = " << 
	double(CounterSTVRESEventsPassedSelection)/ double(CounterSTVEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting DIS events = " << 
	double(CounterSTVDISEventsPassedSelection)/ double(CounterSTVEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;
	
	std::cout << "Success percetage in selecting COH events = " << 
	double(CounterSTVCOHEventsPassedSelection)/ double(CounterSTVEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;	

	// -------------------------------------------------------------------------------------------------------------------------------------------

	// Moving forward towards a cross section extraction
	// Reweight the genie sample with the sum of cv tune weights / number of events
	// Multiply by the flux integrated xsec that has been calculated using flux_averaged_total_xsec

	double ScalingFactor = 1. / (double)(SumCVWeights);
//	double ScalingFactor = 1. / (double)(nentries);

	cout << "# Events = " << nentries << endl;
	cout << "Sum Weights = " << SumCVWeights << endl;

	// -------------------------------------------------------------------------------------------------------------------------------------------

	// Flux integrated cross sections

	if (OutFileName == "Genie_v3_0_6_uB_Tune_1" || OutFileName == "Genie_v3_0_6_Out_Of_The_Box" || OutFileName == "Genie_v3_0_6_Nominal" 
	 || OutFileName == "Genie_v3_0_6_NoFSI"  || OutFileName == "Genie_v3_0_6_hN2018") 
		{ ScalingFactor = ScalingFactor * G18_10a_02_11a_FluxIntegratedXSection ; }
	if (OutFileName == "SuSav2") { ScalingFactor = ScalingFactor * SuSav2FluxIntegratedXSection ; }
	if (OutFileName == "GENIEv2") { ScalingFactor = ScalingFactor * R_2_12_10_FluxIntegratedXSection ; }
	if (OutFileName == "GENIEv3_0_4") { ScalingFactor = ScalingFactor * R_3_0_4_FluxIntegratedXSection ; }
	if (OutFileName == "Genie_v3_0_6_NoRPA") { ScalingFactor = ScalingFactor * R_3_0_6_G18_10a_02_11a_NoRPA_FluxIntegratedXSection ; }
	if (OutFileName == "Genie_v3_0_6_NoCoulomb") { ScalingFactor = ScalingFactor * R_3_0_6_G18_10a_02_11a_NoCoulomb_FluxIntegratedXSection ; }
	if (OutFileName == "Genie_v3_0_6_RFG") { ScalingFactor = ScalingFactor * R_3_0_6_G18_10a_02_11a_RFG_FluxIntegratedXSection ; }
	if (OutFileName == "Genie_v3_0_6_EffSF") { ScalingFactor = ScalingFactor * R_3_0_6_G18_10a_02_11a_EffSF_FluxIntegratedXSection ; }

	file->cd();

	// -------------------------------------------------------------------------------------------------------------------------------------------

	// Division by bin width to get the cross sections in true variables
	
	tools.Reweight(TrueMuonMomentumPlot,ScalingFactor);
	tools.Reweight(TrueMuonPhiPlot,ScalingFactor);
	tools.Reweight(TrueMuonCosThetaPlot,ScalingFactor);
	// Factor of 2 to account for the fact that the bin width is 2, but we want the number of events, as if the bin width is 1
	tools.Reweight(TrueMuonCosThetaSingleBinPlot,2*ScalingFactor);

	tools.Reweight(TrueProtonMomentumPlot,ScalingFactor);
	tools.Reweight(TrueProtonPhiPlot,ScalingFactor);
	tools.Reweight(TrueProtonCosThetaPlot,ScalingFactor);

	tools.Reweight(TrueECalPlot,ScalingFactor);
	tools.Reweight(TrueECalLowPTPlot,ScalingFactor);
	tools.Reweight(TrueECalMidPTPlot,ScalingFactor);	
	tools.Reweight(TrueECalHighPTPlot,ScalingFactor);

	tools.Reweight(TrueEQEPlot,ScalingFactor);	
	tools.Reweight(TrueQ2Plot,ScalingFactor);

	tools.Reweight(TrueCCQEMuonMomentumPlot,ScalingFactor);
	tools.Reweight(TrueCCQEMuonPhiPlot,ScalingFactor);
	tools.Reweight(TrueCCQEMuonCosThetaPlot,ScalingFactor);

	tools.Reweight(TrueCCQEProtonMomentumPlot,ScalingFactor);
	tools.Reweight(TrueCCQEProtonPhiPlot,ScalingFactor);
	tools.Reweight(TrueCCQEProtonCosThetaPlot,ScalingFactor);

	tools.Reweight(TrueCCQEECalPlot,ScalingFactor);
	tools.Reweight(TrueCCQEQ2Plot,ScalingFactor);	
	
	tools.Reweight(TrueDeltaPTPlot,ScalingFactor);
	tools.Reweight(TrueDeltaAlphaTPlot,ScalingFactor);
	tools.Reweight(TrueDeltaPhiTPlot,ScalingFactor);

	tools.Reweight(TrueDeltaPLPlot,ScalingFactor);
	tools.Reweight(TrueDeltaPnPlot,ScalingFactor);
	tools.Reweight(TrueDeltaPtxPlot,ScalingFactor);
	tools.Reweight(TrueDeltaPtyPlot,ScalingFactor);
	tools.Reweight(TrueAPlot,ScalingFactor);

	tools.Reweight(TruekMissPlot,ScalingFactor);
	tools.Reweight(TruePMissPlot,ScalingFactor);
	tools.Reweight(TruePMissMinusPlot,ScalingFactor);

	tools.Reweight2D(TrueCosThetaMuPmuPlot,ScalingFactor);
	tools.Reweight2D(TrueCosThetaPPpPlot,ScalingFactor);

	file->cd();
	file->Write();

	// -------------------------------------------------------------------------------------------------------------------------------------------

//	// Use smearing matrix to get cross sections in reco space
//	// Using the smearing matrices from the available individual runs

//	std::vector<TString> Runs{"Run1"};
//	const int NRuns = Runs.size();

//	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

//		TFile* FileMigrationMatrix = TFile::Open("../mySTVAnalysis/myMigrationMatrices/"+UBCodeVersion+"/FileMigrationMatrices_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root","readonly");

//		TH1D* RecoMuonMomentumPlot = SmearTrueToReco(TrueMuonMomentumPlot,FileMigrationMatrix,"MuonMomentumPlot","_"+Runs[WhichRun]);
//		TH1D* RecoMuonCosThetaPlot = SmearTrueToReco(TrueMuonCosThetaPlot,FileMigrationMatrix,"MuonCosThetaPlot","_"+Runs[WhichRun]);
//		TH1D* RecoMuonPhiPlot = SmearTrueToReco(TrueMuonCosThetaPlot,FileMigrationMatrix,"MuonPhiPlot","_"+Runs[WhichRun]);

//		TH1D* RecoProtonMomentumPlot = SmearTrueToReco(TrueProtonMomentumPlot,FileMigrationMatrix,"ProtonMomentumPlot","_"+Runs[WhichRun]);
//		TH1D* RecoProtonCosThetaPlot = SmearTrueToReco(TrueProtonCosThetaPlot,FileMigrationMatrix,"ProtonCosThetaPlot","_"+Runs[WhichRun]);
//		TH1D* RecoProtonPhiPlot = SmearTrueToReco(TrueProtonCosThetaPlot,FileMigrationMatrix,"ProtonPhiPlot","_"+Runs[WhichRun]);

//		TH1D* RecoDeltaPTPlot = SmearTrueToReco(TrueDeltaPTPlot,FileMigrationMatrix,"DeltaPTPlot","_"+Runs[WhichRun]);
//		TH1D* RecoDeltaAlphaTPlot = SmearTrueToReco(TrueDeltaAlphaTPlot,FileMigrationMatrix,"DeltaAlphaTPlot","_"+Runs[WhichRun]);
//		TH1D* RecoDeltaPhiTPlot = SmearTrueToReco(TrueDeltaPhiTPlot,FileMigrationMatrix,"DeltaPhiTPlot","_"+Runs[WhichRun]);

//		file->cd();

//		RecoMuonMomentumPlot->Write();
//		RecoMuonCosThetaPlot->Write();
//		RecoMuonPhiPlot->Write();

//		RecoProtonMomentumPlot->Write();
//		RecoProtonCosThetaPlot->Write();
//		RecoProtonPhiPlot->Write();

//		RecoDeltaPTPlot->Write();
//		RecoDeltaAlphaTPlot->Write();
//		RecoDeltaPhiTPlot->Write();

//	}
	
	// --------------------------------------------------------------------------------------------------------------------------------------------	

	std::cout << std::endl;
	std::cout << "File " << FileNameAndPath +" has been created created " << std::endl; 
	std::cout << std::endl;	

	// --------------------------------------------------------------------------------------------------------------------------------------------

} // End of the program
