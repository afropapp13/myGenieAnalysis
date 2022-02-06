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

	//----------------------------------------//

	if (fChain == 0) return;
	TH1D::SetDefaultSumw2();
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	//--------------------------------------------------//

	Tools tools;	

	//----------------------------------------//	

	TString FileNameAndPath = "OutputFiles/STVAnalysis_"+fOutFileName+".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");
	std::cout << std::endl << std::endl;

	//----------------------------------------//

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

	// 1D analysis

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

	//--------------------------------------------------//

	// 2D analysis (uncorrelated)

	TH1D* TrueDeltaPT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
	TH1D* TrueDeltaPT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];	
	TH1D* TrueDeltaPT_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];
	TH1D* TrueMuonMomentum_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
	TH1D* TrueProtonMomentum_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];			
	TH1D* TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[TwoDNBinsMuonCosTheta];
	TH1D* TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];		
	TH1D* TrueDeltaAlphaT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];			
	TH1D* TrueProtonCosTheta_InMuonCosThetaTwoDPlot[TwoDNBinsProtonCosTheta];	
	TH1D* TrueDeltaPhiT_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];	
	TH1D* TrueDeltaPn_InDeltaPTTwoDPlot[TwoDNBinsDeltaPT];
	TH1D* TrueDeltaPty_InDeltaPtxTwoDPlot[TwoDNBinsDeltaPtx];
	TH1D* TrueDeltaPtx_InDeltaPtyTwoDPlot[TwoDNBinsDeltaPty];	
	TH1D* TrueECalLowPTPlot = new TH1D("TrueECalLowPTPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueECalMidPTPlot = new TH1D("TrueECalMidPTPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueECalHighPTPlot = new TH1D("TrueECalHighPTPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueECal_InDeltaAlphaTTwoDPlot[TwoDNBinsDeltaAlphaT];	

	//--------------------------------------------------//	

	// 3D analysis (uncorrelated)	

	TH1D* TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];
	TH1D* TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];
	TH1D* TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];

	//--------------------------------------------------//	


	for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

		TString DeltaAlphaTTwoDInDeltaPTLabel = "DeltaAlphaT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
		TrueDeltaAlphaT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("True"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);

		TString DeltaPhiTTwoDInDeltaPTLabel = "DeltaPhiT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
		TrueDeltaPhiT_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("True"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);	

		TString DeltaPnTwoDInDeltaPTLabel = "DeltaPn_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
		TrueDeltaPn_InDeltaPTTwoDPlot[WhichDeltaPT] = new TH1D("True"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);	

		for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {	

			TString ECalTwoDInDeltaPTDeltaAlphaTLabel = "ECal_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";
			TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[WhichDeltaPT][WhichDeltaAlphaT] = new TH1D("True"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);

		}

	}

	for (int WhichMuonCosTheta = 0; WhichMuonCosTheta < TwoDNBinsMuonCosTheta; WhichMuonCosTheta++) {

		TString DeltaAlphaTTwoDInMuonCosThetaLabel = "DeltaAlphaT_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
		TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("True"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);

		TString DeltaPTTwoDInMuonCosThetaLabel = "DeltaPT_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
		TrueDeltaPT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("True"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);		

		TString MuonMomentumTwoDInMuonCosThetaLabel = "MuonMomentum_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
		TrueMuonMomentum_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("True"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);

		TString ProtonCosThetaTwoDInMuonCosThetaLabel = "ProtonCosTheta_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
		TrueProtonCosTheta_InMuonCosThetaTwoDPlot[WhichMuonCosTheta] = new TH1D("True"+ProtonCosThetaTwoDInMuonCosThetaLabel,LabelXAxisProtonCosTheta,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);

		for (int WhichMuonMomentum = 0; WhichMuonMomentum < TwoDNBinsMuonMomentum; WhichMuonMomentum++) {	

			TString ECalTwoDInMuonCosThetaMuonMomentumLabel = "ECal_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"_MuonMomentum_"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum+1])+"Plot";
			TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[WhichMuonCosTheta][WhichMuonMomentum] = new TH1D("True"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);

		}

	}		

	for (int WhichProtonCosTheta = 0; WhichProtonCosTheta < TwoDNBinsProtonCosTheta; WhichProtonCosTheta++) {

		TString ProtonMomentumTwoDInProtonCosThetaLabel = "ProtonMomentum_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
		TrueProtonMomentum_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("True"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);

		TString DeltaAlphaTTwoDInProtonCosThetaLabel = "DeltaAlphaT_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
		TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("True"+DeltaAlphaTTwoDInProtonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);

		TString DeltaPTTwoDInProtonCosThetaLabel = "DeltaPT_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
		TrueDeltaPT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta] = new TH1D("True"+DeltaPTTwoDInProtonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);

		for (int WhichProtonMomentum = 0; WhichProtonMomentum < TwoDNBinsProtonMomentum; WhichProtonMomentum++) {	

			TString ECalTwoDInProtonCosThetaProtonMomentumLabel = "ECal_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"_ProtonMomentum_"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum+1])+"Plot";
			TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[WhichProtonCosTheta][WhichProtonMomentum] = new TH1D("True"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);

		}

	}

	for (int WhichDeltaPtx = 0; WhichDeltaPtx < TwoDNBinsDeltaPtx; WhichDeltaPtx++) {

		TString DeltaPtyTwoDInDeltaPtxLabel = "DeltaPty_DeltaPtx_"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx+1])+"Plot";			
		TrueDeltaPty_InDeltaPtxTwoDPlot[WhichDeltaPtx] = new TH1D("True"+DeltaPtyTwoDInDeltaPtxLabel,LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);	

	}	

	for (int WhichDeltaPty = 0; WhichDeltaPty < TwoDNBinsDeltaPty; WhichDeltaPty++) {

		TString DeltaPtxTwoDInDeltaPtyLabel = "DeltaPtx_DeltaPty_"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty+1])+"Plot";			
		TrueDeltaPtx_InDeltaPtyTwoDPlot[WhichDeltaPty] = new TH1D("True"+DeltaPtxTwoDInDeltaPtyLabel,LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);	

	}		

	for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {

		TString DeltaPTTwoDInDeltaAlphaTLabel = "DeltaPT_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
		TrueDeltaPT_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("True"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

		TString ECalTwoDInDeltaAlphaTLabel = "ECal_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
		TrueECal_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT] = new TH1D("True"+ECalTwoDInDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

	}						

	//--------------------------------------------------//

	// 2D analysis in 1D grid	

	TH1D* SerialTrueDeltaPT_InMuonCosThetaPlot = new TH1D("TrueSerialDeltaPT_MuonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);
	TH1D* SerialTrueDeltaPT_InProtonCosThetaPlot = new TH1D("TrueSerialDeltaPT_ProtonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);	
	TH1D* SerialTrueDeltaPT_InDeltaAlphaTPlot = new TH1D("TrueSerialDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
	TH1D* SerialTrueMuonMomentum_InMuonCosThetaPlot = new TH1D("TrueSerialMuonMomentum_MuonCosThetaPlot",LabelXAxisMuonMomentum,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);
	TH1D* SerialTrueProtonMomentum_InProtonCosThetaPlot = new TH1D("TrueSerialProtonMomentum_ProtonCosThetaPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);	
	TH1D* SerialTrueDeltaAlphaT_InMuonCosThetaPlot = new TH1D("TrueSerialDeltaAlphaT_MuonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);
	TH1D* SerialTrueDeltaAlphaT_InProtonCosThetaPlot = new TH1D("TrueSerialDeltaAlphaT_ProtonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);
	TH1D* SerialTrueDeltaAlphaT_InDeltaPTPlot = new TH1D("TrueSerialDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);		
	TH1D* SerialTrueDeltaPhiT_InDeltaPTPlot = new TH1D("TrueSerialDeltaPhiT_DeltaPTPlot",LabelXAxisDeltaPhiT,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);
	TH1D* SerialTrueDeltaPn_InDeltaPTPlot = new TH1D("TrueSerialDeltaPn_DeltaPTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);
	TH1D* SerialTrueProtonCosTheta_InMuonCosThetaPlot = new TH1D("TrueSerialProtonCosTheta_MuonCosThetaPlot",LabelXAxisProtonCosTheta,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);
	TH1D* SerialTrueDeltaPty_InDeltaPtxPlot = new TH1D("TrueSerialDeltaPty_DeltaPtxPlot",LabelXAxisDeltaPty,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);
	TH1D* SerialTrueDeltaPtx_InDeltaPtyPlot = new TH1D("TrueSerialDeltaPtx_DeltaPtyPlot",LabelXAxisDeltaPtx,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);
	TH1D* SerialTrueECal_InDeltaPTPlot = new TH1D("TrueSerialECal_DeltaPTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);
	TH1D* SerialTrueECal_InDeltaAlphaTPlot = new TH1D("TrueSerialECal_DeltaAlphaTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);

	//--------------------------------------------------//	

	// 3D analysis in 1D grid		
	
	TH1D* SerialTrueECal_InDeltaPTDeltaAlphaTPlot = new TH1D("TrueSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);
	TH1D* SerialTrueECal_InMuonCosThetaMuonMomentumPlot = new TH1D("TrueSerialECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);
	TH1D* SerialTrueECal_InProtonCosThetaProtonMomentumPlot = new TH1D("TrueSerialECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);

	//--------------------------------------------------//

	// Counters STV analysis

	int CounterSTVEventsPassedSelection = 0;
	int CounterSTVQEEventsPassedSelection = 0;
	int CounterSTVMECEventsPassedSelection = 0;
	int CounterSTVRESEventsPassedSelection = 0;
	int CounterSTVDISEventsPassedSelection = 0;
	int CounterSTVCOHEventsPassedSelection = 0;		

	//--------------------------------------------------//

	// T2K Weights to make GENIE consistent with the MicroBooNE tune v2 (or v1)
	
	TFile* fweights = nullptr;
	TTree* tweights = nullptr;
	float cv_weight = -99.;

	if (fOutFileName == "Genie_v3_0_6_uB_Tune_1") {

		fweights = TFile::Open("mySamples/myWeights_uB_Tune_v1.root");
		tweights = (TTree*)fweights->Get("ub_tune_cv");
		tweights->SetBranchAddress("cv_weight", &cv_weight);

	}

	if (fOutFileName == "Genie_v3_0_6_Nominal" || fOutFileName == "Genie_v3_0_6_NoFSI" || fOutFileName == "Genie_v3_0_6_NoRPA" || fOutFileName == "Genie_v3_0_6_NoCoulomb" 
	|| fOutFileName == "Genie_v3_0_6_hN2018" || fOutFileName == "Genie_v3_0_6_RFG" || fOutFileName == "Genie_v3_0_6_EffSF") {

		if (fOutFileName == "Genie_v3_0_6_Nominal" || fOutFileName == "Genie_v3_0_6_NoFSI" || fOutFileName == "Genie_v3_0_6_hN2018") 
			{ fweights = TFile::Open("mySamples/myWeights_uB_Tune_Nominal.root"); }
		if (fOutFileName == "Genie_v3_0_6_NoRPA") { fweights = TFile::Open("mySamples/myWeights_uB_Tune_NoRPA.root"); }
		if (fOutFileName == "Genie_v3_0_6_NoCoulomb") { fweights = TFile::Open("mySamples/myWeights_uB_Tune_NoCoulomb.root"); }
		if (fOutFileName == "Genie_v3_0_6_RFG") { fweights = TFile::Open("mySamples/myWeights_uB_Tune_RFG.root"); }
		if (fOutFileName == "Genie_v3_0_6_EffSF") { fweights = TFile::Open("mySamples/myWeights_uB_Tune_EffSF.root"); }

		tweights = (TTree*)fweights->Get("GenericVectors__VARS");
		tweights->SetBranchAddress("Weight", &cv_weight);

	}

	//--------------------------------------------------//

	int SumCVWeights = 0;
	int SumSelectedCVWeights = 0;	

	//--------------------------------------------------//

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		//--------------------------------------------------//

		Long64_t ientry = LoadTree(jentry); if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

		//--------------------------------------------------//

		double weight = 1.;

		if (fOutFileName == "Genie_v3_0_6_uB_Tune_1" || fOutFileName == "Genie_v3_0_6_Nominal" || 
		    fOutFileName == "Genie_v3_0_6_NoFSI" || fOutFileName == "Genie_v3_0_6_hN2018" || 
			fOutFileName == "Genie_v3_0_6_NoRPA" || fOutFileName == "Genie_v3_0_6_NoCoulomb" || 
			fOutFileName == "Genie_v3_0_6_RFG" || fOutFileName == "Genie_v3_0_6_EffSF") 
			{ tweights->GetEntry(jentry); weight = cv_weight; }

		if (weight <= 0 || weight >= 30) { continue; }
		SumCVWeights += weight;
		
		//--------------------------------------------------//

		double genie_mode = -1.;

		if (qel==1) { genie_mode = 1; }
		if (mec==1) { genie_mode = 2; }
		if (res==1) { genie_mode = 3; }
		if (dis==1) { genie_mode = 4; }
		if (coh==1) { genie_mode = 5; }		

		//--------------------------------------------------//

		if (!cc) { continue; }

		//--------------------------------------------------//

		// Muon

		TVector3 Muon3Vector(pxl,pyl,pzl); // GeV
		TLorentzVector Muon4Vector(pxl,pyl,pzl,El); // GeV
		double MuonMomentum = pl; // GeV / c
		double MuonTheta = TMath::ACos(cthl) * 180. / TMath::Pi(); // deg
		double MuonCosTheta = cthl;
		double MuonPhi = Muon4Vector.Phi() * 180. / TMath::Pi(); // deg		

		//--------------------------------------------------//

		// True struck nucleon momentum (LFG / RFG et al)

		double LFG_pn = TMath::Sqrt(pxn*pxn+pyn*pyn+pzn*pzn);	

		int Index = -1;

		for (int i = 0; i < NRanges; i++) {

			// Don't forget the offste by 1, 0 corresponds to all the events
			if (LFG_pn > range[i] && LFG_pn < range[i+1]) { Index = i+1; }

		}			
		
		//--------------------------------------------------//

		int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0, TrueHeavierMesonCounter = 0;
		vector <int> ProtonID; ProtonID.clear();

		//--------------------------------------------------//

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

		//--------------------------------------------------//

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


		//--------------------------------------------------//

		// Impose signal
		// 1 muon > 0.1 GeV/c
		// 1 proton > 0.3 GeV/c
		// no charged pions > 0.07 GeV/c
		// no neutral pions of any momenta
		// any number of neutrons & photons

		if ( ProtonTagging != 1 || ChargedPionTagging != 0 || NeutralPionTagging != 0 ) { continue; }
		if ( TrueHeavierMesonCounter != 0 ) { continue; }

		//--------------------------------------------------//

		// Proton

		int ProtonIndex = ProtonID.at(0);
		TVector3 Proton3Vector(pxf[ProtonIndex],pyf[ProtonIndex],pzf[ProtonIndex]); // GeV
		TLorentzVector Proton4Vector(pxf[ProtonIndex],pyf[ProtonIndex],pzf[ProtonIndex],Ef[ProtonIndex]); // GeV
		double ProtonMomentum = pf[ProtonIndex]; // GeV / c
//		double ProtonTheta = TMath::ACos(cthf[ProtonIndex]) * 180. / TMath::Pi(); // deg
		double ProtonCosTheta = cthf[ProtonIndex];
		double ProtonPhi = Proton4Vector.Phi() * 180. / TMath::Pi(); // deg
		double ProtonE = Proton4Vector.E(); // GeV

		//--------------------------------------------------//

		// Relative Proton - Muon Angles

		double DeltaThetaProtonMuon = Proton3Vector.Angle(Muon3Vector) * 180. / TMath::Pi();
		if (DeltaThetaProtonMuon < 0.) { DeltaThetaProtonMuon += 180.; }
		if (DeltaThetaProtonMuon > 180.) { DeltaThetaProtonMuon -= 180.; }

		double DeltaPhiProtonMuon = Muon3Vector.DeltaPhi(Proton3Vector)* 180. / TMath::Pi();
		if (DeltaPhiProtonMuon < 0.) { DeltaPhiProtonMuon += 360.; }
		if (DeltaPhiProtonMuon > 360.) { DeltaPhiProtonMuon -= 360.; }

		//--------------------------------------------------//

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

		//--------------------------------------------------//
		//--------------------------------------------------//

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

		//--------------------------------------------------//

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

		//--------------------------------------------------//
		//--------------------------------------------------//

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

				//--------------------------------------------------//

				SumSelectedCVWeights += weight;

				CounterSTVEventsPassedSelection++;
				if (qel==1) { CounterSTVQEEventsPassedSelection++; }
				if (mec==1) { CounterSTVMECEventsPassedSelection++; }
				if (res==1) { CounterSTVRESEventsPassedSelection++; }
				if (dis==1) { CounterSTVDISEventsPassedSelection++; }
				if (coh==1) { CounterSTVCOHEventsPassedSelection++; }	

				//----------------------------------------//

				TrueLFGPnSlicePlot[0]->Fill(LFG_pn,weight);				
				TruekMissSlicePlot[0]->Fill(TruekMiss,weight);
				TruePnProxySlicePlot[0]->Fill(TruePn,weight);

				TrueLFGPnSlicePlot[Index]->Fill(LFG_pn,weight);				
				TruekMissSlicePlot[Index]->Fill(TruekMiss,weight);
				TruePnProxySlicePlot[Index]->Fill(TruePn,weight);							

				//--------------------------------------------------//				

				// 1D analysis

				TrueDeltaPTPlot->Fill(PTmissMomentum,weight);
				TrueDeltaAlphaTPlot->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaPhiTPlot->Fill(TrueDeltaPhiT,weight);
				TrueECalPlot->Fill(ECal,weight);
				TrueEQEPlot->Fill(EQE,weight);
				TrueQ2Plot->Fill(TrueQ2,weight);				
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

				//----------------------------------------//

				// Indices for 2D analysis

				int DeltaPTTwoDIndex = tools.ReturnIndex(PTmissMomentum, TwoDArrayNBinsDeltaPT);
				int DeltaAlphaTTwoDIndex = tools.ReturnIndex(TrueDeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
				int MuonCosThetaTwoDIndex = tools.ReturnIndex(MuonCosTheta, TwoDArrayNBinsMuonCosTheta);
				int ProtonCosThetaTwoDIndex = tools.ReturnIndex(ProtonCosTheta, TwoDArrayNBinsProtonCosTheta);
				int DeltaPtxTwoDIndex = tools.ReturnIndex(TruePtx, TwoDArrayNBinsDeltaPtx);	
				int DeltaPtyTwoDIndex = tools.ReturnIndex(TruePty, TwoDArrayNBinsDeltaPty);	
				int MuonMomentumTwoDIndex = tools.ReturnIndex(MuonMomentum, TwoDArrayNBinsMuonMomentum);
				int ProtonMomentumTwoDIndex = tools.ReturnIndex(ProtonMomentum, TwoDArrayNBinsProtonMomentum);	

				// 2D analysis (uncorrelated)

				if (PTmissMomentum > LowPT[0] && PTmissMomentum < HighPT[0]) { TrueECalLowPTPlot->Fill(ECal,weight); }
				if (PTmissMomentum > LowPT[1] && PTmissMomentum < HighPT[1]) { TrueECalMidPTPlot->Fill(ECal,weight); }
				if (PTmissMomentum > LowPT[2] && PTmissMomentum < HighPT[2]) { TrueECalHighPTPlot->Fill(ECal,weight); }				
				TrueDeltaAlphaT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaPhiT_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(TrueDeltaPhiT,weight);	
				TrueDeltaPn_InDeltaPTTwoDPlot[DeltaPTTwoDIndex]->Fill(TruePn,weight);
				TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);				
				TrueDeltaPT_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(PTmissMomentum,weight);
				TrueDeltaPT_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(PTmissMomentum,weight);				
				TrueMuonMomentum_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(MuonMomentum,weight);
				TrueProtonCosTheta_InMuonCosThetaTwoDPlot[MuonCosThetaTwoDIndex]->Fill(ProtonCosTheta,weight);				
				TrueProtonMomentum_InProtonCosThetaTwoDPlot[ProtonCosThetaTwoDIndex]->Fill(ProtonMomentum,weight);
				TrueDeltaPty_InDeltaPtxTwoDPlot[DeltaPtxTwoDIndex]->Fill(TruePty,weight);
				TrueDeltaPtx_InDeltaPtyTwoDPlot[DeltaPtyTwoDIndex]->Fill(TruePtx,weight);	
				TrueDeltaPT_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(PTmissMomentum,weight);
				TrueECal_InDeltaAlphaTTwoDPlot[DeltaAlphaTTwoDIndex]->Fill(ECal,weight);				

				// 2D analysis in 1D grid

				int SerialMuonMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices,MuonCosThetaTwoDIndex,MuonMomentum);
				int SerialProtonMomentumInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices,ProtonCosThetaTwoDIndex,ProtonMomentum);				
				int SerialDeltaAlphaTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices,MuonCosThetaTwoDIndex,TrueDeltaAlphaT);
				int SerialDeltaAlphaTInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices,ProtonCosThetaTwoDIndex,TrueDeltaAlphaT);
				int SerialDeltaAlphaTInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices,DeltaPTTwoDIndex,TrueDeltaAlphaT);
				int SerialDeltaPhiTInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices,DeltaPTTwoDIndex,TrueDeltaPhiT);
				int SerialDeltaPnInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaPTSlices,DeltaPTTwoDIndex,TruePn);	
				int SerialECalInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsECalInDeltaPTSlices,DeltaPTTwoDIndex,ECal);															
				int SerialECalInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsECalInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,ECal);
				int SerialProtonCosThetaInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices,MuonCosThetaTwoDIndex,ProtonCosTheta);
				int SerialDeltaPtyInDeltaPtxIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices,DeltaPtxTwoDIndex,TruePty);				
				int SerialDeltaPtxInDeltaPtyIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPtxInDeltaPtySlices,DeltaPtyTwoDIndex,TruePtx);
				int SerialDeltaPTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices,MuonCosThetaTwoDIndex,PTmissMomentum);
				int SerialDeltaPTInProtonCosThetaIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices,ProtonCosThetaTwoDIndex,PTmissMomentum);
				int SerialDeltaPTInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,PTmissMomentum);

				SerialTrueMuonMomentum_InMuonCosThetaPlot->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
				SerialTrueProtonMomentum_InProtonCosThetaPlot->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
				SerialTrueDeltaAlphaT_InMuonCosThetaPlot->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
				SerialTrueDeltaAlphaT_InProtonCosThetaPlot->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
				SerialTrueDeltaAlphaT_InDeltaPTPlot->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
				SerialTrueDeltaPhiT_InDeltaPTPlot->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
				SerialTrueDeltaPn_InDeltaPTPlot->Fill(SerialDeltaPnInDeltaPTIndex,weight);	
				SerialTrueECal_InDeltaPTPlot->Fill(SerialECalInDeltaPTIndex,weight);																							
				SerialTrueProtonCosTheta_InMuonCosThetaPlot->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
				SerialTrueDeltaPty_InDeltaPtxPlot->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
				SerialTrueDeltaPtx_InDeltaPtyPlot->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
				SerialTrueDeltaPT_InMuonCosThetaPlot->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
				SerialTrueDeltaPT_InProtonCosThetaPlot->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);
				SerialTrueDeltaPT_InDeltaAlphaTPlot->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);																			
				SerialTrueECal_InDeltaAlphaTPlot->Fill(SerialECalInDeltaAlphaTIndex,weight);

				//----------------------------------------//				

				// 3D analysis (uncorrelated)

				TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
				TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
				TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);								

				// 3D analysis treated in 1D grid

				int SerialECalInDeltaPTDeltaAlphaTIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices,DeltaPTTwoDIndex,DeltaAlphaTTwoDIndex,ECal);
				int SerialECalInMuonCosThetaMuonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices,MuonCosThetaTwoDIndex,MuonMomentumTwoDIndex,ECal);
				int SerialECalInProtonCosThetaProtonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices,ProtonCosThetaTwoDIndex,ProtonMomentumTwoDIndex,ECal);																																																									

				SerialTrueECal_InDeltaPTDeltaAlphaTPlot->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
				SerialTrueECal_InMuonCosThetaMuonMomentumPlot->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
				SerialTrueECal_InProtonCosThetaProtonMomentumPlot->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);				

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

	if (fOutFileName == "Genie_v3_0_6_uB_Tune_1" || fOutFileName == "Genie_v3_0_6_Out_Of_The_Box" || fOutFileName == "Genie_v3_0_6_Nominal" 
	 || fOutFileName == "Genie_v3_0_6_NoFSI"  || fOutFileName == "Genie_v3_0_6_hN2018") 
		{ ScalingFactor = ScalingFactor * G18_10a_02_11a_FluxIntegratedXSection ; }
	if (fOutFileName == "SuSav2") { ScalingFactor = ScalingFactor * SuSav2FluxIntegratedXSection ; }
	if (fOutFileName == "GENIEv2") { ScalingFactor = ScalingFactor * R_2_12_10_FluxIntegratedXSection ; }
	if (fOutFileName == "GENIEv3_0_4") { ScalingFactor = ScalingFactor * R_3_0_4_FluxIntegratedXSection ; }
	if (fOutFileName == "Genie_v3_0_6_NoRPA") { ScalingFactor = ScalingFactor * R_3_0_6_G18_10a_02_11a_NoRPA_FluxIntegratedXSection ; }
	if (fOutFileName == "Genie_v3_0_6_NoCoulomb") { ScalingFactor = ScalingFactor * R_3_0_6_G18_10a_02_11a_NoCoulomb_FluxIntegratedXSection ; }
	if (fOutFileName == "Genie_v3_0_6_RFG") { ScalingFactor = ScalingFactor * R_3_0_6_G18_10a_02_11a_RFG_FluxIntegratedXSection ; }
	if (fOutFileName == "Genie_v3_0_6_EffSF") { ScalingFactor = ScalingFactor * R_3_0_6_G18_10a_02_11a_EffSF_FluxIntegratedXSection ; }

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

	//----------------------------------------//

	// 2D analysis

	for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

		tools.Reweight(TrueDeltaAlphaT_InDeltaPTTwoDPlot[WhichDeltaPT],ScalingFactor);
		tools.Reweight(TrueDeltaPhiT_InDeltaPTTwoDPlot[WhichDeltaPT],ScalingFactor);	
		tools.Reweight(TrueDeltaPn_InDeltaPTTwoDPlot[WhichDeltaPT],ScalingFactor);	

		for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {	

					tools.Reweight(TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[WhichDeltaPT][WhichDeltaAlphaT],ScalingFactor);

		}

	}

	for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {

		tools.Reweight(TrueDeltaPT_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT],ScalingFactor);
		tools.Reweight(TrueECal_InDeltaAlphaTTwoDPlot[WhichDeltaAlphaT],ScalingFactor);		

	}	

	for (int WhichMuonCosTheta = 0; WhichMuonCosTheta < TwoDNBinsMuonCosTheta; WhichMuonCosTheta++) {

		tools.Reweight(TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta],ScalingFactor);
		tools.Reweight(TrueDeltaPT_InMuonCosThetaTwoDPlot[WhichMuonCosTheta],ScalingFactor);		
		tools.Reweight(TrueMuonMomentum_InMuonCosThetaTwoDPlot[WhichMuonCosTheta],ScalingFactor);
		tools.Reweight(TrueProtonCosTheta_InMuonCosThetaTwoDPlot[WhichMuonCosTheta],ScalingFactor);		

		for (int WhichMuonMomentum = 0; WhichMuonMomentum < TwoDNBinsMuonMomentum; WhichMuonMomentum++) {	

					tools.Reweight(TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[WhichMuonCosTheta][WhichMuonMomentum],ScalingFactor);

		}		

	}	

	for (int WhichProtonCosTheta = 0; WhichProtonCosTheta < TwoDNBinsProtonCosTheta; WhichProtonCosTheta++) {

		tools.Reweight(TrueProtonMomentum_InProtonCosThetaTwoDPlot[WhichProtonCosTheta],ScalingFactor);	
		tools.Reweight(TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta],ScalingFactor);
		tools.Reweight(TrueDeltaPT_InProtonCosThetaTwoDPlot[WhichProtonCosTheta],ScalingFactor);				

		for (int WhichProtonMomentum = 0; WhichProtonMomentum < TwoDNBinsProtonMomentum; WhichProtonMomentum++) {	

					tools.Reweight(TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[WhichProtonCosTheta][WhichProtonMomentum],ScalingFactor);

		}		

	}

	for (int WhichDeltaPtx = 0; WhichDeltaPtx < TwoDNBinsDeltaPtx; WhichDeltaPtx++) {

		tools.Reweight(TrueDeltaPty_InDeltaPtxTwoDPlot[WhichDeltaPtx],ScalingFactor);	

	}

	for (int WhichDeltaPty = 0; WhichDeltaPty < TwoDNBinsDeltaPty; WhichDeltaPty++) {

		tools.Reweight(TrueDeltaPtx_InDeltaPtyTwoDPlot[WhichDeltaPty],ScalingFactor);	

	}	

	//----------------------------------------//

	file->cd();
	file->Write();

	//----------------------------------------//

	std::cout << std::endl;
	std::cout << "File " << FileNameAndPath +" has been created created " << std::endl; 
	std::cout << std::endl;	

	fFile->Close();

	//----------------------------------------//

} // End of the program
