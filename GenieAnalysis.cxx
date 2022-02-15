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

	//----------------------------------------//

	// 1D analysis

	TH1D* TrueDeltaPTPlot[NInte];
	TH1D* TrueDeltaAlphaTPlot[NInte];
	TH1D* TrueDeltaPhiTPlot[NInte];
	TH1D* TrueDeltaPLPlot[NInte];
	TH1D* TrueDeltaPnPlot[NInte];
	TH1D* TrueDeltaPtxPlot[NInte];
	TH1D* TrueDeltaPtyPlot[NInte];
	TH1D* TrueAPlot[NInte];
	TH1D* TrueMuonMomentumPlot[NInte];
	TH1D* TrueMuonPhiPlot[NInte];
	TH1D* TrueMuonCosThetaPlot[NInte];
	TH1D* TrueMuonCosThetaSingleBinPlot[NInte];
	TH1D* TrueProtonMomentumPlot[NInte];
	TH1D* TrueProtonPhiPlot[NInte];
	TH1D* TrueProtonCosThetaPlot[NInte];
	TH1D* TrueECalPlot[NInte];
	TH1D* TrueEQEPlot[NInte];
	TH1D* TrueQ2Plot[NInte];
	TH1D* TrueCCQEMuonMomentumPlot[NInte];
	TH1D* TrueCCQEMuonPhiPlot[NInte];
	TH1D* TrueCCQEMuonCosThetaPlot[NInte];	
	TH1D* TrueCCQEProtonMomentumPlot[NInte];
	TH1D* TrueCCQEProtonPhiPlot[NInte];
	TH1D* TrueCCQEProtonCosThetaPlot[NInte];
	TH1D* TrueCCQEECalPlot[NInte];
	TH1D* TrueCCQEQ2Plot[NInte];	
	TH1D* TruekMissPlot[NInte];
	TH1D* TruePMissMinusPlot[NInte];
	TH1D* TruePMissPlot[NInte];

	//--------------------------------------------------//

	// 2D analysis (uncorrelated)

	TH1D* TrueDeltaPT_InMuonCosThetaTwoDPlot[NInte][TwoDNBinsMuonCosTheta];
	TH1D* TrueDeltaPT_InProtonCosThetaTwoDPlot[NInte][TwoDNBinsProtonCosTheta];	
	TH1D* TrueDeltaPT_InDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaAlphaT];
	TH1D* TrueMuonMomentum_InMuonCosThetaTwoDPlot[NInte][TwoDNBinsMuonCosTheta];
	TH1D* TrueProtonMomentum_InProtonCosThetaTwoDPlot[NInte][TwoDNBinsProtonCosTheta];			
	TH1D* TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[NInte][TwoDNBinsMuonCosTheta];
	TH1D* TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[NInte][TwoDNBinsProtonCosTheta];		
	TH1D* TrueDeltaAlphaT_InDeltaPTTwoDPlot[NInte][TwoDNBinsDeltaPT];			
	TH1D* TrueProtonCosTheta_InMuonCosThetaTwoDPlot[NInte][TwoDNBinsProtonCosTheta];	
	TH1D* TrueDeltaPhiT_InDeltaPTTwoDPlot[NInte][TwoDNBinsDeltaPT];	
	TH1D* TrueDeltaPn_InDeltaPTTwoDPlot[NInte][TwoDNBinsDeltaPT];
	TH1D* TrueDeltaPn_InDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaAlphaT];	
	TH1D* TrueDeltaPty_InDeltaPtxTwoDPlot[NInte][TwoDNBinsDeltaPtx];
	TH1D* TrueDeltaPtx_InDeltaPtyTwoDPlot[NInte][TwoDNBinsDeltaPty];	
	TH1D* TrueECal_InDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaAlphaT];	
	TH1D* TrueECal_InDeltaPTTwoDPlot[NInte][TwoDNBinsDeltaPT];

	//--------------------------------------------------//	

	// 3D analysis (uncorrelated)	

	TH1D* TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];
	TH1D* TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[NInte][TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];
	TH1D* TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[NInte][TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];	

	//--------------------------------------------------//

	// 2D analysis in 1D grid	

	TH1D* SerialTrueDeltaPT_InMuonCosThetaPlot[NInte];
	TH1D* SerialTrueDeltaPT_InProtonCosThetaPlot[NInte];
	TH1D* SerialTrueDeltaPT_InDeltaAlphaTPlot[NInte];
	TH1D* SerialTrueMuonMomentum_InMuonCosThetaPlot[NInte];
	TH1D* SerialTrueProtonMomentum_InProtonCosThetaPlot[NInte];	
	TH1D* SerialTrueDeltaAlphaT_InMuonCosThetaPlot[NInte];
	TH1D* SerialTrueDeltaAlphaT_InProtonCosThetaPlot[NInte];
	TH1D* SerialTrueDeltaAlphaT_InDeltaPTPlot[NInte];		
	TH1D* SerialTrueDeltaPhiT_InDeltaPTPlot[NInte];
	TH1D* SerialTrueDeltaPn_InDeltaPTPlot[NInte];
	TH1D* SerialTrueDeltaPn_InDeltaAlphaTPlot[NInte];	
	TH1D* SerialTrueProtonCosTheta_InMuonCosThetaPlot[NInte];
	TH1D* SerialTrueDeltaPty_InDeltaPtxPlot[NInte];
	TH1D* SerialTrueDeltaPtx_InDeltaPtyPlot[NInte];
	TH1D* SerialTrueECal_InDeltaPTPlot[NInte];
	TH1D* SerialTrueECal_InDeltaAlphaTPlot[NInte];

	//--------------------------------------------------//	

	// 3D analysis in 1D grid		
	
	TH1D* SerialTrueECal_InDeltaPTDeltaAlphaTPlot[NInte];
	TH1D* SerialTrueECal_InMuonCosThetaMuonMomentumPlot[NInte];
	TH1D* SerialTrueECal_InProtonCosThetaProtonMomentumPlot[NInte];

	//--------------------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

		//--------------------------------------------------//

		// 1D analysis

		TrueDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
		TrueDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
		TrueDeltaPhiTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
		TrueDeltaPLPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPLPlot",LabelXAxisDeltaPL,NBinsDeltaPL,ArrayNBinsDeltaPL);
		TrueDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
		TrueDeltaPtxPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
		TrueDeltaPtyPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
		TrueAPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueAPlot",LabelXAxisA,NBinsA,ArrayNBinsA);
		TrueMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
		TrueMuonPhiPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
		TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
		TrueMuonCosThetaSingleBinPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaSingleBinPlot",LabelXAxisMuonCosTheta,1,-1.,1.);
		TrueProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
		TrueProtonPhiPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
		TrueProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
		TrueECalPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
		TrueEQEPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
		TrueQ2Plot[inte] = new TH1D(InteractionLabels[inte]+"TrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);
		TrueCCQEMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEMuonMomentumPlot",LabelXAxisMuonMomentum,CCQENBinsMuonMomentum,CCQEArrayNBinsMuonMomentum);
		TrueCCQEMuonPhiPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEMuonPhiPlot",LabelXAxisMuonPhi,CCQENBinsMuonPhi,CCQEArrayNBinsMuonPhi);
		TrueCCQEMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEMuonCosThetaPlot",LabelXAxisMuonCosTheta,CCQENBinsMuonCosTheta,CCQEArrayNBinsMuonCosTheta);	
		TrueCCQEProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEProtonMomentumPlot",LabelXAxisProtonMomentum,CCQENBinsProtonMomentum,CCQEArrayNBinsProtonMomentum);
		TrueCCQEProtonPhiPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEProtonPhiPlot",LabelXAxisProtonPhi,CCQENBinsProtonPhi,CCQEArrayNBinsProtonPhi);
		TrueCCQEProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEProtonCosThetaPlot",LabelXAxisProtonCosTheta,CCQENBinsProtonCosTheta,CCQEArrayNBinsProtonCosTheta);
		TrueCCQEECalPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEECalPlot",LabelXAxisECal,CCQENBinsECal,CCQEArrayNBinsECal);
		TrueCCQEQ2Plot[inte] = new TH1D(InteractionLabels[inte]+"TrueCCQEQ2Plot",LabelXAxisQ2,CCQENBinsQ2,CCQEArrayNBinsQ2);	
		TruekMissPlot[inte] = new TH1D(InteractionLabels[inte]+"TruekMissPlot",LabelXAxiskMiss,NBinskMiss,ArrayNBinskMiss);
		TruePMissMinusPlot[inte] = new TH1D(InteractionLabels[inte]+"TruePMissMinusPlot",LabelXAxisPMissMinus,NBinsPMissMinus,ArrayNBinsPMissMinus);
		TruePMissPlot[inte] = new TH1D(InteractionLabels[inte]+"TruePMissPlot",LabelXAxisPMiss,NBinsPMiss,ArrayNBinsPMiss);

		//--------------------------------------------------//	

		// 2D & 3D analysis (uncorrelated)


		for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

			TString DeltaAlphaTTwoDInDeltaPTLabel = "DeltaAlphaT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
			TrueDeltaAlphaT_InDeltaPTTwoDPlot[inte][WhichDeltaPT] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);

			TString ECalTwoDInDeltaPTLabel = "ECal_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
			TrueECal_InDeltaPTTwoDPlot[inte][WhichDeltaPT] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaPTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsECalInDeltaPTSlices[WhichDeltaPT][0]);		


			TString DeltaPhiTTwoDInDeltaPTLabel = "DeltaPhiT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
			TrueDeltaPhiT_InDeltaPTTwoDPlot[inte][WhichDeltaPT] = new TH1D(InteractionLabels[inte]+"True"+DeltaPhiTTwoDInDeltaPTLabel,LabelXAxisDeltaPhiT,TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPhiTInDeltaPTSlices[WhichDeltaPT][0]);	

			TString DeltaPnTwoDInDeltaPTLabel = "DeltaPn_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";			
			TrueDeltaPn_InDeltaPTTwoDPlot[inte][WhichDeltaPT] = new TH1D(InteractionLabels[inte]+"True"+DeltaPnTwoDInDeltaPTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaPTSlices[WhichDeltaPT][0]);	

			for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {	

				TString ECalTwoDInDeltaPTDeltaAlphaTLabel = "ECal_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";
				TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[inte][WhichDeltaPT][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);

			}

		}

		for (int WhichMuonCosTheta = 0; WhichMuonCosTheta < TwoDNBinsMuonCosTheta; WhichMuonCosTheta++) {

			TString DeltaAlphaTTwoDInMuonCosThetaLabel = "DeltaAlphaT_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
			TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[inte][WhichMuonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlphaTTwoDInMuonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices[WhichMuonCosTheta][0]);

			TString DeltaPTTwoDInMuonCosThetaLabel = "DeltaPT_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
			TrueDeltaPT_InMuonCosThetaTwoDPlot[inte][WhichMuonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+DeltaPTTwoDInMuonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInMuonCosThetaSlices[WhichMuonCosTheta][0]);		

			TString MuonMomentumTwoDInMuonCosThetaLabel = "MuonMomentum_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
			TrueMuonMomentum_InMuonCosThetaTwoDPlot[inte][WhichMuonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+MuonMomentumTwoDInMuonCosThetaLabel,LabelXAxisMuonMomentum,TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices[WhichMuonCosTheta][0]);

			TString ProtonCosThetaTwoDInMuonCosThetaLabel = "ProtonCosTheta_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"Plot";			
			TrueProtonCosTheta_InMuonCosThetaTwoDPlot[inte][WhichMuonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+ProtonCosThetaTwoDInMuonCosThetaLabel,LabelXAxisProtonCosTheta,TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta].size()-1,&TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices[WhichMuonCosTheta][0]);

			for (int WhichMuonMomentum = 0; WhichMuonMomentum < TwoDNBinsMuonMomentum; WhichMuonMomentum++) {	

				TString ECalTwoDInMuonCosThetaMuonMomentumLabel = "ECal_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"_MuonMomentum_"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum+1])+"Plot";
				TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[inte][WhichMuonCosTheta][WhichMuonMomentum] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);

			}

		}		

		for (int WhichProtonCosTheta = 0; WhichProtonCosTheta < TwoDNBinsProtonCosTheta; WhichProtonCosTheta++) {

			TString ProtonMomentumTwoDInProtonCosThetaLabel = "ProtonMomentum_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
			TrueProtonMomentum_InProtonCosThetaTwoDPlot[inte][WhichProtonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+ProtonMomentumTwoDInProtonCosThetaLabel,LabelXAxisProtonMomentum,TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices[WhichProtonCosTheta][0]);

			TString DeltaAlphaTTwoDInProtonCosThetaLabel = "DeltaAlphaT_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
			TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[inte][WhichProtonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlphaTTwoDInProtonCosThetaLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices[WhichProtonCosTheta][0]);

			TString DeltaPTTwoDInProtonCosThetaLabel = "DeltaPT_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"Plot";			
			TrueDeltaPT_InProtonCosThetaTwoDPlot[inte][WhichProtonCosTheta] = new TH1D(InteractionLabels[inte]+"True"+DeltaPTTwoDInProtonCosThetaLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta].size()-1,&TwoDArrayNBinsDeltaPTInProtonCosThetaSlices[WhichProtonCosTheta][0]);

			for (int WhichProtonMomentum = 0; WhichProtonMomentum < TwoDNBinsProtonMomentum; WhichProtonMomentum++) {	

				TString ECalTwoDInProtonCosThetaProtonMomentumLabel = "ECal_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"_ProtonMomentum_"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum+1])+"Plot";
				TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[inte][WhichProtonCosTheta][WhichProtonMomentum] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);

			}

		}

		for (int WhichDeltaPtx = 0; WhichDeltaPtx < TwoDNBinsDeltaPtx; WhichDeltaPtx++) {

			TString DeltaPtyTwoDInDeltaPtxLabel = "DeltaPty_DeltaPtx_"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx+1])+"Plot";			
			TrueDeltaPty_InDeltaPtxTwoDPlot[inte][WhichDeltaPtx] = new TH1D(InteractionLabels[inte]+"True"+DeltaPtyTwoDInDeltaPtxLabel,LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);	

		}	

		for (int WhichDeltaPty = 0; WhichDeltaPty < TwoDNBinsDeltaPty; WhichDeltaPty++) {

			TString DeltaPtxTwoDInDeltaPtyLabel = "DeltaPtx_DeltaPty_"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty+1])+"Plot";			
			TrueDeltaPtx_InDeltaPtyTwoDPlot[inte][WhichDeltaPty] = new TH1D(InteractionLabels[inte]+"True"+DeltaPtxTwoDInDeltaPtyLabel,LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);	

		}		

		for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {

			TString DeltaPTTwoDInDeltaAlphaTLabel = "DeltaPT_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
			TrueDeltaPT_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"True"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

			TString DeltaPnTwoDInDeltaAlphaTLabel = "DeltaPn_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
			TrueDeltaPn_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"True"+DeltaPnTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);			

			TString ECalTwoDInDeltaAlphaTLabel = "ECal_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";			
			TrueECal_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

		}

		//--------------------------------------------------//	

		// 2D analysis in 1D grid	

		SerialTrueDeltaPT_InMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPT_MuonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInMuonCosThetaSlices)[0]);
		SerialTrueDeltaPT_InProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPT_ProtonCosThetaPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInProtonCosThetaSlices)[0]);	
		SerialTrueDeltaPT_InDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
		SerialTrueMuonMomentum_InMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialMuonMomentum_MuonCosThetaPlot",LabelXAxisMuonMomentum,tools.Return2DNBins(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices)[0]);
		SerialTrueProtonMomentum_InProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialProtonMomentum_ProtonCosThetaPlot",LabelXAxisProtonMomentum,tools.Return2DNBins(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices)[0]);	
		SerialTrueDeltaAlphaT_InMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaAlphaT_MuonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]);
		SerialTrueDeltaAlphaT_InProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaAlphaT_ProtonCosThetaPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices)[0]);
		SerialTrueDeltaAlphaT_InDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);		
		SerialTrueDeltaPhiT_InDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPhiT_DeltaPTPlot",LabelXAxisDeltaPhiT,tools.Return2DNBins(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPhiTInDeltaPTSlices)[0]);
		SerialTrueDeltaPn_InDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPn_DeltaPTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaPTSlices)[0]);
		SerialTrueDeltaPn_InDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPn_DeltaAlphaTPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices)[0]);
		SerialTrueProtonCosTheta_InMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialProtonCosTheta_MuonCosThetaPlot",LabelXAxisProtonCosTheta,tools.Return2DNBins(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices),&tools.Return2DBinIndices(TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices)[0]);
		SerialTrueDeltaPty_InDeltaPtxPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPty_DeltaPtxPlot",LabelXAxisDeltaPty,tools.Return2DNBins(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtyInDeltaPtxSlices)[0]);
		SerialTrueDeltaPtx_InDeltaPtyPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialDeltaPtx_DeltaPtyPlot",LabelXAxisDeltaPtx,tools.Return2DNBins(TwoDArrayNBinsDeltaPtxInDeltaPtySlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPtxInDeltaPtySlices)[0]);
		SerialTrueECal_InDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_DeltaPTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaPTSlices)[0]);
		SerialTrueECal_InDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_DeltaAlphaTPlot",LabelXAxisECal,tools.Return2DNBins(TwoDArrayNBinsECalInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsECalInDeltaAlphaTSlices)[0]);

		//--------------------------------------------------//	

		// 3D analysis in 1D grid		
		
		SerialTrueECal_InDeltaPTDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);
		SerialTrueECal_InMuonCosThetaMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);
		SerialTrueECal_InProtonCosThetaProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);

		//--------------------------------------------------//

	} // End of the loop over the interaction processes						

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
				int genie_mode = -1.;
				if (qel==1) { CounterSTVQEEventsPassedSelection++; genie_mode = 1; }
				if (mec==1) { CounterSTVMECEventsPassedSelection++; genie_mode = 2; }
				if (res==1) { CounterSTVRESEventsPassedSelection++; genie_mode = 3; }
				if (dis==1) { CounterSTVDISEventsPassedSelection++; genie_mode = 4; }
				if (coh==1) { CounterSTVCOHEventsPassedSelection++; genie_mode = 5; }	

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
				int SerialDeltaPnInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,TruePn);				

				int SerialECalInDeltaPTDeltaAlphaTIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices,DeltaPTTwoDIndex,DeltaAlphaTTwoDIndex,ECal);
				int SerialECalInMuonCosThetaMuonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices,MuonCosThetaTwoDIndex,MuonMomentumTwoDIndex,ECal);
				int SerialECalInProtonCosThetaProtonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices,ProtonCosThetaTwoDIndex,ProtonMomentumTwoDIndex,ECal);																																																									

				//----------------------------------------//

				TrueLFGPnSlicePlot[0]->Fill(LFG_pn,weight);				
				TruekMissSlicePlot[0]->Fill(TruekMiss,weight);
				TruePnProxySlicePlot[0]->Fill(TruePn,weight);

				TrueLFGPnSlicePlot[Index]->Fill(LFG_pn,weight);				
				TruekMissSlicePlot[Index]->Fill(TruekMiss,weight);
				TruePnProxySlicePlot[Index]->Fill(TruePn,weight);							

				//--------------------------------------------------//				

				// 1D analysis

				TrueDeltaPTPlot[0]->Fill(PTmissMomentum,weight);
				TrueDeltaAlphaTPlot[0]->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaPhiTPlot[0]->Fill(TrueDeltaPhiT,weight);
				TrueECalPlot[0]->Fill(ECal,weight);
				TrueEQEPlot[0]->Fill(EQE,weight);
				TrueQ2Plot[0]->Fill(TrueQ2,weight);				
				TrueMuonMomentumPlot[0]->Fill(MuonMomentum,weight);
				TrueMuonPhiPlot[0]->Fill(MuonPhi,weight);
				TrueMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
				TrueMuonCosThetaSingleBinPlot[0]->Fill(MuonCosTheta,weight);
				TrueProtonMomentumPlot[0]->Fill(ProtonMomentum,weight);
				TrueProtonPhiPlot[0]->Fill(ProtonPhi,weight);
				TrueProtonCosThetaPlot[0]->Fill(ProtonCosTheta,weight);			
				TruekMissPlot[0]->Fill(TruekMiss,weight);
				TruePMissMinusPlot[0]->Fill(TruePMissMinus,weight);
				TruePMissPlot[0]->Fill(TrueMissMomentum,weight);
				TrueDeltaPLPlot[0]->Fill(TruePL,weight);
				TrueDeltaPnPlot[0]->Fill(TruePn,weight);
				TrueDeltaPtxPlot[0]->Fill(TruePtx,weight);
				TrueDeltaPtyPlot[0]->Fill(TruePty,weight);
				TrueAPlot[0]->Fill(TrueA,weight);

				TrueDeltaPTPlot[genie_mode]->Fill(PTmissMomentum,weight);
				TrueDeltaAlphaTPlot[genie_mode]->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaPhiTPlot[genie_mode]->Fill(TrueDeltaPhiT,weight);
				TrueECalPlot[genie_mode]->Fill(ECal,weight);
				TrueEQEPlot[genie_mode]->Fill(EQE,weight);
				TrueQ2Plot[genie_mode]->Fill(TrueQ2,weight);				
				TrueMuonMomentumPlot[genie_mode]->Fill(MuonMomentum,weight);
				TrueMuonPhiPlot[genie_mode]->Fill(MuonPhi,weight);
				TrueMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
				TrueMuonCosThetaSingleBinPlot[genie_mode]->Fill(MuonCosTheta,weight);
				TrueProtonMomentumPlot[genie_mode]->Fill(ProtonMomentum,weight);
				TrueProtonPhiPlot[genie_mode]->Fill(ProtonPhi,weight);
				TrueProtonCosThetaPlot[genie_mode]->Fill(ProtonCosTheta,weight);			
				TruekMissPlot[genie_mode]->Fill(TruekMiss,weight);
				TruePMissMinusPlot[genie_mode]->Fill(TruePMissMinus,weight);
				TruePMissPlot[genie_mode]->Fill(TrueMissMomentum,weight);
				TrueDeltaPLPlot[genie_mode]->Fill(TruePL,weight);
				TrueDeltaPnPlot[genie_mode]->Fill(TruePn,weight);
				TrueDeltaPtxPlot[genie_mode]->Fill(TruePtx,weight);
				TrueDeltaPtyPlot[genie_mode]->Fill(TruePty,weight);
				TrueAPlot[genie_mode]->Fill(TrueA,weight);					

				//----------------------------------------//

				// 2D analysis (uncorrelated)
			
				TrueDeltaAlphaT_InDeltaPTTwoDPlot[0][DeltaPTTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaPhiT_InDeltaPTTwoDPlot[0][DeltaPTTwoDIndex]->Fill(TrueDeltaPhiT,weight);	
				TrueDeltaPn_InDeltaPTTwoDPlot[0][DeltaPTTwoDIndex]->Fill(TruePn,weight);
				TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[0][MuonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[0][ProtonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);				
				TrueDeltaPT_InMuonCosThetaTwoDPlot[0][MuonCosThetaTwoDIndex]->Fill(PTmissMomentum,weight);
				TrueDeltaPT_InProtonCosThetaTwoDPlot[0][ProtonCosThetaTwoDIndex]->Fill(PTmissMomentum,weight);				
				TrueMuonMomentum_InMuonCosThetaTwoDPlot[0][MuonCosThetaTwoDIndex]->Fill(MuonMomentum,weight);
				TrueProtonCosTheta_InMuonCosThetaTwoDPlot[0][MuonCosThetaTwoDIndex]->Fill(ProtonCosTheta,weight);				
				TrueProtonMomentum_InProtonCosThetaTwoDPlot[0][ProtonCosThetaTwoDIndex]->Fill(ProtonMomentum,weight);
				TrueDeltaPty_InDeltaPtxTwoDPlot[0][DeltaPtxTwoDIndex]->Fill(TruePty,weight);
				TrueDeltaPtx_InDeltaPtyTwoDPlot[0][DeltaPtyTwoDIndex]->Fill(TruePtx,weight);	
				TrueDeltaPT_InDeltaAlphaTTwoDPlot[0][DeltaAlphaTTwoDIndex]->Fill(PTmissMomentum,weight);
				TrueDeltaPn_InDeltaAlphaTTwoDPlot[0][DeltaAlphaTTwoDIndex]->Fill(TruePn,weight);				
				TrueECal_InDeltaAlphaTTwoDPlot[0][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
				TrueECal_InDeltaPTTwoDPlot[0][DeltaPTTwoDIndex]->Fill(ECal,weight);

				TrueDeltaAlphaT_InDeltaPTTwoDPlot[genie_mode][DeltaPTTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaPhiT_InDeltaPTTwoDPlot[genie_mode][DeltaPTTwoDIndex]->Fill(TrueDeltaPhiT,weight);	
				TrueDeltaPn_InDeltaPTTwoDPlot[genie_mode][DeltaPTTwoDIndex]->Fill(TruePn,weight);
				TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[genie_mode][MuonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[genie_mode][ProtonCosThetaTwoDIndex]->Fill(TrueDeltaAlphaT,weight);				
				TrueDeltaPT_InMuonCosThetaTwoDPlot[genie_mode][MuonCosThetaTwoDIndex]->Fill(PTmissMomentum,weight);
				TrueDeltaPT_InProtonCosThetaTwoDPlot[genie_mode][ProtonCosThetaTwoDIndex]->Fill(PTmissMomentum,weight);				
				TrueMuonMomentum_InMuonCosThetaTwoDPlot[genie_mode][MuonCosThetaTwoDIndex]->Fill(MuonMomentum,weight);
				TrueProtonCosTheta_InMuonCosThetaTwoDPlot[genie_mode][MuonCosThetaTwoDIndex]->Fill(ProtonCosTheta,weight);				
				TrueProtonMomentum_InProtonCosThetaTwoDPlot[genie_mode][ProtonCosThetaTwoDIndex]->Fill(ProtonMomentum,weight);
				TrueDeltaPty_InDeltaPtxTwoDPlot[genie_mode][DeltaPtxTwoDIndex]->Fill(TruePty,weight);
				TrueDeltaPtx_InDeltaPtyTwoDPlot[genie_mode][DeltaPtyTwoDIndex]->Fill(TruePtx,weight);	
				TrueDeltaPT_InDeltaAlphaTTwoDPlot[genie_mode][DeltaAlphaTTwoDIndex]->Fill(PTmissMomentum,weight);
				TrueDeltaPn_InDeltaAlphaTTwoDPlot[genie_mode][DeltaAlphaTTwoDIndex]->Fill(TruePn,weight);				
				TrueECal_InDeltaAlphaTTwoDPlot[genie_mode][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
				TrueECal_InDeltaPTTwoDPlot[genie_mode][DeltaPTTwoDIndex]->Fill(ECal,weight);

				//----------------------------------------//

				// 3D analysis (uncorrelated)

				TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[0][DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
				TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[0][MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
				TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[0][ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);

				TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[genie_mode][DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
				TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[genie_mode][MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
				TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[genie_mode][ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);									

				//----------------------------------------//															

				// 2D analysis in 1D grid

				SerialTrueMuonMomentum_InMuonCosThetaPlot[0]->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
				SerialTrueProtonMomentum_InProtonCosThetaPlot[0]->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
				SerialTrueDeltaAlphaT_InMuonCosThetaPlot[0]->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
				SerialTrueDeltaAlphaT_InProtonCosThetaPlot[0]->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
				SerialTrueDeltaAlphaT_InDeltaPTPlot[0]->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
				SerialTrueDeltaPhiT_InDeltaPTPlot[0]->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
				SerialTrueDeltaPn_InDeltaPTPlot[0]->Fill(SerialDeltaPnInDeltaPTIndex,weight);	
				SerialTrueECal_InDeltaPTPlot[0]->Fill(SerialECalInDeltaPTIndex,weight);																							
				SerialTrueProtonCosTheta_InMuonCosThetaPlot[0]->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
				SerialTrueDeltaPty_InDeltaPtxPlot[0]->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
				SerialTrueDeltaPtx_InDeltaPtyPlot[0]->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
				SerialTrueDeltaPT_InMuonCosThetaPlot[0]->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
				SerialTrueDeltaPT_InProtonCosThetaPlot[0]->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);
				SerialTrueDeltaPT_InDeltaAlphaTPlot[0]->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
				SerialTrueDeltaPn_InDeltaAlphaTPlot[0]->Fill(SerialDeltaPnInDeltaAlphaTIndex,weight);																							
				SerialTrueECal_InDeltaAlphaTPlot[0]->Fill(SerialECalInDeltaAlphaTIndex,weight);

				SerialTrueMuonMomentum_InMuonCosThetaPlot[genie_mode]->Fill(SerialMuonMomentumInMuonCosThetaIndex,weight);			
				SerialTrueProtonMomentum_InProtonCosThetaPlot[genie_mode]->Fill(SerialProtonMomentumInProtonCosThetaIndex,weight);
				SerialTrueDeltaAlphaT_InMuonCosThetaPlot[genie_mode]->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
				SerialTrueDeltaAlphaT_InProtonCosThetaPlot[genie_mode]->Fill(SerialDeltaAlphaTInProtonCosThetaIndex,weight);
				SerialTrueDeltaAlphaT_InDeltaPTPlot[genie_mode]->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
				SerialTrueDeltaPhiT_InDeltaPTPlot[genie_mode]->Fill(SerialDeltaPhiTInDeltaPTIndex,weight);
				SerialTrueDeltaPn_InDeltaPTPlot[genie_mode]->Fill(SerialDeltaPnInDeltaPTIndex,weight);	
				SerialTrueECal_InDeltaPTPlot[genie_mode]->Fill(SerialECalInDeltaPTIndex,weight);																							
				SerialTrueProtonCosTheta_InMuonCosThetaPlot[genie_mode]->Fill(SerialProtonCosThetaInMuonCosThetaIndex,weight);
				SerialTrueDeltaPty_InDeltaPtxPlot[genie_mode]->Fill(SerialDeltaPtyInDeltaPtxIndex,weight);
				SerialTrueDeltaPtx_InDeltaPtyPlot[genie_mode]->Fill(SerialDeltaPtxInDeltaPtyIndex,weight);	
				SerialTrueDeltaPT_InMuonCosThetaPlot[genie_mode]->Fill(SerialDeltaPTInMuonCosThetaIndex,weight);
				SerialTrueDeltaPT_InProtonCosThetaPlot[genie_mode]->Fill(SerialDeltaPTInProtonCosThetaIndex,weight);
				SerialTrueDeltaPT_InDeltaAlphaTPlot[genie_mode]->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
				SerialTrueDeltaPn_InDeltaAlphaTPlot[genie_mode]->Fill(SerialDeltaPnInDeltaAlphaTIndex,weight);																							
				SerialTrueECal_InDeltaAlphaTPlot[genie_mode]->Fill(SerialECalInDeltaAlphaTIndex,weight);				

				//----------------------------------------//												

				// 3D analysis treated in 1D grid

				SerialTrueECal_InDeltaPTDeltaAlphaTPlot[0]->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
				SerialTrueECal_InMuonCosThetaMuonMomentumPlot[0]->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
				SerialTrueECal_InProtonCosThetaProtonMomentumPlot[0]->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);

				SerialTrueECal_InDeltaPTDeltaAlphaTPlot[genie_mode]->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
				SerialTrueECal_InMuonCosThetaMuonMomentumPlot[genie_mode]->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
				SerialTrueECal_InProtonCosThetaProtonMomentumPlot[genie_mode]->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);								

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

					TrueCCQEMuonMomentumPlot[0]->Fill(MuonMomentum,weight);
					TrueCCQEMuonPhiPlot[0]->Fill(MuonPhi,weight);
					TrueCCQEMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
					TrueCCQEProtonMomentumPlot[0]->Fill(ProtonMomentum,weight);
					TrueCCQEProtonPhiPlot[0]->Fill(ProtonPhi,weight);
					TrueCCQEProtonCosThetaPlot[0]->Fill(ProtonCosTheta,weight);	
					TrueCCQEECalPlot[0]->Fill(ECal,weight);
					TrueCCQEQ2Plot[0]->Fill(TrueQ2,weight);

					TrueCCQEMuonMomentumPlot[genie_mode]->Fill(MuonMomentum,weight);
					TrueCCQEMuonPhiPlot[genie_mode]->Fill(MuonPhi,weight);
					TrueCCQEMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
					TrueCCQEProtonMomentumPlot[genie_mode]->Fill(ProtonMomentum,weight);
					TrueCCQEProtonPhiPlot[genie_mode]->Fill(ProtonPhi,weight);
					TrueCCQEProtonCosThetaPlot[genie_mode]->Fill(ProtonCosTheta,weight);	
					TrueCCQEECalPlot[genie_mode]->Fill(ECal,weight);
					TrueCCQEQ2Plot[genie_mode]->Fill(TrueQ2,weight);					
													
				}	

			}

		} // End of the selection cuts && the demand that we fill the plots with the same events

		//----------------------------------------//

	} // End of the loop over the events

	//----------------------------------------//

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

	//----------------------------------------//

	// Moving forward towards a cross section extraction
	// Reweight the genie sample with the sum of cv tune weights / number of events
	// Multiply by the flux integrated xsec that has been calculated using flux_averaged_total_xsec

	double ScalingFactor = 1. / (double)(SumCVWeights);
//	double ScalingFactor = 1. / (double)(nentries);

	cout << "# Events = " << nentries << endl;
	cout << "Sum Weights = " << SumCVWeights << endl;

	//----------------------------------------//

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

	//----------------------------------------//

	// Division by bin width to get the cross sections in true variables

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

		//----------------------------------------//
	
		// 1D analysis

		tools.Reweight(TrueMuonMomentumPlot[inte],ScalingFactor);
		tools.Reweight(TrueMuonPhiPlot[inte],ScalingFactor);
		tools.Reweight(TrueMuonCosThetaPlot[inte],ScalingFactor);
		tools.Reweight(TrueMuonCosThetaSingleBinPlot[inte],2*ScalingFactor); // Factor of 2 to account for the fact that the bin width is 2, but we want the number of events, as if the bin width is 1
		tools.Reweight(TrueProtonMomentumPlot[inte],ScalingFactor);
		tools.Reweight(TrueProtonPhiPlot[inte],ScalingFactor);
		tools.Reweight(TrueProtonCosThetaPlot[inte],ScalingFactor);
		tools.Reweight(TrueECalPlot[inte],ScalingFactor);
		tools.Reweight(TrueEQEPlot[inte],ScalingFactor);	
		tools.Reweight(TrueQ2Plot[inte],ScalingFactor);
		tools.Reweight(TrueCCQEMuonMomentumPlot[inte],ScalingFactor);
		tools.Reweight(TrueCCQEMuonPhiPlot[inte],ScalingFactor);
		tools.Reweight(TrueCCQEMuonCosThetaPlot[inte],ScalingFactor);
		tools.Reweight(TrueCCQEProtonMomentumPlot[inte],ScalingFactor);
		tools.Reweight(TrueCCQEProtonPhiPlot[inte],ScalingFactor);
		tools.Reweight(TrueCCQEProtonCosThetaPlot[inte],ScalingFactor);
		tools.Reweight(TrueCCQEECalPlot[inte],ScalingFactor);
		tools.Reweight(TrueCCQEQ2Plot[inte],ScalingFactor);	
		tools.Reweight(TrueDeltaPTPlot[inte],ScalingFactor);
		tools.Reweight(TrueDeltaAlphaTPlot[inte],ScalingFactor);
		tools.Reweight(TrueDeltaPhiTPlot[inte],ScalingFactor);
		tools.Reweight(TrueDeltaPLPlot[inte],ScalingFactor);
		tools.Reweight(TrueDeltaPnPlot[inte],ScalingFactor);
		tools.Reweight(TrueDeltaPtxPlot[inte],ScalingFactor);
		tools.Reweight(TrueDeltaPtyPlot[inte],ScalingFactor);
		tools.Reweight(TrueAPlot[inte],ScalingFactor);
		tools.Reweight(TruekMissPlot[inte],ScalingFactor);
		tools.Reweight(TruePMissPlot[inte],ScalingFactor);
		tools.Reweight(TruePMissMinusPlot[inte],ScalingFactor);

		//----------------------------------------//

		// 2D & 3D analysis (uncorrelated)

		for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

			tools.Reweight(TrueDeltaAlphaT_InDeltaPTTwoDPlot[inte][WhichDeltaPT],ScalingFactor);
			tools.Reweight(TrueDeltaPhiT_InDeltaPTTwoDPlot[inte][WhichDeltaPT],ScalingFactor);	
			tools.Reweight(TrueDeltaPn_InDeltaPTTwoDPlot[inte][WhichDeltaPT],ScalingFactor);	
			tools.Reweight(TrueECal_InDeltaPTTwoDPlot[inte][WhichDeltaPT],ScalingFactor);		

			for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {	

						tools.Reweight(TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[inte][WhichDeltaPT][WhichDeltaAlphaT],ScalingFactor);

			}

		}

		for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {

			tools.Reweight(TrueDeltaPT_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT],ScalingFactor);
			tools.Reweight(TrueDeltaPn_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT],ScalingFactor);			
			tools.Reweight(TrueECal_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT],ScalingFactor);		

		}	

		for (int WhichMuonCosTheta = 0; WhichMuonCosTheta < TwoDNBinsMuonCosTheta; WhichMuonCosTheta++) {

			tools.Reweight(TrueDeltaAlphaT_InMuonCosThetaTwoDPlot[inte][WhichMuonCosTheta],ScalingFactor);
			tools.Reweight(TrueDeltaPT_InMuonCosThetaTwoDPlot[inte][WhichMuonCosTheta],ScalingFactor);		
			tools.Reweight(TrueMuonMomentum_InMuonCosThetaTwoDPlot[inte][WhichMuonCosTheta],ScalingFactor);
			tools.Reweight(TrueProtonCosTheta_InMuonCosThetaTwoDPlot[inte][WhichMuonCosTheta],ScalingFactor);		

			for (int WhichMuonMomentum = 0; WhichMuonMomentum < TwoDNBinsMuonMomentum; WhichMuonMomentum++) {	

						tools.Reweight(TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[inte][WhichMuonCosTheta][WhichMuonMomentum],ScalingFactor);

			}		

		}	

		for (int WhichProtonCosTheta = 0; WhichProtonCosTheta < TwoDNBinsProtonCosTheta; WhichProtonCosTheta++) {

			tools.Reweight(TrueProtonMomentum_InProtonCosThetaTwoDPlot[inte][WhichProtonCosTheta],ScalingFactor);	
			tools.Reweight(TrueDeltaAlphaT_InProtonCosThetaTwoDPlot[inte][WhichProtonCosTheta],ScalingFactor);
			tools.Reweight(TrueDeltaPT_InProtonCosThetaTwoDPlot[inte][WhichProtonCosTheta],ScalingFactor);				

			for (int WhichProtonMomentum = 0; WhichProtonMomentum < TwoDNBinsProtonMomentum; WhichProtonMomentum++) {	

						tools.Reweight(TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[inte][WhichProtonCosTheta][WhichProtonMomentum],ScalingFactor);

			}		

		}

		for (int WhichDeltaPtx = 0; WhichDeltaPtx < TwoDNBinsDeltaPtx; WhichDeltaPtx++) {

			tools.Reweight(TrueDeltaPty_InDeltaPtxTwoDPlot[inte][WhichDeltaPtx],ScalingFactor);	

		}

		for (int WhichDeltaPty = 0; WhichDeltaPty < TwoDNBinsDeltaPty; WhichDeltaPty++) {

			tools.Reweight(TrueDeltaPtx_InDeltaPtyTwoDPlot[inte][WhichDeltaPty],ScalingFactor);	

		}		

		//----------------------------------------//

		// 2D & 3D analysis (correlated)			

		tools.Reweight(SerialTrueDeltaPT_InMuonCosThetaPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueDeltaPT_InProtonCosThetaPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueDeltaPT_InDeltaAlphaTPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueDeltaPn_InDeltaAlphaTPlot[inte],ScalingFactor);		
		tools.Reweight(SerialTrueMuonMomentum_InMuonCosThetaPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueProtonMomentum_InProtonCosThetaPlot[inte],ScalingFactor);	
		tools.Reweight(SerialTrueDeltaAlphaT_InMuonCosThetaPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueDeltaAlphaT_InProtonCosThetaPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueDeltaAlphaT_InDeltaPTPlot[inte],ScalingFactor);		
		tools.Reweight(SerialTrueDeltaPhiT_InDeltaPTPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueDeltaPn_InDeltaPTPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueProtonCosTheta_InMuonCosThetaPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueDeltaPty_InDeltaPtxPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueDeltaPtx_InDeltaPtyPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueECal_InDeltaPTPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueECal_InDeltaAlphaTPlot[inte],ScalingFactor);	
			
		tools.Reweight(SerialTrueECal_InDeltaPTDeltaAlphaTPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueECal_InMuonCosThetaMuonMomentumPlot[inte],ScalingFactor);
		tools.Reweight(SerialTrueECal_InProtonCosThetaProtonMomentumPlot[inte],ScalingFactor);

		//----------------------------------------//


	} // End of the loop over the interaction processes	

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
