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

using namespace std;
using namespace Constants;

void GenieAnalysis::Loop() {

	if (fChain == 0) return;
	TH1D::SetDefaultSumw2();
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	// ---------------------------------------------------------------------------------------------------------------------------------	

	TString FileNameAndPath = "../mySTVAnalysis/myFiles/"+UBCodeVersion+"/CCQEAnalysis_Genie_"+UBCodeVersion+".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");
	std::cout << std::endl << std::endl;

	// ---------------------------------------------------------------------------------------------------------------------------------	

	// 1D Reco Transverse Variables

	TH1D* RecoDeltaPTPlot = new TH1D("RecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* RecoDeltaAlphaTPlot = new TH1D("RecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* RecoDeltaPhiTPlot = new TH1D("RecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* RecoMuonMomentumPlot = new TH1D("RecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* RecoMuonPhiPlot = new TH1D("RecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* RecoMuonCosThetaPlot = new TH1D("RecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH1D* RecoProtonMomentumPlot = new TH1D("RecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* RecoProtonPhiPlot = new TH1D("RecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* RecoProtonCosThetaPlot = new TH1D("RecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH1D* RecoECalPlot = new TH1D("RecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* RecoEQEPlot = new TH1D("RecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
	TH1D* RecoQ2Plot = new TH1D("RecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	// ---------------------------------------------------------------------------------------------------------------------------------	

	// CC1p 1D Reco True Transverse Variables

	TH1D* CC1pRecoDeltaPTPlot = new TH1D("CC1pRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* CC1pRecoDeltaAlphaTPlot = new TH1D("CC1pRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* CC1pRecoDeltaPhiTPlot = new TH1D("CC1pRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* CC1pRecoMuonMomentumPlot = new TH1D("CC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* CC1pRecoMuonPhiPlot = new TH1D("CC1pRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* CC1pRecoMuonCosThetaPlot = new TH1D("CC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH1D* CC1pRecoProtonMomentumPlot = new TH1D("CC1pRecoProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* CC1pRecoProtonPhiPlot = new TH1D("CC1pRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* CC1pRecoProtonCosThetaPlot = new TH1D("CC1pRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH1D* CC1pRecoECalPlot = new TH1D("CC1pRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* CC1pRecoEQEPlot = new TH1D("CC1pRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
	TH1D* CC1pRecoQ2Plot = new TH1D("CC1pRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	TH1D* CC1pTrueDeltaPTPlot = new TH1D("CC1pTrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* CC1pTrueDeltaAlphaTPlot = new TH1D("CC1pTrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* CC1pTrueDeltaPhiTPlot = new TH1D("CC1pTrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* CC1pTrueMuonMomentumPlot = new TH1D("CC1pTrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* CC1pTrueMuonPhiPlot = new TH1D("CC1pTrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* CC1pTrueMuonCosThetaPlot = new TH1D("CC1pTrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH1D* CC1pTrueProtonMomentumPlot = new TH1D("CC1pTrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* CC1pTrueProtonPhiPlot = new TH1D("CC1pTrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* CC1pTrueProtonCosThetaPlot = new TH1D("CC1pTrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH1D* CC1pTrueECalPlot = new TH1D("CC1pTrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* CC1pTrueEQEPlot = new TH1D("CC1pTrueEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
	TH1D* CC1pTrueQ2Plot = new TH1D("CC1pTrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	// 2D Reco True Transverse Variables

	TH2D* CC1pRecoTrueDeltaPTPlot2D = new TH2D("CC1pRecoTrueDeltaPTPlot2D",
		LabelXAxisDeltaPT2D,NBinsDeltaPT,ArrayNBinsDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH2D* CC1pRecoTrueDeltaAlphaTPlot2D = new TH2D("CC1pRecoTrueDeltaAlphaTPlot2D",
		LabelXAxisDeltaAlphaT2D,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH2D* CC1pRecoTrueDeltaPhiTPlot2D = new TH2D("CC1pRecoTrueDeltaPhiTPlot2D",
		LabelXAxisDeltaPhiT2D,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH2D* CC1pRecoTrueMuonMomentumPlot2D = new TH2D("CC1pRecoTrueMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,NBinsMuonMomentum,
		ArrayNBinsMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH2D* CC1pRecoTrueMuonPhiPlot2D = new TH2D("CC1pRecoTrueMuonPhiPlot2D",LabelXAxisMuonPhi2D,NBinsMuonPhi,
		ArrayNBinsMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH2D* CC1pRecoTrueMuonCosThetaPlot2D = new TH2D("CC1pRecoTrueMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,
		ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH2D* CC1pRecoTrueProtonMomentumPlot2D = new TH2D("CC1pRecoTrueProtonMomentumPlot2D",
		LabelXAxisProtonMomentum2D,NBinsProtonMomentum,ArrayNBinsProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH2D* CC1pRecoTrueProtonPhiPlot2D = new TH2D("CC1pRecoTrueProtonPhiPlot2D",LabelXAxisProtonPhi2D,NBinsProtonPhi,
		ArrayNBinsProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH2D* CC1pRecoTrueProtonCosThetaPlot2D = new TH2D("CC1pRecoTrueProtonCosThetaPlot2D",
		LabelXAxisProtonCosTheta2D,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH2D* CC1pRecoTrueECalPlot2D = new TH2D("CC1pRecoTrueECalPlot2D",LabelXAxisECal2D,NBinsECal,ArrayNBinsECal,NBinsECal,ArrayNBinsECal);
	TH2D* CC1pRecoTrueEQEPlot2D = new TH2D("CC1pRecoTrueEQEPlot2D",LabelXAxisEQE2D,NBinsEQE,ArrayNBinsEQE,NBinsEQE,ArrayNBinsEQE);
	TH2D* CC1pRecoTrueQ2Plot2D = new TH2D("CC1pRecoTrueQ2Plot2D",LabelXAxisQ22D,NBinsQ2,ArrayNBinsQ2,NBinsQ2,ArrayNBinsQ2);

	// --------------------------------------------------------------------------------------------------------------------------------------------

	// Bkg 1D Reco True Transverse Variables

	TH1D* BkgRecoTrueDeltaPTPlot = new TH1D("NonCC1pRecoDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* BkgRecoTrueDeltaAlphaTPlot = new TH1D("NonCC1pRecoDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* BkgRecoTrueDeltaPhiTPlot = new TH1D("NonCC1pRecoDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* BkgRecoTrueMuonMomentumPlot = new TH1D("NonCC1pRecoMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* BkgRecoTrueMuonPhiPlot = new TH1D("NonCC1pRecoMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* BkgRecoTrueMuonCosThetaPlot = new TH1D("NonCC1pRecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH1D* BkgRecoTrueProtonMomentumPlot = new TH1D("NonCC1pRecoProtonMomentumPlot",LabelXAxisProtonMomentum,
		NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* BkgRecoTrueProtonPhiPlot = new TH1D("NonCC1pRecoProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* BkgRecoTrueProtonCosThetaPlot = new TH1D("NonCC1pRecoProtonCosThetaPlot",LabelXAxisProtonCosTheta,
		NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH1D* BkgRecoTrueECalPlot = new TH1D("NonCC1pRecoECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* BkgRecoTrueEQEPlot = new TH1D("NonCC1pRecoEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
	TH1D* BkgRecoTrueQ2Plot = new TH1D("NonCC1pRecoQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	// 2D Reco True Transverse Variables

	TH2D* BkgRecoTrueDeltaPTPlot2D = new TH2D("NonCC1pRecoDeltaPTPlot2D",
		LabelXAxisDeltaPT2D,NBinsDeltaPT,ArrayNBinsDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH2D* BkgRecoTrueDeltaAlphaTPlot2D = new TH2D("NonCC1pRecoDeltaAlphaTPlot2D",
		LabelXAxisDeltaAlphaT2D,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH2D* BkgRecoTrueDeltaPhiTPlot2D = new TH2D("NonCC1pRecoDeltaPhiTPlot2D",
		LabelXAxisDeltaPhiT2D,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH2D* BkgRecoMuonMomentumPlot2D = new TH2D("BkgRecoMuonMomentumPlot2D",LabelXAxisMuonMomentum2D,NBinsMuonMomentum,
		ArrayNBinsMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH2D* BkgRecoMuonPhiPlot2D = new TH2D("BkgRecoMuonPhiPlot2D",LabelXAxisMuonPhi2D,NBinsMuonPhi,ArrayNBinsMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH2D* BkgRecoMuonCosThetaPlot2D = new TH2D("BkgRecoMuonCosThetaPlot2D",LabelXAxisMuonCosTheta2D,NBinsMuonCosTheta,
		ArrayNBinsMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH2D* BkgRecoProtonMomentumPlot2D = new TH2D("BkgRecoProtonMomentumPlot2D",
		LabelXAxisProtonMomentum2D,NBinsProtonMomentum,ArrayNBinsProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH2D* BkgRecoProtonPhiPlot2D = new TH2D("BkgRecoProtonPhiPlot2D",LabelXAxisProtonPhi2D,
		NBinsProtonPhi,ArrayNBinsProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH2D* BkgRecoProtonCosThetaPlot2D = new TH2D("BkgRecoProtonCosThetaPlot2D",
		LabelXAxisProtonCosTheta2D,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH2D* BkgRecoTrueECalPlot2D = new TH2D("NonCC1pRecoECalPlot2D",LabelXAxisECal2D,NBinsECal,ArrayNBinsECal,NBinsECal,ArrayNBinsECal);
	TH2D* BkgRecoTrueEQEPlot2D = new TH2D("NonCC1pRecoEQEPlot2D",LabelXAxisEQE2D,NBinsEQE,ArrayNBinsEQE,NBinsEQE,ArrayNBinsEQE);
	TH2D* BkgRecoTrueQ2Plot2D = new TH2D("NonCC1pRecoQ2Plot2D",LabelXAxisQ22D,NBinsQ2,ArrayNBinsQ2,NBinsQ2,ArrayNBinsQ2);

	// -----------------------------------------------------------------------------------------------------------------------------------------------

	// 1D True Transverse Variables

	TH1D* TrueDeltaPTPlot = new TH1D("TrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	TH1D* TrueDeltaAlphaTPlot = new TH1D("TrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	TH1D* TrueDeltaPhiTPlot = new TH1D("TrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);

	TH1D* TrueMuonMomentumPlot = new TH1D("TrueMuonMomentumPlot",LabelXAxisMuonMomentum,NBinsMuonMomentum,ArrayNBinsMuonMomentum);
	TH1D* TrueMuonPhiPlot = new TH1D("TrueMuonPhiPlot",LabelXAxisMuonPhi,NBinsMuonPhi,ArrayNBinsMuonPhi);
	TH1D* TrueMuonCosThetaPlot = new TH1D("TrueMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

	TH1D* TrueProtonMomentumPlot = new TH1D("TrueProtonMomentumPlot",LabelXAxisProtonMomentum,NBinsProtonMomentum,ArrayNBinsProtonMomentum);
	TH1D* TrueProtonPhiPlot = new TH1D("TrueProtonPhiPlot",LabelXAxisProtonPhi,NBinsProtonPhi,ArrayNBinsProtonPhi);
	TH1D* TrueProtonCosThetaPlot = new TH1D("TrueProtonCosThetaPlot",LabelXAxisProtonCosTheta,NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);

	TH1D* TrueECalPlot = new TH1D("TrueECalPlot",LabelXAxisECal,NBinsECal,ArrayNBinsECal);
	TH1D* TrueEQEPlot = new TH1D("TrueEQEPlot",LabelXAxisEQE,NBinsEQE,ArrayNBinsEQE);
	TH1D* TrueQ2Plot = new TH1D("TrueQ2Plot",LabelXAxisQ2,NBinsQ2,ArrayNBinsQ2);

	// ------------------------------------------------------------------------------------------------------------------------------------------------

	TH2D* CVWeightsVsEvPlot = new TH2D("CVWeightsVsEvPlot",";Ev [GeV];cv weights",200,0.,2.,300,0.,3.);
	TH2D* CVWeightsVsInteractionModePlot = new TH2D("CVWeightsVsInteractionModePlot",";CVWeightsVsInteractionModePlot;cv weights",4,-0.5,3.5,300,0.,3.);

	// ------------------------------------------------------------------------------------------------------------------------------------------------

	SetLArTools();
	Tools tools;

	// ---------------------------------------------------------------------------------------------------------------------------------

	// Counters STV analysis

	int CounterSTVEventsPassedSelection = 0;
	int CounterSTVQEEventsPassedSelection = 0;
	int CounterSTVMECEventsPassedSelection = 0;
	int CounterSTVRESEventsPassedSelection = 0;
	int CounterSTVDISEventsPassedSelection = 0;	

	// ---------------------------------------------------------------------------------------------------------------------------------	

	// T2K Weights to make GENIE consistent with the MicroBooNE tune
	
	TFile* fweights = TFile::Open("mySamples/myWeights.root");

	TTree *tweights = (TTree*)fweights->Get("ub_tune_cv");

	float cv_weight;
	tweights->SetBranchAddress("cv_weight", &cv_weight);

	// ---------------------------------------------------------------------------------------------------------------------------------	

	int SumCVWeights = 0;
	int SumSelectedCVWeights = 0;	

	// ---------------------------------------------------------------------------------------------------------------------------------	

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry); if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

		// ---------------------------------------------------------------------------------------------------------------------------------

		tweights->GetEntry(jentry);
		double weight = cv_weight;
//		double weight = 1.;		
		if (weight <= 0 || weight >= 10) { continue; }
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
		
		// ---------------------------------------------------------------------------------------------------------------------------------	

		int ProtonTagging = 0, ChargedPionTagging = 0;
		vector <int> ProtonID; ProtonID.clear();

		for (int i = 0; i < nf; i++) {

			if (pdgf[i] == ProtonPdg && pf[i] > ArrayNBinsProtonMomentum[0]) {

				ProtonTagging ++;
				ProtonID.push_back(i);

			}

			if (fabs(pdgf[i]) == AbsChargedPionPdg && pf[i] > ChargedPionMomentumThres)  {

				ChargedPionTagging ++;

			}

		} // End of the loop over the final state particles

		if ( ProtonTagging != 1 || ChargedPionTagging != 0 ) { continue; }

		// -----------------------------------------------------------------------------------------------------------------------------------

		// Muon

		TVector3 Muon3Vector(pxl,pyl,pzl); // GeV
		TLorentzVector Muon4Vector(pxl,pyl,pzl,El); // GeV
		double MuonMomentum = pl; // GeV / c
		double MuonTheta = TMath::ACos(cthl) * 180. / TMath::Pi(); // deg
		double MuonCosTheta = cthl;
		double MuonPhi = Muon4Vector.Phi() * 180. / TMath::Pi(); // deg

		// ------------------------------------------------------------------------------------------------------------------------------------

		// Proton

		int ProtonIndex = ProtonID.at(0);
		TVector3 Proton3Vector(pxf[ProtonIndex],pyf[ProtonIndex],pzf[ProtonIndex]); // GeV
		TLorentzVector Proton4Vector(pxf[ProtonIndex],pyf[ProtonIndex],pzf[ProtonIndex],Ef[ProtonIndex]); // GeV
		double ProtonMomentum = pf[ProtonIndex]; // GeV / c
		double ProtonTheta = TMath::ACos(cthf[ProtonIndex]) * 180. / TMath::Pi(); // deg
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

		// ----------------------------------------------------------------------------------------------------------------------------------

		if (
			MuonMomentum > ArrayNBinsMuonMomentum[0]
		     && ProtonMomentum > ArrayNBinsProtonMomentum[0]
		) {

			if (		
			    // Same events fill all the STV plots
			        
			    PTmissMomentum > ArrayNBinsDeltaPT[0] && PTmissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
			    && TrueDeltaAlphaT > ArrayNBinsDeltaAlphaT[0] && TrueDeltaAlphaT < ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]
		 	    && TrueDeltaPhiT > ArrayNBinsDeltaPhiT[0] && TrueDeltaPhiT < ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]
		 	    
			    && MuonMomentum < ArrayNBinsMuonMomentum[NBinsMuonMomentum]  
			    && ProtonMomentum < ArrayNBinsProtonMomentum[NBinsProtonMomentum]
			    
			    && MuonCosTheta > ArrayNBinsMuonCosTheta[0]
			    && MuonCosTheta < ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]
			    && ProtonCosTheta > ArrayNBinsProtonCosTheta[0]
			    && ProtonCosTheta < ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]
			    
			    && ECal > ArrayNBinsECal[0] && ECal < ArrayNBinsECal[NBinsECal]
			    && EQE > ArrayNBinsEQE[0] && EQE < ArrayNBinsEQE[NBinsEQE]
			    && TrueQ2 > ArrayNBinsQ2[0] && TrueQ2 < ArrayNBinsQ2[NBinsQ2]
			    
			) {
			
				SumSelectedCVWeights += weight;

				TrueDeltaPTPlot->Fill(PTmissMomentum,weight);
				TrueDeltaAlphaTPlot->Fill(TrueDeltaAlphaT,weight);
				TrueDeltaPhiTPlot->Fill(TrueDeltaPhiT,weight);

				TrueECalPlot->Fill(ECal,weight);
				TrueEQEPlot->Fill(EQE,weight);
				TrueQ2Plot->Fill(TrueQ2,weight);

				CounterSTVEventsPassedSelection++;
				if (qel==1) { CounterSTVQEEventsPassedSelection++; }
				if (mec==1) { CounterSTVMECEventsPassedSelection++; }
				if (res==1) { CounterSTVRESEventsPassedSelection++; }
				if (dis==1) { CounterSTVDISEventsPassedSelection++; }

				TrueMuonMomentumPlot->Fill(MuonMomentum,weight);
				TrueMuonPhiPlot->Fill(MuonPhi,weight);
				TrueMuonCosThetaPlot->Fill(MuonCosTheta,weight);

				TrueProtonMomentumPlot->Fill(ProtonMomentum,weight);
				TrueProtonPhiPlot->Fill(ProtonPhi,weight);
				TrueProtonCosThetaPlot->Fill(ProtonCosTheta,weight);			

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
	double(CounterSTVDISEventsPassedSelection)/ double(CounterSTVEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;	

	// -------------------------------------------------------------------------------------------------------------------------------------------

	// Reweight the genie sample with the cv tune weights to NEvents in the nominal samples

//	double ScalingFactor = (double)(CounterSTVEventsPassedSelection) / (double)(SumSelectedCVWeights);
	double ScalingFactor = 1. / (double)(SumCVWeights);
//	double ScalingFactor = (double)(SumSelectedCVWeights) / (double)(SumCVWeights);
//	double ScalingFactor = 1.;	
//	double ScalingFactor = 1. / (double)(nentries);	
	
	TrueMuonMomentumPlot->Scale(ScalingFactor);
	TrueMuonPhiPlot->Scale(ScalingFactor);
	TrueMuonCosThetaPlot->Scale(ScalingFactor);

	TrueProtonMomentumPlot->Scale(ScalingFactor);
	TrueProtonPhiPlot->Scale(ScalingFactor);
	TrueProtonCosThetaPlot->Scale(ScalingFactor);

	TrueECalPlot->Scale(ScalingFactor);
	TrueEQEPlot->Scale(ScalingFactor);	
	TrueQ2Plot->Scale(ScalingFactor);
	
	TrueDeltaPTPlot->Scale(ScalingFactor);
	TrueDeltaAlphaTPlot->Scale(ScalingFactor);
	TrueDeltaPhiTPlot->Scale(ScalingFactor);
	
	// --------------------------------------------------------------------------------------------------------------------------------------------	
	
	file->cd();
	file->Write();

	std::cout << std::endl;
	std::cout << "File " << FileNameAndPath +" has been created created " << std::endl; 
	std::cout << std::endl;	

	// --------------------------------------------------------------------------------------------------------------------------------------------

} // End of the program
