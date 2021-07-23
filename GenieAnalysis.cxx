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
#include "../../../../../Secondary_Code/mySimFunctions.cpp"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------------------

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

	// ---------------------------------------------------------------------------------------------------------------------------------

	// Counters STV analysis

	int CounterSTVEventsPassedSelection = 0;
	int CounterSTVQEEventsPassedSelection = 0;
	int CounterSTVMECEventsPassedSelection = 0;
	int CounterSTVRESEventsPassedSelection = 0;
	int CounterSTVDISEventsPassedSelection = 0;
	int CounterSTVCOHEventsPassedSelection = 0;		

	// ---------------------------------------------------------------------------------------------------------------------------------	

	// T2K Weights to make GENIE consistent with the MicroBooNE tune v1
	
	TFile* fweights_uB_Tune_v1 = TFile::Open("mySamples/myWeights_uB_Tune_v1.root");
	TTree *tweights_uB_Tune_v1 = (TTree*)fweights_uB_Tune_v1->Get("ub_tune_cv");
	float cv_weight_uB_Tune_v1;
	tweights_uB_Tune_v1->SetBranchAddress("cv_weight", &cv_weight_uB_Tune_v1);

	// ---------------------------------------------------------------------------------------------------------------------------------	

	int SumCVWeights = 0;
	int SumSelectedCVWeights = 0;	

	// ---------------------------------------------------------------------------------------------------------------------------------	

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry); if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

		// ---------------------------------------------------------------------------------------------------------------------------------

		tweights_uB_Tune_v1->GetEntry(jentry);
		double weight = 1.;
		double cv_weight = 1.;

		if (OutFileName == "Genie_v3_0_6_uB_Tune_1") { 
		
			cv_weight = cv_weight_uB_Tune_v1;	
			weight = cv_weight; 

		}

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
		
		// ---------------------------------------------------------------------------------------------------------------------------------	

		int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0, TrueHeavierMesonCounter = 0;
		vector <int> ProtonID; ProtonID.clear();

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

		// ----------------------------------------------------------------------------------------------------------------------------------

		if (
			MuonMomentum > ArrayNBinsMuonMomentum[0]
		     && ProtonMomentum > ArrayNBinsProtonMomentum[0]
		) {

			if (		
			    // Same events fill all the plots
			        
			    PTmissMomentum > ArrayNBinsDeltaPT[0] && PTmissMomentum < ArrayNBinsDeltaPT[NBinsDeltaPT]
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

				// 2D Analysis
		
				TrueCosThetaMuPmuPlot->Fill(MuonCosTheta,MuonMomentum,weight);
				TrueCosThetaPPpPlot->Fill(ProtonCosTheta,ProtonMomentum,weight);	

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

//	double ScalingFactor = 1. / (double)(SumCVWeights);

	double ScalingFactor = 1. / (double)(nentries);

	// -------------------------------------------------------------------------------------------------------------------------------------------

	// Flux integrated cross sections

	if (OutFileName == "Genie_v3_0_6_uB_Tune_1" || OutFileName == "Genie_v3_0_6_Out_Of_The_Box") { ScalingFactor = ScalingFactor * G18_10a_02_11a_FluxIntegratedXSection ; }
	if (OutFileName == "SuSav2") { ScalingFactor = ScalingFactor * SuSav2FluxIntegratedXSection ; }
	if (OutFileName == "GENIEv2") { ScalingFactor = ScalingFactor * R_2_12_10_FluxIntegratedXSection ; }
	if (OutFileName == "GENIEv3_0_4") { ScalingFactor = ScalingFactor * R_3_0_4_FluxIntegratedXSection ; }

	file->cd();

	// -------------------------------------------------------------------------------------------------------------------------------------------

	// Division by bin width to get the cross sections in true variables
	
	Reweight(TrueMuonMomentumPlot,ScalingFactor);
	Reweight(TrueMuonPhiPlot,ScalingFactor);
	Reweight(TrueMuonCosThetaPlot,ScalingFactor);

	Reweight(TrueProtonMomentumPlot,ScalingFactor);
	Reweight(TrueProtonPhiPlot,ScalingFactor);
	Reweight(TrueProtonCosThetaPlot,ScalingFactor);

	Reweight(TrueECalPlot,ScalingFactor);
	Reweight(TrueEQEPlot,ScalingFactor);	
	Reweight(TrueQ2Plot,ScalingFactor);

	Reweight(TrueCCQEMuonMomentumPlot,ScalingFactor);
	Reweight(TrueCCQEMuonPhiPlot,ScalingFactor);
	Reweight(TrueCCQEMuonCosThetaPlot,ScalingFactor);

	Reweight(TrueCCQEProtonMomentumPlot,ScalingFactor);
	Reweight(TrueCCQEProtonPhiPlot,ScalingFactor);
	Reweight(TrueCCQEProtonCosThetaPlot,ScalingFactor);

	Reweight(TrueCCQEECalPlot,ScalingFactor);
	Reweight(TrueCCQEQ2Plot,ScalingFactor);	
	
	Reweight(TrueDeltaPTPlot,ScalingFactor);
	Reweight(TrueDeltaAlphaTPlot,ScalingFactor);
	Reweight(TrueDeltaPhiTPlot,ScalingFactor);

	// Reweight(TruekMissPlot,ScalingFactor);
	// Reweight(TruePMissPlot,ScalingFactor);
	// Reweight(TruePMissMinusPlot,ScalingFactor);

	Reweight2D(TrueCosThetaMuPmuPlot,ScalingFactor);
	Reweight2D(TrueCosThetaPPpPlot,ScalingFactor);

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
