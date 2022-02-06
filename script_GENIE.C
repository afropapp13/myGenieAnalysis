{

	vector<TString> WhichSample; vector<TString> WhichName;

	//----------------------------------------//

	WhichSample.push_back("argon40_CCinclMEC_BNBFlux_R-3_0_6_Nominal"); WhichName.push_back("Genie_v3_0_6_Nominal"); // Nominal + ub Tune v2
	WhichSample.push_back("argon40_CCinclMEC_BNBFlux_R-3_0_6_Nominal"); WhichName.push_back("Genie_v3_0_6_NoFSI");
	WhichSample.push_back("argon40_CCinclMEC_BNBFlux_R-3_0_6_hN2018"); WhichName.push_back("Genie_v3_0_6_hN2018");
	WhichSample.push_back("argon40_CCinclMEC_BNBFlux_R-3_0_6_NoRPA"); WhichName.push_back("Genie_v3_0_6_NoRPA");
	WhichSample.push_back("argon40_CCinclMEC_BNBFlux_R-3_0_6_NoCoulomb"); WhichName.push_back("Genie_v3_0_6_NoCoulomb");
	WhichSample.push_back("argon40_CCinclMEC_BNBFlux_R-3_0_6_RFG"); WhichName.push_back("Genie_v3_0_6_RFG");
	WhichSample.push_back("argon40_CCinclMEC_BNBFlux"); WhichName.push_back("Genie_v3_0_6_uB_Tune_1");
	WhichSample.push_back("argon40_CCinclMEC_BNBFlux"); WhichName.push_back("Genie_v3_0_6_Out_Of_The_Box");
	WhichSample.push_back("SuSav2_argon40_CCinclMEC_BNBFlux"); WhichName.push_back("SuSav2");
	WhichSample.push_back("argon40_CCinclMEC_BNBFlux_R-2_12_10"); WhichName.push_back("GENIEv2");
	WhichSample.push_back("argon40_CCinclMEC_BNBFlux_R-3_0_4"); WhichName.push_back("GENIEv3_0_4");				

	//----------------------------------------//

	gROOT->ProcessLine(".L ../myClasses/Util.C+");
	gROOT->ProcessLine(".L ../myClasses/STV_Tools.cxx+");
	gROOT->ProcessLine(".L ../myClasses/Tools.cxx+");	

	gROOT->ProcessLine(".L GenieAnalysis.cxx+");

	for (int i =0;i < (int)(WhichSample.size()); i++) {

		gROOT->ProcessLine("GenieAnalysis(\""+WhichSample[i]+"\",\""+WhichName[i]+"\").Loop()");

	}
	//gROOT->ProcessLine(".q");
};
