{

	gROOT->ProcessLine(".L ../myClasses/Util.C+");
	gROOT->ProcessLine(".L ../myClasses/STV_Tools.cxx+");
	gROOT->ProcessLine(".L ../myClasses/Tools.cxx+");	

	gROOT->ProcessLine(".L GenieAnalysis.cxx+");
	gROOT->ProcessLine("GenieAnalysis().Loop()");
	//gROOT->ProcessLine(".q");
};
