{
//gROOT->ProcessLine(".L Secondary_Code/ToString.cpp");
//gROOT->ProcessLine(".L acceptance_c.cpp");

gROOT->ProcessLine(".L GenieAnalysis.cxx+");
gROOT->ProcessLine("GenieAnalysis().Loop()");
//gROOT->ProcessLine(".q");
};
