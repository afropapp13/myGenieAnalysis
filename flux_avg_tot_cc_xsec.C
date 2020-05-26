#include <fstream>
#include <sstream>

constexpr double FLUX_INTEGRAL_RELATIVE_TOLERANCE = 1e-5;

// Computes the flux-averaged total cross section given a histogram
// of the neutrino flux and a total cross section TGraph (both with
// arbitrary units)
double flux_averaged_total_xsec(TH1& flux_hist,
  const TGraph& xsec_spline)
{
  double flux_integral = flux_hist.Integral( "width" );

  // Integrate up to the smaller of either the final point on the
  // total cross section spline or the left edge of the overflow bin
  // of the flux histogram
  // TODO: Add more careful checking here. This assumes that you're doing
  // something reasonable.
  double spline_max_energy = xsec_spline.GetX()[ xsec_spline.GetN() - 1 ];
  double flux_max_energy = flux_hist.GetBinLowEdge( flux_hist.GetNbinsX() + 1 );
  double max_energy = std::min( spline_max_energy, flux_max_energy );

  // Function to use for numerical integration. Element x[0] is the
  // neutrino energy
  TF1 xsec_weighted_flux_func("temp_func", [&](double* x, double*)
    {
      int flux_bin = flux_hist.FindBin( x[0] );
      double flux = flux_hist.GetBinContent( flux_bin );

      double total_xsec = xsec_spline.Eval( x[0] );

      return flux * total_xsec / flux_integral;
    }, 0., max_energy, 0);

  return xsec_weighted_flux_func.Integral(0., max_energy,
    FLUX_INTEGRAL_RELATIVE_TOLERANCE);
}

void flux_avg_tot_cc_xsec(const std::string& flux_filename,
  const std::string& flux_hist_name, const std::string& spline_filename)
{
  TFile flux_file( flux_filename.c_str(), "read" );
  TH1D* flux_hist = nullptr;
  flux_file.GetObject( flux_hist_name.c_str(), flux_hist );

  TFile spline_file( spline_filename.c_str(), "read" );
  TGraph* spline_graph = nullptr;
  spline_file.GetObject("nu_mu_Ar40/tot_cc", spline_graph); // CC inclusive

  double totalXSecAvg = flux_averaged_total_xsec(*flux_hist, *spline_graph);

  // The gspl2root files use units of 1e-38 cm^2 for the cross section
  // splines
  std::cout << "Flux-averaged total CC cross section: " << totalXSecAvg
    << " * 1e-38 cm^2\n";
}
