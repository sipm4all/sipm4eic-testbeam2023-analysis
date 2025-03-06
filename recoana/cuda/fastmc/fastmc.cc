/** fastmc main.cc **/

#include <boost/program_options.hpp>
#include <iostream>
//#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include <cmath>
#include "mapping.h"

#include "common.h"

extern void fastmc_init(float *x_coords, float *y_coords);
extern void fastmc_process(int nhits, float *x_hits, float *y_hits, int *channels);
extern void fastmc_free();

struct program_options_t {
  int seed;
  std::string recodata;
  int Nevents;
  float Nsig, Nbkg, X0, Y0, R, X0sigma, Y0sigma, Rsigma;
  bool full, ideal;
};

void
process_program_options(int argc, char *argv[], program_options_t &opt)
{
  /** process arguments **/
  namespace po = boost::program_options;
  po::options_description desc("Options");
  try {
    desc.add_options()
      ("help"             , "Print help messages")
      ("seed"             , po::value<int>(&opt.seed)->default_value(123456789), "Random seed")
      ("recodata"         , po::value<std::string>(&opt.recodata)->required(), "Reconstructed data output filename")
      ("Nevents"          , po::value<int>(&opt.Nevents)->required(), "Number of events")
      ("Nsig"             , po::value<float>(&opt.Nsig)->default_value(20.), "Average number of signal hits")
      ("Nbkg"             , po::value<float>(&opt.Nbkg)->default_value(12.), "Average number of background hits")
      ("X0"               , po::value<float>(&opt.X0)->default_value(0.), "X position of ring centre")
      ("X0sigma"          , po::value<float>(&opt.X0sigma)->default_value(0.), "Spread of X position of ring centre")
      ("Y0"               , po::value<float>(&opt.Y0)->default_value(0.), "Y position of ring centre")
      ("Y0sigma"          , po::value<float>(&opt.Y0sigma)->default_value(0.), "Spread of Y position of ring centre")
      ("R"                , po::value<float>(&opt.R)->default_value(70.), "Radius of the ring")
      ("Rsigma"           , po::value<float>(&opt.Rsigma)->default_value(0.), "Spread of radius of the ring")
      ("full"             , po::bool_switch(&opt.full), "Full acceptance")
      ("ideal"            , po::bool_switch(&opt.ideal), "Ideal detector")
      ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
      std::cout << desc << std::endl;
      exit(1);
    }
  }
  catch(std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    std::cout << desc << std::endl;
    exit(1);
  }
}

bool
is_efficient(std::array<int, 3> geo)
{
  if (geo[0] == -1 || geo[1] == -1 || geo[2] == -1) return false;
  
  if (geo[0] == 6) {
    if (geo[1] == 2) {
      if (geo[2] >= 8 && geo[2] < 16) return false;
     else return true;
    }
    else return false;
  }
  
  if (geo[0] == 1) {
    if (geo[1] == 1 && (geo[2] >= 0 && geo[2] < 32)) return false;
    if (geo[1] == 2 && (geo[2] >= 32 && geo[2] < 64)) return false;
    if (geo[1] == 3 && (geo[2] >= 40 && geo[2] < 48)) return false;
    return true;
  }
  
  if (geo[0] == 7) {
    if (geo[1] == 4) return true;
    return false;
  }
  
  if (geo[0] == 8) {
    if (geo[1] == 4) return true;
    return false;
  }
  
  if (geo[0] == 5) {
    if (geo[1] == 2) return true;
    return false;
  }
  
  if (geo[0] == 3) {
    if (geo[1] == 2 && (geo[2] >= 8 && geo[2] < 16)) return false;
    return true;
  }
  
  if (geo[0] == 2) {
    if (geo[1] == 1 && (geo[2] >= 40 && geo[2] < 48)) return false;
    if (geo[1] == 4 && (geo[2] >= 48 && geo[2] < 56)) return false;
    return true;
  }
  
  if (geo[0] == 4) {
    if (geo[1] == 2 && (geo[2] >= 32 && geo[2] < 40)) return false;
    if (geo[1] == 2 && (geo[2] >= 40 && geo[2] < 48)) return false;
    if (geo[1] == 2 && (geo[2] >= 48 && geo[2] < 56)) return false;
    if (geo[1] == 4 && (geo[2] >= 16 && geo[2] < 24)) return false;
    return true;
  }
  
  return true;
}

void
init_map(float *x_coords, float *y_coords, bool *eff)
{
  std::cout << " --- initialise map " << std::endl;
  for (int ich = 0; ich < N_CHANNELS; ++ich) {
    auto pdu = 1 + ich / 256;
    auto matrix = 1 + (ich % 256) / 64;
    auto doch = ich % 64;
    auto geo = sipm4eic::get_geo(pdu, matrix, doch);
    auto pos = sipm4eic::get_position(geo);
    x_coords[ich] = pos[0];
    y_coords[ich] = pos[1];
    eff[ich] = is_efficient({pdu, matrix, doch});
  }
}

int
main(int argc, char *argv[])
{

  program_options_t opt;
  process_program_options(argc, argv, opt);

  /** create output reconstructed data tree **/
  unsigned short N;
  float X0[65535];
  float Y0[65535];
  float R[65535];
  unsigned short n;
  float x[65535];
  float y[65535];
  float t[65535];
  unsigned short id[65535];
  auto fout = TFile::Open(opt.recodata.c_str(), "RECREATE");
  auto tout = new TTree("recodata", "recodata");
  tout->Branch("N", &N, "N/s");
  tout->Branch("X0", &X0, "X0[N]/F");
  tout->Branch("Y0", &Y0, "Y0[N]/F");
  tout->Branch("R", &R, "R[N]/F");
  tout->Branch("n", &n, "n/s");
  tout->Branch("x", &x, "x[n]/F");
  tout->Branch("y", &y, "y[n]/F");
  tout->Branch("t", &t, "t[n]/F");
  tout->Branch("id", &id, "id[n]/s");

  float x_coords[N_CHANNELS] = {0.}, y_coords[N_CHANNELS] = {0.};
  bool eff[N_CHANNELS] = {false}, hit[N_CHANNELS] = {false};
  int n_hits[MAX_EVENTS];
  float x_hits[MAX_HITS], y_hits[MAX_HITS];
  int channels[MAX_HITS];
  unsigned short id_hits[MAX_HITS];
  
  /** init **/  
  init_map(x_coords, y_coords, eff);
  fastmc_init(x_coords, y_coords);
  gRandom->SetSeed(opt.seed);

  /** this is fixed until we simulate more than one ring **/
  N = 0;
  X0[N] = opt.X0;
  Y0[N] = opt.Y0;
  R[N] = opt.R;
  ++N;
  
  /** loop over events **/
  int bufevents = 0;
  int bufhits = 0;
  for (int iev = 0; iev < opt.Nevents; ++iev) {

    auto Nsig = gRandom->Poisson(opt.Nsig);
    auto Nbkg = gRandom->Poisson(opt.Nbkg);

    /** process buffered events if we cannot store more **/
    if (bufevents >= MAX_EVENTS ||
	bufhits + Nsig + Nbkg >= MAX_HITS) {
      std::cout << " --- process " << bufevents << " buffered events: " << bufhits << " hits " << std::endl;
      int offhit = 0;
      fastmc_process(bufhits, x_hits, y_hits, channels);
      for (int iiev = 0; iiev < bufevents; ++iiev) {
	n = 0;
	for (int ihit = 0; ihit < n_hits[iiev]; ++ihit) {
	  if (opt.ideal) {
	    auto ch = channels[offhit + ihit];
	    x[n] = x_hits[offhit + ihit];
	    y[n] = y_hits[offhit + ihit];
	    t[n] = 0.;
	    if (ch == -1) t[n] = -1.;
	  } else {
	    auto ch = channels[offhit + ihit];
	    if (ch == -1) continue;
	    if (hit[ch]) continue;
	    if (!opt.full && !eff[ch]) continue;
	    x[n] = x_coords[ch];
	    y[n] = y_coords[ch];
	    t[n] = 0.;
	    hit[ch] = true;
	  }
	  id[n] = id_hits[offhit + ihit];
	  ++n;
	}
	/** reset hit channels **/
	for (int ihit = 0; ihit < n_hits[iiev]; ++ihit) {
	  auto ch = channels[offhit + ihit];
	  if (ch == -1) continue;
	  hit[ch] = false;
	}	
	tout->Fill();
	offhit += n_hits[iiev];
      }
      
      bufhits = 0;
      bufevents = 0;
    }

    /** generate event hits **/
    
    n_hits[bufevents] = 0;

    /** generate signal hits **/
    auto X0 = gRandom->Gaus(opt.X0, opt.X0sigma);
    auto Y0 = gRandom->Gaus(opt.Y0, opt.Y0sigma);
    for (int isig = 0; isig < Nsig; ++isig) {
      auto R = gRandom->Gaus(opt.R, opt.Rsigma);
      auto phi = gRandom->Uniform(0., 2. * M_PI);
      x_hits[bufhits] = X0 + R * std::cos(phi);
      y_hits[bufhits] = Y0 + R * std::sin(phi);
      id_hits[bufhits] = 0;
      ++bufhits;
      ++n_hits[bufevents];
    }
    
    /** generate background hits **/
    for (int ibkg = 0; ibkg < Nbkg; ++ibkg) {
      x_hits[bufhits] = gRandom->Uniform(-100., 100.);
      y_hits[bufhits] = gRandom->Uniform(-100., 100.);
      id_hits[bufhits] = 65535;
      ++bufhits;
      ++n_hits[bufevents];
    }

    ++bufevents;
  } /** end of loop over events **/
  
  /** process remaining buffered events **/
  std::cout << " --- process " << bufevents << " buffered events: " << bufhits << " hits " << std::endl;
  int offhit = 0;
  fastmc_process(bufhits, x_hits, y_hits, channels);
  for (int iiev = 0; iiev < bufevents; ++iiev) {
    n = 0;
    for (int ihit = 0; ihit < n_hits[iiev]; ++ihit) {
      if (opt.ideal) {
	x[n] = x_hits[offhit + ihit];
	y[n] = y_hits[offhit + ihit];
	t[n] = 0.;
      } else {
	auto ch = channels[offhit + ihit];
	if (ch == -1) continue;
	if (hit[ch]) continue;
	if (!opt.full && !eff[ch]) continue;
	x[n] = x_coords[ch];
	y[n] = y_coords[ch];
	t[n] = 0.;
	hit[ch] = true;
      }
      id[n] = id_hits[offhit + ihit];
      ++n;
    }
    /** reset hit channels **/
    for (int ihit = 0; ihit < n_hits[iiev]; ++ihit) {
      auto ch = channels[offhit + ihit];
      if (ch == -1) continue;
      hit[ch] = false;
    }	
    tout->Fill();
    offhit += n_hits[iiev];
  }
  
  /** write output and close **/
  fout->cd();
  tout->Write();
  fout->Close();

  return 0;
}
