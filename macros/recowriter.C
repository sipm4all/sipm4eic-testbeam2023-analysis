#include "../lib/lightio.h"
#include "../lib/mapping.h"

void
recowriter(std::string lightdata_infilename, std::string recodata_outfilename)
{

  /** read input data **/
  sipm4eic::lightio io;
  io.read_from_tree(lightdata_infilename);

  /** prepare output data **/
  unsigned short n;
  float x[65534];
  float y[65534];
  float t[65534];
  auto fout = TFile::Open(recodata_outfilename.c_str(), "RECREATE");
  auto tout = new TTree("recodata", "recodata");
  tout->Branch("n", &n, "n/s");
  tout->Branch("x", &x, "x[n]/F");
  tout->Branch("y", &y, "y[n]/F");
  tout->Branch("t", &t, "t[n]/F");

  int n_spills = 0;
  while (io.next_spill()) {
    std::cout << " --- processing spill: " << n_spills << std::endl;
                 
    while (io.next_frame()) {

      /** reset event **/
      n = 0;

      /** define reference time **/
      auto trigger0_vector = io.get_trigger0_vector();
      auto ref = trigger0_vector[0].coarse;

      /** loop over cherenkov hits **/
      auto cherenkov_map = io.get_cherenkov_map();
      for (auto &[index, hits] : cherenkov_map) {
        std::sort(hits.begin(), hits.end());
        auto hit = hits[0];
	auto coarse = hit.coarse;
	auto delta = coarse - ref;
        
	if (fabs(delta) > 25.) continue;
     
	auto geo = sipm4eic::get_geo(hit);
	auto pos = sipm4eic::get_position(geo);
        
        x[n] = pos[0];
        y[n] = pos[1];
        t[n] = delta * sipm4eic::lightdata::coarse_to_ns;
        ++n;
      }

      tout->Fill();
      
    }
    ++n_spills;
  }

  std::cout << " --- collected " << tout->GetEntries() << " events, " << n_spills << " spills " << std::endl;
  std::cout << " --- output written: " << recodata_outfilename << std::endl;
  
  fout->cd();
  tout->Write();
  fout->Close();
  
}
