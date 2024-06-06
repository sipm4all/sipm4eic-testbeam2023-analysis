#include "../lib/lightio.h"
#include "../lib/mapping.h"

enum EReferenceTime_t {
  kTrigger,
  kTiming
};

float
reference_time(sipm4eic::lightio &io, EReferenceTime_t method = kTrigger)
{
  
  /** time from trigger **/
  if (method == kTrigger) {
    auto trigger0_vector = io.get_trigger0_vector();
    auto Ttrg = trigger0_vector[0].coarse;
    return Ttrg;
  }

  /** reference time from timing system **/
  if (method == kTiming) {
    auto timing_vector = io.get_timing_vector();
    
    /** collect timing hits **/
    std::map<int, sipm4eic::lightdata> timing_hits;
    for (auto &hit : timing_vector) {
      auto index = hit.index;
      if (timing_hits.count(index) && timing_hits[index].time() < hit.time())
	continue;
      timing_hits[index] = hit;
    }
    
    /** compute timing time **/
    int Nref = timing_hits.size();
    float Tref = 0.;
    for (auto &[index, hit] : timing_hits)
      Tref += hit.time();
    if (Nref > 0)
      Tref /= Nref;
    
    return Tref;
  }

  return 0.;
}
  
void
recowriter(std::string lightdata_infilename,
	   std::string recodata_outfilename,
	   std::string finecalib_infilename = "",
	   EReferenceTime_t reference_method = kTrigger)
{

  /** read input data **/
  sipm4eic::lightio io;
  io.read_from_tree(lightdata_infilename);

  /** read calib data **/
  sipm4eic::lightdata::load_fine_calibration(finecalib_infilename);

  /** prepare output data **/
  unsigned short n;
  unsigned short ch[65543];
  float x[65534];
  float y[65534];
  float t[65534];
  auto fout = TFile::Open(recodata_outfilename.c_str(), "RECREATE");
  auto tout = new TTree("recodata", "recodata");
  tout->Branch("n", &n, "n/s");
  tout->Branch("ch", &ch, "ch[n]/s");
  tout->Branch("x", &x, "x[n]/F");
  tout->Branch("y", &y, "y[n]/F");
  tout->Branch("t", &t, "t[n]/F");

  int n_spills = 0;
  while (io.next_spill()) {
    std::cout << " --- processing spill: " << n_spills << std::endl;
                 
    while (io.next_frame()) {

      /** reset event **/
      n = 0;

      /** time from scintillators **/
      auto Tref = reference_time(io, reference_method);

      /** loop over cherenkov hits **/
      auto cherenkov_map = io.get_cherenkov_map();
      for (auto &[idx, hits] : cherenkov_map) {
        std::sort(hits.begin(), hits.end());
        auto hit = hits[0];
	auto device = hit.device;
	auto index = hit.index;
	auto chip = index / 32;
	auto time = hit.time();
	auto delta = time - Tref;
        
	if (fabs(delta) > 25.) continue;
     
	auto geo = sipm4eic::get_geo(hit);
	auto pos = sipm4eic::get_position(geo);

	ch[n] = (device - 192) * 256 + index;
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
