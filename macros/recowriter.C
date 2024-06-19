#include "../lib/lightio.h"
#include "../lib/mapping.h"

#define TREFCUT
float Tref_min = 32;
float Tref_max = 224;

enum EReferenceTime_t {
  kTrigger,
  kTiming
};

unsigned short
trigger_mask(sipm4eic::lightio *io, float Tref)
{
  unsigned short trgmask = 0;

  auto trigger0_vector = io->get_trigger0_vector();
  for (auto &hit : trigger0_vector) {
    if (std::fabs(hit.coarse - Tref) < 25)
      trgmask |= (1 << 0);
  }
  
  auto trigger1_vector = io->get_trigger1_vector();
  for (auto &hit : trigger1_vector) {
    if (std::fabs(hit.coarse - Tref) < 25)
      trgmask |= (1 << 1);
  }

  auto trigger2_vector = io->get_trigger2_vector();
  for (auto &hit : trigger2_vector) {
    if (std::fabs(hit.coarse - Tref) < 25)
      trgmask |= (1 << 2);
  }

  auto trigger3_vector = io->get_trigger3_vector();
  for (auto &hit : trigger3_vector) {
    if (std::fabs(hit.coarse - Tref) < 25)
      trgmask |= (1 << 3);
  }

  return trgmask;
}

float
reference_time(sipm4eic::lightio *io, EReferenceTime_t method = kTrigger)
{

  /** time from trigger **/
  if (method == kTrigger) {
    auto trigger0_vector = io->get_trigger0_vector();
    auto Ttrg = trigger0_vector.size() > 0 ? trigger0_vector[0].coarse : -666.;
    return Ttrg;
  }

  /** reference time from timing system **/
  if (method == kTiming) {
    auto timing_vector = io->get_timing_vector();
    
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
  auto io = new sipm4eic::lightio;
  io->read_from_tree(lightdata_infilename);

  /** read calib data **/
  //  sipm4eic::lightdata::load_fine_calibration(finecalib_infilename);

  /** read calib data **/
  sipm4eic::lightdata::load_fine_calibration(finecalib_infilename);

  /** prepare output data **/
  unsigned short spill;
  unsigned int frame;
  float tref;
  unsigned short trgmask;
  unsigned short n;
  unsigned short ch[65543];
  float x[65534];
  float y[65534];
  float t[65534];
  auto fout = TFile::Open(recodata_outfilename.c_str(), "RECREATE");
  auto tout = new TTree("recodata", "recodata");
  tout->Branch("spill", &spill, "spill/s");
  tout->Branch("frame", &frame, "frame/i");
  tout->Branch("tref", &tref, "tref/F");
  tout->Branch("trgmask", &trgmask, "trgmask/s");
  tout->Branch("n", &n, "n/s");
  tout->Branch("ch", &ch, "ch[n]/s");
  tout->Branch("x", &x, "x[n]/F");
  tout->Branch("y", &y, "y[n]/F");
  tout->Branch("t", &t, "t[n]/F");

  int n_spills = 0;
  while (io->next_spill()) {
    std::cout << " --- processing spill: " << n_spills << std::endl;
    spill = n_spills;
    
    while (io->next_frame()) {
      frame = io->get_current_frame_id();
      
      /** reset event **/
      n = 0;

      /** reference time information **/
      auto Tref = reference_time(io, reference_method);
      if (Tref == -666.) continue;
#ifdef TREFCUT
      if (Tref < Tref_min || Tref > Tref_max) continue;
#endif
      tref = Tref * sipm4eic::lightdata::coarse_to_ns;

      /** trigger information **/
      trgmask = trigger_mask(io, Tref);
      
      /** loop over cherenkov hits **/
      auto cherenkov_map = io->get_cherenkov_map();
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
