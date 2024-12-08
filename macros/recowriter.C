#include "../lib/lightio.h"
#include "../lib/mapping.h"

/*
  COMMENTS BY ANDREA N. ON THE ADDITIONS

  The recowriter macro now accepts a boolean flag called "triggerless", this specifies
  the modality through which the input ROOT file's lightdata was selected: either by triggerless
  or trigger selection. It is important to specify the correct mode because otherwise the
  macro may crash if the flag is wrong (for example many frames selected by the triggerless
  algorithm are not selected by the trigger, so there are no trigger hits for those frames and 
  the recowriter macro in trigger mode cannot work without them)

  The selection of the hits in trigger mode works by cutting off all the hits with a temporal
  distance from the first trigger hit higher than a certain threshold. In triggerless mode the
  macro keeps all the hits who are on the same temporal window as trigger mode, but from half 
  of one of the selected subframes (instead of the first trigger hit)
*/

void recowriter(std::string lightdata_infilename, std::string recodata_outfilename, bool triggerless = false)
{

  /* hardcoded subframe size */
  const int subframe_size = 8;

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

  int n_spills = 0, n_frames = 0, n_hits = 0;
  while (io.next_spill()) {
    std::cout << " --- processing spill: " << n_spills << std::endl;
                 
    while (io.next_frame()) {

      /** reset event **/
      n = 0;

      /** define reference time if trigger mode **/
      int ref = 0;
      if (!triggerless) {
        auto trigger0_vector = io.get_trigger0_vector();
        ref = trigger0_vector[0].coarse;
      }

      /* fetch selected subframes indexes if triggerless mode */
      std::vector<unsigned int> subframes_vector;
      if (triggerless) {
        subframes_vector = io.get_subframe_vector();
      }

      /** loop over cherenkov hits **/
      auto cherenkov_map = io.get_cherenkov_map();
      for (auto &[index, hits] : cherenkov_map) {
        std::sort(hits.begin(), hits.end());
        auto hit = hits[0];
	      auto coarse = hit.coarse;
	      int delta = 0;
        
        /* hit selection */
        if (triggerless) {
          bool select_hit = false;
          for (int isubframe : subframes_vector) {
            int subframe_coarse = isubframe * subframe_size + subframe_size / 2; // coarse time is taken in half of subframe (better than the start)
            delta = subframe_coarse - coarse;
            if (std::fabs(delta) <= 25.) {
              select_hit = true;
              break;
            }
          }
          if (!select_hit) continue;
        
        } else {
          delta = coarse - ref;
	        if (std::fabs(delta) > 25.) continue;
        }

        auto geo = sipm4eic::get_geo(hit);
        auto pos = sipm4eic::get_position(geo);
              
        x[n] = pos[0];
        y[n] = pos[1];
        t[n] = delta * sipm4eic::lightdata::coarse_to_ns;
        ++n;
      }

      tout->Fill();
      n_hits += n;
      ++n_frames;
    }
    ++n_spills;
  }

  std::cout << "Analyzed " << n_frames << " frames in " << n_spills << " spills with " << n_hits
            << " hits which passed the time selection\n";
  std::cout << " --- output written: " << recodata_outfilename << std::endl;

  fout->cd();
  tout->Write();
  fout->Close();
  
}
