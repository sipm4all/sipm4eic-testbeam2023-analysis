#include "../lib/framer.h"
#include "../lib/lightio.h"

/*
  COMMENTS BY ANDREA N. ON THE ADDITIONS

  The lightwriter macro now accepts a boolean flag called 
  "triggerless" which lets a user choose if it prefers to perform a 
  trigger or triggerless selection.

  The triggerless selection is performed by the "triggerless_selection" function
  which divides a given frame in subframes of a give size and selects those who 
  have a number of hits higher than a chosen threshold, then it returns the indexes
  of those subframe relative to the given frame. The macro then writes to a tree only
  the frames with at least one selected subframe and takes care of recording the indexes
  of all selected subframes for each frame, for subsequent analysis

*/

const int frame_size = 256;

std::vector<std::string> devices = {
  "kc705-192",
  "kc705-193",
  "kc705-194",
  "kc705-195",
  "kc705-196",
  "kc705-197",
  "kc705-198",
  "kc705-207"
};

/*
  checks for subframes inside the given frame that have a number of hits higher than the given threshold
  returns a vector with the indexes of subframes (relative to the frame) which satisfy the selection, if 
  no acceptable subframe is found, an empty vector is returned.
  Note: subframe_size should be a divisor of frame_size otherwise some rounding errors may happen
*/
std::vector<unsigned int> triggerless_selection(std::pair<const int, sipm4eic::framer::frame_t> &frame, int frame_size,
                                                  int subframe_size, int min_subframe_hits) {

  int n_subframes = frame_size / subframe_size;
  std::vector<int> subframe_hits(n_subframes, 0);

  for (auto &device : frame.second) { // the devices are the kc705 boards
    auto idevice = device.first;
    auto adevice = device.second;
    if (idevice == 207) continue;
    
    for (auto &chip : adevice.hits) { // the chips are the ALCOR chips
      auto achip = chip.second;
    
      for (auto &channel : achip) { // the channels are the SiPM sensors
        auto hits = channel.second;
    
        for (auto &hit : hits) {
          auto coarse = hit.coarse_time_clock() - frame.first * frame_size;
          int isubframe = coarse / subframe_size;
          ++subframe_hits[isubframe];
        }
      }
    } 
  }

  std::vector<unsigned int> subframe_index;
  for (unsigned int i = 0; i < n_subframes; ++i) {
    if (subframe_hits[i] >= min_subframe_hits) {
      subframe_index.push_back(i);
    }
  }
  
  return subframe_index;
}

void lightwriter(std::vector<std::string> filenames, std::string outfilename, std::string fineoutfilename, bool triggerless = false, unsigned int max_spill = kMaxUInt, bool verbose = false)
{
  // hardcoded values, we keep a standard frame size of 256 clock hits
  const int frame_size = 256;
  const int subframe_size = 8;
  const int min_subframe_hits = 5;

  /**
   ** CREATE OUTPUT TREE
   **/

  sipm4eic::lightio io;
  io.write_to_tree(outfilename);

  /** 
   ** FINE OUTPUT 
   **/ 

  std::map<int, TH2F *> h_fine_device;

  
  /** 
   ** INITIALIZE FRAMER AND PROCESS
   **/

  // Note: the framer class is agnostic to the subframe functionality
  std::cout << " --- initialize framer: frame size = " << frame_size << std::endl;
  sipm4eic::framer framer(filenames, frame_size);
  framer.verbose(verbose);
  framer.set_trigger_coarse_offset(192, 112);
  
  /** loop over spills **/
  int n_spills = 0, n_frames = 0, n_subframes = 0, n_hits = 0;
  for (int ispill = 0; ispill < max_spill && framer.next_spill(); ++ispill) {

/** keep this, but it is not needed **/

    /**
     ** FINE FILL
     **/

    /** loop over frames **/
    for (auto &frame : framer.frames()) {
      auto iframe = frame.first;
      auto aframe = frame.second;

      /** fill hits **/
      for (auto &device : aframe) {
        auto idevice = device.first;
        auto adevice = device.second;
        for (auto &chip : adevice.hits) {
          auto ichip = chip.first;
          auto achip = chip.second;
          for (auto &channel : achip) {
            auto ichannel = channel.first;
            auto hits = channel.second;
            for (auto &hit : hits) {
                  auto device = idevice;
                  if (!h_fine_device.count(device))
                    h_fine_device[device] = new TH2F(Form("hFine_%d", device), "hFine", 768, 0, 768, 256, 0, 256);

                  auto fine = hit.fine;
                  auto index = hit.device_index();
                  auto tdc = hit.tdc;
                  auto cindex = tdc + 4 * index;
                  h_fine_device[device]->Fill(cindex, fine);
            }
          }
        }   
      } /** end of loop over devices and hits **/
    } /** end of loop over frames **/

/** end of keep this, but is it not needed **/

    /**
     ** LIGHT DATA
     **/
    
    io.new_spill(ispill);

/** keep this, but it it not needed **/

    for (auto &part : framer.part_mask()) {
      auto idevice = part.first;
      auto amask = part.second;
      io.add_part(idevice, amask);
    }
    for (auto &dead : framer.dead_mask()) {
      auto idevice = dead.first;
      auto amask = dead.second;
      io.add_dead(idevice, amask);
    }

/** end of keep this, but is it not needed **/

    /** loop over frames **/
    for (auto &frame : framer.frames()) {
      auto iframe = frame.first;
      auto aframe = frame.second;

      io.new_frame(iframe);

      std::vector<unsigned int> subframes_index;
      if (triggerless) {
        /* triggerless selection */
        subframes_index = triggerless_selection(frame, frame_size, subframe_size, min_subframe_hits);
        if (subframes_index.size() == 0) continue;
      } else {
        /* trigger selection */

        /** selection on Luca's trigger, device 192 **/
        if (aframe[192].triggers.size() != 1) continue;
        
        /** selection on timing scintillators, device 207 **/
        auto nsipm4 = aframe[207].hits[4].size();
        auto nsipm5 = aframe[207].hits[5].size();
        if (nsipm4 == 0 && nsipm5 == 0) continue;  
      }

/** keep this, but it is not needed **/

      /** fill trigger0 hits **/
      auto trigger0 = aframe[192].triggers;
      for (auto &trigger : trigger0)
      	io.add_trigger0(trigger.coarse_time_clock() - iframe * frame_size);

      /** fill timing hits **/
      for (auto &chip : aframe[207].hits) {
        auto ichip = chip.first;
        auto achip = chip.second;
        for (auto &channel : achip) {
          auto ichannel = channel.first;
          auto hits = channel.second;
          for (auto &hit : hits) {
            auto coarse = hit.coarse_time_clock() - iframe * frame_size;
            io.add_timing(207, hit.device_index(), coarse, hit.fine, hit.tdc);
	        }
        }    
      }

/** end of keep this, but it is not needed **/

      /** fill cherenkov hits **/
      for (auto &device : aframe) { // the devices are the kc705 boards
        auto idevice = device.first;
        auto adevice = device.second;
        if (idevice == 207) continue; // skip scintillators

        for (auto &chip : adevice.hits) { // the chips are the ALCOR chips
          auto ichip = chip.first;
          auto achip = chip.second;

          for (auto &channel : achip) { // the channels are the SiPM sensors
            auto ichannel = channel.first;
            auto hits = channel.second;

            for (auto &hit : hits) {
              auto coarse = hit.coarse_time_clock() - iframe * frame_size;
              io.add_cherenkov(idevice, hit.device_index(), coarse, hit.fine, hit.tdc);
            }
          }
        }
      } /** end of loop over devices and hits **/

      /** fill selected subframes **/
      if (triggerless) {
        for (auto index : subframes_index){
          io.add_subframe(index);
          ++n_subframes;
        }
      }

      io.add_frame();
      
    } /** end of loop over frames **/

    io.fill();
    ++n_spills;

  } /** end of loop over spills **/

  /** 
   ** WRITE OUTPUT TO FILE
   **/

  std::cout << " --- writing light data output file: " << outfilename << std::endl;
  io.write_and_close();

  if (!fineoutfilename.empty()) {
    std::cout << " --- writing fine data output file: " << fineoutfilename << std::endl;
    auto fout = TFile::Open(fineoutfilename.c_str(), "RECREATE");
    for (auto &h : h_fine_device)
      h.second->Write();
    fout->Close();
  }

  std::cout << " --- completed: " << n_spills << " spills " << std::endl;

}

void
lightwriter(std::string dirname, std::string outfilename, std::string fineoutfilename, bool triggerless = false, unsigned int max_spill = kMaxUInt, bool verbose = false)
{

  /** 
   ** BUILD INPUT FILE LIST
   **/

  std::vector<std::string> filenames;
  for (auto device : devices) {
    for (int ififo = 0; ififo < 25; ++ififo) {
      std::string filename = dirname + "/" + device + "/decoded/alcdaq.fifo_" + std::to_string(ififo) + ".root";
      filenames.push_back(filename);
    }
  }

  lightwriter(filenames, outfilename, fineoutfilename, triggerless, max_spill, verbose);
}

