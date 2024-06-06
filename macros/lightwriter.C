#include "../lib/framer.h"
#include "../lib/lightio.h"

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

void
lightwriter(std::vector<std::string> filenames, std::string outfilename, std::string fineoutfilename, unsigned int max_spill = kMaxUInt, bool verbose = false)
{

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

  std::cout << " --- initialize framer: frame size = " << frame_size << std::endl;
  sipm4eic::framer framer(filenames, frame_size);
  framer.verbose(verbose);
  framer.set_trigger_coarse_offset(192, 112);
  
  /** loop over spills **/
  int n_spills = 0, n_frames = 0;
  for (int ispill = 0; ispill < max_spill && framer.next_spill(); ++ispill) {

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

	    }}}
	
      } /** end of loop over devices and hits **/
    } /** end of loop over frames **/

    /**
     ** LIGHT DATA
     **/
    
    io.new_spill(ispill);

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

    /** loop over frames **/
    for (auto &frame : framer.frames()) {
      auto iframe = frame.first;
      auto aframe = frame.second;

      io.new_frame(iframe);
      
      /** selection on Luca's trigger, device 192 **/
      if (aframe[192].triggers.size() != 1) continue;
      
      /** selection on timing scintillators, device 207 **/
      auto nsipm4 = aframe[207].hits[4].size();
      auto nsipm5 = aframe[207].hits[5].size();
      if (nsipm4 == 0 && nsipm5 == 0) continue;

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
	  }}}
		
      /** fill cherenkov hits **/
      for (auto &device : aframe) {
	auto idevice = device.first;
	auto adevice = device.second;
	if (idevice == 207) continue;
	for (auto &chip : adevice.hits) {
	  auto ichip = chip.first;
	  auto achip = chip.second;
	  for (auto &channel : achip) {
	    auto ichannel = channel.first;
	    auto hits = channel.second;
	    for (auto &hit : hits) {
	      auto coarse = hit.coarse_time_clock() - iframe * frame_size;
	      io.add_cherenkov(idevice, hit.device_index(), coarse, hit.fine, hit.tdc);
	    }}}
	
      } /** end of loop over devices and hits **/

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
lightwriter(std::string dirname, std::string outfilename, std::string fineoutfilename, unsigned int max_spill = kMaxUInt, bool verbose = false)
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

  lightwriter(filenames, outfilename, fineoutfilename, max_spill, verbose);
}

