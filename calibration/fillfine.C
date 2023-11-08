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
fillfine(std::string dirname, std::string outfilename = "finedata.root", unsigned int max_spill = kMaxUInt)
{

  std::map<int, TH2F *> h_fine_device;

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

  /** 
   ** INITIALIZE FRAMER AND PROCESS
   **/

  std::cout << " --- initialize framer: frame size = " << frame_size << std::endl;
  sipm4eic::framer framer(filenames, frame_size);
  framer.set_trigger_coarse_offset(192, 112);
  
  /** loop over spills **/
  int n_spills = 0, n_frames = 0;
  for (int ispill = 0; ispill < max_spill && framer.next_spill(); ++ispill) {
    std::cout << " --- new spill: " << ispill << std::endl;

    /** loop over frames **/
    for (auto &[iframe, aframe] : framer.frames()) {

      /** fill hits **/
      for (auto &[idevice, adevice] : aframe) {
	for (auto &[ichip, achip] : adevice.hits) {
	  for (auto &[ichannel, hits] : achip) {
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
    
    std::cout << "     spill completed " << std::endl;
  } /** end of loop over spills **/

  auto fout = TFile::Open("finedata.root", "RECREATE");
  for (auto &h : h_fine_device)
    h.second->Write();
  fout->Close();  

}
