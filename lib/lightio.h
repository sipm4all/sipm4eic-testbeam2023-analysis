#pragma once

#include "lightdata.h"

/*
  COMMENTS BY ANDREA N. ON THE ADDITIONS

  Each frame of standard frame size 256 (clock intervals) is subdivided in 
  subframes of size 8, if a subframe is selected (i.e. the number of hits
  contained is over a chosen threshold) then its index relative to the frame
  (i.e. how many subframes have been already counted after the start of the 
  frame) is written in the lightio::subframe_index vector, while the subframe_n
  vector records the numer of saved subframe indexes for each frame, then to 
  recover the position of a subframe inside a spill (in terms of clock hits) you
  just need to apply the formula
  
      subframe_pos = frame_index * frame_size + subframe_index * subframe_size
  
  note that, as already stated, we hardcoded frame_size to 256 and subframe_size to 8
  in the code, they are not yet completely generalizable
   
*/

namespace sipm4eic {
  
class lightio {

 public:

  static const int frame_size = 256;
  static const int max_devices = 256;
  static const int max_frames = 65534;   // maximum number of frames in a spill
  static const int max_triggers = 65534; // maximum number of triggers in a spill
  static const int max_hits = 262144;   // maximum number of hits in a spill
  static const int max_subframes = 262144; // maximum number of subframes in a spill

/* 
  Note: the correct number of max_subframes should be 32 (num. of size 8 subframes in a size 256 frame) * max_frames = 2097088
        but for some reasons it is too big for the ROOT compiler and it does not run on my machine, so we use max_hits since 
        a subframe must have at least 1 hit which means there are at most max_hits total subframes
*/

  unsigned char part_n;
  unsigned char part_device[max_devices];
  unsigned int part_mask[max_devices];
  //
  unsigned char dead_n;
  unsigned char dead_device[max_devices];
  unsigned int dead_mask[max_devices];
  //  
  unsigned short frame_n;         // number of frames
  unsigned int frame[max_frames]; // frame index relative to a spill
  //
  unsigned int subframe_size;
  unsigned int subframe_n[max_frames]; // number of selected subframes in each frame
  unsigned int subframe_index[max_subframes]; // index of subframe relative to a frame
  //
  unsigned int trigger0_size;           // trigger hits in spill
  unsigned char trigger0_n[max_frames]; // trigger hits in frame
  unsigned char trigger0_coarse[max_triggers];
  //
  unsigned int timing_size;            // timing hits in spill
  unsigned short timing_n[max_frames]; // timing hits in frame
  unsigned char timing_device[max_hits];
  unsigned char timing_index[max_hits];
  unsigned char timing_coarse[max_hits];
  unsigned char timing_fine[max_hits];
  unsigned char timing_tdc[max_hits];
  //
  unsigned int cherenkov_size;            // cherenkov hits in spill
  unsigned short cherenkov_n[max_frames]; // cherenkov hits in frame
  unsigned char cherenkov_device[max_hits];
  unsigned char cherenkov_index[max_hits];
  unsigned char cherenkov_coarse[max_hits];
  unsigned char cherenkov_fine[max_hits];
  unsigned char cherenkov_tdc[max_hits];
  
  lightio() = default;
  
  void new_spill(unsigned int ispill);
  void new_frame(unsigned int iframe);
  void add_part(unsigned char device, unsigned int mask);
  void add_dead(unsigned char device, unsigned int mask);
  void add_trigger0(unsigned char coarse);
  void add_timing(unsigned char device, unsigned char index, unsigned char coarse, unsigned char fine, unsigned char tdc);
  void add_cherenkov(unsigned char device, unsigned char index, unsigned char coarse, unsigned char fine, unsigned char tdc);
  void add_frame() { ++frame_n; };
  void add_subframe(unsigned int index);
  void fill();
  void write_and_close();  
  void write_to_tree(std::string filename, std::string treename = "lightdata");
  void write_to_tree(TTree *t);

  void read_from_tree(std::string filename, std::string treename = "lightdata");
  void read_from_tree(TTree *t);
  bool next_spill();
  bool next_frame();
  void reset() { spill_current = frame_current = 0; };
  
  std::vector<lightdata> &get_trigger0_vector() { return trigger0_vector; };
  std::vector<lightdata> &get_timing_vector() { return timing_vector; };
  std::vector<lightdata> &get_cherenkov_vector() { return cherenkov_vector; };
  std::vector<unsigned int> &get_subframe_vector() { return subframe_vector; };

  std::map<std::array<unsigned char, 2>, std::vector<lightdata>> &get_timing_map() { return timing_map; };
  std::map<std::array<unsigned char, 2>, std::vector<lightdata>> &get_cherenkov_map() { return cherenkov_map; };

  TTree *get_tree() { return tree; };
  
 private:

  TFile *file = nullptr;
  TTree *tree = nullptr;

  int spill_current = 0;
  int frame_current = 0;
  
  int trigger0_offset = 0;
  int timing_offset = 0;
  int cherenkov_offset = 0;
  int subframe_offset = 0;

  std::vector<lightdata> trigger0_vector;
  std::vector<lightdata> timing_vector;
  std::vector<lightdata> cherenkov_vector;
  std::vector<unsigned int> subframe_vector;

  std::map<std::array<unsigned char, 2>, std::vector<lightdata>> timing_map;
  std::map<std::array<unsigned char, 2>, std::vector<lightdata>> cherenkov_map;
  
};

void
lightio::new_spill(unsigned int ispill)
{
  std::cout << " --- new spill: " << ispill << std::endl;
  part_n = 0;
  dead_n = 0;
  frame_n = 0;
  trigger0_size = 0;
  timing_size = 0;
  cherenkov_size = 0;
  subframe_size = 0;
};

void
lightio::new_frame(unsigned int iframe) {
  frame[frame_n] = iframe;
  trigger0_n[frame_n] = 0;
  timing_n[frame_n] = 0;
  cherenkov_n[frame_n] = 0;
  subframe_n[frame_n] = 0;
}

void
lightio::add_part(unsigned char device, unsigned int mask) {
  part_device[part_n] = device;
  part_mask[part_n] = mask;
  ++part_n;
}
 
void
lightio::add_dead(unsigned char device, unsigned int mask) {
  dead_device[dead_n] = device;
  dead_mask[dead_n] = mask;
  ++dead_n;
}

void
lightio::add_trigger0(unsigned char coarse) {
  trigger0_coarse[trigger0_size] = coarse;
  ++trigger0_n[frame_n];
  ++trigger0_size;
}

void
lightio::add_timing(unsigned char device, unsigned char index, unsigned char coarse, unsigned char fine, unsigned char tdc) {
  timing_device[timing_size] = device;
  timing_index[timing_size] = index;
  timing_coarse[timing_size] = coarse;
  timing_fine[timing_size] = fine;
  timing_tdc[timing_size] = tdc;
  ++timing_n[frame_n];
  ++timing_size;
}

void
lightio::add_cherenkov(unsigned char device, unsigned char index, unsigned char coarse, unsigned char fine, unsigned char tdc) {
  cherenkov_device[cherenkov_size] = device;
  cherenkov_index[cherenkov_size] = index;
  cherenkov_coarse[cherenkov_size] = coarse;
  cherenkov_fine[cherenkov_size] = fine;
  cherenkov_tdc[cherenkov_size] = tdc;
  ++cherenkov_n[frame_n];
  ++cherenkov_size;
}

void lightio::add_subframe(unsigned int index) {
  subframe_index[subframe_size] = index;
  ++subframe_n[frame_n];
  ++subframe_size;
}

void
lightio::fill() {
  std::cout << " --- fill tree: trigger0_size = " << trigger0_size << std::endl;
  std::cout << "                  timing_size = " << timing_size << std::endl;
  std::cout << "                     cherenkov_size = " << cherenkov_size << std::endl;
  std::cout << "                        subframe_size = " << subframe_size << std::endl;
  tree->Fill();
};
 
void
lightio::write_and_close() {
#if 0
  auto n_spills = tree->GetEntries();
  auto n_frames = 0;
  for (int ispill = 0; ispill < n_spills; ++ispill) {
    tree->GetEvent(ispill);
    n_frames += frame_n;
  }
  std::cout << " --- write and close " << std::endl;
  std::cout << " --- collected " << n_spills << " spills " << std::endl;
  std::cout << "               " << n_frames << " frames " << std::endl;
#endif
  file->cd();
  tree->Write();
  file->Close();
}

void
lightio::write_to_tree(std::string filename, std::string treename)
{
  file = TFile::Open(filename.c_str(), "RECREATE");
  tree = new TTree(treename.c_str(), treename.c_str());
  write_to_tree(tree);
}

void
lightio::write_to_tree(TTree *t)
{
  t->Branch("part_n", &part_n, "part_n/b");
  t->Branch("part_device", &part_device, "part_device[part_n]/b");
  t->Branch("part_mask", &part_mask, "part_mask[part_n]/i");
  t->Branch("dead_n", &dead_n, "dead_n/b");
  t->Branch("dead_device", &dead_device, "dead_device[dead_n]/b");
  t->Branch("dead_mask", &dead_mask, "dead_mask[dead_n]/i");
  t->Branch("frame_n", &frame_n, "frame_n/s");
  t->Branch("frame", &frame, "frame[frame_n]/i");
  t->Branch("trigger0_size", &trigger0_size, "trigger0_size/i");
  t->Branch("trigger0_n", &trigger0_n, "trigger0_n[frame_n]/b");
  t->Branch("trigger0_coarse", &trigger0_coarse, "trigger0_coarse[trigger0_size]/b");
  t->Branch("timing_size", &timing_size, "timing_size/i");
  t->Branch("timing_n", &timing_n, "timing_n[frame_n]/s");
  t->Branch("timing_device", &timing_device, "timing_device[timing_size]/b");
  t->Branch("timing_index", &timing_index, "timing_index[timing_size]/b");
  t->Branch("timing_coarse", &timing_coarse, "timing_coarse[timing_size]/b");
  t->Branch("timing_fine", &timing_fine, "timing_fine[timing_size]/b");
  t->Branch("timing_tdc", &timing_tdc, "timing_tdc[timing_size]/b");
  t->Branch("cherenkov_size", &cherenkov_size, "cherenkov_size/i");
  t->Branch("cherenkov_n", &cherenkov_n, "cherenkov_n[frame_n]/s");
  t->Branch("cherenkov_device", &cherenkov_device, "cherenkov_device[cherenkov_size]/b");
  t->Branch("cherenkov_index", &cherenkov_index, "cherenkov_index[cherenkov_size]/b");
  t->Branch("cherenkov_coarse", &cherenkov_coarse, "cherenkov_coarse[cherenkov_size]/b");
  t->Branch("cherenkov_fine", &cherenkov_fine, "cherenkov_fine[cherenkov_size]/b");
  t->Branch("cherenkov_tdc", &cherenkov_tdc, "cherenkov_tdc[cherenkov_size]/b");
  t->Branch("subframe_size", &subframe_size, "subframe_size/i");
  t->Branch("subframe_n", &subframe_n, "subframe_n[frame_n]/i");
  t->Branch("subframe_index", &subframe_index, "subframe_index[subframe_size]/i");
}
  
void
lightio::read_from_tree(std::string filename, std::string treename)
{
  file = TFile::Open(filename.c_str());
  tree = (TTree *)file->Get(treename.c_str());
  read_from_tree(tree);
}

void
lightio::read_from_tree(TTree *t)
{
  t->SetBranchAddress("part_n", &part_n);
  t->SetBranchAddress("part_device", &part_device);
  t->SetBranchAddress("part_mask", &part_mask);
  t->SetBranchAddress("dead_n", &dead_n);
  t->SetBranchAddress("dead_device", &dead_device);
  t->SetBranchAddress("dead_mask", &dead_mask);
  t->SetBranchAddress("frame_n", &frame_n);
  t->SetBranchAddress("frame", &frame);
  t->SetBranchAddress("trigger0_size", &trigger0_size);
  t->SetBranchAddress("trigger0_n", &trigger0_n);
  t->SetBranchAddress("trigger0_coarse", &trigger0_coarse);
  t->SetBranchAddress("timing_size", &timing_size);
  t->SetBranchAddress("timing_n", &timing_n);
  t->SetBranchAddress("timing_device", &timing_device);
  t->SetBranchAddress("timing_index", &timing_index);
  t->SetBranchAddress("timing_coarse", &timing_coarse);
  t->SetBranchAddress("timing_fine", &timing_fine);
  t->SetBranchAddress("timing_tdc", &timing_tdc);
  t->SetBranchAddress("cherenkov_size", &cherenkov_size);
  t->SetBranchAddress("cherenkov_n", &cherenkov_n);
  t->SetBranchAddress("cherenkov_device", &cherenkov_device);
  t->SetBranchAddress("cherenkov_index", &cherenkov_index);
  t->SetBranchAddress("cherenkov_coarse", &cherenkov_coarse);
  t->SetBranchAddress("cherenkov_fine", &cherenkov_fine);
  t->SetBranchAddress("cherenkov_tdc", &cherenkov_tdc);
  t->SetBranchAddress("subframe_size", &subframe_size);
  t->SetBranchAddress("subframe_n", &subframe_n);
  t->SetBranchAddress("subframe_index", &subframe_index);
}
  
bool
lightio::next_spill()
{
  if (spill_current >= tree->GetEntries())
    return false;
  
  tree->GetEntry(spill_current);
  frame_current = 0;
  trigger0_offset = 0;
  timing_offset = 0;
  cherenkov_offset = 0;
  subframe_offset = 0;
  
  ++spill_current;
  return true;
}

bool
lightio::next_frame()
{
  if (frame_current >= frame_n)
    return false;

// fill subframe vector
  subframe_vector.clear();
  for (int i = 0; i < subframe_n[frame_current]; ++i) {
    auto ii = subframe_offset + i;
    subframe_vector.push_back(subframe_index[ii]);
  }
  subframe_offset += subframe_n[frame_current];
  
  // fill trigger0 vector
  trigger0_vector.clear();
  for (int i = 0; i < trigger0_n[frame_current]; ++i) {
    auto ii = trigger0_offset + i;
    trigger0_vector.push_back(lightdata(0, 0, trigger0_coarse[ii], 0, 0));
  }
  trigger0_offset += trigger0_n[frame_current];
  
  // fill timing vector and map
  timing_vector.clear();
  timing_map.clear();
  for (int i = 0; i < timing_n[frame_current]; ++i) {
    auto ii = timing_offset + i;
    timing_vector.push_back(lightdata(timing_device[ii], timing_index[ii], timing_coarse[ii], timing_fine[ii], timing_tdc[ii]));
    timing_map[{timing_device[ii], timing_index[ii]}].push_back(lightdata(timing_device[ii], timing_index[ii], timing_coarse[ii], timing_fine[ii], timing_tdc[ii]));
  }
  timing_offset += timing_n[frame_current];
  
  // fill cherenkov vector and map
  cherenkov_vector.clear();
  cherenkov_map.clear();
  for (int i = 0; i < cherenkov_n[frame_current]; ++i) {
    auto ii = cherenkov_offset + i;
    cherenkov_vector.push_back(lightdata(cherenkov_device[ii], cherenkov_index[ii], cherenkov_coarse[ii], cherenkov_fine[ii], cherenkov_tdc[ii]));
    cherenkov_map[{cherenkov_device[ii], cherenkov_index[ii]}].push_back(lightdata(cherenkov_device[ii], cherenkov_index[ii], cherenkov_coarse[ii], cherenkov_fine[ii], cherenkov_tdc[ii]));
  }
  cherenkov_offset += cherenkov_n[frame_current];

  ++frame_current;
  return true;
}

} /** namespace sipm4eic **/

