#pragma once

namespace sipm4eic
{

  class data {

    /** the actual data **/
    
  public:
    
    int device;
    int fifo;
    int type;
    int counter;
    int column;
    int pixel;
    int tdc;
    int rollover;
    int coarse;
    int fine;

    bool operator<(const data &rhs) const { return fine_time_clock() < rhs.fine_time_clock(); };

    /** calibration **/

    static double fine_min[768];
    static double fine_max[768];
    static double fine_off[768];
    static bool load_fine_calibration(std::string filename);

    static bool close_to_cut(int _index, int _fine, int dist = 1) {
      if (_fine == 0.) return false;      
      double min = fine_min[_index];
      double max = fine_max[_index];
      if (min == 0. || max == 0.) return false;
      double cut = 0.5 * (max + min) - 0.5;
      if (std::fabs(_fine - cut) < dist) return true;
      return false;
    }
    
    static double fine_phase(int _index, int _fine, bool _nice = true) {
      if (_fine == 0.) return 0.;
      double min = fine_min[_index];
      double max = fine_max[_index];
      if (min == 0. || max == 0.) return 0.;
      double phase = (_fine - min) / (max - min);
      if (_nice && _fine >= 0.5 * (max + min)) phase -= 1.;
      return phase;
    }

    static double fine_offset(int index) {
      return fine_off[index];
    }

    bool close_to_cut() const {
      auto index = calib_index();
      return close_to_cut(index, fine);
    }
    
    double fine_phase() const {
      auto index = calib_index();
      return fine_phase(index, fine);
    }

    double fine_offset() const {
      auto index = calib_index();
      return fine_offset(index);
    }
    
    /** indices **/

    //    static const int eo2do[32];

    int chip() const { return fifo / 4; };
    int eo_channel() const { return pixel + 4 * column; };
    //    int do_channel() const { return eo2do[eo_channel()]; };
    int calib_index() const { return tdc + 4 * pixel + 16 * column + 128 * chip(); };
    int device_index() const { return eo_channel() + 32 * chip(); };
    int global_index() const { return (device - 192) * 256 + device_index(); };
    
    /** conversion **/

    static const int rollover_to_clock;
    static const double coarse_to_ns;
    static const double rollover_to_ns;

    int coarse_time_clock() const { return coarse + rollover * rollover_to_clock; };
    double coarse_time_ns() const { return coarse * coarse_to_ns + rollover * rollover_to_ns; };

    double fine_time_clock() const { return coarse_time_clock() - fine_phase(); };
    double fine_time_ns() const { return coarse_time_ns() - fine_phase() * coarse_to_ns - fine_offset(); };
    
    
    /** hit type **/
    
    enum type {
      alcor_hit   = 1,
      trigger_tag = 9,
      start_spill = 7,
      end_spill   = 15
    };

    bool is_alcor_hit()   const { return type == alcor_hit; };
    bool is_trigger_tag() const { return type == trigger_tag; };
    bool is_start_spill() const { return type == start_spill; };
    bool is_end_spill()   const { return type == end_spill; };

    /** tree utils **/

    void link_to_tree(TTree *t);
    
  };
  
  //  const int data::eo2do[32] = {22, 20, 18, 16, 24, 26, 28, 30, 25, 27, 29, 31, 23, 21, 19, 17, 9, 11, 13, 15, 7, 5, 3, 1, 6, 4, 2, 0, 8, 10, 12, 14};
  
  const int data::rollover_to_clock = 32768;
  const double data::coarse_to_ns = 3.125;
  const double data::rollover_to_ns = 102400.;

  double data::fine_min[768] = {0.};
  double data::fine_max[768] = {0.};
  double data::fine_off[768] = {0.};

  void data::link_to_tree(TTree *t)
  {
    if (!t) return;
    t->SetBranchAddress("device", &device);
    t->SetBranchAddress("fifo", &fifo);
    t->SetBranchAddress("type", &type);
    t->SetBranchAddress("counter", &counter);
    t->SetBranchAddress("column", &column);
    t->SetBranchAddress("pixel", &pixel);
    t->SetBranchAddress("tdc", &tdc);
    t->SetBranchAddress("rollover", &rollover);
    t->SetBranchAddress("coarse", &coarse);
    t->SetBranchAddress("fine", &fine);
  }

  bool data::load_fine_calibration(std::string filename)
  {
    std::cout << " --- loading fine calibration: " << filename << std::endl;
    auto fin = TFile::Open(filename.c_str());
    auto hFine_min = (TH1 *)fin->Get("hFine_min");
    auto hFine_max = (TH1 *)fin->Get("hFine_max");
    auto hFine_off = (TH1 *)fin->Get("hFine_off");
    int found = 0;
    for (int i = 0; i < 768; ++i) {
      if (hFine_min->GetBinError(i + 1) <= 0. ||
          hFine_max->GetBinError(i + 1) <= 0.) continue;
      fine_min[i] = hFine_min ? hFine_min->GetBinContent(i + 1) : 0.;
      fine_max[i] = hFine_max ? hFine_max->GetBinContent(i + 1) : 0.;
      fine_off[i] = hFine_off ? hFine_off->GetBinContent(i + 1) : 0.;
    found++;
    }
    std::cout << " --- loaded fine calibration: found " << found << " channels " << std::endl;
    return true;
  }
  
} /** namespace sipm4eic **/
