#pragma once
//  --- Includes  ---
#include "framer.h"
//
//  Addition of fine tune class
namespace sipm4eic
{
  std::vector<std::string> timing_device = {"kc705-207"};
  //
  std::map<int, std::string> device_int2string = {
      {192, "kc705-192"},
      {193, "kc705-193"},
      {194, "kc705-194"},
      {195, "kc705-195"},
      {196, "kc705-196"},
      {197, "kc705-197"},
      {198, "kc705-198"},
      {207, "kc705-207"}};
  //
  std::vector<std::string> devices = {
      "kc705-192",
      "kc705-193",
      "kc705-194",
      "kc705-195",
      "kc705-196",
      "kc705-197",
      "kc705-198",
      "kc705-207"};
  //
  //  CSV calibration file
  const std::vector<std::string> calibration_csv_fields = {"device", "cindex", "IF", "CUT", "OFFSET"};
  //
  //! --- Convertion
  //! ---    --- Coarse
  const double coarse_to_s = 3.1250000e-09;
  const double coarse_to_ms = 3.1250000e-06;
  const double coarse_to_us = 3.1250000e-03;
  const double coarse_to_ns = 3.1250000;
  //! --- --- Rollover
  const int rollover_to_coarse = 32768;
  const double rollover_to_s = 0.0001024;
  const double rollover_to_ms = 0.1024;
  const double rollover_to_us = 102.4;
  const double rollover_to_ns = 102400.;
  //
  //  --- [0] IF    --- [1] CUT   --- [2] OFFSET
  typedef std::tuple<float, float, float> calibration_unit;
  //  --- Fit function on the fine distribution
  TF1 *fine_analysis_fit_function = new TF1("fine_analysis_fit_function", "[0]*(1./(exp((x-[1])/[2])+1))*(1./(exp(([3]-x)/[4])+1))", 10., 200.);
  //
  class finecalibration
  {
  public:
    //  --- Variables ---
    //  calibration_scheme[device][cindex]  [0] IF
    //                                      [1] CUT
    //                                      [2] OFFSET
    std::map<std::string, std::map<int, calibration_unit>> calibration_scheme;

    //  --- Methods ---
    //  --- --- Calculate phase based on current calibration scheme
    float get_phase(std::string device, int cindex, float fine);
    float calculate_phase_(sipm4eic::calibration_unit current_calibration_unit, float current_fine);

    //  --- --- I/O calibration scheme from file
    void load_calibration(std::string filename);
    void dump_calibration(std::string filename);
    void read_TDC_dist_file(std::string filename, std::vector<std::string> target_devices);

    //  --- --- Calibration utilities
    std::map<std::string, TH2F *> make_fine_TDC_dist_(std::string dirname, std::vector<std::string> target_devices, std::string outfilename = "finedata.root", unsigned int max_spill = kMaxUInt, int frame_size = 256);
    void calib_over_TDC_dist_(std::map<std::string, TH2F *> TDC_dist, std::vector<std::string> target_devices = devices);

    //  --- --- Getters
    calibration_unit get_calibration_unit(std::string device, int cindex) { return calibration_scheme[device][cindex]; }
    float get_if(std::string device, int cindex) { return get<0>(get_calibration_unit(device, cindex)); }
    float get_cut(std::string device, int cindex) { return get<1>(get_calibration_unit(device, cindex)); }
    float get_offset(std::string device, int cindex) { return get<2>(get_calibration_unit(device, cindex)); }

    //  --- --- Setters
    void set_calibration_unit(std::string device, int cindex, calibration_unit current_calibration = {0, 0, 0}) { calibration_scheme[device][cindex] = current_calibration; }
    void set_if(std::string device, int cindex, float current_if) { get<0>(calibration_scheme[device][cindex]) = current_if; }
    void set_cut(std::string device, int cindex, float current_cut) { get<1>(calibration_scheme[device][cindex]) = current_cut; }
    void set_offset(std::string device, int cindex, float current_offset) { get<2>(calibration_scheme[device][cindex]) = current_offset; }

    //  --- Generators ---
    finecalibration(){};

    //  --- --- --- Ongoing R&D DO NOT USE
    std::map<std::string, TH3F *> make_fine_DeltaT_dist_(std::string dirname, std::vector<std::string> target_devices, std::string outfilename = "finedata.root", unsigned int max_spill = kMaxUInt, int frame_size = 256);
    void calib_over_DeltaT_dist_(std::map<std::string, TH3F *> DeltaT_dist, std::vector<std::string> target_devices);
    finecalibration(std::string filename) { load_calibration(filename); };

  private:
  };
}
//
//  --- Definitions ---
//  --- --- Calculate phase based on current calibration scheme
float sipm4eic::finecalibration::get_phase(std::string device, int cindex, float fine)
{
  return calculate_phase_(this->get_calibration_unit(device, cindex), fine);
}
float sipm4eic::finecalibration::calculate_phase_(sipm4eic::calibration_unit current_calibration_unit, float current_fine)
{
  auto IF = get<0>(current_calibration_unit);
  auto CUT = get<1>(current_calibration_unit);
  auto OFFSET = get<2>(current_calibration_unit);
  auto MIN = CUT - 0.5 * IF;
  auto MAX = CUT + 0.5 * IF;
  auto result = -10.;
  if ((current_fine < MIN) || (current_fine > MAX))
  {
    //  VERBOSE
    //  cout << << endl;
    return -100.;
  }
  result = (current_fine - MIN) / IF;
  if (current_fine > CUT)
    result -= 1;
  result += OFFSET;
  return result;
}
//
//  --- --- I/O calibration scheme from file
void sipm4eic::finecalibration::load_calibration(std::string filename)
{
  std::ifstream data_stream(filename);
  std::string current_line;
  while (std::getline(data_stream, current_line))
  {
    //  Skip comment characters
    if (current_line[0] == '#' || current_line[0] == ' ')
      continue;
    //  Read database
    std::stringstream string_in_stream(current_line);
    std::map<std::string, std::string> data_by_field;
    std::string current_data;
    for (auto current_field : sipm4eic::calibration_csv_fields)
    {
      string_in_stream >> current_data;
      data_by_field[current_field] = current_data;
    }
    this->set_calibration_unit(data_by_field["device"], std::stoi(data_by_field["cindex"]), {std::stof(data_by_field["IF"]), std::stof(data_by_field["CUT"]), std::stof(data_by_field["OFFSET"])});
  }
};
void sipm4eic::finecalibration::dump_calibration(std::string filename)
{
  std::ofstream data_stream(filename);
  data_stream << "#\tCalibration file for Fine Tune\n";
  auto current_time = std::chrono::system_clock::now();
  std::time_t current_time_stamp = std::chrono::system_clock::to_time_t(current_time);
  data_stream << "#\tGenerated on " << std::ctime(&current_time_stamp);
  data_stream << "#\t## Device ## Calibration index ## IF ## CUT ## OFFSET ##\n";
  for (auto [device, calibration_of_device] : calibration_scheme)
  {
    for (auto [cindex, calibration_unit_] : calibration_of_device)
    {
      data_stream << device << "\t" << cindex << "\t" << get<0>(calibration_unit_) << "\t" << get<1>(calibration_unit_) << "\t" << get<2>(calibration_unit_) << "\n";
    }
  }
};
void sipm4eic::finecalibration::read_TDC_dist_file(std::string filename, std::vector<std::string> target_devices)
{
  //  Open file with TDC fine distributions
  TFile *input_file = new TFile(filename.c_str());
  std::map<std::string, TH2F *> input_TDC_dist;
  for (auto [current_device_index, current_device] : sipm4eic::device_int2string)
  {
    //  Take 2D TDC distribution for fit
    TH2F *current_calib_histo = (TH2F *)(input_file->Get(Form("hFine_%i", current_device_index)));
    if (!current_calib_histo)
      continue;

    input_TDC_dist[current_device] = (TH2F *)current_calib_histo->Clone();
  }
  calib_over_TDC_dist_(input_TDC_dist);
};
//
//  --- --- Calibration utilities
std::map<std::string, TH2F *> sipm4eic::finecalibration::make_fine_TDC_dist_(std::string dirname, std::vector<std::string> target_devices, std::string outfilename = "finedata.root", unsigned int max_spill = kMaxUInt, int frame_size = 256)
{
  //  Output file
  std::map<std::string, TH2F *> h_fine_device;
  std::vector<std::string> filenames;
  for (auto device : target_devices)
  {
    for (int ififo = 0; ififo < 25; ++ififo)
    {
      std::string filename = dirname + "/" + device + "/decoded/alcdaq.fifo_" + std::to_string(ififo) + ".root";
      filenames.push_back(filename);
    }
  }
  //  Framer initialisation
  std::cout << " --- initialize framer: frame size = " << frame_size << std::endl;
  sipm4eic::framer framer(filenames, frame_size);
  framer.set_trigger_coarse_offset(192, 112);
  //  Loop over spills
  int n_spills = 0, n_frames = 0;
  for (int ispill = 0; ispill < max_spill && framer.next_spill(); ++ispill)
  {
    std::cout << " --- new spill: " << ispill << std::endl;
    //  Loop over spills
    for (auto &[iframe, aframe] : framer.frames())
    {
      //  Loop over hits
      for (auto &[idevice, adevice] : aframe)
      {
        for (auto &[ichip, achip] : adevice.hits)
        {
          for (auto &[ichannel, hits] : achip)
          {
            for (auto &hit : hits)
            {
              auto device = idevice;
              if (!h_fine_device.count(device_int2string[device]))
                h_fine_device[device_int2string[device]] = new TH2F(Form("hFine_%d", device), "hFine", 768, 0, 768, 256, 0, 256);
              auto fine = hit.fine;
              auto index = hit.device_index();
              auto tdc = hit.tdc;
              auto cindex = tdc + 4 * index;
              h_fine_device[device_int2string[device]]->Fill(cindex, fine);
            } //  End of Hits loop
          }   //  End of Channel loop
        }     //  End of Chip loop
      }       //  End of Device loop
    }         //  End of Frame loop
    std::cout << "     spill completed " << std::endl;
  } //  End of Spill loop
  // Output
  auto fout = TFile::Open(outfilename.c_str(), "RECREATE");
  for (auto &h : h_fine_device)
    h.second->Write();
  fout->Close();
  return h_fine_device;
}
void sipm4eic::finecalibration::calib_over_TDC_dist_(std::map<std::string, TH2F *> TDC_dist, std::vector<std::string> target_devices = devices)
{
  //  Loop over available TDC distributions
  for (auto [current_device, current_calib_histo] : TDC_dist)
  {
    //  Only consider target devices
    if (std::find(target_devices.begin(), target_devices.end(), current_device) == target_devices.end())
      continue;
    //  Take 2D TDC distribution for fit
    for (auto iBin = 1; iBin <= current_calib_histo->GetNbinsX(); iBin++)
    {
      //  Take slice to act on single TDC
      auto current_histo = current_calib_histo->ProjectionY("tmp", iBin, iBin);
      //  Book results variable
      auto IF = 0.;
      auto CUT = 0.;
      auto OFFSET = 0.;
      //  If
      if (current_histo->GetEntries() < 1)
      {
        //  Assign dummy values in the calibration object
        this->set_calibration_unit(current_device, iBin, {IF, CUT, OFFSET});
        continue;
      }
      //  Setup fit function
      fine_analysis_fit_function->SetParameter(0, current_histo->GetBinContent(70));
      fine_analysis_fit_function->SetParameter(1, 100);
      fine_analysis_fit_function->SetParLimits(1, 80, 120);
      fine_analysis_fit_function->SetParameter(2, 0.5);
      fine_analysis_fit_function->SetParLimits(2, 0.1, 1.);
      fine_analysis_fit_function->SetParameter(3, 35);
      fine_analysis_fit_function->SetParLimits(3, 15, 55);
      fine_analysis_fit_function->SetParameter(4, 0.5);
      fine_analysis_fit_function->SetParLimits(4, 0.1, 1.);
      //  Fit Fine distribution
      current_histo->Fit(fine_analysis_fit_function, "IMRESQ", "", 20, 120);
      //  Recover MIN and MAX from fit
      auto maximum = fine_analysis_fit_function->GetParameter(1);
      auto minimum = fine_analysis_fit_function->GetParameter(3);
      //  Calculate IF and CUT
      IF = maximum - minimum;
      CUT = 0.5 * (minimum + maximum);
      //  Assign the values in the calibration object
      this->set_calibration_unit(current_device, iBin, {IF, CUT, 0});
      //  VEROBSE
    }
  }
}

//
//  --- --- --- Ongoing R&D DO NOT USE
/*
std::map<std::string, TH3F *> sipm4eic::finecalibration::make_fine_DeltaT_dist_(std::string dirname, std::vector<std::string> target_devices, std::string outfilename = "finedata.root", unsigned int max_spill = kMaxUInt, int frame_size = 256)
{
  //  Output file
  std::map<std::string, TH3F *> h_fine_device;
  std::vector<std::string> filenames;
  for (auto device : target_devices)
  {
    for (int ififo = 0; ififo < 25; ++ififo)
    {
      std::string filename = dirname + "/" + device + "/decoded/alcdaq.fifo_" + std::to_string(ififo) + ".root";
      filenames.push_back(filename);
    }
  }
  //  Framer initialisation
  std::cout << " --- initialize framer: frame size = " << frame_size << std::endl;
  sipm4eic::framer framer(filenames, frame_size);
  framer.set_trigger_coarse_offset(192, 112);
  //  Loop over spills
  int n_spills = 0, n_frames = 0;
  for (int ispill = 0; ispill < max_spill && framer.next_spill(); ++ispill)
  {
    //  Loop over frames
    for (auto &[iframe, aframe] : framer.frames())
    {
      //  Loop over devices
      for (auto &[idevice, adevice] : aframe)
      {
        if (!h_fine_device.count(device_int2string[idevice]))
          h_fine_device[device_int2string[idevice]] = new TH3F(Form("hFine_%d", idevice), "hFine", 768, 0, 768, 256, 0, 256, 1000, -10, 10);
        //  Loop over chips
        for (auto &[ichip, achip] : adevice.hits)
        {
          //  Loop over channels
          for (auto &[ichannel, hits] : achip)
          {
            auto hit = hits[0];
            //  Loop over hits
            //  for (auto &hit : hits)
            //{
            auto current_index = hit.device_index();
            auto current_tdc = hit.tdc;
            auto current_cindex = current_tdc + 4 * current_index;
            auto current_rollover = hit.rollover;
            auto current_coarse = hit.coarse;
            auto current_fine = hit.fine;
            auto current_time = current_rollover * sipm4eic::rollover_to_ns + current_coarse * sipm4eic::coarse_to_ns;
            // auto current_deltaT = current_time - timing_signal;
            // h_fine_device[device_int2string[idevice]]->Fill(current_cindex, current_fine, current_deltaT);
            // cout << "timing_signal: " << timing_signal << "current_time: " << current_time << " current_cindex: " << current_cindex << " current_fine: " << current_fine << " current_deltaT: " << current_deltaT << endl;
            // } //  End of Hits loop
          } //  End of Channel loop
        }   //  End of Chip loop
      }     //  End of Device loop
    }       //  End of Frame loop
    std::cout << "     spill completed " << std::endl;
  } //  End of Spill loop
  // Output
  auto fout = TFile::Open("finedata.root", "RECREATE");
  for (auto &h : h_fine_device)
    h.second->Write();
  fout->Close();
  return h_fine_device;
}

std::map<std::string, TH3F *> sipm4eic::finecalibration::make_fine_DeltaT_dist_(std::string dirname, std::vector<std::string> target_devices, std::string outfilename = "finedata.root", unsigned int max_spill = kMaxUInt, int frame_size = 256)
{
  //  Output file
  std::map<std::string, TH3F *> h_fine_device;
  std::vector<std::string> filenames;
  for (auto device : target_devices)
  {
    for (int ififo = 0; ififo < 25; ++ififo)
    {
      std::string filename = dirname + "/" + device + "/decoded/alcdaq.fifo_" + std::to_string(ififo) + ".root";
      filenames.push_back(filename);
    }
  }
  //  Framer initialisation
  std::cout << " --- initialize framer: frame size = " << frame_size << std::endl;
  sipm4eic::framer framer(filenames, frame_size);
  framer.set_trigger_coarse_offset(192, 112);
  //  Loop over spills
  int n_spills = 0, n_frames = 0;
  for (int ispill = 0; ispill < max_spill && framer.next_spill(); ++ispill)
  {
    //  Loop over frames
    for (auto &[iframe, aframe] : framer.frames())
    {
      bool skipevent_bad_timing = true;
      std::map<int, float> chip_timing_signal;
      float timing_signal = 0.;
      //  Loop on timing frames
      for (auto &[ichip, achip] : aframe[207].hits)
      {
        //  Timing information per chip
        int channel_check = 0;
        chip_timing_signal[ichip] = 0;
        for (auto &[ichannel, hits] : achip)
        {
          auto hit = hits[0];
          channel_check++;
          auto current_index = hit.device_index();
          auto current_tdc = hit.tdc;
          auto current_cindex = current_tdc + 4 * current_index;
          auto current_rollover = hit.rollover;
          auto current_coarse = hit.coarse;
          auto current_fine = hit.fine;
          auto ref_time = current_rollover * sipm4eic::rollover_to_ns + (current_coarse - this->get_phase("kc705-207", current_cindex, current_fine)) * sipm4eic::coarse_to_ns;
          chip_timing_signal[ichip] += ref_time;
        } //  End of Channel loop
        skipevent_bad_timing *= (channel_check == 32);
        //  Calculate average time for timing chips
        chip_timing_signal[ichip] /= 32;
      }
      if (skipevent_bad_timing)
        continue;
      auto number_of_chips = 0;
      auto average_timing_of_chips = 0;
      for (auto [current_chip, current_timing] : chip_timing_signal)
      {
        number_of_chips++;
        average_timing_of_chips += current_timing;
      }
      timing_signal = average_timing_of_chips / number_of_chips;
      //  Loop over devices
      for (auto &[idevice, adevice] : aframe)
      {
        if (idevice == 207)
          continue;
        if (!h_fine_device.count(device_int2string[idevice]))
          h_fine_device[device_int2string[idevice]] = new TH3F(Form("hFine_%d", idevice), "hFine", 768, 0, 768, 256, 0, 256, 1000, -10, 10);
        for (auto &[ichip, achip] : adevice.hits)
        {
          for (auto &[ichannel, hits] : achip)
          {
            auto hit = hits[0];
            // for (auto &hit : hits)
            //{
            auto current_index = hit.device_index();
            auto current_tdc = hit.tdc;
            auto current_cindex = current_tdc + 4 * current_index;
            auto current_rollover = hit.rollover;
            auto current_coarse = hit.coarse;
            auto current_fine = hit.fine;
            auto current_time = current_rollover * sipm4eic::rollover_to_ns + current_coarse * sipm4eic::coarse_to_ns;
            auto current_deltaT = current_time - timing_signal;
            h_fine_device[device_int2string[idevice]]->Fill(current_cindex, current_fine, current_deltaT);
            cout << "timing_signal: " << timing_signal << "current_time: " << current_time << " current_cindex: " << current_cindex << " current_fine: " << current_fine << " current_deltaT: " << current_deltaT << endl;
            //} //  End of Hits loop
          } //  End of Channel loop
        }   //  End of Chip loop
      }     //  End of Device loop
    }       //  End of Frame loop
    std::cout << "     spill completed " << std::endl;
  } //  End of Spill loop
  // Output
  auto fout = TFile::Open("finedata.root", "RECREATE");
  for (auto &h : h_fine_device)
    h.second->Write();
  fout->Close();
  return h_fine_device;
}
void sipm4eic::finecalibration::calib_over_DeltaT_dist_(std::map<std::string, TH3F *> DeltaT_dist, std::vector<std::string> target_devices)
{
  //  Loop over available TDC distributions
  for (auto [current_device, current_calib_histo] : DeltaT_dist)
  {
    //  Only consider target devices
    if (std::find(target_devices.begin(), target_devices.end(), current_device) == target_devices.end())
      continue;
    //  Take 2D TDC distribution for fit
    for (auto iBin = 1; iBin <= current_calib_histo->GetNbinsX(); iBin++)
    {
      //  Take slice to act on single TDC
      auto current_histo = current_calib_histo->ProjectionY("tmp", iBin, iBin);
      //  Book results variable
      auto IF = 0.;
      auto CUT = 0.;
      auto OFFSET = 0.;
      //  If
      if (current_histo->GetEntries() < 1)
      {
        //  Assign dummy values in the calibration object
        this->set_calibration_unit(current_device, iBin, {IF, CUT, OFFSET});
        continue;
      }
      //  Setup fit function
      fine_analysis_fit_function->SetParameter(0, current_histo->GetBinContent(70));
      fine_analysis_fit_function->SetParameter(1, 100);
      fine_analysis_fit_function->SetParLimits(1, 80, 120);
      fine_analysis_fit_function->SetParameter(2, 0.5);
      fine_analysis_fit_function->SetParLimits(2, 0.1, 1.);
      fine_analysis_fit_function->SetParameter(3, 35);
      fine_analysis_fit_function->SetParLimits(3, 15, 55);
      fine_analysis_fit_function->SetParameter(4, 0.5);
      fine_analysis_fit_function->SetParLimits(4, 0.1, 1.);
      //  Fit Fine distribution
      current_histo->Fit(fine_analysis_fit_function, "IMRESQ", "", 20, 120);
      //  Recover MIN and MAX from fit
      auto maximum = fine_analysis_fit_function->GetParameter(1);
      auto minimum = fine_analysis_fit_function->GetParameter(3);
      //  Calculate IF and CUT
      IF = maximum - minimum;
      CUT = 0.5 * (minimum + maximum);
      //  Assign the values in the calibration object
      this->set_calibration_unit(current_device, iBin, {IF, CUT, 0});
      //  VEROBSE
      //  cout << "Device: " << current_device.c_str() << " - iBin: " << iBin << " - minimum: " << minimum << " - maximum: " << maximum << " - IF: " << IF << " - CUT: " << CUT << endl;
    }
  }
}
//
*/