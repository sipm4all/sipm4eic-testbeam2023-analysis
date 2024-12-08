#include "../lib/framer.h"

/*
    COMMENTS BY ANDREA N.

    In this file all the macros used to generate the figures included in
    my thesis are present. The data used is a combination of both data with
    exclusively Cherenkov events, and data which includes also background.
    Since generating the data in the necessary format from raw data (as is 
    performed by the lightwriter.C macro) is heavy, some macros were added 
    with the goal of generating the necessary data only once and save it in
    ROOT files. Inside those macros you will find the naming convention of each
    data object saved to the ROOT file.

    The macros which generate the necessary data expect the raw data to be in a 
    folder named "data" relative to the place where you run ROOT, otherwise you
    can specify the path.

    The generated ROOT files are saved in an output folder present in the project
    path (if it's not there you should create it).
    
    Note that many of the generation macros create an image of the generated 
    histograms themeselves, but the images used in the thesis are created by 
    the respective macros. All images are saved in the "/output/images" folder
    in the project path (you should create those folders beforehand as well).

*/

/* UTILITIES */

void print_hist(TH1 *hist) {
  std::cout << "\n## " << hist->GetName() << ": " << hist->GetTitle() << " ##";
  std::cout << "\n Mean: " << hist->GetMean() << " +/- " << hist->GetMeanError();
  std::cout << "\nRMS: " << hist->GetRMS();
  std::cout << "\nEntries: " << hist->GetEntries() << "\n";
}
void print_hist(std::string name) {
  TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects-F.root", "READ");  
  TH1I *hist = file->Get<TH1I>(name.c_str());
  print_hist(hist);
}

std::vector<std::string> build_filenames(std::string data_dir = "") {

  if (data_dir.empty()) {
    data_dir = "./data"; // default raw data folder
  }

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

  std::vector<std::string> filenames;
  for (auto device : devices) {
    for (int ififo = 0; ififo < 25; ++ififo) {
      std::string filename = data_dir + "/" + device + "/decoded/alcdaq.fifo_" + std::to_string(ififo) + ".root";
      filenames.push_back(filename);
    }
  }

  return filenames;
}

// this is the muber of different frame sizes considered in the analysis
const int N_SIZES = 13;

// this are the bin numbers for each frame size used in the respective histograms
const int N_BINS_VEC[N_SIZES] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 200, 200, 400 };

// these are the frame sizes currently used
const int FRAME_SIZE_VEC[N_SIZES] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096 };

// checks for subframes inside the given frame that have a number of hits higher than the given threshold
bool triggerless_selection(std::pair<const int, sipm4eic::framer::frame_t> &frame, int frame_size,
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

  for (int n_hits : subframe_hits) {
    if (n_hits >= min_subframe_hits) return true;
  }

  return false;
}
// checks for subframes inside the given frame that have a number of hits higher than the threshold
//    returns a vector with the indexes of subframes which satisfy the selection, if vector of 
//    subframe indexes returned is empty it means there was no triggerless selected subframe
std::vector<unsigned int> triggerless_selection_indexes(std::pair<const int, sipm4eic::framer::frame_t> &frame, int frame_size,
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

/* end of UTILITIES */

/* FIGURES GENERATION */

/*
    2D histograms showing the number of hits vs the time (for a single spill) and hits vs frame number for all the spills
    for frame size 512
    --- Thesis' figure [3.1]
*/
void generate_hits_vs_time_frame(std::string data_dir = "./data", unsigned int frame_size = 512, unsigned int max_spill = kMaxUInt) {

  // initialize framer
  std::cout << " --- initialize framer: frame size = " << frame_size << std::endl;
  sipm4eic::framer framer(build_filenames(), frame_size);
  framer.verbose(false);
  framer.set_trigger_coarse_offset(192, 112);

  // histogram for given frame size
  TH2I *hist = new TH2I("hits-vs-framenum", "", 7900, 0, 7900000, 100, 0, 100); // in 25 spills there are 7812500 frames of size 512

  TH2D *hist_1spill = new TH2D("hits-vs-time", "", 650, 0, 650, 100, 0, 100); // in 1 spill there are 600ms
  const double num_to_ms = 512 * 0.000003125;
  
  // max_spill = 1;

  // loop over spills
  int n_spills = 0, n_frames = 0;
  for (int ispill = 0;  ispill < max_spill && framer.next_spill(); ++ispill) {
    std::cout << "Analyzing spill " << ispill << std::endl;

    // loop over frames         
    for (auto &frame : framer.frames()) {
      auto iframe = frame.first;
      auto aframe = frame.second;

      // loop over cherenkov hits
      int n_hits = 0;
      for (auto &device : aframe) {
        auto idevice = device.first;
        auto adevice = device.second;

        if (idevice == 207) continue;
        
        for (auto &chip : adevice.hits) {
          auto achip = chip.second;
        
          for (auto &channel : achip) {
            auto hits = channel.second;

            for (auto &hit : hits) {
              ++n_hits;
            }
          }
        }
      }
      if (n_spills == 0) {
        hist_1spill->Fill(n_frames * num_to_ms, n_hits);
      }

      ++n_frames;
      hist->Fill(n_frames, n_hits);
    }
    ++n_spills;
  }
  std::cout << "\nFound " << n_frames << " frames in " << n_spills << " spills\n";

  hist->SetStats(0);
  hist->SetXTitle("Numero di frame");
  hist->SetYTitle("N_{hits}");
  hist->SetZTitle("N_{frames}"); 

  hist_1spill->SetStats(0);
  hist_1spill->SetXTitle("t (ms)");
  hist_1spill->SetYTitle("N_{hits}");
  hist_1spill->SetZTitle("N_{frames}"); 

  TCanvas *canvas = new TCanvas("canvas", nullptr);
  canvas->SetRightMargin(0.18); // otherwise palette axis title is not shown properly
  hist->Draw("COLZ");

  TCanvas *canvas2 = new TCanvas("canvas2", nullptr);
  canvas2->SetRightMargin(0.18); // otherwise palette axis title is not shown properly
  hist_1spill->Draw("COLZ");

  canvas2->SaveAs(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", hist_1spill->GetName()));


  // save histograms to data file
  TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "UPDATE");
  hist_1spill->Write(nullptr, TObject::kOverwrite);
  hist->Write(nullptr, TObject::kOverwrite);
  file->Close();
}
void hits_vs_time_frame() {
  // gStyle->SetPalette(kRainBow);

  // open the histograms data file
  TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "READ");

  // fetch histogram
  TH2I *hist_num = file->Get<TH2I>("hits-vs-framenum");
  TH2D *hist_time = file->Get<TH2D>("hits-vs-time");
  
  hist_num->SetXTitle("numero del frame");
  hist_num->SetYTitle("N_{hit}");
  hist_num->SetZTitle("N_{frame}");
  hist_num->SetStats(0);

  TGaxis::SetMaxDigits(4);
  TGaxis::SetExponentOffset(-0.03, -0.036, "x");

  hist_time->SetXTitle("t (ms)");
  hist_time->SetYTitle("N_{hit}");
  hist_time->SetZTitle("N_{frame}");
  hist_time->SetStats(0);

  TCanvas *canvas_num = new TCanvas("canvas-num", nullptr);
  canvas_num->SetRightMargin(0.2);
  hist_num->Draw("COLZ");

  TCanvas *canvas_time = new TCanvas("canvas-time", nullptr);
  canvas_time->SetRightMargin(0.2);
  hist_time->Draw("COLZ");

  // gPad->Update();
  canvas_num->SaveAs(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", hist_num->GetName()));
  canvas_time->SaveAs(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", hist_time->GetName()));
}

/*
    Canvas with 4 histograms showing the frames per hit number distributions for frame sizes 8, 256, 1024, 4096 
    Note: to run this macro first run "generate_mean_vs_size" which generates the hits distributions
    --- Thesis' figure [3.3] (used also for figure [3.6])
*/
void frames_vs_hitnum() {

  // in stats box print entries, mean, RMS
  gStyle->SetOptStat("emr");
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);

  // change stats box and font size
  gStyle->SetStatW(0.4);
  gStyle->SetStatH(0.15);

  // histograms data file
  TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "READ");
  
  // fetch the histograms for the given frame size
  TH1I *hist_8 = file->Get<TH1I>("frames-vs-hit-tot-size8-bin100");
  hist_8->SetTitle("Dimensione frame: 8;N_{hit};N_{frame}");
  hist_8->SetLineColor(kOrange);

  TH1I *hist_256 = file->Get<TH1I>("frames-vs-hit-tot-size256-bin100");
  hist_256->SetTitle("Dimensione frame: 256;N_{hit};N_{frame}");
  hist_256->SetLineColor(kMagenta+4);
  
  TH1I *hist_1024 = file->Get<TH1I>("frames-vs-hit-tot-size1024-bin200");
  hist_1024->SetTitle("Dimensione frame: 1024;N_{hit};N_{frame}");
  hist_1024->SetLineColor(kOrange);

  TH1I *hist_4096 = file->Get<TH1I>("frames-vs-hit-tot-size4096-bin400");
  hist_4096->SetTitle("Dimensione frame: 4096;N_{hit};N_{frame}");
  hist_4096->SetLineColor(kOrange);

  // draw
  TCanvas *canvas = new TCanvas("canvas", "");
  canvas->Divide(2,2);
  canvas->cd(1);
  gPad->SetLogy();
  hist_8->Draw();

  canvas->cd(2);
  gPad->SetLogy();
  hist_256->Draw();

  canvas->cd(3);
  gPad->SetLogy();
  hist_1024->Draw();

  canvas->cd(4);
  gPad->SetLogy();
  hist_4096->Draw();

  canvas->Print(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", "frames-vs-hitnum"));
}

/*
    Canvas with 2 histograms, one of the total distribution the other of the bkg distribution of frames per hit num., for size 256
    --- Thesis' figure [3.4]
*/
void frames_vs_hitnum_tot_bkg() {
    
  // in stats box print entries, mean, RMS
  gStyle->SetOptStat("emr");
  gStyle->SetPadLeftMargin(0.15);

  // change stats box and font size
  gStyle->SetStatW(0.4);
  gStyle->SetStatH(0.15);

  // histograms data file
  TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "READ");
  
  // fetch the histograms for the given frame size
  TH1I *hist_tot = file->Get<TH1I>("frames-vs-hit-tot-size256-bin100");
  hist_tot->SetNameTitle("frames-vs-hitnum-tot-256", "segnale + rumore;N_{hit};N_{frame}");
  hist_tot->SetLineColor(kOrange);

  TH1I *hist_bkg = file->Get<TH1I>("frames-vs-hit-bkg-size256-bin100");
  hist_bkg->SetNameTitle("frames-vs-hitnum-bkg-256", "rumore;N_{hit};N_{frame}");
  hist_bkg->SetLineColor(kRed);

  // draw
  TCanvas *canvas_tot = new TCanvas("canvas-tot", "");
  canvas_tot->SetLogy();
  hist_tot->Draw();

  TCanvas *canvas_bkg = new TCanvas("canvas-bkg", "");
  canvas_bkg->SetLogy();
  hist_bkg->Draw();

  print_hist(hist_tot);
  print_hist(hist_bkg);

  canvas_tot->Print(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", hist_tot->GetName()));
  canvas_bkg->Print(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", hist_bkg->GetName()));
}

/*
    Graph showing on x-axis the frame size and and on y-axis the hits distribution mean for the give frame size
  --- Thesis' figure [3.7-a]
*/
void generate_mean_vs_size(std::string data_dir = "./data") {
  
  // HITS DISTRIBUTION GEN

  for (int isize = 0; isize < N_SIZES; ++isize) {
    const int frame_size = FRAME_SIZE_VEC[isize];
    const int n_bins = N_BINS_VEC[isize];

    // initialize framer
    std::cout << " --- initialize framer: frame size = " << frame_size << std::endl;
    sipm4eic::framer framer(build_filenames(), frame_size);
    framer.verbose(false);
    framer.set_trigger_coarse_offset(192, 112);

    // background and signal histograms
    std::unique_ptr<TH1I> tot_hist(new TH1I(Form("frames-vs-hit-tot-size%i-bin%i", frame_size, n_bins), "",n_bins, 0, n_bins));
    std::unique_ptr<TH1I> bkg_hist(new TH1I(Form("frames-vs-hit-bkg-size%i-bin%i", frame_size, n_bins), "",n_bins, 0, n_bins));
    std::unique_ptr<TH1I> signal_hist(new TH1I(Form("frames-vs-hit-signal-size%i-bin%i", frame_size, n_bins), "",n_bins, 0, n_bins));

    // loop over spills
    int n_spills = 0, n_frames = 0, frames_acc = 0;
    const int spill_length = 0.5 * 320e6;
    int max_frames = spill_length / frame_size;
    std::cout << "\nMax frames: " << max_frames << "\n";
    for (int ispill = 0; framer.next_spill(); ++ispill) {
      std::cout << "Analyzing spill " << ispill << " - non zero frame count: "
                << framer.frames().size() <<  " - total frame count: " << max_frames << std::endl;

      frames_acc += framer.frames().size();

      // loop over frames         
      bool event_found = true;
      for (int iframe = 0; iframe < max_frames; ++iframe) {
        if (!framer.frames().count(iframe)) {
          tot_hist->Fill(0);
          bkg_hist->Fill(0);
        } else {
          auto aframe = framer.frames()[iframe];

          // filtering
          auto nsipm4 = aframe[207].hits[4].size();
          auto nsipm5 = aframe[207].hits[5].size();
          if ((aframe[192].triggers.size() != 1) || (nsipm4 == 0 && nsipm5 == 0)) {
            event_found = false;
          }

          // loop over cherenkov hits
          int n_hits = 0;
          for (auto &device : aframe) {
            auto idevice = device.first;
            auto adevice = device.second;

            if (idevice == 207) continue;

            for (auto &chip : adevice.hits) {
              auto achip = chip.second;

              for (auto &channel : achip) {
                auto hits = channel.second;

                for (auto &hit : hits) {
                  ++n_hits;
                }
              }
            }
          }

          if(event_found) {
            signal_hist->Fill(n_hits);
          }
          else {
            bkg_hist->Fill(n_hits);
          }
          tot_hist->Fill(n_hits);
        }
        event_found = true;
        ++n_frames;
      }
      ++n_spills;
    }
    std::cout << "\nFound " << n_frames << " frames of which " << frames_acc 
              << " are non zero, in "<< n_spills << " spills\n";
    std::cout << "Hits number mean of frame size " << frame_size << " is " << bkg_hist->GetMean() 
              << " for background and " << signal_hist->GetMean() << " for signal \n\n";

    print_hist(tot_hist.get());
    print_hist(signal_hist.get());
    print_hist(bkg_hist.get());

    // save histograms to data file
    TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "UPDATE");
    tot_hist->Write(nullptr, TObject::kOverwrite);
    bkg_hist->Write(nullptr, TObject::kOverwrite);
    signal_hist->Write(nullptr, TObject::kOverwrite);
    file->Close();
  }

  // GRAPH GENERATION

  TFile *charts_file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "UPDATE");

  // vectors to save data and error points (NB: TGraphError::AddPointError is not available in the latest ROOT release)
  std::vector<double> tot_vx;
  std::vector<double> tot_vy;
  std::vector<double> tot_vey;
  std::vector<double> tot_vey_RMS;

  std::vector<double> bkg_vx;
  std::vector<double> bkg_vy;
  std::vector<double> bkg_vey;
  std::vector<double> bkg_vey_RMS;

  std::vector<double> signal_vx;
  std::vector<double> signal_vy;
  std::vector<double> signal_vey;
  std::vector<double> signal_vey_RMS;

  // loop over frame sizes (we assume the histograms for each frame size are already generated)
  for (int isize = 0; isize < N_SIZES; ++isize) {
    const int frame_size = FRAME_SIZE_VEC[isize];
    const int n_bins = N_BINS_VEC[isize];
    std::unique_ptr<TH1I> tot_hist(charts_file->Get<TH1I>(Form("frames-vs-hit-tot-size%i-bin%i", frame_size, n_bins)));  
    std::unique_ptr<TH1I> bkg_hist(charts_file->Get<TH1I>(Form("frames-vs-hit-bkg-size%i-bin%i", frame_size, n_bins)));  
    std::unique_ptr<TH1I> signal_hist(charts_file->Get<TH1I>(Form("frames-vs-hit-signal-size%i-bin%i", frame_size, n_bins)));  

    // save mean hits as a graph point
    tot_vx.push_back(frame_size);
    tot_vy.push_back(tot_hist->GetMean());
    tot_vey.push_back(tot_hist->GetMeanError());
    tot_vey_RMS.push_back(tot_hist->GetRMS());

    bkg_vx.push_back(frame_size);
    bkg_vy.push_back(bkg_hist->GetMean());
    bkg_vey.push_back(bkg_hist->GetMeanError());
    bkg_vey_RMS.push_back(bkg_hist->GetRMS());

    signal_vx.push_back(frame_size);
    signal_vy.push_back(signal_hist->GetMean());
    signal_vey.push_back(signal_hist->GetMeanError());
    signal_vey_RMS.push_back(signal_hist->GetRMS());
  }

  // hits mean vs frame size error graphs
  TGraphErrors *tot_graph = new TGraphErrors(tot_vx.size(), tot_vx.data(), tot_vy.data(), 
                                              std::vector(tot_vx.size(), 0.).data(), tot_vey.data());
  tot_graph->SetNameTitle("mean-vs-size-tot", ""); 
  TGraphErrors *tot_graph_RMS = new TGraphErrors(tot_vx.size(), tot_vx.data(), tot_vy.data(), 
                                              std::vector(tot_vx.size(), 0.).data(), tot_vey_RMS.data());
  tot_graph_RMS->SetNameTitle("mean-vs-size-tot-RMS",""); 

  TGraphErrors *bkg_graph = new TGraphErrors(bkg_vx.size(), bkg_vx.data(), bkg_vy.data(), 
                                              std::vector(bkg_vx.size(), 0.).data(), bkg_vey.data());
  bkg_graph->SetNameTitle("mean-vs-size-bkg", ""); 
  TGraphErrors *bkg_graph_RMS = new TGraphErrors(bkg_vx.size(), bkg_vx.data(), bkg_vy.data(), 
                                              std::vector(bkg_vx.size(), 0.).data(), bkg_vey_RMS.data());
  bkg_graph_RMS->SetNameTitle("mean-vs-size-bkg-RMS",""); 

  TGraphErrors *signal_graph = new TGraphErrors(signal_vx.size(), signal_vx.data(), signal_vy.data(), 
                                              std::vector(signal_vx.size(), 0.).data(), signal_vey.data());
  signal_graph->SetNameTitle("mean-vs-size-signal", ""); 
  TGraphErrors *signal_graph_RMS = new TGraphErrors(signal_vx.size(), signal_vx.data(), signal_vy.data(), 
                                              std::vector(signal_vx.size(), 0.).data(), signal_vey_RMS.data());
  signal_graph_RMS->SetNameTitle("mean-vs-size-signal-RMS",""); 


  tot_graph->Write(nullptr, TObject::kOverwrite);
  tot_graph_RMS->Write(nullptr, TObject::kOverwrite);
  
  bkg_graph->Write(nullptr, TObject::kOverwrite);
  bkg_graph_RMS->Write(nullptr, TObject::kOverwrite);
  
  signal_graph->Write(nullptr, TObject::kOverwrite);
  signal_graph_RMS->Write(nullptr, TObject::kOverwrite);
  
  charts_file->Close();
}
void mean_vs_size() {

  // fetch hits mean vs frame size error graphs
  TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "READ");

  TGraphErrors *bkg_graph = file->Get<TGraphErrors>("mean-vs-size-bkg");
  bkg_graph->SetMarkerStyle(kCircle);
  bkg_graph->SetMarkerColor(kRed);
  bkg_graph->SetLineColor(kRed);
  bkg_graph->SetMarkerSize(0.5);
  
  TGraphErrors *signal_graph = file->Get<TGraphErrors>("mean-vs-size-signal");
  signal_graph->SetMarkerStyle(kCircle);
  signal_graph->SetMarkerColor(kBlue);
  signal_graph->SetLineColor(kBlue);
  signal_graph->SetMarkerSize(0.5);

  // create total graph
  TMultiGraph *multi_graph = new TMultiGraph("mean-vs-size-bkg-signal", "");
  multi_graph->Add(bkg_graph);
  multi_graph->Add(signal_graph);
  multi_graph->GetXaxis()->SetTitle("dimensione frame");
  multi_graph->GetYaxis()->SetTitle("<N_{hits}>");
  
  // create graph Legend
  TLegend *leg = new TLegend(0.15, 0.73, 0.4, 0.87);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(bkg_graph, "rumore", "P");
  leg->AddEntry(signal_graph, "segnale + rumore", "P");

  // draw
  TCanvas* canvas = new TCanvas("canvas", nullptr);
  // canvas->SetLogx();
  multi_graph->Draw("APZ");
  leg->Draw();

  // save graph as image
  canvas->Print(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", multi_graph->GetName()));

  // // RMS data visualization
  TGraphErrors *tot_graph_RMS = file->Get<TGraphErrors>("mean-vs-size-tot-RMS");
  tot_graph_RMS->SetMarkerStyle(kFullCircle);
  tot_graph_RMS->SetMarkerSize(0.5);
  tot_graph_RMS->SetMarkerColor(kOrange);
  tot_graph_RMS->SetLineColor(kOrange);
  
  TGraphErrors *bkg_graph_RMS = file->Get<TGraphErrors>("mean-vs-size-bkg-RMS");
  bkg_graph_RMS->SetMarkerStyle(kFullCircle);
  bkg_graph_RMS->SetMarkerSize(0.5);
  bkg_graph_RMS->SetMarkerColor(kRed);
  bkg_graph_RMS->SetLineColor(kRed);
  
  TGraphErrors *signal_graph_RMS = file->Get<TGraphErrors>("mean-vs-size-signal-RMS");
  signal_graph_RMS->SetMarkerStyle(kFullCircle);
  signal_graph_RMS->SetMarkerSize(0.5);
  signal_graph_RMS->SetMarkerColor(kBlue);
  signal_graph_RMS->SetLineColor(kBlue);

  // create total graph
  TMultiGraph *multi_graph_RMS = new TMultiGraph("mean-vs-size-bkg-signal-RMS","");
  multi_graph_RMS->Add(bkg_graph_RMS);
  multi_graph_RMS->Add(signal_graph_RMS);
  multi_graph_RMS->GetXaxis()->SetTitle("Dimensione frame");
  multi_graph_RMS->GetYaxis()->SetTitle("<N_{hit}>");
  
  // create graph Legend
  TLegend *leg_RMS = new TLegend(0.15, 0.73, 0.4, 0.87);
  leg_RMS->SetFillStyle(0);
  leg_RMS->SetBorderSize(0);
  leg_RMS->AddEntry(bkg_graph_RMS, "rumore", "P");
  leg_RMS->AddEntry(signal_graph_RMS, "segnale + rumore", "P");

  // draw
  TCanvas* canvas_RMS = new TCanvas("canvas-RMS", nullptr);
  // canvas_RMS->SetLogx();
  multi_graph_RMS->Draw("APZ");
  leg_RMS->Draw();

  // save graph as image
  canvas_RMS->Print(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", multi_graph_RMS->GetName()));
}

/*
    Graph showing on x-axis the frame size and on y-axis the difference between the hits means of signal and background/DCR model
    Note: to run "signal_DCR_diff" you first need to compute the fit model by running the macro "fit_bkg"
    --- Thesis' figure [3.7-b]
*/
void signal_bkg_diff() {
  
  // open the output ROOT file in READ mode
  TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "READ");

  // from the output ROOT file retrieve the background and signal graphs
  TGraphErrors *bkg_graph = file->Get<TGraphErrors>("mean-vs-size-bkg");
  TGraphErrors *signal_graph = file->Get<TGraphErrors>("mean-vs-size-signal");

  // vectors to save data and error points 
  std::vector<double> diff_vx;
  std::vector<double> diff_vy;
  std::vector<double> diff_vey;
  for (int i = 0; i < bkg_graph->GetN(); i++) {
    // compute error using error propagation (the two variables are considered statistically independent)
    double i_error = std::sqrt(std::pow(bkg_graph->GetErrorY(i), 2) + std::pow(signal_graph->GetErrorY(i), 2));

    diff_vx.push_back(bkg_graph->GetPointX(i));
    diff_vy.push_back(signal_graph->GetPointY(i) - bkg_graph->GetPointY(i));
    diff_vey.push_back(i_error);
  }

  // setup difference graph
  TGraphErrors *diff_graph = new TGraphErrors(diff_vx.size(), diff_vx.data(), diff_vy.data(),
                                                std::vector(diff_vx.size(), 0.).data(), diff_vey.data());

  diff_graph->SetNameTitle("mean-vs-size-diff", "");
  diff_graph->GetXaxis()->SetTitle("dimensione frame");
  diff_graph->GetYaxis()->SetTitle("diff. <N_{hit}>"); // TODO: Change the axis title
  diff_graph->SetMarkerStyle(kCircle);
  diff_graph->SetLineColor(kGreen);
  // diff_graph->SetMarkerSize(0.5);
  diff_graph->SetMarkerColor(kGreen);

  // draw
  TCanvas* canvas = new TCanvas("canvas", nullptr);
  canvas->SetLogx();
  diff_graph->Draw("APZ");

  canvas->Print(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", diff_graph->GetName()));
  // diff_graph->Write(nullptr, TObject::kOverwrite);
  // file->Close();
}
void signal_DCR_diff() {
  // NOTE: the DCR model is implemented using the slope from the fit performed on the background data points

  // retrieve data files
  TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "UPDATE");

  // fetch signal graph for 13 repetitions
  TGraphErrors *signal_graph = file->Get<TGraphErrors>("mean-vs-size-signal");

  // fetch background fit function and create DCR model
  TF1 *bkg_fit_func = file->Get<TF1>("mean-vs-size-bkg-fitfunc");
  TF1 *dcr_func = new TF1("mean-vs-size-bkg-DCRfunc", "[0]*x", 0., 4096.);
  dcr_func->SetParameter(0, bkg_fit_func->GetParameter(1));

  // compute the difference points
  std::vector<double> diff_vx;
  std::vector<double> diff_vy;
  std::vector<double> diff_vey;
  for (int i = 0; i < signal_graph->GetN(); i++) {

    // as error consider the fit parameter a error times the x (error propagation)
    double i_error = std::sqrt(std::pow(bkg_fit_func->GetParError(1) * signal_graph->GetPointX(i), 2) 
                                          + std::pow(signal_graph->GetErrorY(i), 2));

    diff_vx.push_back(signal_graph->GetPointX(i));
    diff_vy.push_back(signal_graph->GetPointY(i) - dcr_func->Eval(signal_graph->GetPointX(i)));
    diff_vey.push_back(i_error);
  }

  // setup difference graph
  // TGraphErrors *diff_graph = new TGraphErrors(diff_vx.size(), diff_vx.data(), diff_vy.data(),
  //                                               std::vector(diff_vx.size(), 0.).data(), diff_vey.data());
  TGraphErrors *diff_graph = new TGraphErrors(diff_vx.size(), diff_vx.data(), diff_vy.data());
  diff_graph->SetNameTitle("signal-DCR-diff", "");
  diff_graph->GetYaxis()->SetTitle("diff. <N_{hit}>"); // TODO: Change the axis title
  diff_graph->GetXaxis()->SetTitle("dimensione del frame");
  diff_graph->SetMarkerStyle(kCircle);
  diff_graph->SetLineColor(kBlue);
  diff_graph->SetMarkerSize(0.5);
  diff_graph->SetMarkerColor(kBlue);

  // draw
  TCanvas *canvas = new TCanvas("canvas", nullptr);
  canvas->SetLogx();
  diff_graph->Draw("APZ");

  canvas->Print(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", diff_graph->GetName()));
  diff_graph->Write(nullptr, TObject::kOverwrite);
}

/*
    Graph showing the background hits means per frame size with a linear fit
    --- Thesis' figure [3.5]
*/
void fit_bkg() {

  // on stats box print fit chi square and fit parameters with errors
  gStyle->SetOptFit(1);
  // change stats box and font size
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.15);

  // graphs data file
  TFile *objects_file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "UPDATE");

  // from the output ROOT file retrieve the background graph
  TGraphErrors *bkg_graph = objects_file->Get<TGraphErrors>("mean-vs-size-bkg");
  bkg_graph->SetName("mean-vs-size-bkg-fit");
  bkg_graph->SetMarkerColor(kRed);
  bkg_graph->SetMarkerStyle(kFullCircle);
  bkg_graph->SetLineColor(kRed);
  // bkg_graph->GetXaxis()->SetRangeUser(0, 8000);
  bkg_graph->GetXaxis()->SetTitle("dimensione del frame");
  bkg_graph->GetYaxis()->SetTitle("<N_{hit}>");
  // bkg_graph->GetYaxis()->SetRangeUser(0, 300);

  // fit graph with a linear polynomial: p0 + p1*x. The range includes points from 256 to 8192
  // NB: the NDF are calculated as num. points - num. params
  bkg_graph->Fit("pol1", "", "", 256, 4096);
  
  TF1* fit_func = bkg_graph->GetFunction("pol1");
  fit_func->SetLineColor(kViolet);
  fit_func->SetParName(0, "b");
  fit_func->SetParName(1, "a");
  
  // create graph Legend
  TLegend *leg = new TLegend(0.15, 0.77, 0.3, 0.87);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(bkg_graph, "rumore", "P");
  leg->AddEntry(fit_func, "fit", "L");

  // create "y= ax + b" label for fit line
  TLatex* fit_label = new TLatex();
  // fit_label->SetTextAngle(35);
  fit_label->SetTextSize(0.035);

  bkg_graph->GetXaxis()->SetRangeUser(0, 7000);

  // draw
  TCanvas* canvas = new TCanvas("canvas", nullptr);
  canvas->SetLogx();
  bkg_graph->Draw("APZ");
  leg->Draw();
  fit_label->DrawLatexNDC(0.7, 0.5, "y = ax + b");

  canvas->SaveAs(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", bkg_graph->GetName()));

  // save it to file with a different name and its associated fit function
  bkg_graph->Write("mean-vs-size-bkg-fit", TObject::kOverwrite);
  fit_func->Write("mean-vs-size-bkg-fitfunc", TObject::kOverwrite);
  
  objects_file->Close();
}

/*
    Graph showing on x-axis the frame size and on y-axis the SNR: the ratio between effective signal and background hits means
    Note: to compute the SNR you first need to create the difference graph by running the macro "signal_DCR_diff"
    --- Thesis' figure [3.8]
*/
void signal_to_noise_ratio() {
    // retrieve graphs files
  TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "READ");

  // fetch diff and background graph
  TGraphErrors *diff_graph = file->Get<TGraphErrors>("signal-DCR-diff");
  TGraphErrors *bkg_graph = file->Get<TGraphErrors>("mean-vs-size-bkg");

  // compute ratio and fill SNR graph
  TGraph * snr_graph = new TGraph();
  snr_graph->SetNameTitle("snr-graph", "");
  snr_graph->GetXaxis()->SetTitle("dimensione del frame");
  snr_graph->GetYaxis()->SetTitle("SNR"); // TODO: change axis title
  snr_graph->SetMarkerColor(kMagenta+4);
  snr_graph->SetMarkerStyle(kFullTriangleUp);
  // snr_graph->SetMarkerSize(0.8);

  for (int i = 0; i < diff_graph->GetN(); ++i) {
    snr_graph->AddPoint(diff_graph->GetPointX(i), 
                          diff_graph->GetPointY(i) / bkg_graph->GetPointY(i));
  }
  
  // draw
  TCanvas *canvas = new TCanvas("canvas", nullptr);
  canvas->SetLogx();
  snr_graph->Draw("APZ");

  canvas->Print(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", snr_graph->GetName()));
  file->Close();
}

/*
  Poissonian bkg fit on the frame size 8 hits distribution
  --- Thesis' figure [3.9]
*/
void optimal_poisson_fit() {
  // on stats box print fit chi square and fit parameters with errors
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  // change stats box and font size
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);


  TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "READ"); 
  TH1I *hist = file->Get<TH1I>("frames-vs-hit-tot-size8-bin100");
  hist->SetXTitle("N_{hit}");
  hist->SetYTitle("N_{frame}");
  hist->SetName("optimal-poisson-fit");
  // hist->Sumw2();

  TF1 *fit_func = new TF1("fit", "[0] * TMath::Poisson(x, [1])");
  fit_func->SetParameters(hist->GetMaximum(), hist->GetMean());
  fit_func->SetLineColor(kViolet);
  fit_func->SetParName(0, "A");
  fit_func->SetParName(1, "#lambda");

  hist->Fit(fit_func, "0");

  TH1D *dist_hist = new TH1D("dist-hist", "", 100, 0, 100);
  dist_hist->FillRandom("fit", hist->GetEntries());
  dist_hist->SetLineColor(kViolet);

  TLegend *leg = new TLegend(0.6, 0.6, 0.88, 0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hist, "segnale + rumore", "L");
  leg->AddEntry(dist_hist, "fit", "L");

  // create "y= ax + b" label for fit line
  TLatex* fit_label = new TLatex();
  // fit_label->SetTextAngle(35);
  fit_label->SetTextSize(0.05 );

  // draw
  TCanvas* canvas = new TCanvas("canvas", nullptr);
  canvas->SetLogy();
  hist->Draw();
  dist_hist->Draw("SAME");
  leg->Draw();
  fit_label->DrawLatexNDC(0.16, 0.7, "y = A#frac{#lambda^{n}}{n!}e^{-#lambda}");

  canvas->SaveAs(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", hist->GetName()));
}

/*
    histogram showing on x-axis the time in a selected frame of size 256 and on the y-axis the hits. num, with 256 bins
    --- Thesis' figure [3.2]
*/
void hits_in_frame() {
  gStyle->SetOptStat(0);

  const int frame_size = 256;
  const int subframe_size = 8;
  const int min_subframe_size = 5;

  const double coarse_to_ns = 3.125;

  // initialize framer
  std::cout << " --- initialize framer: frame size = " << frame_size << std::endl;
  sipm4eic::framer framer(build_filenames(), frame_size);
  framer.verbose(false);
  framer.set_trigger_coarse_offset(192, 112);

  // Histogram for the single frame
  TH1D *hist = new TH1D("hits-in-frame", "", 400, 0, 800);
  hist->SetLineColor(kBlue);
  hist->SetXTitle("t (ns)");
  hist->SetYTitle("N_{hit}");

  std::cout << "Analyzing first spill " << std::endl;

  int n_selected_frames = 0;
  framer.next_spill();
  for (auto &frame : framer.frames()) {

  // for the chosen frame select the 200th selected frame
    auto iframe = frame.first;
    auto aframe = frame.second;
    double frame_time = (iframe * frame_size / 320.e6);

    // filter selection
    bool filter_found = (triggerless_selection(frame, frame_size, subframe_size, min_subframe_size));
    if ((n_selected_frames == 205) && filter_found) {

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
              hist->Fill(coarse * coarse_to_ns);
            }
          }
        } 
      }

      break;
    }
    
    if (filter_found) {
      ++n_selected_frames;
    }
  }

  TCanvas* canvas = new TCanvas("canvas", "");
  hist->Draw();

  canvas->Print(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", hist->GetName()));
}

/*
    histograms comparing the number of selected frames with trigger and triggerless mode vs time
    Note: this uses the developed triggerless selection algorithm 
    --- Thesis' figure [3.10]
*/
void generate_selected_vs_time(unsigned int max_spill = kMaxUInt) {

  const int frame_size = 256;
  const int subframe_size = 8;
  const int min_subframe_size = 5;

  gStyle->SetOptStat("emr");

  // build decoded data filenames
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
  std::vector<std::string> filenames;
  for (auto device : devices) {
    for (int ififo = 0; ififo < 25; ++ififo) {
      std::string filename = "./data/" + device + "/decoded/alcdaq.fifo_" + std::to_string(ififo) + ".root";
      filenames.push_back(filename);
    }
  }

  // initialize framer
  std::cout << " --- initialize framer: frame size = " << frame_size << std::endl;
  sipm4eic::framer framer(filenames, frame_size);
  framer.verbose(false);
  framer.set_trigger_coarse_offset(192, 112);

  // histograms with number of selected frames per frame number
  TH1I *triggerless_hist = new TH1I("selected-vs-number-triggerless-size256-spills25-I", "triggerless frame selection;t (s);N_{frames}", 60, 0, 0.6);
  triggerless_hist->SetLineColor(kBlue);

  TH1I *trigger_hist = new TH1I("selected-vs-number-trigger-size256-spills25-I", "trigger frame selection;t (s);N_{frames}", 60, 0, 0.6);
  trigger_hist->SetLineColor(kOrange);

  // loop over spills
  int n_spills = 0, n_frames = 0;
  for (int ispill = 0;  ispill < max_spill && framer.next_spill(); ++ispill) {
    std::cout << "Analyzing spill " << ispill << std::endl;

    // loop over frames
    for (auto &frame : framer.frames()) {
      auto iframe = frame.first;
      auto aframe = frame.second;
      double frame_time = (iframe * frame_size / 320.e6);

      // filter selection
      if(triggerless_selection_indexes(frame, frame_size, subframe_size, min_subframe_size).size() != 0) {
        triggerless_hist->Fill(frame_time);
      }

      // trigger selection
      auto nsipm4 = aframe[207].hits[4].size();
      auto nsipm5 = aframe[207].hits[5].size();
      if (!((aframe[192].triggers.size() != 1) || (nsipm4 == 0 && nsipm5 == 0))) {
        trigger_hist->Fill(frame_time);
      }

      ++n_frames;
    }
    ++n_spills;
  }

  std::cout << "\nFound " << n_frames << " in " << n_spills << " spills";
  print_hist(triggerless_hist);
  print_hist(trigger_hist);
  
  // save histograms
  TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "UPDATE");
  trigger_hist->Write(nullptr, TObject::kOverwrite);
  triggerless_hist->Write(nullptr, TObject::kOverwrite);
  file->Close();
}
void selected_vs_time() {

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.1);
  // gStyle->SetOptStat("emr");
  // gStyle->SetStatW(0.4);
  // gStyle->SetStatH(0.15);

  // fetch histograms
  TFile *file = TFile::Open("./sipm4eic-testbeam2023-analysis/output/objects.root", "READ");
  
  TH1I *triggerless_hist = file->Get<TH1I>("selected-vs-number-triggerless-size256-spills25-I");
  triggerless_hist->SetTitle("");  
  triggerless_hist->Sumw2();

  TH1I *trigger_hist = file->Get<TH1I>("selected-vs-number-trigger-size256-spills25-I");
  trigger_hist->Sumw2();
 
  TH1D *triggerless_hist_d = new TH1D;
  TH1D *trigger_hist_d = new TH1D;
  triggerless_hist->Copy(*triggerless_hist_d);
  trigger_hist->Copy(*trigger_hist_d);

  triggerless_hist_d->SetYTitle("N_{frames}/N_{tot}");
  trigger_hist_d->SetMarkerColor(kMagenta);
  trigger_hist_d->SetLineColor(kMagenta);

  TLegend *leg = new TLegend(0.72, 0.75, 0.92, 0.83);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(triggerless_hist_d, "triggerless", "L");
  leg->AddEntry(trigger_hist_d, "trigger", "L");

  // draw
  TCanvas *canvas = new TCanvas("canvas", "frame selection overlap - log");
  canvas->SetLogy();
  triggerless_hist_d->DrawNormalized(); // do not draw error bars
  trigger_hist_d->DrawNormalized("SAME");
  leg->Draw();

  canvas->Print(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png","selected-vs-time-overlap-log"));
}

/*
    histogram showing the distribution of the radius of each cherenkov hit for trigger and triggerless mode
    Note: the two expected parameters are the filenames of the recodata root files generated using the recowriter.C
          macro in triggerless and trigger modes
    --- Thesis' figure [3.13]
*/
void selected_hits_radius(std::string triggerless_recodata, std::string trigger_recodata) {
  gStyle->SetOptStat(0);

  // fetch reconstructed hits Tree  
  unsigned short n;
  float x[65534];
  float y[65534];
  float t[65534];
  auto triggerless_fin = TFile::Open(triggerless_recodata.c_str());
  auto triggerless_tin = triggerless_fin->Get<TTree>("recodata");
  auto triggerless_nev = triggerless_tin->GetEntries();
  triggerless_tin->SetBranchAddress("n", &n);
  triggerless_tin->SetBranchAddress("x", &x);
  triggerless_tin->SetBranchAddress("y", &y);
  triggerless_tin->SetBranchAddress("t", &t);

  TH1F *triggerless_hist = new TH1F("triggerless-radius-dist-normalized", ";raggio (cm);N_{hit} / Tot_{hit}", 100, 0, 100);
  triggerless_hist->SetLineColor(kBlue);

  int n_hits_triggerless = 0;
  for (int iev = 0; iev < triggerless_nev; ++iev) {
    triggerless_tin->GetEntry(iev);
    for (int i = 0 ; i < n; ++i) {
      int ix = gRandom->Uniform(x[i] - 1.5, x[i] + 1.5);
      int iy = gRandom->Uniform(y[i] - 1.5, y[i] + 1.5);
      int iradius = std::sqrt(ix * ix + iy * iy);
      triggerless_hist->Fill(iradius);
    }
    n_hits_triggerless += n;
  }
  // "./sipm4eic-testbeam2023-analysis/output/trigger_recodata_I.root"
  auto trigger_fin = TFile::Open(trigger_recodata.c_str());
  auto trigger_tin = trigger_fin->Get<TTree>("recodata");
  auto trigger_nev = trigger_tin->GetEntries();
  trigger_tin->SetBranchAddress("n", &n);
  trigger_tin->SetBranchAddress("x", &x);
  trigger_tin->SetBranchAddress("y", &y);
  trigger_tin->SetBranchAddress("t", &t);

  TH1F *trigger_hist = new TH1F("trigger-radius-dist-normalized", "cherenkov hits radius (normalized);radius;N_{hits} / tot hits", 100, 0, 100);
  trigger_hist->SetLineColor(kRed);

  int n_hits_trigger = 0;
  for (int iev = 0; iev < trigger_nev; ++iev) {
    triggerless_tin->GetEntry(iev);
    for (int i = 0 ; i < n; ++i) {
      int ix = gRandom->Uniform(x[i] - 1.5, x[i] + 1.5);
      int iy = gRandom->Uniform(y[i] - 1.5, y[i] + 1.5);
      int iradius = std::sqrt(ix * ix + iy * iy);
      trigger_hist->Fill(iradius);
    }
    n_hits_trigger += n;
  }

  // normalize the histograms
  triggerless_hist->Sumw2();
  trigger_hist->Sumw2();
  triggerless_hist->Scale(1. / triggerless_hist->Integral());
  trigger_hist->Scale(1./ trigger_hist->Integral());

  TLegend *leg = new TLegend(0.72, 0.75, 0.92, 0.83);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(triggerless_hist, "triggerless", "L");
  leg->AddEntry(trigger_hist, "trigger", "L");

  // draw
  TCanvas *canvas = new TCanvas("canvas", nullptr);
  canvas->SetLogy();
  triggerless_hist->Draw();
  trigger_hist->Draw("SAME");
  leg->Draw();

  canvas->Print(Form("./sipm4eic-testbeam2023-analysis/output/images/%s.png", trigger_hist->GetTitle()));

}

/* end of FIGURES GENERATION */