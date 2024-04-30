#include "../lib/lightio.h"
#include "../lib/data.h"

TCanvas *get_std_canvas()
{
  TCanvas *result = new TCanvas("", "", 1000, 1000);
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  return result;
}

void lightQA(std::string input_file = "lightdata.root", std::string output_file = "out.root")
{
  //  Define output objects
  //  === Trigger
  auto hTriggerHitsTimeInSpill = new TH2F("hTriggerHitsTimeInSpill", ";trigger channels time;entries", 1050, -0.1, 1, 99, 1, 100);
  auto hTriggerChannelInFrame = new TH2F("hTriggerChannelInFrame", ";trigger channels fired;entries", 200, 0, 200, 99, 1, 100);
  //  === Timing
  auto hTimingChannelInFrame = new TH2F("hTimingChannelInFrame", ";timing channels fired;entries", 64, 0, 64, 99, 1, 100);
  auto hTimingChannelMap = new TH2F("hTimingChannelMap", ";timing channels on chip 0 fired;timing channels on chip 1 fired", 33, -0.5, 32.5, 33, -0.5, 32.5);
  //  === Cherenkov
  auto hCherenkovChannelInFrame = new TH2F("hCherenkovChannelInFrame", ";cherenkov channels fired;entries", 200, 0, 200, 99, 1, 100);

  cout << "[INFO] Starting the light QA: reading " << input_file << endl;
  sipm4eic::lightio *io = new sipm4eic::lightio();
  io->read_from_tree(input_file);
  cout << "[INFO] Finished reading from tree in file " << input_file << endl;

  //  Loop on spills
  auto full_frames = 0;
  auto current_spill = 0;
  while (io->next_spill())
  {
    current_spill++;
    cout << "[INFO] Start spill: " << current_spill << endl;

    //  Loop on frames
    auto current_frame = 0;
    while (io->next_frame())
    {
      //  === Frame info
      current_frame++;
      auto frame_id = io->frame[current_frame];

      // === Trigger
      auto trigger0_vector = io->get_trigger0_vector();
      hTriggerChannelInFrame->Fill(trigger0_vector.size(), current_spill);
      auto trigger_time = -99999999.;
      if (trigger0_vector.size() != 0)
        trigger_time = trigger0_vector[0].coarse;
      hTriggerHitsTimeInSpill->Fill((trigger_time + 256 * (frame_id)) * sipm4eic::data::coarse_to_ns * 1.e-9, current_spill);

      // === Timing
      auto timing_vector = io->get_timing_vector();
      std::map<int, std::vector<float>> timing_channels_times;
      auto contributors_first_chip = 0;
      auto contributors_second_chip = 0;
      for (auto &timing : timing_vector)
      {
        auto timing_coarse = timing.coarse;
        auto timing_channel = timing.eoch();
        timing_channels_times[timing_channel].push_back(timing_coarse);
        if (timing_channel > 31)
          contributors_first_chip++;
        else
          contributors_second_chip++;
      }
      hTimingChannelInFrame->Fill(timing_channels_times.size(), current_spill);
      hTimingChannelMap->Fill(contributors_first_chip, contributors_second_chip);

      // === Cherenkov
      auto cherenkov_vector = io->get_cherenkov_vector();
      std::map<int, std::vector<float>> cherenkov_channels_times;
      auto average_cherenkov_time = 0.;
      for (auto &cherenkov : cherenkov_vector)
      {
        auto cherenkov_coarse = cherenkov.coarse;
        auto cherenkov_channel = cherenkov.eoch();
        cherenkov_channels_times[cherenkov_channel].push_back(cherenkov_coarse);
      }
      hCherenkovChannelInFrame->Fill(cherenkov_channels_times.size(), current_spill);
    }
  }

  gStyle->SetPalette(kInvertedDarkBodyRadiator);

  auto current_canvas = get_std_canvas();
  gPad->SetLogz();
  hTimingChannelMap->Draw();
  current_canvas->SaveAs("hTimingChannelMap.png");

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  auto hTriggerChannelTimeInSpillIntegrated = hTriggerHitsTimeInSpill->ProjectionX("hTriggerChannelTimeInSpillIntegrated", 0, 99);
  hTriggerChannelTimeInSpillIntegrated->Draw();
  current_canvas->SaveAs("hTriggerChannelTimeInSpillIntegrated.png");

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  auto hTriggerChannelInFrameIntegrated = hTriggerChannelInFrame->ProjectionX("hTriggerChannelInFrameIntegrated", 0, 99);
  hTriggerChannelInFrameIntegrated->Draw();
  current_canvas->SaveAs("hTriggerChannelInFrameIntegrated.png");

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  auto hTimingChannelInFrameIntegrated = hTimingChannelInFrame->ProjectionX("hTimingChannelInFrameIntegrated", 0, 99);
  hTimingChannelInFrameIntegrated->Draw();
  current_canvas->SaveAs("hTimingChannelInFrameIntegrated.png");

  current_canvas = get_std_canvas();
  gPad->SetLogy();
  auto hCherenkovChannelInFrameIntegrated = hCherenkovChannelInFrame->ProjectionX("hCherenkovChannelInFrameIntegrated", 0, 99);
  hCherenkovChannelInFrameIntegrated->Draw();
  current_canvas->SaveAs("hCherenkovChannelInFrameIntegrated.png");

  TFile *out = new TFile("out.root", "RECREATE");
  hTriggerHitsTimeInSpill->Write();
  hTriggerChannelInFrame->Write();
  hTimingChannelInFrame->Write();
  hCherenkovChannelInFrame->Write();
  hTimingChannelMap->Write();
  hTriggerChannelTimeInSpillIntegrated->Write();
  hTriggerChannelInFrameIntegrated->Write();
  hTimingChannelInFrameIntegrated->Write();
  hCherenkovChannelInFrameIntegrated->Write();
  out->Close();
}
