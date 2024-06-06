#include "../lib/lightio.h"
#include "../lib/data.h"
#include "../lib/mapping.h"
#include "../lib/utility.h"

void recoQA(std::string input_file = "recodata.root", std::string output_file = "out.root", std::string save_dir = "./images/")
{
    //  Output
    auto hPersistance2D = new TH2F("hPersistance2D", ";X (mm);Y (mm); t (ns)", 396, -99, 99, 396, -99, 99);
    auto hMap_fullsetup_SiPM = new TH2F("hMap_fullsetup_SiPM", ";X (mm);Y (mm)", 4000, -100, 100, 4000, -100, 100);
    auto hMap_availsetup_SiPM = new TH2F("hMap_availsetup_SiPM", ";X (mm);Y (mm)", 4000, -100, 100, 4000, -100, 100);

    //  Link TTree to local data instance
    recodata reco_data;
    auto reco_tree = load_data(input_file, reco_data);

    //  Store available SiPMs
    std::vector<std::array<float, 2>> list_of_available_SiPMs;

    //  First loop on events
    cout << "[INFO] Start of preliminary loop for X_{0}, Y_{0} and R_{0}" << endl;
    for (int iEv = 0; iEv < reco_tree->GetEntries(); iEv++)
    {
        //  Recover recodata entry form tree
        reco_tree->GetEntry(iEv);

        if (iEv % 1000 == 0)
            cout << "[INFO] event: " << iEv << endl;

        //  Persistance plot
        fill_persistance(hPersistance2D, reco_data);

        //  Loop on hits
        for (int iPnt = 0; iPnt < reco_data.n; iPnt++)
        {
            if (!count(list_of_available_SiPMs.begin(), list_of_available_SiPMs.end(), std::array<float, 2>({reco_data.x[iPnt], reco_data.y[iPnt]})))
                //  Add sensor sensitive area
                list_of_available_SiPMs.push_back({reco_data.x[iPnt], reco_data.y[iPnt]});
        }
    }
    for (auto iPDU = 0; iPDU < 8; iPDU++)
        for (auto iCol = 0; iCol < 16; iCol++)
            for (auto iRow = 0; iRow < 16; iRow++)
                fill_with_SiPM_coverage(hMap_fullsetup_SiPM, sipm4eic::get_position({iPDU, iCol, iRow}));
    for (auto current_position : list_of_available_SiPMs)
        fill_with_SiPM_coverage(hMap_availsetup_SiPM, current_position);

    //  === Graphics
    gROOT->SetBatch();
    system(Form("mkdir -p %s/", save_dir.c_str()));

    gStyle->SetPalette(kInvertedDarkBodyRadiator);

    //  === === Persistance
    auto current_canvas = get_std_canvas();
    hPersistance2D->Draw();
    current_canvas->SaveAs(Form("%s/hPersistance2D.png", save_dir.c_str()));

    //  === === SiPM and ring coverage
    //  === === === All SiPMs
    current_canvas = get_std_canvas();
    hMap_fullsetup_SiPM->Draw();
    current_canvas->SaveAs(Form("%s/hMap_fullsetup_SiPM.png", save_dir.c_str()));
    //  === === === Available SiPMs
    current_canvas = get_std_canvas();
    hMap_availsetup_SiPM->Draw();
    current_canvas->SaveAs(Form("%s/hMap_availsetup_SiPM.png", save_dir.c_str()));

    TFile *out = new TFile(output_file.c_str(), "RECREATE");
    hPersistance2D->Write();
    hMap_fullsetup_SiPM->Write();
    hMap_availsetup_SiPM->Write();
    out->Close();

    gROOT->SetBatch(kFALSE);
}
