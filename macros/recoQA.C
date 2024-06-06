#include "../lib/lightio.h"
#include "../lib/data.h"
#include "../lib/utility.h"

void recoQA(std::string input_file = "recodata.root", std::string output_file = "out.root", std::string save_dir = "./images/")
{
    //  Output
    auto hPersistance2D = new TH2F("hPersistance2D", ";X (mm);Y (mm); t (ns)", 396, -99, 99, 396, -99, 99);

    //  Link TTree to local data instance
    recodata reco_data;
    auto reco_tree = load_data(input_file, reco_data);

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
    }

    //  === Graphics
    gROOT->SetBatch();
    system(Form("mkdir -p %s/", save_dir.c_str()));

    gStyle->SetPalette(kInvertedDarkBodyRadiator);

    //  === === Trigger
    auto current_canvas = get_std_canvas();
    hPersistance2D->Draw();
    current_canvas->SaveAs(Form("%s/hPersistance2D.png", save_dir.c_str()));

    TFile *out = new TFile(output_file.c_str(), "RECREATE");
    hPersistance2D->Write();
    out->Close();

    gROOT->SetBatch(kFALSE);
}
