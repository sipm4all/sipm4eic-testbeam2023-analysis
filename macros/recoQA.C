#include "../lib/lightio.h"
#include "../lib/data.h"

struct recodata
{
    unsigned short n = 0;
    float x[1024];
    float y[1024];
    float t[1024];
};

void load_data(TTree *reco_data_tree, recodata &target_data_struct)
{
    reco_data_tree->SetBranchAddress("n", &target_data_struct.n);
    reco_data_tree->SetBranchAddress("x", &target_data_struct.x);
    reco_data_tree->SetBranchAddress("y", &target_data_struct.y);
    reco_data_tree->SetBranchAddress("t", &target_data_struct.t);
}

TTree *load_data(std::string filename, recodata &target_data_struct)
{
    auto input_file = TFile::Open(filename.c_str());
    auto input_tree = (TTree *)input_file->Get("recodata");
    load_data(input_tree, target_data_struct);
    return input_tree;
}

void recoQA(std::string input_file = "recodata.root", std::string output_file = "out.root", std::string save_dir = "./images/")
{
    
}
