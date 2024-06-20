#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMultiGraph.h>

using namespace std;

// Function to compute efficiencies for a given run
unordered_map<int, pair<double, double>> computeEfficiency(const string& filename) {
    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << filename << endl;
        return {};
    }

    TTree *tree = (TTree*)file->Get("raw");
    if (!tree) {
        cerr << "Tree 'raw' not found in file: " << filename << endl;
        file->Close();
        return {};
    }

    ROOT::RDataFrame df(*tree);
    auto hasDetector13 = [](const vector<unsigned int>& ids) {
        return find(ids.begin(), ids.end(), 13) != ids.end();
    };
    auto filtered = df.Filter(hasDetector13, {"apv_id"});

    unordered_map<int, int> counts;
    filtered.Foreach([&](const vector<unsigned int>& ids) {
        unordered_set<int> unique_ids(ids.begin(), ids.end());
        for (int id : {8,9,10,11,12,13}) {
            if (unique_ids.count(id)) counts[id]++;
        }
    }, {"apv_id"});

    double N = filtered.Count().GetValue();
    unordered_map<int, pair<double, double>> efficiencies;

    for (int id : {8,9,10,11,12,13}) {
        double n = counts[id];
        double efficiency = N > 0 ? (n / N) * 100 : 0;
        double error = N > 0 ? (1 / sqrt(N)) * sqrt(n / N * (1 - n / N)) * 100 : 0;
        efficiencies[id] = make_pair(efficiency, error);
    }

    file->Close();
    return efficiencies;
}

// // Main function to process multiple runs
// int efficiency() {
//     vector<string> files = {"Data/run6578.root", "Data/run6583.root", "Data/run6586.root", "Data/run6592.root", "Data/run6596.root", "Data/run6599.root"};
//     vector<int> hvs = {500, 520, 520, 540, 560, 580}; // HV for detectors 8-12, detector 13 always at 600
    
//     // Store efficiency data
//     unordered_map<int, vector<pair<double, double>>> data; // ID -> List of (efficiency, error)
//     vector<double> hv_levels = {500, 520, 520, 540, 560, 580};  // HV levels without 600 since it's constant for ID 13

//     for (size_t i = 0; i < files.size(); ++i) {
//         auto eff = computeEfficiency(files[i]);
//         for (int id : {8,9,10,11,12}) {
//             data[id].push_back(eff[id]);
//         }
//         // Handle detector 13 separately if needed
//     }

//     // Create graphs for each detector
//     for (int id : {8,9,10,11,12}) {
//         vector<double> effs, errs;.
//         for (auto& e : data[id]) {
//             effs.push_back(e.first);
//             errs.push_back(e.second);
//         }

//         TGraphErrors *graph = new TGraphErrors(hv_levels.size(), hv_levels.data(), effs.data(), nullptr, errs.data());
//         graph->SetTitle(Form("Efficiency vs HV for Detector %d;HV;Efficiency (%)", id));
//         TCanvas *c = new TCanvas(Form("c%d", id), Form("Canvas for Detector %d", id), 600, 400);
//         graph->Draw("AP");
//         c->SaveAs(Form("Efficiency_HV_Detector%d.png", id));
//     }

//     return 0;
// }


// Main function to process multiple runs
int efficiency() {
    vector<string> files = {"Data/run6578.root", "Data/run6583.root", "Data/run6586.root", "Data/run6592.root", "Data/run6596.root", "Data/run6599.root"};
    vector<int> hvs = {500, 520, 520, 540, 560, 580}; // HV for detectors 8-12, detector 13 always at 600
    
    map<int, string> idToName = {{8, "2"}, {9, "3"}, {10, "4"}, {11, "5"}, {12, "6"}, {13, "7"}};
    unordered_map<int, vector<pair<double, double>>> data; // ID -> List of (efficiency, error)
    vector<double> hv_levels = {500, 520, 520, 540, 560, 580};  // HV levels without 600 since it's constant for ID 13

    for (size_t i = 0; i < files.size(); ++i) {
        auto eff = computeEfficiency(files[i]);
        for (int id : {8,9,10,11,12}) {
            data[id].push_back(eff[id]);
        }
    }

    TMultiGraph *mg = new TMultiGraph();
    TLegend *legend = new TLegend(0.7, 0.5, 0.9, 0.7);

    int colors[] = {kRed, kBlue, kGreen+1, kMagenta, kOrange};

    for (int id = 8, ci = 0; id <= 12; id++, ci++) {
        vector<double> effs, errs;
        for (auto& e : data[id]) {
            effs.push_back(e.first);
            errs.push_back(e.second);
        }

        TGraphErrors *graph = new TGraphErrors(hv_levels.size(), hv_levels.data(), effs.data(), nullptr, errs.data());
        graph->SetLineColor(colors[ci]);
        graph->SetMarkerColor(colors[ci]);
        graph->SetMarkerStyle(21 + ci);

        mg->Add(graph);
        legend->AddEntry(graph, Form("Detector %s", idToName[id].c_str()), "lp");
    }

    TCanvas *c = new TCanvas("c", "Efficiencies vs HV", 800, 600);
    mg->SetTitle("Efficiency vs HV for Detectors;HV;Efficiency (%)");
    c->SetGridx();
    c->SetGridy();
    mg->Draw("AP");
    legend->Draw();

    c->SaveAs("Efficiency_HV_AllDetectors.png");

    return 0;
}