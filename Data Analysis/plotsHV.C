#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>
#include <ROOT/RDataFrame.hxx>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/TThreadExecutor.hxx>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraph.h>

using namespace std;

using StatsMap = tuple<map<int, pair<double, double>>, map<int, pair<double, double>>, map<int, pair<double, double>>>;

StatsMap stats(const string& filePath)
{
    ROOT::EnableImplicitMT(8);

    TFile *file = TFile::Open(filePath.c_str());
    if (!file || file->IsZombie()) {
        cerr << "Error opening file or file not found: " << filePath << endl;
        return {};
    }

    TTree *tree = (TTree*)file->Get("raw");
    if (!tree) {
        cerr << "Tree 'raw' not found in file: " << filePath << endl;
        file->Close();
        return {};
    }

    ROOT::RDataFrame df(*tree);

    vector<int> ids = {8, 9, 10, 11, 12, 13};
    map<int, string> idToName = {{8, "2"}, {9, "3"}, {10, "4"}, {11, "5"}, {12, "6"}, {13, "7"}};
    
    map<int, TH1I*> histogramsHits;
    map<int, TH1F*> histogramsQ;
    map<int, TH1F*> histogramsTQ;
    
    for (int id : ids) {
        string name = idToName[id];
        
        histogramsHits[id] = new TH1I(Form("hCounts_%s", name.c_str()),
                                      Form("Hits of Detector in plane %s per Event; Number of hits; Number of Events", name.c_str()),
                                      150, 0, 150);

        histogramsQ[id] = new TH1F(Form("hMaxQ_%s", name.c_str()),
                                   Form("Charge of Detector %s;Charge (ADC);Number of Events", name.c_str()),
                                   200, 0, 2000);

        histogramsTQ[id] = new TH1F(Form("hSumMaxCharge_%s", name.c_str()), 
                                    Form("Sum of Maximum Charges for Detector %s;Charge;Number of Events", name.c_str()), 
                                    200, 0, 4000);
    }

    auto hasDetector13 = [](const vector<unsigned int>& ids) {
        return find(ids.begin(), ids.end(), 13) != ids.end();
    };

    auto filtered = df.Filter(hasDetector13, {"apv_id"});

    filtered.Foreach([&](const vector<unsigned int>& apv_id) {
        map<int, int> hitsCount;
        for (int id : ids) {
            hitsCount[id] = count(apv_id.begin(), apv_id.end(), id);
            if (hitsCount[id] != 0 && hitsCount[id] <= 10) {
                histogramsHits[id]->Fill(hitsCount[id]);
            }  
        }
    }, {"apv_id"});

    filtered.Foreach([&](const vector<unsigned int>& apv_id, const vector<vector<short>>& apv_q) {
        for (size_t i = 0; i < apv_id.size(); i++) {
            int current_id = apv_id[i];
            if (histogramsQ.count(current_id) > 0 && histogramsQ[current_id] != nullptr) {
                short maxQ = *max_element(apv_q[i].begin(), apv_q[i].end());
                histogramsQ[current_id]->Fill(maxQ);
            }
        }   
    }, {"apv_id", "apv_q"});

    filtered.Foreach([&](const vector<unsigned int>& apv_id, const vector<vector<short>>& apv_q) {
        map<int, double> sumMaxCharges;

        for (size_t i = 0; i < apv_id.size(); ++i) {
            int current_id = apv_id[i];
            if (histogramsTQ.find(current_id) != histogramsTQ.end()) {
                const auto& charges = apv_q[i];
                if (!charges.empty()) {
                    short maxCharge = *max_element(charges.begin(), charges.end());
                    sumMaxCharges[current_id] += maxCharge;
                }
            }
        }

        for (const auto& sum : sumMaxCharges) {
            histogramsTQ[sum.first]->Fill(sum.second);
        }
    }, {"apv_id", "apv_q"});

    map<int, pair<double, double>> statsHits; // mean, mean_error
    map<int, pair<double, double>> statsQ;    // mean, mean_error
    map<int, pair<double, double>> statsTQ;   // mean, mean_error

    for (int id : ids) {
        double mean = histogramsHits[id]->GetMean();
        int N = histogramsHits[id]->GetEntries();
        double stddev = histogramsHits[id]->GetStdDev();
        double mean_error = stddev / sqrt(N);
        statsHits[id] = make_pair(mean, mean_error);
    }

    for (int id : ids) {
        double mean = histogramsQ[id]->GetMean();
        int N = histogramsQ[id]->GetEntries();
        double stddev = histogramsQ[id]->GetStdDev();
        double mean_error = stddev / sqrt(N);
        statsQ[id] = make_pair(mean, mean_error);
    }

    for (int id : ids) {
        double mean = histogramsTQ[id]->GetMean();
        int N = histogramsTQ[id]->GetEntries();
        double stddev = histogramsTQ[id]->GetStdDev();
        double mean_error = stddev / sqrt(N);
        statsTQ[id] = make_pair(mean, mean_error);
    }

    file->Close();
//     for (auto& hist : histogramsHits) delete hist.second;
//     for (auto& hist : histogramsQ) delete hist.second;
//     for (auto& hist : histogramsTQ) delete hist.second;

    return {statsHits, statsQ, statsTQ};
}

void plotsHV()
{
    vector<string> files = {
        "Data/run6578.root",
        "Data/run6583.root",
        "Data/run6586.root",
        "Data/run6592.root",
        "Data/run6596.root",
        "Data/run6599.root",
        "Data/run6602.root",
        "Data/run6610.root",
        "Data/run6612.root",
        "Data/run6614.root",
        "Data/run6616.root"
    };

    vector<double> hv_levels = {500, 520, 520, 540, 560, 580, 600, 480, 460, 440, 420};
    map<int, vector<double>> hitsMeans;
    map<int, vector<double>> hitsStdDevs;
    map<int, vector<double>> tqMeans;
    map<int, vector<double>> tqMeanErrors;

    vector<int> ids = {8, 9, 10, 11, 12, 13};
    int numFiles = files.size();

    for (int id : ids) {
        hitsMeans[id] = vector<double>(numFiles, 0);
        hitsStdDevs[id] = vector<double>(numFiles, 0);
        tqMeans[id] = vector<double>(numFiles, 0);
        tqMeanErrors[id] = vector<double>(numFiles, 0);
    }

    for (size_t i = 0; i < files.size(); ++i) {
        auto statsMap = stats(files[i]);
        auto& statsHits = get<0>(statsMap);
        auto& statsTQ = get<2>(statsMap);

        for (int id : ids) {
            hitsMeans[id][i] = statsHits[id].first; // Store the mean of hits
            hitsStdDevs[id][i] = statsHits[id].second; // Store the standard deviation of hits
            tqMeans[id][i] = statsTQ[id].first;
            tqMeanErrors[id][i] = statsTQ[id].second;
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Cluster size vs High Voltage", 850, 600);
    c1->SetGrid();

    TLegend *legendHits = new TLegend(0.1, 0.7, 0.25, 0.9);
//     legendHits->SetHeader("Detectors", "C");

    vector<int> colors = {kRed, kBlue, kGreen+2, kMagenta+2, kCyan+2, kOrange-3};
    vector<int> markers = {20, 21, 22, 23, 24, 25};
    map<int, string> idToName = {{8, "2"}, {9, "3"}, {10, "4"}, {11, "5"}, {12, "6"}, {13, "7"}};  

    TMultiGraph *mgHits = new TMultiGraph();

    for (size_t i = 0; i < ids.size(); ++i) {
        int id = ids[i];
        string name = idToName[id];
        TGraphErrors *graph = new TGraphErrors(numFiles, hv_levels.data(), hitsMeans[id].data(), nullptr, hitsStdDevs[id].data());
        graph->SetTitle(Form("Detector %s", name.c_str()));
        graph->SetMarkerColor(colors[i]);
        graph->SetMarkerStyle(markers[i]);
        graph->SetLineColor(colors[i]);
        graph->SetLineWidth(1);
        mgHits->Add(graph, "P"); // Add points only with error bars
        legendHits->AddEntry(graph, Form("Detector %s", name.c_str()), "p");
    }

    mgHits->SetTitle("Cluster size vs High Voltage;HV (V);Cluster size");
    mgHits->Draw("A");
    legendHits->Draw();
    c1->SaveAs("HV_vs_Clustersize.png");



    TCanvas *c2 = new TCanvas("c2", "Total Charge vs High Voltage", 850, 600);
    c2->SetGrid();

    TLegend *legendTQ = new TLegend(0.1, 0.7, 0.25, 0.9);
//     legendTQ->SetHeader("Detectors", "C"); 

    TMultiGraph *mgTQ = new TMultiGraph();

    for (size_t i = 0; i < ids.size(); ++i) {
        int id = ids[i];
        string name = idToName[id];
        TGraphErrors *graph = new TGraphErrors(numFiles, hv_levels.data(), tqMeans[id].data(), nullptr, tqMeanErrors[id].data());
        graph->SetTitle(Form("Detector %s", name.c_str()));
        graph->SetMarkerColor(colors[i]);
        graph->SetMarkerStyle(markers[i]);
        graph->SetLineColor(colors[i]);
        graph->SetLineWidth(1);
        mgTQ->Add(graph, "P"); // Add points only with error bars
        legendTQ->AddEntry(graph, Form("Detector %s", name.c_str()), "p");
    }

    mgTQ->SetTitle("Total Charge vs High Voltage;HV (V);Q (ADC)");
    mgTQ->Draw("A");
    legendTQ->Draw();
    c2->SaveAs("HV_vs_Totalcharge.png");



    delete c1;
    delete legendHits;
    delete c2;
    delete legendTQ;
}

// int main() {
//     plotsHV();
//     return 0;
// }


