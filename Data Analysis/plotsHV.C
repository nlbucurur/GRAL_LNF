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
#include <tuple>

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
            if (hitsCount[id] != 0) {
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

    map<int, pair<double, double>> statsHits;
    map<int, pair<double, double>> statsQ;
    map<int, pair<double, double>> statsTQ;

    for (int id : ids) {
        double mean = histogramsHits[id]->GetMean();
        double stddev = histogramsHits[id]->GetStdDev();
        statsHits[id] = make_pair(mean, stddev);
    }

    for (int id : ids) {
        double mean = histogramsQ[id]->GetMean();
        double stddev = histogramsQ[id]->GetStdDev();
        statsQ[id] = make_pair(mean, stddev);
    }

    for (int id : ids) {
        double mean = histogramsTQ[id]->GetMean();
        double stddev = histogramsTQ[id]->GetStdDev();
        statsTQ[id] = make_pair(mean, stddev);
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

    for (const auto& file : files) {
        auto statsMap = stats(file);
        auto& statsHits = get<0>(statsMap);
        auto& statsQ = get<1>(statsMap);
        auto& statsTQ = get<2>(statsMap);

        cout << "File: " << file << endl;

        cout << "Hits Statistics:" << endl;
        for (const auto& stat : statsHits) {
            cout << "ID: " << stat.first << " Mean: " << stat.second.first << " StdDev: " << stat.second.second << endl;
        }

        cout << "Charge Statistics:" << endl;
        for (const auto& stat : statsQ) {
            cout << "ID: " << stat.first << " Mean: " << stat.second.first << " StdDev: " << stat.second.second << endl;
        }

        cout << "Total Charge Statistics:" << endl;
        for (const auto& stat : statsTQ) {
            cout << "ID: " << stat.first << " Mean: " << stat.second.first << " StdDev: " << stat.second.second << endl;
        }

        cout << endl;
    }
}




