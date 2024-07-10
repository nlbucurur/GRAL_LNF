#include <iostream>
#include <vector>
#include <algorithm>           // Include for std::max_element
#include <cmath>               // Include for sqrt function
#include <map>
#include <ROOT/RDataFrame.hxx> // For modern data frame tool
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/TThreadExecutor.hxx> // Include for setting number of threads

using namespace std;

void histograms()
{

    // Enable multi-threading and set the number of threads
    ROOT::EnableImplicitMT(8); // Replace 4 with the desired number of threads

    // Open the ROOT file containing the tree.
    TFile *file = TFile::Open("Data/run6599.root");
    if (!file || file->IsZombie()) { // Check if the file is opened successfully.
        cerr << "Error opening file or file not found." << endl;
        return;
    }

    // Retrieve the tree named "raw" from the file.
    TTree *tree = (TTree*)file->Get("raw");
    if (!tree) { // Check if the tree exists.
        cerr << "Tree 'raw' not found in file." << endl;
        file->Close(); // Close the file if the tree is not found.
        return; // Exit the function.
    }

    // Create a data frame from the TTree for easier processing.
    ROOT::RDataFrame df(*tree);

    vector<int> ids = {8, 9, 10, 11, 12, 13};
    map<int, string> idToName = {{8, "2"}, {9, "3"}, {10, "4"}, {11, "5"}, {12, "6"}, {13, "7"}};
    
    // Map to store a histogram for each detector ID
    map<int, TH1I*> histogramsHits;
    map<int, TH1F*> histogramsQ;
    map<int, TH1F*> histogramsTQ;
    
    for (int id : ids) {
        string name = idToName[id];
        
        histogramsHits[id] = new TH1I(Form("hCounts_%s", name.c_str()),
                                      Form("Hits of Detector in plane %s per Event; Number of hits; Number of Events", name.c_str()),
                                      150, 0, 150);

        histogramsQ[id] = new TH1F(Form("hMaxQ_%s", name.c_str()),
                                   Form("Maximum Charge of Detector %s;Charge (ADC);Number of Events", name.c_str()),
                                   200, 0, 2000);

        histogramsTQ[id] = new TH1F(Form("hSumMaxCharge_%s", name.c_str()), 
                                    Form("Sum of Maximum Charges for Detector %s;Charge;Number of Events", name.c_str()), 
                                    200, 0, 4000);
    }

    // Lambda function to check if detector with id 13 is activated in the event.
    auto hasDetector13 = [](const vector<unsigned int>& ids) {
        return find(ids.begin(), ids.end(), 13) != ids.end();
    };

    // Filter the data frame to include only those events where detector 13 is activated.
    auto filtered = df.Filter(hasDetector13, {"apv_id"});

    /**********/
    /** Hits **/
    /**********/
    filtered.Foreach([&](const vector<unsigned int>& apv_id) {
        map<int, int> hitsCount;
        for (int id : ids) {
            hitsCount[id] = count(apv_id.begin(), apv_id.end(), id);
            // histogramsHits[id]->Fill(hitsCount[id]);
            if (hitsCount[id] != 0)
            {
                histogramsHits[id]->Fill(hitsCount[id]);
            }  
        }
    }, {"apv_id"});

    /********************/
    /** Charge Maximum **/
    /********************/
    filtered.Foreach([&](const vector<unsigned int>& apv_id, const vector<vector<short>>& apv_q) {
        for (size_t i = 0; i < apv_id.size(); i++) {
            int current_id = apv_id[i];
            if (histogramsQ.count(current_id) > 0 && histogramsQ[current_id] != nullptr) {
                short maxQ = *max_element(apv_q[i].begin(), apv_q[i].end());
                histogramsQ[current_id]->Fill(maxQ);
            }
        }   
    }, {"apv_id", "apv_q"});

    /************************/
    /** Maximum Charge Sum **/
    /************************/
    filtered.Foreach([&](const vector<unsigned int>& apv_id, const vector<vector<short>>& apv_q) {
        map<int, double> sumMaxCharges; // Map to store the sum of max charges for each ID in this event

        for (size_t i = 0; i < apv_id.size(); ++i) {
            int current_id = apv_id[i];
            if (histogramsTQ.find(current_id) != histogramsTQ.end()) { // Check if we have a histogram for this ID
                const auto& charges = apv_q[i];
                if (!charges.empty()) {
                    short maxCharge = *max_element(charges.begin(), charges.end());
                    sumMaxCharges[current_id] += maxCharge; // Add max charge to sum for this detector
                }
            }
        }

        // Fill histograms with the sum of maximum charges for each relevant detector
        for (const auto& sum : sumMaxCharges) {
            histogramsTQ[sum.first]->Fill(sum.second);
        }
    }, {"apv_id", "apv_q"});

    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "Histograms of Hits for Different Detectors", 1200, 1200);
    canvas->Divide(2, 3); // Divide the canvas into a 2x3 grid

    // Define colors for each histogram
    vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange};

    // Draw each histogram on the canvas
    for (size_t i = 0; i < ids.size(); ++i) {
        int id = ids[i];
        string name = idToName[id];
        
        canvas->cd(i + 1); // Move to the next pad
        histogramsHits[id]->SetLineColor(colors[i]); // Set the color of the histogram

        // Set the font size for the axes titles and labels
        float sizelabel = 0.035;
        histogramsHits[id]->GetXaxis()->SetTitleSize(sizelabel);
        histogramsHits[id]->GetXaxis()->SetLabelSize(sizelabel);
        histogramsHits[id]->GetYaxis()->SetTitleSize(sizelabel);
        histogramsHits[id]->GetYaxis()->SetLabelSize(sizelabel);

        histogramsHits[id]->Draw(); // Draw the histogram

        histogramsHits[id]->GetXaxis()->SetRangeUser(1, 10); // Set x-axis range

        // Move the statistics box
        gPad->Update(); // Update the pad to ensure the stats box is created
        TPaveStats *stats = (TPaveStats*)histogramsHits[id]->FindObject("stats");
        if (stats) {
            stats->SetX1NDC(0.7); // New X1 position
            stats->SetX2NDC(0.9); // New X2 position
            stats->SetY1NDC(0.7); // New Y1 position
            stats->SetY2NDC(0.9); // New Y2 position
            stats->Draw(); // Draw the stats box on the pad
        }
        
    }

    canvas->SaveAs("Figures/histograms_all_detectors.png"); // Save the canvas as an image file


    // Create a canvas for the overlaid histograms in log scale
    TCanvas *logCanvas = new TCanvas("logCanvas", "Log-Scale Histograms of Hits for Different Detectors", 1200, 800);
    logCanvas->SetLogy(); // Set the y-axis to log scale

    // Set the font size for the entire canvas
    gStyle->SetTitleSize(0.05, "xyz"); // Set title size for x, y, and z axes
    gStyle->SetLabelSize(0.05, "xyz"); // Set label size for x, y, and z axes
    gStyle->SetOptStat(0); // Hide statistics boxes for the overlaid histograms

    // Draw the overlaid histograms
    for (size_t i = 0; i < ids.size(); ++i) {
        int id = ids[i];
        string name = idToName[id];

        histogramsHits[id]->SetLineColor(colors[i]); // Set the color of the histogram
        histogramsHits[id]->SetLineWidth(1 + i); // Set varying line width for distinction
        histogramsHits[id]->SetStats(0); // Hide the individual stats boxes for the overlaid histograms
        histogramsHits[id]->GetXaxis()->SetRangeUser(1, 121); // Set x-axis range for log-scale plot

        // Set the font size for the axes titles and labels
        float sizelabel = 0.035;
        histogramsHits[id]->GetXaxis()->SetTitleSize(sizelabel);
        histogramsHits[id]->GetXaxis()->SetLabelSize(sizelabel);
        histogramsHits[id]->GetYaxis()->SetTitleSize(sizelabel);
        histogramsHits[id]->GetYaxis()->SetLabelSize(sizelabel);

        if (i == 0) {
            histogramsHits[id]->Draw(); // Draw the first histogram
        } else {
            histogramsHits[id]->Draw("same"); // Overlay the subsequent histograms
        }
    }

    // Create a legend for the overlaid histograms
    TLegend *legend = new TLegend(0.75, 0.7, 0.9, 0.9);
    for (size_t i = 0; i < ids.size(); ++i) {
        legend->AddEntry(histogramsHits[ids[i]], Form("Detector %s", idToName[ids[i]].c_str()), "l");
        float sizelabel = 0.03;
        legend->SetTextSize(sizelabel);
    }
    legend->Draw();

    logCanvas->SaveAs("Figures/logscale_histograms_all_detectors.png"); // Save the log-scale canvas as an image file



    // Clean up
    file->Close();
    delete file;
    delete canvas;
    delete logCanvas;
    delete legend;
    // for (auto& hist : histogramsHits) delete hist.second;
    // for (auto& hist : histogramsQ) delete hist.second;
    // for (auto& hist : histogramsTQ) delete hist.second;
}
