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
#include <TStyle.h>
#include <TLegend.h>

using namespace std;

void histograms()
{

    // Enable multi-threading and set the number of threads
    ROOT::EnableImplicitMT(8); // Replace 4 with the desired number of threads

    // Open the ROOT file containing the tree.
    TFile *file = TFile::Open("Data/run6578.root");
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
                                   Form("Charge of Detector %s;Charge (ADC);Number of Events", name.c_str()),
                                   200, 0, 2000);

        histogramsTQ[id] = new TH1F(Form("hSumMaxCharge_%s", name.c_str()), 
                                    Form("Total Charge for Detector %s;Charge;Number of Events", name.c_str()), 
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
            if (hitsCount[id] != 0 && hitsCount[id] <= 10){
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


    /****************/
    /** Total Hits **/
    /****************/

    // Create a canvas
    TCanvas *canvasHits = new TCanvas("canvasHits", "Histograms of Hits for Different Detectors", 1200, 1200);
    canvasHits->Divide(2, 3); // Divide the canvas into a 2x3 grid
    
    // Grid for both axes
    canvasHits->SetGridx();
    canvasHits->SetGridy();

    // Define colors for each histogram
    vector<int> colors = {kRed, kBlue, kGreen+2, kMagenta+2, kCyan+2, kOrange-3};

    // Draw each histogram on the canvas
    for (size_t i = 0; i < ids.size(); ++i) {
        int id = ids[i];
        string name = idToName[id];

        // double mean = histogramsHits[id]->GetMean();
        // cout << "ID: " << id << "\t" << mean << endl;
        
        canvasHits->cd(i + 1); // Move to the next pad
        histogramsHits[id]->SetLineColor(colors[i]); // Set the color of the histogram
        histogramsHits[id]->SetLineWidth(2); // Set varying line width for distinction
        // Set the font size for the axes titles and labels
        float sizelabel = 0.035;
        histogramsHits[id]->GetXaxis()->SetTitleSize(sizelabel);
        histogramsHits[id]->GetXaxis()->SetLabelSize(sizelabel);
        histogramsHits[id]->GetYaxis()->SetTitleSize(sizelabel);
        histogramsHits[id]->GetYaxis()->SetLabelSize(sizelabel);

        histogramsHits[id]->Draw(); // Draw the histogram
        gPad->SetGrid(); // Add grid to the current pad
        histogramsHits[id]->GetXaxis()->SetRangeUser(1, 10); // Set x-axis range
        

        // Move the statistics box
        gPad->Update(); // Update the pad to ensure the stats box is created
        TPaveStats *statsHits = (TPaveStats*)histogramsHits[id]->FindObject("stats");
        if (statsHits) {
            statsHits->SetX1NDC(0.7); // New X1 position
            statsHits->SetX2NDC(0.9); // New X2 position
            statsHits->SetY1NDC(0.7); // New Y1 position
            statsHits->SetY2NDC(0.9); // New Y2 position
            statsHits->Draw(); // Draw the stats box on the pad
        }

        histogramsHits[id]->ResetStats();
        
    }

    canvasHits->SetGrid(); // Enable grid
    canvasHits->SaveAs("Figures/histograms_all_detectors.png"); // Save the canvas as an image file


    /********************/
    /** Total Hits Log **/
    /********************/

    // Create a canvas for the overlaid histograms in log scale
    TCanvas *logcanvasHits = new TCanvas("logcanvasHits", "Log-Scale Histograms of Hits for Different Detectors", 1200, 800);
    logcanvasHits->SetTitle("Hits on detectors per event"); // Set title for log-scale canvas
    logcanvasHits->SetLogy(); // Set the y-axis to log scale
    logcanvasHits->SetGrid(); // Enable grid

    // Set the font size for the entire canvas
    gStyle->SetTitleSize(0.05, "xyz"); // Set title size for x, y, and z axes
    gStyle->SetLabelSize(0.05, "xyz"); // Set label size for x, y, and z axes
    gStyle->SetOptStat(0); // Hide statistics boxes for the overlaid histograms
    gStyle->SetOptTitle(0);

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
    TLegend *legendHits = new TLegend(0.78, 0.7, 0.9, 0.9);
    for (size_t i = 0; i < ids.size(); ++i) {
        legendHits->AddEntry(histogramsHits[ids[i]], Form("Detector %s", idToName[ids[i]].c_str()), "l");
        float sizelabel = 0.03;
        legendHits->SetTextSize(sizelabel);
    }
    legendHits->Draw();

    // logcanvasHits->SaveAs("Figures/logscale_histograms_all_detectors.png"); // Save the log-scale canvas as an image file

    gStyle->SetOptStat(1111); // Hide statistics boxes for the overlaid histograms
    gStyle->SetOptTitle(1);


    /********************/
    /** Maximum Charge **/
    /********************/

    // Create a canvas for individual total charge histograms
    TCanvas *canvasQ = new TCanvas("canvasQ", "Histograms of Maximum Charge for Different Detectors", 1200, 1200);
    canvasQ->Divide(2, 3); // Divide the canvas into a 2x3 grid

    // Grid for both axes
    canvasQ->SetGridx();
    canvasQ->SetGridy();

    // Draw each total charge histogram on the canvas
    for (size_t i = 0; i < ids.size(); ++i) {
        int id = ids[i];
        string name = idToName[id];
        
        canvasQ->cd(i + 1); // Move to the next pad
        histogramsQ[id]->SetLineColor(colors[i]); // Set the color of the histogram

        // Set the font size for the axes titles and labels
        float sizelabel = 0.035;
        histogramsQ[id]->GetXaxis()->SetTitleSize(sizelabel);
        histogramsQ[id]->GetXaxis()->SetLabelSize(sizelabel);
        histogramsQ[id]->GetYaxis()->SetTitleSize(sizelabel);
        histogramsQ[id]->GetYaxis()->SetLabelSize(sizelabel);

        histogramsQ[id]->Draw(); // Draw the histogram
        gPad->SetGrid(); // Add grid to the current pad
        histogramsQ[id]->GetXaxis()->SetRangeUser(0, 2000); // Set x-axis range

        // Move the statistics box
        gPad->Update(); // Update the pad to ensure the stats box is created
        TPaveStats *statsQ = (TPaveStats*)histogramsQ[id]->FindObject("stats");
        if (statsQ) {
            statsQ->SetX1NDC(0.7); // New X1 position
            statsQ->SetX2NDC(0.9); // New X2 position
            statsQ->SetY1NDC(0.7); // New Y1 position
            statsQ->SetY2NDC(0.9); // New Y2 position
            statsQ->Draw(); // Draw the stats box on the pad
        }
        gPad->Update();
    }

    canvasQ->SetGrid(); // Enable grid
    canvasQ->SaveAs("Figures/histograms_maximum_charge.png"); // Save the canvas as an image file

    


    /******************/
    /** Total Charge **/
    /******************/

    
    // Create a canvas for individual total charge histograms
    TCanvas *canvasTQ = new TCanvas("canvasTQ", "Histograms of Total Charge for Different Detectors", 1200, 1200);
    canvasTQ->Divide(2, 3); // Divide the canvas into a 2x3 grid

    // Grid for both axes
    canvasTQ->SetGridx();
    canvasTQ->SetGridy();

    // Draw each total charge histogram on the canvas
    for (size_t i = 0; i < ids.size(); ++i) {
        int id = ids[i];
        string name = idToName[id];
        
        canvasTQ->cd(i + 1); // Move to the next pad
        histogramsTQ[id]->SetLineColor(colors[i]); // Set the color of the histogram

        // Set the font size for the axes titles and labels
        float sizelabel = 0.035;
        histogramsTQ[id]->GetXaxis()->SetTitleSize(sizelabel);
        histogramsTQ[id]->GetXaxis()->SetLabelSize(sizelabel);
        histogramsTQ[id]->GetYaxis()->SetTitleSize(sizelabel);
        histogramsTQ[id]->GetYaxis()->SetLabelSize(sizelabel);

        histogramsTQ[id]->Draw(); // Draw the histogram
        gPad->SetGrid(); // Add grid to the current pad
        histogramsTQ[id]->GetXaxis()->SetRangeUser(0, 2000); // Set x-axis range

        // Move the statistics box
        gPad->Update(); // Update the pad to ensure the stats box is created
        TPaveStats *statsTQ = (TPaveStats*)histogramsTQ[id]->FindObject("stats");
        if (statsTQ) {
            statsTQ->SetX1NDC(0.7); // New X1 position
            statsTQ->SetX2NDC(0.9); // New X2 position
            statsTQ->SetY1NDC(0.7); // New Y1 position
            statsTQ->SetY2NDC(0.9); // New Y2 position
            statsTQ->Draw(); // Draw the stats box on the pad
        }
        gPad->Update();
    }

    canvasTQ->SetGrid(); // Enable grid
    canvasTQ->SaveAs("Figures/histograms_total_charge.png"); // Save the canvas as an image file




    // Clean up
    file->Close();
    delete file;
    delete canvasHits;
    delete logcanvasHits;
    delete canvasQ;
    delete canvasTQ;
    delete legendHits;
    // for (auto& hist : histogramsHits) delete hist.second;
    // for (auto& hist : histogramsQ) delete hist.second;
    // for (auto& hist : histogramsTQ) delete hist.second;
}
