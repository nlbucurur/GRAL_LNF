#include <iostream>
#include <vector>
#include <algorithm>           // Include for std::max_element
#include <cmath>               // Include for sqrt function
#include <map>
#include <ROOT/RDataFrame.hxx> // For modern data frame tool
#include <TH1F.h>
#include <TCanvas.h>

using namespace std;

// Main function to compute the efficiency of detectors, particularly focusing on detector 13.
void histograms()
{
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

    // Create a histogram to store the counts of detector 13 in each event
    TH1I *hCounts13 = new TH1I("hCounts13", "Counts of Detector 13 per Event;Number of hits;Number of Events", 10, 0, 10);
    TH1F *hMaxSum = new TH1F("hMaxSum13", "Sum of Maximum Charges for Detector 13;Charge;Number of Events", 200, 0, 2000);
    
    std::vector<int> ids = {8,9,10,11,12,13};
    
    // Map to store a histogram for each detector ID
    std::map<int, TH1F*> histogramsQ;
    std::map<int, TH1F*> histogramsTQ;

    for (int id : ids) {
        histogramsQ[id] = new TH1F(Form("hMaxQ_%d", id),
                                   Form("Maximum Charge of Detector %d;Charge (ADC);Number of Events", id),
                                   200, 0, 2000);

        histogramsTQ[id] = new TH1F(Form("hSumMaxCharge_%d", id), 
                                    Form("Sum of Maximum Charges for Detector %d;Charge;Number of Events", id), 
                                    200, 0, 2000);
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

    filtered.Foreach([&](const vector<unsigned int>& ids) {
        int count13 = count(ids.begin(), ids.end(), 13);
        hCounts13->Fill(count13);
    }, {"apv_id"});
    
    TCanvas *canvas = new TCanvas("Canvas hits", "Hits in Detector", 800, 600);
    hCounts13->Draw();
    canvas->SaveAs("Figures/hits_detector_13.png");



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


    /****************/
    /** Charge Sum **/
    /****************/

    filtered.Foreach([&](const vector<unsigned int>& apv_id, const vector<vector<short>>& apv_q) {
        std::map<int, double> sumMaxCharges; // Map to store the sum of max charges for each ID in this event

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

    
    
    // Draw and save each histogram
    for (int id : ids) {
        TCanvas *canvas = new TCanvas(Form("canvas_%d", id), Form("Histogram of Maximum Charge for Detector %d", id), 800, 600);
        histogramsQ[id]->Draw();
        canvas->SaveAs(Form("Figures/maxQ_detector_%d.png", id));

        TCanvas *canvas1 = new TCanvas(Form("canvas_%d", id), Form("Histogram of Sum of Maximum Charges for Detector %d", id), 800, 600);
        histogramsTQ[id]->Draw();
        canvas1->SaveAs(Form("Figures/TotalQ_detector_%d.png", id));
    }

    // Process the filtered data
    df.Foreach([&](const std::vector<unsigned int>& apv_id, const std::vector<std::vector<short>>& apv_q) {
        int sumMax = 0;
        for (size_t i = 0; i < apv_id.size(); ++i) {
            if (apv_id[i] == 13) {  // Check for detector 13 at index i
                if (i < apv_q.size() && !apv_q[i].empty()) {  // Ensure index is valid and vector is not empty
                    sumMax += *max_element(apv_q[i].begin(), apv_q[i].end()); // Find max and add to sum
                }
            }
        }

        if (sumMax > 0) {  // Only fill the histogram if the sum is greater than zero
            hMaxSum->Fill(sumMax);
        }
    }, {"apv_id", "apv_q"});

    TCanvas *canvas2 = new TCanvas("canvas_13", "Histogram of Maximum Charge for Detector 13", 800, 600);
    hMaxSum->Draw();
    canvas2->SaveAs("Figures/maxQ_detector_13Proof.png");

    cout << filtered.Describe();
}