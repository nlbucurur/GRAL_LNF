#include <iostream>
#include <vector>
#include <algorithm>           // Include for std::max_element
#include <cmath>               // Include for sqrt function
#include <map>
#include <ROOT/RDataFrame.hxx> // For modern data frame tool
#include <TH1F.h>
#include <TCanvas.h>

using namespace std;


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

    vector<int> ids = {8,9,10,11,12,13};
    
    // Map to store a histogram for each detector ID
    map<int, TH1F*> histogramsQ;
    map<int, TH1F*> histogramsTQ;
    map<int, TH1I*> histogramsHits;

    for (int id : ids) {
        histogramsHits[id] = new TH1I(Form("hCounts_%d", id),
                                      Form("Counts of Detector %d per Event;Number of hits;Number of Events", id),
                                      10, 0, 10);

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

    filtered.Foreach([&](const vector<unsigned int>& apv_id) {
        map<int, int> hitsCount;
        for (int id : ids) {
            hitsCount[id] = count(apv_id.begin(), apv_id.end(), id);
            histogramsHits[id]->Fill(hitsCount[id]);
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

    
    
    // Draw and save each histogram
    for (int id : ids) {
        TCanvas *canvas0 = new TCanvas(Form("canvas_%d", id), Form("Hits in Detector %d", id), 900, 600);
        histogramsHits[id]->Draw();
        canvas0->SaveAs(Form("Figures/hits_detector_%d.png", id)); 

        TCanvas *canvas1 = new TCanvas(Form("canvas_%d", id), Form("Histogram of Maximum Charge for Detector %d", id), 900, 600);
        histogramsQ[id]->Draw();
        canvas1->SaveAs(Form("Figures/maxQ_detector_%d.png", id));

        TCanvas *canvas2 = new TCanvas(Form("canvas_%d", id), Form("Histogram of Sum of Maximum Charges for Detector %d", id), 900, 600);
        histogramsTQ[id]->Draw();
        canvas2->SaveAs(Form("Figures/TotalQ_detector_%d.png", id));
    }
}