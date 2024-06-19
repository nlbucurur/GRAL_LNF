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

    vector<int> ids = {8,9,10,11,12,13};

    map<int, string> idToName = {{8, "2"}, {9, "3"}, {10, "4"}, {11, "5"}, {12, "6"}, {13, "7"}};
    
    // Map to store a histogram for each detector ID
    map<int, TH1F*> histogramsQ;
    map<int, TH1F*> histogramsTQ;
    map<int, TH1I*> histogramsHits;


    
    for (int id : ids) {
        string name = idToName[id];
        
        histogramsHits[id] = new TH1I(Form("hCounts_%s", name.c_str()),
                                      Form("Counts of Detector in plane %s per Event; Number of hits; Number of Events", name.c_str()),
                                      121, 0, 121);

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

    
    TCanvas *canvas0 = new TCanvas("", "", 900, 600);
    //canvas0->SetLogy();
    
    // Draw and save each histogram
    for (int id : ids) {
        //string name = idToName[id];
        
        //TCanvas *canvas0 = new TCanvas(Form("canvas_%s", name.c_str()), Form("Hits in Detector in plane %s", name.c_str()), 900, 600);
        //canvas0->SetLogy();
        
        histogramsHits[id]->Draw();
        // //canvas0->SaveAs(Form("Figures/hits_detector_plane_%s.png", name.c_str()));
        // Int_t bin_under = histogramsHits[id]->GetBinContent(0); 
        // Int_t bin_over = histogramsHits[id]->GetBinContent(11);
        // cout << bin_under << "\t" << bin_over << endl;

        // // TCanvas *canvas1 = new TCanvas(Form("canvas_%s", name.c_str()), Form("Histogram of Maximum Charge for Detector in plane %s", name.c_str()), 900, 600);
        //histogramsQ[id]->Draw("same");
        // // canvas1->SaveAs(Form("Figures/maxQ_detector_plane_%s.png", name.c_str()));

        // // TCanvas *canvas2 = new TCanvas(Form("canvas_%s", name.c_str()), Form("Histogram of Sum of Maximum Charges for Detector in plane %s", name.c_str()), 900, 600);
        // histogramsTQ[id]->Draw("same");
        // // canvas2->SaveAs(Form("Figures/TotalQ_detector_plane_%s.png", name.c_str()));
    }
    canvas0->Draw();
}