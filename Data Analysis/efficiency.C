#include <iostream>
#include <unordered_set>
#include <ROOT/RDataFrame.hxx> // For modern data frame tool
#include <cmath>              // Include for sqrt function

// Main function to compute the efficiency of detectors, particularly focusing on detector 13.
void efficiency()
{	
	// Open the ROOT file containing the tree.
    TFile *file = TFile::Open("run6578.root");
    if (!file || file->IsZombie()) { // Check if the file is opened successfully.
        std::cerr << "Error opening file or file not found." << std::endl;
        return;
    }

	// Retrieve the tree named "raw" from the file.
    TTree *tree = (TTree*)file->Get("raw");
    if (!tree) { // Check if the tree exists.
        std::cerr << "Tree 'raw' not found in file." << std::endl;
        file->Close(); // Close the file if the tree is not found.
        return; // Exit the function.
    }

	// Create a data frame from the TTree for easier processing.
    ROOT::RDataFrame df(*tree);

	// Lambda function to check if detector with id 13 is activated in the event.
    auto hasDetector13 = [](const std::vector<unsigned int>& ids) {
        return std::find(ids.begin(), ids.end(), 13) != ids.end();
    };

	// Filter the data frame to include only those events where detector 13 is activated.
    auto filtered = df.Filter(hasDetector13, {"apv_id"});

	// Map to store counts of activations for each detector.
    std::unordered_map<int, int> counts;

	// Process each entry in the filtered data frame.
    filtered.Foreach([&](const vector<vector<short>& ids) {
        std::unordered_set<int> unique_ids(ids.begin(), ids.end()); // Use a set to avoid counting duplicates.
        for (int id : {8,9,10,11,12,13}) {
            if (unique_ids.count(id)) {
                counts[id]++; // Increment count for each detector found.
            }
        }
    }, {"apv_id"});

	// Get the number of entries where detector 13 was hit to calculate efficiencies.
    double n = filtered.Count().GetValue();

	// Output the efficiency results and their errors for each detector.
    for (int id : {8,9,10,11,12,13}) {
        int n = counts[id];
        double efficiency = N > 0 ? (n / N) * 100 : 0; // Calculate efficiency as a percentage.
        double error = N > 0 ? (1 / sqrt(N)) * sqrt(n / N * (1 - n / N)) * 100 : 0; // Calculate error as a percentage.

        // Print counts as integers without decimal places
        std::cout << "Detector " << id << "\t" << "Counts " << n << "\t"; // Display count as an integer

        // Set fixed point and two decimal places for efficiencies and errors
        std::cout << std::fixed << std::setprecision(2) 
                << "Err Counts " << sqrt(n) << "\t" << "Efficiencies " << efficiency << "%\t" << "Error Eff " << error << "%" << std::endl;
    }

}
