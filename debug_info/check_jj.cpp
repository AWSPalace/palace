#include <iostream>
#include <fstream>
#include <string>

// This simple diagnostic tool checks if a configuration file contains Josephson junctions
// and reports basic information that can help debug EPR analysis issues

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config.json>" << std::endl;
        return 1;
    }

    std::string config_file = argv[1];
    std::ifstream file(config_file);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << config_file << std::endl;
        return 1;
    }

    std::string line;
    bool has_josephson = false;
    int jj_count = 0;
    
    std::cout << "Scanning " << config_file << " for Josephson junction definitions..." << std::endl;
    
    while (std::getline(file, line)) {
        // Look for Josephson junction type in the config
        if (line.find("\"type\"") != std::string::npos && 
            line.find("\"josephson\"", line.find("\"type\"")) != std::string::npos) {
            has_josephson = true;
            jj_count++;
            std::cout << "Found Josephson junction: " << line << std::endl;
        }
    }
    
    if (has_josephson) {
        std::cout << "Found " << jj_count << " Josephson junctions in the configuration." << std::endl;
        std::cout << "EPR analysis with custom convergence should be active." << std::endl;
    } else {
        std::cout << "No Josephson junctions found in the configuration." << std::endl;
        std::cout << "EPR analysis with custom convergence won't be triggered." << std::endl;
    }
    
    return 0;
}