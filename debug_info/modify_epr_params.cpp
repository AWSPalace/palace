#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <regex>

// This tool allows modifying the EPR convergence parameters in errorindicator.hpp

void modifyParameters(const std::string& filepath, double global_tol, double relative_tol, 
                     int required_consecutive, double jj_weight) {
    // Read the file content
    std::ifstream inFile(filepath);
    if (!inFile.is_open()) {
        std::cerr << "Error: Cannot open file " << filepath << std::endl;
        return;
    }
    
    std::string line;
    std::vector<std::string> content;
    
    while (std::getline(inFile, line)) {
        content.push_back(line);
    }
    inFile.close();
    
    // Create regex patterns for the parameters we want to modify
    std::regex global_tol_re("double global_tol\\{([^}]*)\\}");
    std::regex relative_tol_re("double relative_tol\\{([^}]*)\\}");
    std::regex required_consecutive_re("int required_consecutive\\{([^}]*)\\}");
    std::regex jj_weight_re("double jj_weight\\{([^}]*)\\}");
    
    // Process each line
    for (auto& line : content) {
        // Update global_tol
        line = std::regex_replace(line, global_tol_re, 
                                 "double global_tol{" + std::to_string(global_tol) + "}");
        
        // Update relative_tol
        line = std::regex_replace(line, relative_tol_re, 
                                 "double relative_tol{" + std::to_string(relative_tol) + "}");
        
        // Update required_consecutive
        line = std::regex_replace(line, required_consecutive_re, 
                                 "int required_consecutive{" + std::to_string(required_consecutive) + "}");
        
        // Update jj_weight
        line = std::regex_replace(line, jj_weight_re, 
                                 "double jj_weight{" + std::to_string(jj_weight) + "}");
    }
    
    // Write back to the file
    std::ofstream outFile(filepath);
    if (!outFile.is_open()) {
        std::cerr << "Error: Cannot write to file " << filepath << std::endl;
        return;
    }
    
    for (const auto& line : content) {
        outFile << line << std::endl;
    }
    
    std::cout << "Successfully updated EPR convergence parameters in " << filepath << std::endl;
    std::cout << "New values:" << std::endl;
    std::cout << "  global_tol: " << global_tol << std::endl;
    std::cout << "  relative_tol: " << relative_tol << std::endl;
    std::cout << "  required_consecutive: " << required_consecutive << std::endl;
    std::cout << "  jj_weight: " << jj_weight << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cout << "Usage: " << argv[0] << " <errorindicator.hpp> <global_tol> <relative_tol> <required_consecutive> <jj_weight>" << std::endl;
        std::cout << "Example: " << argv[0] << " /path/to/errorindicator.hpp 1e-3 1e-4 2 1.5" << std::endl;
        std::cout << std::endl;
        std::cout << "Parameters explanation:" << std::endl;
        std::cout << "  global_tol: Global error tolerance (default: 1e-2)" << std::endl;
        std::cout << "  relative_tol: Relative tolerance for consecutive iterations (default: 1e-3)" << std::endl;
        std::cout << "  required_consecutive: Number of consecutive iterations needed for convergence (default: 3)" << std::endl;
        std::cout << "  jj_weight: Weight factor for Josephson Junction elements (default: 3.0)" << std::endl;
        std::cout << std::endl;
        std::cout << "Recommended values for better performance:" << std::endl;
        std::cout << "  For faster convergence: 5e-3 5e-4 2 1.5" << std::endl;
        std::cout << "  For higher accuracy: 1e-3 1e-4 3 3.0" << std::endl;
        return 1;
    }
    
    std::string filepath = argv[1];
    double global_tol = std::stod(argv[2]);
    double relative_tol = std::stod(argv[3]);
    int required_consecutive = std::stoi(argv[4]);
    double jj_weight = std::stod(argv[5]);
    
    modifyParameters(filepath, global_tol, relative_tol, required_consecutive, jj_weight);
    
    return 0;
}