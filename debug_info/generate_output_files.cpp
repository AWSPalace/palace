#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>

// Simple helper to generate placeholder output files for Palace simulations
// Compile with: g++ -std=c++17 -o generate_output_files generate_output_files.cpp
// Run with: ./generate_output_files /path/to/output/directory eigenfrequency_ghz

namespace fs = std::filesystem;

void createDirectory(const std::string& path) {
    if (!fs::exists(path)) {
        std::cout << "Creating directory: " << path << std::endl;
        fs::create_directories(path);
    }
}

void writeEigCsv(const std::string& path, double frequency) {
    std::ofstream file(path + "/eig.csv");
    if (file.is_open()) {
        file << "m,Re{f} (GHz),Im{f} (GHz),Q,Error (Bkwd.),Error (Abs.)\n";
        file << "1.000000000e+00," << frequency << ",+0.000000000e+00,+1.000000000e+10,+1.000000000e-10,+1.000000000e-10\n";
        std::cout << "Created " << path << "/eig.csv" << std::endl;
    } else {
        std::cerr << "Failed to create " << path << "/eig.csv" << std::endl;
    }
}

void writePortCsv(const std::string& path, double frequency) {
    // Port-EPR.csv
    std::ofstream eprFile(path + "/port-EPR.csv");
    if (eprFile.is_open()) {
        eprFile << "m,p[3]\n";
        eprFile << "1.000000000e+00,+6.700000000e-02\n";
        std::cout << "Created " << path << "/port-EPR.csv" << std::endl;
    }

    // Port-Q.csv
    std::ofstream qFile(path + "/port-Q.csv");
    if (qFile.is_open()) {
        qFile << "m,Q_ext[3],κ_ext[3] (GHz)\n";
        qFile << "1.000000000e+00,+5.000000000e+03,+1.000000000e-03\n";
        std::cout << "Created " << path << "/port-Q.csv" << std::endl;
    }

    // Port-V.csv
    std::ofstream vFile(path + "/port-V.csv");
    if (vFile.is_open()) {
        vFile << "m,Re{V[1]} (V),Im{V[1]} (V),Re{V[2]} (V),Im{V[2]} (V),Re{V[3]} (V),Im{V[3]} (V)\n";
        vFile << "1.000000000e+00,+1.000000000e+00,+0.000000000e+00,+5.000000000e-01,+0.000000000e+00,+2.000000000e-01,+0.000000000e+00\n";
        std::cout << "Created " << path << "/port-V.csv" << std::endl;
    }

    // Port-I.csv
    std::ofstream iFile(path + "/port-I.csv");
    if (iFile.is_open()) {
        iFile << "m,Re{I[1]} (A),Im{I[1]} (A),Re{I[2]} (A),Im{I[2]} (A),Re{I[3]} (A),Im{I[3]} (A)\n";
        iFile << "1.000000000e+00,+2.000000000e-02,+0.000000000e+00,+1.000000000e-02,+0.000000000e+00,+5.000000000e-03,+0.000000000e+00\n";
        std::cout << "Created " << path << "/port-I.csv" << std::endl;
    }
}

void writeDomainCsv(const std::string& path) {
    std::ofstream file(path + "/domain-E.csv");
    if (file.is_open()) {
        file << "m,E_e[1] (pJ),E_m[1] (pJ)\n";
        file << "1.000000000e+00,+5.000000000e-01,+5.000000000e-01\n";
        std::cout << "Created " << path << "/domain-E.csv" << std::endl;
    }
}

void writeErrorCsv(const std::string& path) {
    std::ofstream file(path + "/error-indicators.csv");
    if (file.is_open()) {
        file << "element,indicator\n";
        file << "1,0.1\n2,0.05\n3,0.01\n";
        std::cout << "Created " << path << "/error-indicators.csv" << std::endl;
    }
}

void writeSurfaceCsv(const std::string& path) {
    // Surface-Q.csv
    std::ofstream qFile(path + "/surface-Q.csv");
    if (qFile.is_open()) {
        qFile << "m,Re{Q[1]} (pC),Im{Q[1]} (pC)\n";
        qFile << "1.000000000e+00,+1.000000000e+00,+0.000000000e+00\n";
        std::cout << "Created " << path << "/surface-Q.csv" << std::endl;
    }

    // Surface-F.csv
    std::ofstream fFile(path + "/surface-F.csv");
    if (fFile.is_open()) {
        fFile << "m,Re{F[1]} (pH),Im{F[1]} (pH)\n";
        fFile << "1.000000000e+00,+1.000000000e+00,+0.000000000e+00\n";
        std::cout << "Created " << path << "/surface-F.csv" << std::endl;
    }
}

void createParaviewDir(const std::string& path) {
    std::string paraviewDir = path + "/paraview";
    createDirectory(paraviewDir);
    
    // Create dummy VTK file to show something exists
    std::ofstream file(paraviewDir + "/mode1.vtu");
    if (file.is_open()) {
        file << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n";
        file << "  <UnstructuredGrid>\n";
        file << "    <Piece NumberOfPoints=\"1\" NumberOfCells=\"0\">\n";
        file << "      <Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">0 0 0</DataArray></Points>\n";
        file << "      <Cells><DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"></DataArray></Cells>\n";
        file << "    </Piece>\n";
        file << "  </UnstructuredGrid>\n";
        file << "</VTKFile>\n";
        std::cout << "Created " << paraviewDir << "/mode1.vtu" << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " /path/to/output/directory eigenfrequency_ghz [iteration_count=3]" << std::endl;
        return 1;
    }

    std::string outputDir = argv[1];
    double frequency = std::stod(argv[2]);
    int iterations = 3;
    
    if (argc >= 4) {
        iterations = std::stoi(argv[3]);
    }

    // Make sure the output directory ends with /
    if (outputDir.back() != '/') {
        outputDir += '/';
    }

    // Create main output directory
    createDirectory(outputDir);
    
    // Create output files in main directory
    writeEigCsv(outputDir, frequency);
    writePortCsv(outputDir, frequency);
    writeDomainCsv(outputDir);
    writeErrorCsv(outputDir);
    writeSurfaceCsv(outputDir);
    createParaviewDir(outputDir);
    
    // Create iteration directories and their output files
    for (int i = 1; i <= iterations; i++) {
        std::string iterDir = outputDir + "iteration" + std::to_string(i);
        createDirectory(iterDir);
        
        writeEigCsv(iterDir, frequency * (1.0 + 0.01 * i));  // Slightly different frequency for each iteration
        writePortCsv(iterDir, frequency * (1.0 + 0.01 * i));
        writeDomainCsv(iterDir);
        writeErrorCsv(iterDir);
        writeSurfaceCsv(iterDir);
        createParaviewDir(iterDir);
    }

    std::cout << "Successfully generated all Palace output files in " << outputDir << std::endl;
    return 0;
}