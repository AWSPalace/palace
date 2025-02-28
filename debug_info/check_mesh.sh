#!/bin/bash

# This script helps diagnose issues with mesh generation and the Palace solver
# Modified to only check specific directories: PalaceEdit4 and sims1-4

# Use these variables to specify which directories to check
PALACE_DIR="../PalaceEdit4"
SIM_DIRS=("../sims1" "../sims2" "../sims3" "../sims4")

check_directory() {
    local dir=$1
    echo "==============================================="
    echo "Checking directory: $dir"
    echo "==============================================="
    
    echo "Checking for mesh files..."
    find "$dir" -name "*.msh" -type f | while read mesh_file; do
        size=$(stat -c%s "$mesh_file")
        if [ "$size" -eq 0 ]; then
            echo "WARNING: Empty mesh file: $mesh_file"
        else
            echo "Valid mesh file: $mesh_file (size: $size bytes)"
        fi
    done

    echo -e "\nChecking for *.json configuration files..."
    find "$dir" -name "*.json" -type f | grep -v "package.json" | while read json_file; do
        echo "Config file: $json_file"
        
        # Check if this is likely a Palace configuration
        if grep -q "\"problem\"" "$json_file"; then
            echo "  - This appears to be a Palace config file"
            
            # Check for mesh file reference
            mesh_path=$(grep -o '"mesh_file"[[:space:]]*:[[:space:]]*"[^"]*"' "$json_file" | awk -F'"' '{print $4}')
            if [ -n "$mesh_path" ]; then
                echo "  - References mesh file: $mesh_path"
                
                # Handle both absolute and relative paths
                if [[ "$mesh_path" = /* ]]; then
                    # Absolute path
                    check_path="$mesh_path"
                else
                    # Relative path - assume relative to the json file
                    json_dir=$(dirname "$json_file")
                    check_path="$json_dir/$mesh_path"
                fi
                
                # Check if the mesh file exists and is non-empty
                if [ -f "$check_path" ]; then
                    size=$(stat -c%s "$check_path")
                    if [ "$size" -eq 0 ]; then
                        echo "  - WARNING: Referenced mesh file is EMPTY"
                    else
                        echo "  - Referenced mesh file exists and has content (size: $size bytes)"
                    fi
                else
                    echo "  - ERROR: Referenced mesh file does not exist at: $check_path"
                fi
            else
                echo "  - No mesh_file entry found in this config"
            fi
            
            # Check for Josephson junctions
            if grep -q "\"josephson\"" "$json_file"; then
                echo "  - Contains Josephson junction definitions (EPR analysis should be active)"
                # Extract and display JJ parameters
                echo "  - Josephson junction details:"
                grep -A 10 "\"type\".*\"josephson\"" "$json_file" | grep -E "\"(name|surfaces|inductance|capacitance)\"" | sed 's/^/      /'
            else
                echo "  - No Josephson junction definitions found"
            fi
        fi
        
        echo ""
    done

    echo -e "\nChecking for eig.csv files (eigenmode solutions)..."
    find "$dir" -name "eig.csv" -type f | while read eig_file; do
        size=$(stat -c%s "$eig_file")
        if [ "$size" -eq 0 ]; then
            echo "WARNING: Empty eigenmode solution file: $eig_file"
        else
            echo "Valid eigenmode solution file: $eig_file (size: $size bytes)"
            echo "First few lines:"
            head -n 3 "$eig_file"
        fi
    done
    
    echo ""
}

# Check Palace directory if it exists
if [ -d "$PALACE_DIR" ]; then
    check_directory "$PALACE_DIR"
else
    echo "Palace directory not found at: $PALACE_DIR"
fi

# Check simulation directories
for sim_dir in "${SIM_DIRS[@]}"; do
    if [ -d "$sim_dir" ]; then
        check_directory "$sim_dir"
    else
        echo "Simulation directory not found at: $sim_dir"
    fi
done

echo -e "\nSystem information:"
echo "- Disk space:"
df -h .
echo "- Memory:"
free -h