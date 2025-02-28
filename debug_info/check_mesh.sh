#!/bin/bash

# This script helps diagnose issues with mesh generation and the Palace solver

echo "Checking for mesh files..."
find . -name "*.msh" -type f | while read mesh_file; do
    size=$(stat -c%s "$mesh_file")
    if [ "$size" -eq 0 ]; then
        echo "WARNING: Empty mesh file: $mesh_file"
    else
        echo "Valid mesh file: $mesh_file (size: $size bytes)"
    fi
done

echo -e "\nChecking for *.json configuration files..."
find . -name "*.json" -type f | grep -v "package.json" | while read json_file; do
    echo "Config file: $json_file"
    
    # Check if this is likely a Palace configuration
    if grep -q "\"problem\"" "$json_file"; then
        echo "  - This appears to be a Palace config file"
        
        # Check for mesh file reference
        mesh_path=$(grep -o '"mesh_file"[[:space:]]*:[[:space:]]*"[^"]*"' "$json_file" | awk -F'"' '{print $4}')
        if [ -n "$mesh_path" ]; then
            echo "  - References mesh file: $mesh_path"
            
            # Check if the mesh file exists and is non-empty
            if [ -f "$mesh_path" ]; then
                size=$(stat -c%s "$mesh_path")
                if [ "$size" -eq 0 ]; then
                    echo "  - WARNING: Referenced mesh file is EMPTY"
                else
                    echo "  - Referenced mesh file exists and has content (size: $size bytes)"
                fi
            else
                echo "  - ERROR: Referenced mesh file does not exist!"
            fi
        else
            echo "  - No mesh_file entry found in this config"
        fi
        
        # Check for Josephson junctions
        if grep -q "\"josephson\"" "$json_file"; then
            echo "  - Contains Josephson junction definitions (EPR analysis should be active)"
        else
            echo "  - No Josephson junction definitions found"
        fi
    fi
    
    echo ""
done

echo -e "\nChecking for eig.csv files (eigenmode solutions)..."
find . -name "eig.csv" -type f | while read eig_file; do
    size=$(stat -c%s "$eig_file")
    if [ "$size" -eq 0 ]; then
        echo "WARNING: Empty eigenmode solution file: $eig_file"
    else
        echo "Valid eigenmode solution file: $eig_file (size: $size bytes)"
        echo "First few lines:"
        head -n 3 "$eig_file"
    fi
done

echo -e "\nSystem information:"
echo "- Disk space:"
df -h .
echo "- Memory:"
free -h