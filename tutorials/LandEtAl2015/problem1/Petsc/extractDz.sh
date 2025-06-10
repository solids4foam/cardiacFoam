#!/bin/bash

# Output file
output="dz_values.dat"
echo "# mesh dz" > "$output"

# Loop through mesh2 to mesh5
for mesh in mesh2 mesh3 mesh4 mesh5; do
    file="$mesh/postProcessing/0/solidPointDisplacement_pointDisp.dat"
    
    if [[ -f "$file" ]]; then
        # Get the dz value: 4th entry of 2nd line
	dz=$(awk 'END { print $4 }' "$file")
        #dz=$#(awk 'NR==2 { print $4 }' "$file")
        echo "$mesh $dz" >> "$output"
    else
        echo "Warning: $file not found"
    fi
done
