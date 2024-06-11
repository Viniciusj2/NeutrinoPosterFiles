#!/usr/bin/env python
# coding: utf-8

import ROOT as rt

# Open the existing ROOT file
input_file_path = "10_Marley_upto200MeV.root"
input_file = rt.TFile(input_file_path, "READ")

# Retrieve the histograms
energy_hist = input_file.Get("energy_hist")
integral_hist = input_file.Get("integral_hist")

# Check if histograms are successfully retrieved
if not energy_hist or not integral_hist:
    print("Error: Histograms not found in the input file.")
    input_file.Close()
    sys.exit(1)

# Scale the histograms by dividing the y-axis by livetime
scale_factor = 30.04
energy_hist.Scale(1.0 / scale_factor)
integral_hist.Scale(1.0 / scale_factor)

# Create a new output file
output_file_path = "10_Scaled_Marley_upto200MeV.root"
output_file = rt.TFile(output_file_path, "RECREATE")

# Write the scaled histograms to the output file
energy_hist.Write()
integral_hist.Write()

# Close the files
input_file.Close()
output_file.Close()

print(f"Scaled histograms saved to {output_file_path}")
