#!/usr/bin/env python

import ROOT as rt
import sys

# Open the existing ROOT file
input_file_path = "Combined_cosmic_comp1thru4_vetoed.root"
input_file = rt.TFile(input_file_path, "READ")

# Retrieve the histograms
#integral_hist = input_file.Get("integral_hist")
integral_hist_MeV = input_file.Get("integral_hist_MeV")

# Define the target integral value (expected number of Cosmic events per year)
#target_integral = 230170
#4000 case target
target_integral = 212760
#target_integral = 230170

#Calculate the scaling factors needed to reach the target integral
current_integral_integral = integral_hist_MeV.Integral()
scaling_factor_integral = target_integral / current_integral_integral

# Scale the histograms by the calculated factors
integral_hist_MeV.Scale(scaling_factor_integral)

# Create and Write to output File
output_file_path = "scaled_Combined_cosmic_comp1thru4_vetoed.root"
output_file = rt.TFile(output_file_path, "RECREATE")

integral_hist_MeV.Write()
final_integral= integral_hist_MeV.Integral()

input_file.Close()
output_file.Close()

print("final integral", final_integral)
print("current integral",current_integral_integral)
print("Scaling factor =", scaling_factor_integral)
print(f"Scaled histograms saved to {output_file_path}")
