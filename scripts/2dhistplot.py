import sys
import os
import numpy as np
import ROOT as rt
from ROOT import cenns as cenns

#modified on 7/16/24 To Serve as an Example on How to Access Waveform info ~Vinicius Da Silva 

cennsly_path = os.environ['CENNS10_BASEDIR'] + '/vistools/cennsly'
if cennsly_path not in sys.path:
    if os.path.exists(cennsly_path):
        print("Adding cennsly path: ", cennsly_path)
        sys.path.insert(0, cennsly_path)
    else:
        print("Warning did not find cennsly path: ", cennsly_path)
import cennsly
cenns.pytools.CENNSPyTools

# Define the ROOT file path and the TTree name
rootfile_path = "/home/vdasil01/Research/Coherent/g4-coh-ar-750/marley_output_10000.root"
rootfile = rt.TFile(rootfile_path)
CENNS = rootfile.Get("CENNS")

if not CENNS:
    print(f"Error: 'CENNS' TTree not found in {rootfile_path}")
    sys.exit(1)

nentries = CENNS.GetEntries()

# Create a 2D histogram
n_bins = 75  # Define the number of bins for both true and visible energy
comparison_hist = rt.TH2F("comparison_hist", "True vs Visible Energy;True Energy (MeV);Visible Energy (MeV)",
                          n_bins, 0, 75, n_bins, 0, 75)
integral_hist_MeV = rt.TH1D("integral_hist_MeV", "Integral Values of Waveforms;Integral Waveform Value (MeV);Entries", 75, 0, 75)

# Loop over the entries and fill the 2D histogram
for ientry in range(nentries):
    CENNS.GetEntry(ientry)

    # Extract the true energy
    true_energy = CENNS.Event.fEventEnergy
    
    # Extract the visible energy from waveform (In this script for simplicity we are just summing over the waveform for each entry instead of running the peak identification script)
    daq_branch_v = CENNS.DAQ
    wfms = None
    for idaq in range(daq_branch_v.size()):
        daq = daq_branch_v.at(idaq)
        if daq.name == "caen_vx2740":
            wfms = cenns.pytools.CENNSPyTools.to_numpy(daq)["waveforms"]
            break

    if wfms is None:
        continue

    # Sum the waveforms to get the visible energy
    summed_waveform = np.sum(wfms, axis=0)
    summed_integral = np.sum(summed_waveform)
    visible_energy = summed_integral / 60  # Convert  to MeV

    # Fill the 2D histogram
    
    comparison_hist.Fill(true_energy, visible_energy)
    integral_hist_MeV.Fill(visible_energy)
# Create a canvas to draw the histogram
canvas = rt.TCanvas("canvas", "True vs Visible Energy", 800, 600)
comparison_hist.Draw("COLZ")
canvas.SaveAs("comparison_hist2.png")  # Save as an image file

c3 = rt.TCanvas("c3", "Integral Values of Waveforms MeV", 800, 600)
integral_hist_MeV.Draw()
c3.Draw()

# Save the histogram to a ROOT file
output_file = rt.TFile("comparison_hist2.root", "RECREATE")
comparison_hist.Write()
integral_hist_MeV.Write()
output_file.Close()

# Close the ROOT file
rootfile.Close()










