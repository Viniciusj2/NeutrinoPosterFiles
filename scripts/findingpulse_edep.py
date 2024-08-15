#!/usr/bin/env python
# coding: utf-8

import sys
import os
import numpy as np
import ROOT as rt

rootfile_path = "combined_cosmic_comp1thru4.root"
rootfile = rt.TFile(rootfile_path)

edeptree = rootfile.Get("EDepSimEvents")
edepentries = edeptree.GetEntries()

# Define the cryostat boundaries
cryostat_x_min, cryostat_x_max = -400, 400
cryostat_y_min, cryostat_y_max = -400, 400
cryostat_z_min, cryostat_z_max = -521, 521

# Initialize histograms
total_hist = rt.TH1D("total_hist", "Total Energy Deposition;Energy (MeV);Entries", 100, 0, 500)
inside_hist = rt.TH1D("inside_hist", "Energy Deposition Inside Cryostat;Energy (MeV);Entries", 100, 0, 500)
outside_hist = rt.TH1D("outside_hist", "Energy Deposition Outside Cryostat;Energy (MeV);Entries", 100, 0, 500)

# Loop over all entries
for ientry in range(edepentries):
    edeptree.GetEntry(ientry)
    
    seghits = edeptree.Event.SegmentDetectors["edepseg"]
    nseghits = seghits.size()
    
    total_edep = 0.0
    inside_edep = 0.0
    outside_edep = 0.0
    
    for ihit in range(nseghits):
        hit = seghits.at(ihit)
        energy = hit.GetEnergyDeposit()
        total_edep += energy

        x, y, z = hit.GetStart()[0], hit.GetStart()[1], hit.GetStart()[2]
        
        # Check if the hit is inside the cryostat
        if (cryostat_x_min <= x <= cryostat_x_max) and (cryostat_y_min <= y <= cryostat_y_max) and (cryostat_z_min <= z <= cryostat_z_max):
            inside_edep += energy
        else:
            outside_edep += energy

    # Fill histograms
    total_hist.Fill(total_edep)

    # Debugging: Check if energy is being added to both inside and outside
    if inside_edep > 0 and outside_edep > 0:
        print(f"Warning: Event {ientry} has energy deposition both inside and outside the cryostat.")
        print(f"    Inside energy: {inside_edep} MeV")
        print(f"    Outside energy: {outside_edep} MeV")
        print(f"    Total energy: {total_edep} MeV")
    
    if inside_edep > 0:
        inside_hist.Fill(inside_edep)
    if outside_edep > 0:
        outside_hist.Fill(outside_edep)

# Draw histograms (optional, requires a ROOT environment)
canvas = rt.TCanvas("canvas", "Energy Deposition", 800, 600)
canvas.Divide(2, 2)

canvas.cd(1)
total_hist.Draw()

canvas.cd(2)
inside_hist.Draw()

canvas.cd(3)
outside_hist.Draw()

canvas.SaveAs("energy_deposition_histograms.png")

# Save histograms to a ROOT file
output_file = rt.TFile("energy_deposition_histograms.root", "RECREATE")
total_hist.Write()
inside_hist.Write()
outside_hist.Write()
output_file.Close()


