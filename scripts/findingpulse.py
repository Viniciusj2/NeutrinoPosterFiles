#!/usr/bin/env python
# coding: utf-8

# This was converted in .py from a jupyter notebook file -- 5/30/2024 Vinicius Da Silva
import sys
import os
import numpy as np
import ROOT as rt
from scipy.signal import find_peaks
from math import sqrt
from ROOT import cenns as cenns
import plotly.graph_objects as go

cennsly_path = os.environ['CENNS10_BASEDIR'] + '/vistools/cennsly'
if cennsly_path not in sys.path:
    if os.path.exists(cennsly_path):
        print("Adding cennsly path: ", cennsly_path)
        sys.path.insert(0, cennsly_path)
    else:
        print("Warning did not find cennsly path: ", cennsly_path)
import cennsly
cenns.pytools.CENNSPyTools
from vis_seghits_root import vis_seghits_root

PLOT_HITWFM_AND_PULSES = True
VIS_SEGHITS = True

#rootfile_path = "/home/vdasil01/Research/Coherent/g4-coh-ar-750/output_cosmic_test/CosmicCrywDAQ_combined_27k.root"
#rootfile_path = "combined_cosmic_comp1thru4.root" 
#rootfile_path = "/home/vdasil01/Research/Coherent/g4-coh-ar-750/marley_newgeo_test_300.root"
rootfile_path = "/home/vdasil01/Research/Coherent/g4-coh-ar-750/marley_output_1000_test.root"
#rootfile_path = "/home/vdasil01/Research/Coherent/g4-coh-ar-750/output_cosmic_test/10MeV_Electron_1000.root"

rootfile = rt.TFile(rootfile_path)

CENNS = rootfile.Get("CENNS")
edeptree = rootfile.Get("EDepSimEvents")
nentries = CENNS.GetEntries()
edepentries = edeptree.GetEntries()

# Define the cryostat boundaries
cryostat_x_min, cryostat_x_max = -400, 400
cryostat_y_min, cryostat_y_max = -400, 400
cryostat_z_min, cryostat_z_max = -521, 521

# Sampling rate in MHz (from /dat/daq_configs/caen_vx2470_digitizer.json)
sampling_rate_MHz = 1.25

# Calculate the time per sample in microseconds
time_per_sample_us = 1 / sampling_rate_MHz  # This will be 0.8 microseconds per sample

# Offset in nanoseconds converted to microseconds
t0_ns = -100.0
t0_us = t0_ns / 1000.0
window_pretrig_samples = 2
#window_posttrig_samples = 18

window_posttrig_samples = 6

# Initialize  
waveform_integrals = []
peak_info = {}
pulse_counts = []
energy_integrals = []

energy_hist = rt.TH1D("energy_hist", "Total Energy Deposition in Peak Windows;Energy (MeV);Entries", 40, 0, 200)
integral_hist = rt.TH1D("integral_hist", "Integral Values of Waveforms;Integral Waveform Value;Entries", 20, 0, 6000)
integral_hist_MeV = rt.TH1D("integral_hist_MeV", "Integral Values of Waveforms;Integral Waveform Value (MeV);Entries", 75, 0, 75)
hnpulses_per_entry = rt.TH1D("hnpulses_per_entry", "Pulses found per entry", 10, 0, 10)

total_hist = rt.TH1D("total_hist", "Total Energy Deposition;Energy (MeV);Entries", 100, 0, 500)
inside_hist = rt.TH1D("inside_hist", "Energy Deposition Inside Cryostat;Energy (MeV);Entries", 100, 0, 500)
outside_hist = rt.TH1D("outside_hist", "Energy Deposition Outside Cryostat;Energy (MeV);Entries", 100, 0, 500)
# shorten run to debug
#nentries = 10

# Process waveform data to detect peaks
for ientry in range(nentries):  # nentries
    if ientry % 100 == 0:
        print(f"Processing entry {ientry}/{nentries}")

    CENNS.GetEntry(ientry)
    edeptree.GetEntry(ientry)

    # DAQ ANALYSIS ==============================
    daq_branch_v = CENNS.DAQ
    daq_data = {}
    for idaq in range(daq_branch_v.size()):
        daq = daq_branch_v.at(idaq)
        daqname = daq.name
        daqdict = cenns.pytools.CENNSPyTools.to_numpy(daq)
        daq_data[daqname] = daqdict

    # Plot data setup
    DAQ = "caen_vx2740"
    wfms = daq_data[DAQ]["waveforms"]

    # Sum the waveforms
    summed_waveform = np.sum(wfms, axis=0)

    summed_integral = np.sum( summed_waveform )
    if summed_integral==0.0:
        print("No pulses in waveform, skip entry")
        continue
    else:
        print("Entry[%d] summed_integral=%.2f"%(ientry,summed_integral))
    
    nbins = summed_waveform.shape[0]

    hwfm = rt.TH1D("hwfm","",nbins,0,summed_waveform.shape[0]*time_per_sample_us)
    for ii in range(nbins):
        hwfm.SetBinContent( ii+1, summed_waveform[ii] )

    # Initialize a list to track the indices of identified peaks
    identified_peaks = []
    peak_ranges = []

    # Process the summed waveform
    max_value = np.max(summed_waveform)
    for sample_idx, value in enumerate(summed_waveform):
        if value > 10:
            print(summed_waveform[sample_idx])
            print(sample_idx)
            print(ientry)
            # Define the time window around the first value greater than threshold 
            start_index = max(sample_idx - window_pretrig_samples, 0)
            end_index = min(sample_idx + window_posttrig_samples, summed_waveform.shape[0] - 1)
                
            # Check if the current peak overlaps with other peak ranges
            if not any(start_index <= r[1] and end_index >= r[0] for r in peak_ranges):
                identified_peaks.append(sample_idx)
                peak_ranges.append((start_index, end_index))

                # Convert the sample indices to time in microseconds and apply the offset
                start_time_us = start_index * time_per_sample_us + t0_us
                end_time_us = end_index * time_per_sample_us + t0_us

                # Sum the value of the waveform inside the time window
                integral_value = np.sum(summed_waveform[start_index:end_index + 1])
                waveform_integrals.append(integral_value)
                #integral_hist.Fill(integral_value)
                #integral_hist_MeV.Fill(integral_value / 60)
                # Zero out the part of the waveform used for integration
                summed_waveform[start_index:end_index + 1] = 0

                # Store information in a dictionary
                if f"Entry {ientry}" not in peak_info:
                    peak_info[f"Entry {ientry}"] = []

                peak_info[f"Entry {ientry}"].append({
                    "Peak value": max_value,
                    "Peak index": sample_idx,
                    "Integral value": integral_value,
                    "Start time (us)": start_time_us,
                    "End time (us)": end_time_us
                })

    # Track the number of pulses in this entry
    num_pulses_in_entry = len(identified_peaks)
    pulse_counts.append(num_pulses_in_entry)
    hnpulses_per_entry.Fill( num_pulses_in_entry )

    # EDEP ANALYSIS ==============================
    ndetectors = edeptree.Event.SegmentDetectors.size()
    print("NDETECTORS: ",ndetectors)    
    if ndetectors != 1:
        continue

    seghits = edeptree.Event.SegmentDetectors["edepseg"]
    nseghits = seghits.size()
    event_edep = 0.0

    # Check for peaks in the current entry
    if f"Entry {ientry}" in peak_info:
        entry_peak_info = peak_info[f"Entry {ientry}"]
        entry_npeak = len(entry_peak_info)
        print(f"Peaks found for EDep entry {ientry}: NUM={entry_npeak}")
        


        for ipeak,peak in enumerate(entry_peak_info):
            # reset the energy integral for the peak, not the entire event
            peak_edep = 0.0
            peintegral = peak["Integral value"]
            total_edep = 0.0
            inside_edep = 0.0
            outside_edep = 0.0

            for ihit in range(nseghits):
                hit = seghits.at(ihit)

                            # Get hit position in cm
                x, y, z = hit.GetStart()[0], hit.GetStart()[1], hit.GetStart()[2]
                
                # Print the x, y, z coordinates for debugging
                print(f"Hit[{ihit}] position: x = {x:.2f} cm, y = {y:.2f} cm, z = {z:.2f} cm")
                
                # Use the time component of the start and stop positions
                hit_start_time_us = hit.GetStart()[3] / 1e3  # Convert time to microseconds
                hit_stop_time_us = hit.GetStop()[3] / 1e3

                if (peak["Start time (us)"] <= hit_start_time_us <= peak["End time (us)"]) or (peak["Start time (us)"] <= hit_stop_time_us <= peak["End time (us)"]):
                    #print("Peak contains seghit[",ihit,"] pos=(%.2f,%.2f,%.2f)"%( hit.GetStart()[0], hit.GetStart()[1], hit.GetStart()[2]))
                    event_edep += hit.GetEnergyDeposit() 
                    peak_edep += hit.GetEnergyDeposit() 
                     # Check if the hit is inside the cryostat
                    if (cryostat_x_min <= x <= cryostat_x_max) and (cryostat_y_min <= y <= cryostat_y_max) and (cryostat_z_min <= z <= cryostat_z_max):
                        inside_edep += hit.GetEnergyDeposit() 
                    else:
                        outside_edep += hit.GetEnergyDeposit() 
            print("PEAK[",ipeak,"] Edep=",peak_edep," MeV twindow=[",hit_start_time_us,",",hit_stop_time_us,"] PMTintegral(",peintegral,")/6.0peMeVee=",peintegral/6.0," MeV")

        # Store the total energy deposition for this entry if there was any peak detected
        if outside_edep < 15:
            if inside_edep > 0:
                inside_hist.Fill(inside_edep)
            if outside_edep > 0:
                outside_hist.Fill(outside_edep)
            if event_edep > 0:
                energy_integrals.append(event_edep)
                energy_hist.Fill(event_edep)  
                integral_hist.Fill(integral_value)
                integral_hist_MeV.Fill(integral_value / 60)
        

    #if PLOT_HITWFM_AND_PULSES or VIS_SEGHITS:

        #if PLOT_HITWFM_AND_PULSES:
            # visualize pulses
            #print("Visualize waveform and pulses: integral=%.2f npulses=%d"%(hwfm.Integral(),num_pulses_in_entry))
            #c = rt.TCanvas("c","",1600,400)
            #c.Draw()
            #hwfm.Draw("hist")
    
            #hmax = hwfm.GetMaximum()
            #line_v = []
            #if f"Entry {ientry}" in peak_info:
            #    for pulse in peak_info[f"Entry {ientry}"]:
            #        start_us = pulse['Start time (us)']
            #        end_us = pulse['End time (us)']
            #        l1 = rt.TBox( start_us, 0.0, end_us, hmax )
            #        l1.SetFillStyle(0)
            #        l1.SetLineColor( rt.kRed )
            #        l1.Draw()
            #        line_v.append(l1)
            #c.Update()


        #temp = None
        #if VIS_SEGHITS:
            #temp = vis_seghits_root( seghits, entry="ENTRY %d"%(ientry))

        #print("[ENTER] to continue")
        #input()

# for ientry in range(edepentries):  # edepentries
#     #if ientry % 100 == 0:
#         #print(f"Processing EDep entry {ientry}/{edepentries}")

#     if ientry>=nentries:
#         break

#     edeptree.GetEntry(ientry)
#     ndetectors = edeptree.Event.SegmentDetectors.size()
#     print("NDETECTORS: ",ndetectors)    
#     if ndetectors != 1:
#         continue

#     seghits = edeptree.Event.SegmentDetectors["edepseg"]
#     event_edep = 0.0
#     nseghits = seghits.size()

#     # Check for peaks in the current entry
#     if f"Entry {ientry}" in peak_info:
#         entry_peak_info = peak_info[f"Entry {ientry}"]

#         print(f"Peaks found for EDep entry {ientry}: {entry_peak_info}")

#         for ihit in range(nseghits):
#             hit = seghits.at(ihit)

#             print("seghit[",ihit,"] pos=(%.2f,%.2f,%.2f)"%( hit.GetStart()[0], hit.GetStart()[1], hit.GetStart()[2]))
            
#             # Use the time component of the start and stop positions
#             hit_start_time_us = hit.GetStart()[3] / 1e3  # Convert time to microseconds
#             hit_stop_time_us = hit.GetStop()[3] / 1e3

#             # Check if the hit time is within any of the peak time windows
#             for peak in entry_peak_info:
#                 if (peak["Start time (us)"] <= hit_start_time_us <= peak["End time (us)"]) or (peak["Start time (us)"] <= hit_stop_time_us <= peak["End time (us)"]):
#                     event_edep += hit.GetEnergyDeposit()
#                     break 

#         # Store the total energy deposition for this entry if there was any peak detected
#         if event_edep > 0:
#             energy_integrals.append(event_edep)
#             energy_hist.Fill(event_edep)

#         if VIS_SEGHITS:
#             _ = vis_seghits_root( seghits, entry="ENTRY %d"%(ientry))

# Number of pulses
#num_pulses = len(waveform_integrals)
#print(f"Number of pulses: {num_pulses}")

# Calculate the mean number of pulses
#mean_pulses = np.mean(pulse_counts)
#print(f"Mean number of pulses per entry: {mean_pulses}")
#num_0_pulses = pulse_counts.count(0) 
#num_1_pulse = pulse_counts.count(1) 
#print(f"Num of Entries with 0 Pulses: {num_0_pulses}")
#print(f"Num of Entries with 1 Pulse: {num_1_pulse}")

# Print the peak information dictionary
for key, value in peak_info.items():
    print(f"{key}: {value}")

# Plot the histograms
# Plot the histograms
c1 = rt.TCanvas("c1", "Total Energy Deposition in Peak Windows", 800, 600)
energy_hist.Draw()
c1.Draw()

c2 = rt.TCanvas("c2", "Integral Values of Waveforms", 800, 600)
integral_hist.Draw()
c2.Draw()

c3 = rt.TCanvas("c3", "Integral Values of Waveforms MeV", 800, 600)
integral_hist_MeV.Draw()
c3.Draw()

canvas = rt.TCanvas("canvas", "Energy Deposition", 800, 600)
canvas.Divide(2, 2)

canvas.cd(1)
total_hist.Draw()

canvas.cd(2)
inside_hist.Draw()

canvas.cd(3)
outside_hist.Draw()

#canvas.SaveAs("energy_deposition_histograms.png")

# Save into output file
output_file = rt.TFile("Combined_marley_vetoed15.root", "RECREATE")
total_hist.Write()
inside_hist.Write()
outside_hist.Write()
energy_hist.Write()
integral_hist.Write()
integral_hist_MeV.Write()
output_file.Close()

