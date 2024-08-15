#!/usr/bin/env python
# coding: utf-8

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

PLOT_HITWFM_AND_PULSES = False
VIS_SEGHITS = False

#rootfile_path = "/home/vdasil01/Research/Coherent/g4-coh-ar-750/output_cosmic_test/CosmicCrywDAQ_combined_27k.root"
#rootfile_path = "../CosmicCrywDAQ_combined_27k.root"
rootfile_path = "/home/vdasil01/Research/Coherent/g4-coh-ar-750/marley_output_10000.root"
#rootfile_path = "/home/vdasil01/Research/Coherent/g4-coh-ar-750/output_cosmic_test/10MeV_Electron_1000.root"

rootfile = rt.TFile(rootfile_path)

CENNS = rootfile.Get("CENNS")
edeptree = rootfile.Get("EDepSimEvents")
nentries = CENNS.GetEntries()
edepentries = edeptree.GetEntries()

# Sampling rate in MHz (from /dat/daq_configs/caen_vx2470_digitizer.json)
sampling_rate_MHz = 1.25

# Calculate the time per sample in microseconds
time_per_sample_us = 1 / sampling_rate_MHz  # This will be 0.8 microseconds per sample

# Offset in nanoseconds converted to microseconds
t0_ns = -100.0
t0_us = t0_ns / 1000.0
window_pretrig_samples = 2
window_posttrig_samples = 18

# Initialize  
waveform_integrals = []
peak_info = {}
pulse_counts = []
energy_integrals = []

energy_hist = rt.TH1D("energy_hist", "Total Energy Deposition in Peak Windows;Energy (MeV);Entries", 20, 0, 100)
integral_hist = rt.TH1D("integral_hist", "Integral Values of Waveforms;Integral Waveform Value;Entries", 20, 0, 6000)
integral_hist_veto = rt.TH1D("integral_hist_veto", "Integral Values of Waveforms;Integral Waveform Value (MeV);Entries", 20, 0, 100)
integral_hist_MeV = rt.TH1D("integral_hist_MeV", "Integral Values of Waveforms;Integral Waveform Value (MeV);Entries", 75, 0, 75)
hnpulses_per_entry = rt.TH1D("hnpulses_per_entry", "Pulses found per entry", 10, 0, 10)

# shorten run to debug
# nentries = 10

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

    summed_integral = np.sum(summed_waveform)
    if summed_integral == 0.0:
        print("No pulses in waveform, skip entry")
        continue
    else:
        print("Entry[%d] summed_integral=%.2f" % (ientry, summed_integral))

    nbins = summed_waveform.shape[0]

    hwfm = rt.TH1D("hwfm", "", nbins, 0, summed_waveform.shape[0] * time_per_sample_us)
    for ii in range(nbins):
        hwfm.SetBinContent(ii + 1, summed_waveform[ii])

    # Initialize a list to track the indices of identified peaks
    identified_peaks = []
    peak_ranges = []

    # Process the summed waveform
    max_value = np.max(summed_waveform)
    for sample_idx, value in enumerate(summed_waveform):
        if value > 10:
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
                integral_hist.Fill(integral_value)
                integral_hist_MeV.Fill(integral_value / 60)  # Conversion from PE to MeV
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
    hnpulses_per_entry.Fill(num_pulses_in_entry)

    # EDEP ANALYSIS ==============================
    ndetectors = edeptree.Event.SegmentDetectors.size()
    print("NDETECTORS: ", ndetectors)
    if ndetectors != 1:
        continue

    seghits = edeptree.Event.SegmentDetectors["edepseg"]
    nseghits = seghits.size()
    event_edep = 0.0

    # Check for peaks in the current entry
    if f"Entry {ientry}" in peak_info:
        entry_peak_info = peak_info[f"Entry {ientry}"]
        entry_npeak = len(entry_peak_info)

        veto_hit = False  # Veto flag

        for ipeak, peak in enumerate(entry_peak_info):
            # reset the energy integral for the peak, not the entire event
            peak_edep = 0.0
            peintegral = peak["Integral value"]

            for ihit in range(nseghits):
                hit = seghits.at(ihit)

                # Use the time component of the start and stop positions
                hit_start_time_us = hit.GetStart()[3] / 1e3  # Convert time to microseconds
                hit_stop_time_us = hit.GetStop()[3] / 1e3

                if (peak["Start time (us)"] <= hit_start_time_us <= peak["End time (us)"]) or (
                        peak["Start time (us)"] <= hit_stop_time_us <= peak["End time (us)"]):
                    event_edep += hit.GetEnergyDeposit()
                    peak_edep += hit.GetEnergyDeposit()

                    # Check if the energy deposition is in the veto system
                    if hit.GetStart()[2] > 1110 or hit.GetStart()[2] < -1110 or hit.GetStart()[1] > 300 or hit.GetStart()[1] < -300 or hit.GetStart()[0] > 300 or hit.GetStart()[0] < -300:  # Assuming z > veto_threshold indicates veto
                        veto_hit = True

            if veto_hit:
                break  # Exit early if veto hit is detected

        if not veto_hit:
            # Store the total energy deposition for this entry if there was any peak detected
            if event_edep > 0:
                energy_integrals.append(event_edep)
                energy_hist.Fill(event_edep)

            if event_edep < 1000:
                integral_hist_veto.Fill(integral_value / 60)

    if PLOT_HITWFM_AND_PULSES or VIS_SEGHITS:

        if PLOT_HITWFM_AND_PULSES:
            # visualize pulses
            print("Visualize waveform and pulses: integral=%.2f npulses=%d" % (hwfm.Integral(), num_pulses_in_entry))
            c = rt.TCanvas("c", "", 1600, 400)
            c.Draw()
            hwfm.Draw("hist")

            hmax = hwfm.GetMaximum()
            line_v = []
            if f"Entry {ientry}" in peak_info:
                for pulse in peak_info[f"Entry {ientry}"]:
                    start_us = pulse['Start time (us)']
                    end_us = pulse['End time (us)']
                    l1 = rt.TBox(start_us, 0.0, end_us, hmax)
                    l1.SetFillStyle(0)
                    l1.SetLineColor(rt.kRed)
                    l1.Draw()
                    line_v.append(l1)
            c.Update()

        temp = None
        if VIS_SEGHITS:
            temp = vis_seghits_root(seghits, entry="ENTRY %d" % (ientry))

        print("[ENTER] to continue")
        input()

# Number of pulses
num_pulses = len(waveform_integrals)
print(f"Number of pulses: {num_pulses}")

# Calculate the mean number of pulses
mean_pulses = np.mean(pulse_counts)
print(f"Mean number of pulses per entry: {mean_pulses}")
num_0_pulses = pulse_counts.count(0)
num_1_pulse = pulse_counts.count(1)
print(f"Num of Entries with 0 Pulses: {num_0_pulses}")
print(f"Num of Entries with 1 Pulse: {num_1_pulse}")

# Print the peak information dictionary
for key, value in peak_info.items():
    print(f"{key}: {value}")

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

c4 = rt.TCanvas("c4", "Integral Values of Waveforms MeV (Vetoed)", 800, 600)
integral_hist_veto.Draw()
c4.Draw()

# Save into output file
output_file = rt.TFile("NEW_Cosmic_wvetotest.root", "RECREATE")
energy_hist.Write()
integral_hist.Write()
integral_hist_MeV.Write()
integral_hist_veto.Write()
output_file.Close()

