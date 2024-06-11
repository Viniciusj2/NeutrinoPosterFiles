import ROOT

# File and histogram names
file1_name = "10_Scaled_CRY_upto150MeV.root"
file2_name = "10_Scaled_Marley_upto150MeV.root"
hist_names = ["energy_hist", "integral_hist"]
output_file_name = "10_upto150MeV_cosmicvsmarley_comp.root"

# Open the ROOT files
file1 = ROOT.TFile.Open(file1_name)
file2 = ROOT.TFile.Open(file2_name)

# Create a new ROOT file to save the stacked histograms
output_file = ROOT.TFile(output_file_name, "RECREATE")

# Loop over each histogram name
for hist_name in hist_names:
    # Get the histograms from the files
    hist1 = file1.Get(hist_name)
    hist2 = file2.Get(hist_name)
    
    # Check if the histograms exist
    if not hist1 or not hist2:
        print(f"Histogram {hist_name} not found in one of the files.")
        continue
    
    # Set different line colors and styles
    hist1.SetLineColor(ROOT.kRed)
    hist2.SetLineColor(ROOT.kBlue)
    hist1.SetLineWidth(1)
    hist2.SetLineWidth(1)
    
    # Remove fill color if set (optional, for clarity)
    hist1.SetFillColor(0)
    hist2.SetFillColor(0)

    # Create a THStack for the current histogram
    stack = ROOT.THStack(hist_name + "_stack", f"Stacked {hist_name}")
    
    # Add the histograms to the stack
    stack.Add(hist1)
    stack.Add(hist2)
    
    # Write the stack to the output file
    stack.Write()
    
    # Draw the stack and save as image/PDF
    #canvas = ROOT.TCanvas(hist_name + "_canvas", f"Stacked {hist_name}", 800, 600)
    #stack.Draw("Hist")  # Use "HIST" to draw the histograms as line plots
    

    # Add a legend
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # Adjust the position as needed
    legend.AddEntry(hist1, "Cosmic Background", "l")
    legend.AddEntry(hist2, "Marley Signal", "l")
    legend.Draw()

    # Save the canvas
    #canvas.SaveAs(hist_name + "_stacked.png")
    #canvas.SaveAs(hist_name + "_stacked.pdf")

# Close the output file
output_file.Close()

# Close the input files
file1.Close()
file2.Close()


