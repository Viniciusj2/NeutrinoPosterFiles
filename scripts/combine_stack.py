import ROOT

# File and histogram names
file1_name = "scaled_4us_4000_New_Cosmic_hist.root"
file2_name = "scaled_1000_New_Marley_hist.root"
hist_names = ["integral_hist_MeV"]
output_file_name = "Cosmicvsmarley_4000.root"

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
    stack = ROOT.THStack(hist_name + "_stack", "Estimated Signal and Background Per Year")
    
    # Add the histograms to the stack
    stack.Add(hist1)
    stack.Add(hist2)
    
    # Create a canvas and draw the stack to initialize axis properties
    canvas = ROOT.TCanvas(hist_name + "_canvas", "Estimated Signal and Background Per Year", 800, 600)
    stack.Draw("HIST")  # Use "HIST" to draw the histograms as line plots

    # Set the axis ranges after drawing the stack Based on old Plot made for this comparison
    stack.GetXaxis().SetRangeUser(5, 55)
    stack.SetMinimum(0)
    stack.SetMaximum(7000)
    
    # Set the y-axis title
    stack.GetYaxis().SetTitle("Events/year")
    stack.GetXaxis().SetTitle("Energy (MeV)") 
    # Add a legend
    
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # Adjust the position as needed
    legend.AddEntry(hist1, "Cosmic Background", "l")
    legend.AddEntry(hist2, "Marley Signal", "l")
    legend.Draw()

    # Save the canvas
    canvas.SaveAs(hist_name + "_stacked.png")
    canvas.SaveAs(hist_name + "_stacked.pdf")

    # Write the entire canvas to the output file
    canvas.Write()

# Close the output file
output_file.Close()

# Close the input files
file1.Close()
file2.Close()


