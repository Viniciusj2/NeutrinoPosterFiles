import ROOT

# File and histogram names
file1_name = "scaled_Combined_cosmic_comp1thru4_vetoed.root"
file2_name = "scaled_1000_New_Marley_hist.root"
hist_names = ["integral_hist_MeV"]
output_file_name = "Cosmicvsmarley_4000_vetoed.root"

# Desired bin width
desired_bin_width = 5  # MeV

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
    
    # Get current bin width and range
    current_bin_width = hist1.GetBinWidth(1)
    x_min = hist1.GetXaxis().GetXmin()
    x_max = hist1.GetXaxis().GetXmax()
    
    # Calculate the number of bins with the desired width
    n_bins = int((x_max - x_min) / desired_bin_width)

    # Create new histograms with the desired binning
    hist1_rebinned = ROOT.TH1F(hist1.GetName() + "_rebinned", hist1.GetTitle(), n_bins, x_min, x_max)
    hist2_rebinned = ROOT.TH1F(hist2.GetName() + "_rebinned", hist2.GetTitle(), n_bins, x_min, x_max)

    # Fill the new histograms with the data from the original histograms
    for i in range(1, hist1.GetNbinsX() + 1):
        content1 = hist1.GetBinContent(i)
        content2 = hist2.GetBinContent(i)
        bin_center = hist1.GetBinCenter(i)
        hist1_rebinned.Fill(bin_center, content1)
        hist2_rebinned.Fill(bin_center, content2)

    # Set different line colors and styles
    hist1_rebinned.SetLineColor(ROOT.kRed)
    hist2_rebinned.SetLineColor(ROOT.kBlue)
    hist1_rebinned.SetLineWidth(1)
    hist2_rebinned.SetLineWidth(1)
    
    # Remove fill color if set (optional, for clarity)
    hist1_rebinned.SetFillColor(0)
    hist2_rebinned.SetFillColor(0)

    # Create a THStack for the current histogram
    stack = ROOT.THStack(hist_name + "_stack", "Estimated Signal and Background Per Year")
    
    # Add the histograms to the stack
    stack.Add(hist1_rebinned)
    stack.Add(hist2_rebinned)
    
    # Create a canvas and draw the stack to initialize axis properties
    canvas = ROOT.TCanvas(hist_name + "_canvas", "Estimated Signal and Background Per Year", 800, 600)
    stack.Draw("HIST")  # Use "HIST" to draw the histograms as line plots

    # Set the axis ranges after drawing the stack Based on old Plot made for this comparison
    stack.GetXaxis().SetRangeUser(0, 75)
    stack.SetMinimum(0)
    stack.SetMaximum(100)
    
    # Set the y-axis title
    stack.GetYaxis().SetTitle("Events/year")
    stack.GetXaxis().SetTitle("Energy (MeV)") 
    
    # Add a legend
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # Adjust the position as needed
    legend.AddEntry(hist1_rebinned, "Cosmic Background", "l")
    legend.AddEntry(hist2_rebinned, "Marley Signal", "l")
    legend.Draw()

    # Save the canvas
    #canvas.SaveAs(hist_name + "_stacked.png")
    #canvas.SaveAs(hist_name + "_stacked.pdf")

    # Write the entire canvas to the output file
    canvas.Write()

# Close the output file
output_file.Close()

# Close the input files
file1.Close()
file2.Close()

