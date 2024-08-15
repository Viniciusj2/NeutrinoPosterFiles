import ROOT

def change_histogram_name(filename, old_hist_name, new_hist_name):
    # Open the ROOT file in update mode
    file = ROOT.TFile(filename, "UPDATE")
    
    # Get the histogram
    histogram = file.Get(old_hist_name)
    
    if histogram:
        # Set a new name for the histogram
        histogram.SetName(new_hist_name)
        
        # Write the modified histogram back to the file
        file.Write("", ROOT.TObject.kOverwrite)
        print(f"Changed histogram name from '{old_hist_name}' to '{new_hist_name}' in file '{filename}'.")
    else:
        print(f"Histogram '{old_hist_name}' not found in file '{filename}'.")
    
    # Close the file
    file.Close()

# Example usage
filename = "ScaledNEW_10CRYwveto.root"
old_hist_name = "integral_hist"
new_hist_name = "integral_hist_MeV"

change_histogram_name(filename, old_hist_name, new_hist_name)

