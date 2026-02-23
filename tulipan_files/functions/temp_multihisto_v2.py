import sys
import numpy as np

def doane_bins(data):
    outparam = "OUTPATH/hist_param_do.txt"  # Output file for fitted parameters    
    """Compute the number of bins using Doane's formula."""
    n = len(data)
    if n < 2:
        raise ValueError("Not enough data points!")

    mean_val = np.mean(data)
    stddev_val = np.std(data, ddof=0)
    skew_val = np.mean(((data - mean_val) / stddev_val) ** 3)

    g1_std = np.sqrt(6 * (n - 2) / ((n + 1) * (n + 3)))  # Standard deviation of skewness
    bins = int(1 + np.log2(n) + np.log2(1 + abs(skew_val) / g1_std))
    bin_width = (np.max(data)-np.min(data))/bins
    with open(outparam, "w") as f:
        f.write(f"Doane Bins: {bins}\n")
        f.write(f"Doane BWidth: {bin_width:.5f}\n")
    
    return max(bins, 1),bin_width # Ensure at least 1 bin
    
def sturges_bins(data):
    outparam = "OUTPATH/hist_param_st.txt"  # Output file for fitted parameters  
    """Compute the number of bins using Sturges' rule."""
    n = len(data)
    if n < 2:
        raise ValueError("Not enough data points!")

    bins = int(1 + np.log2(n))
    bin_width = (np.max(data)-np.min(data))/bins
    
    with open(outparam, "w") as f:
        f.write(f"Sturges Bins: {bins}\n")
        f.write(f"Sturges BWidth: {bin_width:.5f}\n")
    
    return max(bins, 1),bin_width # Ensure at least 1 bin

def freedman_diaconis_bins(data):
    outparam = "OUTPATH/hist_param_fd.txt"  # Output file for fitted parameters    
    """Compute the number of bins using Freedman-Diaconis rule."""
    n = len(data)
    if n < 2:
        raise ValueError("Not enough data points!")

    sorted = np.sort(data)
    q1, q3 = np.percentile(sorted, [25, 75],method='closest_observation')  # Compute IQR
    iqr = q3 - q1
    bin_width = 2 * iqr / n**(1/3)
    
    if bin_width == 0:
        bins = int(np.sqrt(n))  # Fallback: use sqrt(n) rule if IQR is 0
    else:
        bins = int((max(data) - min(data)) / bin_width) + 1
    
    with open(outparam, "w") as f:
        f.write(f"FD Bins: {bins}\n")
        f.write(f"FD BWidth: {bin_width:.5f}\n")
    
    return max(bins, 1),bin_width # Ensure at least 1 bin

def absolute_to_relative(counts):
    """Convert an array of absolute counts to a relative population distribution."""
    total = np.sum(counts)  # Sum of all counts
    if total == 0:
        return np.zeros_like(counts)  # Avoid division by zero
    return counts / total  # Normalize
    

def generate_histogram(data, bins, out_file):
    """Generate histogram bins and save to file."""
    hist, bin_edges = np.histogram(data, bins=bins)
    rhist = absolute_to_relative(hist)
    with open(out_file, "w") as f:
        for i in range(len(hist)):
            bin_center = (bin_edges[i] + bin_edges[i+1]) / 2  # Compute bin center
            f.write(f"{bin_center:.2f} {rhist[i]}\n")

def generate_histogram_dns(data, bins, out_file):
    """Generate histogram bins and save to file."""
    hist, bin_edges = np.histogram(data, bins=bins, density=True)
    with open(out_file, "w") as f:
        for i in range(len(hist)):
            bin_center = (bin_edges[i] + bin_edges[i+1]) / 2  # Compute bin center
            f.write(f"{bin_center:.2f} {hist[i]}\n")


def main():

    """Main function to read data, compute histograms, and save to file."""
    if len(sys.argv) != 3:
        print("Usage: python multihisto_v2.py <datafile> <column> ")
        print("Methods included: doane | freedman | sturgers")
        sys.exit(1)
    outrep = "OUTPATH/report_histo.txt"  # Output file for fitted parameters 
    outrepdns = "OUTPATH/report_dns_histo.txt"  # Output file for fitted parameters 
    filename = sys.argv[1]
    column_index = int(sys.argv[2]) - 1  # Convert to zero-based index

    try:
        with open(filename, "r") as f:
            data = []
            for line in f:
                cols = line.strip().split()
                if len(cols) > column_index:
                    try:
                        data.append(float(cols[column_index]))
                    except ValueError:
                        continue  # Skip invalid data
        
        if not data:
            raise ValueError("No valid data found in the specified column.")
        #Relative population
        #Doane
        bins_d,bw_d = doane_bins(data)
        out_file = "OUTPATH/histo_doane.txt"
        generate_histogram(data, bins_d, out_file)
        #FD
        bins_fd,bw_fd = freedman_diaconis_bins(data)
        out_file = "OUTPATH/histo_fd.txt"
        generate_histogram(data, bins_fd, out_file)
        #Sturgers
        bins_stu,bw_stu = sturges_bins(data)
        out_file = "OUTPATH/histo_stu.txt"
        generate_histogram(data, bins_stu, out_file)
        with open(outrep, "w") as f:
            f.write(f"#SturgersBins #SturgerBinW #DoaneBins #DoaneBinW #FreedmanBins #FreedmanBinW\n")     
            f.write(f"{bins_stu:d} {bw_stu:.2f} {bins_d:d} {bw_d:.2f} {bins_fd:d} {bw_fd:.2f}\n")
        #Density Histogram    
        #Doane
        bins_d,bw_d = doane_bins(data)
        out_file = "OUTPATH/histo_dns_doane.txt"
        generate_histogram_dns(data, bins_d, out_file)
        #FD
        bins_fd,bw_fd = freedman_diaconis_bins(data)
        out_file = "OUTPATH/histo_dns_fd.txt"
        generate_histogram_dns(data, bins_fd, out_file)
        #Sturgers
        bins_stu,bw_stu = sturges_bins(data)
        out_file = "OUTPATH/histo_dns_stu.txt"
        generate_histogram_dns(data, bins_stu, out_file)
        with open(outrepdns, "w") as f:
            f.write(f"#SturgersBins #SturgerBinW #DoaneBins #DoaneBinW #FreedmanBins #FreedmanBinW\n")     
            f.write(f"{bins_stu:d} {bw_stu:.2f} {bins_d:d} {bw_d:.2f} {bins_fd:d} {bw_fd:.2f}\n")     
        
        #print(f"Histogram saved to {out_file}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
