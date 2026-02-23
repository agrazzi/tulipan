import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load data
file_path = "OUTPATH/HISTO"  # Update if needed
output_file = "OUTPATH/FITPRM"  # Output file for fitted parameters
output_plot="OUTPATH/PLOT"
data = np.loadtxt(file_path)

# Extract columns (x: bin centers, y: frequency counts)
x = data[:, 0]
y = data[:, 1]

# Define the exponential function to fit
def exponential(x, A, tau_):
    return A * np.exp(- x / tau_ )

# Perform curve fitting with error handling
try:
    popt, pcov = curve_fit(exponential, x, y, p0=(np.max(y), np.mean(x)))  
    A_fit, tau_fit = popt # Best-fit values
    A_err, tau_err = np.sqrt(np.diag(pcov)) # Errors (square root of diagonal of covariance matrix)
except RuntimeError:
    print("Error: Curve fitting failed. Could not find optimal parameters.")
    exit(1)  # Exit script if fitting fails
except Exception as e:
    print(f"Unexpected error during fitting: {e}")
    exit(1)


# Compute fitted values
y_fitted = exponential(x, A_fit, tau_fit)

# Compute RMSE
rmse = np.sqrt(np.mean((y - y_fitted) ** 2))

# Save results to a file
with open(output_file, "w") as f:
    f.write(f"Fitted Exponential Function: y = {A_fit:.6f} * exp(-{tau_fit:.6f} * x)\n")
    f.write(f"Estimated Parameters:\n")
    f.write(f"A       = {A_fit:.6f} ± {A_err:.6f}\n")
    f.write(f"Tau  = {tau_fit:.6f} ± {tau_err:.6f}\n")
    f.write(f"RMSE    = {rmse:.6f}\n")  # Save RMSE

# Generate smooth x values for plotting the fitted curve
x_fit = np.linspace(min(x), max(x), 100)
y_fit = exponential(x_fit, A_fit, tau_fit)

# Plot histogram and fitted curve
plt.scatter(x, y, label="Histogram Data", color="blue")
plt.plot(x_fit, y_fit, 
         label=fr"Fit: $y = {A_fit:.4f} e^{{-\frac{{t}}{{{tau_fit:.0f}}}}}$", 
         color="red")
plt.xlabel("Event Duration [ns]")
plt.xlim(0,(max(x)*1.05))
plt.ylabel("Relative Probability")
plt.title("Residence time Exponential Fit")
plt.legend()
plt.grid()

# Save the plot as PNG
plt.savefig(output_plot, dpi=240)  # High-quality image
plt.close()  # Close the plot to avoid inline display



