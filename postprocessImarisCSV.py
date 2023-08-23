import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import debugpy

# Read the CSV file
filename = "/Users/AndrewTon/Downloads/MSD_wt.csv"
T = pd.read_csv(filename)

# Set plot settings
plt.rcParams["figure.figsize"] = (10, 6)
plt.rcParams["font.size"] = 24

# Initialize variables
time_spacing = 3  # minutes

plt.figure()

# Loop over unique samples
unique_samples = T["sample"].unique()
for sample_num in range(len(unique_samples)):
    # Filter data for the current sample
    sample_data = T[T["sample"] == unique_samples[sample_num]]
    unique_track_ids = sample_data["TrackID"].unique()

    msd = np.full((len(unique_track_ids), sample_data["normTime"].max()), np.nan)

    # Loop over unique TrackIDs
    for ii in range(len(unique_track_ids)):
        track_data = sample_data[sample_data["TrackID"] == unique_track_ids[ii]]
        norm_time = track_data["normTime"].values
        msd_cell = track_data["MSD"].values
        if len(msd[ii, :]) > len(msd_cell):
            msd[ii, :] = np.pad(
                msd_cell,
                (0, len(msd[ii, :]) - len(msd_cell)),
                "constant",
                constant_values=np.nan,
            )

    msd = np.nanmean(msd, axis=0) / (np.pi * 25)
    msd = msd[~np.isnan(msd)]

    print("plotting msd")
    plt.plot(np.arange(len(msd)) * time_spacing, msd, linewidth=2)
    plt.xlabel("Time (min)")
    plt.ylabel("MSD/a0")

    # Fit a power law model
    def power_law_model(x, A, B):
        return A * x**B

    start_ind = round(len(msd) / 4)
    fit_params, _ = curve_fit(
        power_law_model, np.arange(start_ind, len(msd)) * time_spacing, msd[start_ind:]
    )
    fit_A, fit_B = fit_params
    print(f"Fit: A={fit_A:.2f}, B={fit_B:.2f}")

plt.xscale("log")
plt.yscale("log")

xref = np.logspace(0.2, 2.5, num=100)
plt.plot(xref, 10**-2 * xref, "k--", linewidth=2)
plt.plot(xref, 10**-2 * xref**2, "r--", linewidth=2)
debugpy.breakpoint()
print("closing program")

plt.show()
