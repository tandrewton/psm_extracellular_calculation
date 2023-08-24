import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.optimize import curve_fit
import debugpy


def calculate_msd(sample_data, unique_track_ids):
    max_norm_time = sample_data["normTime"].max()
    msd = np.full((len(unique_track_ids), max_norm_time), np.nan)

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
    return msd[~np.isnan(msd)]


def plot_power_law_fits(msd, time_spacing):
    def power_law_model(x, A, B):
        return A * x**B

    start_ind = round(len(msd) / 4)
    fit_params, _ = curve_fit(
        power_law_model, np.arange(start_ind, len(msd)) * time_spacing, msd[start_ind:]
    )
    fit_A, fit_B = fit_params
    print(f"Fit: A={fit_A:.2f}, B={fit_B:.2f}")

    xref = np.logspace(0.2, 2.5, num=100)
    plt.plot(xref, 10**-2 * xref, "k--", linewidth=2)
    plt.plot(xref, 10**-2 * xref**2, "r--", linewidth=2)


def plot_3d_movie(T):
    times = T["Time"].unique()
    # Calculate global axis ranges from the entire dataset
    x_range_global = [min(T["posx"]), max(T["posx"])]
    y_range_global = [min(T["posy"]), max(T["posy"])]
    z_range_global = [min(T["posz"]), max(T["posz"])]

    # Calculate aspect ratios based on global axis ranges
    x_range_length = x_range_global[1] - x_range_global[0]
    y_range_length = y_range_global[1] - y_range_global[0]
    z_range_length = z_range_global[1] - z_range_global[0]

    aspect_ratio_x = x_range_length / z_range_length
    aspect_ratio_y = y_range_length / z_range_length
    aspect_ratio_z = 1.0

    # Create a figure with a 3D scatter plot
    fig = go.Figure(
        data=[
            go.Scatter3d(
                x=[], y=[], z=[], mode="markers", marker=dict(color="red", size=1)
            )
        ]
    )
    # Define frames for the animation
    frames = []
    for t in times:
        filtered_data = T[T["Time"] == t]
        scatter = go.Scatter3d(
            x=filtered_data["posx"],
            y=filtered_data["posy"],
            z=filtered_data["posz"],
            mode="markers",
            name=f"Time {t}",
        )
        frame = go.Frame(data=[scatter], name=f"frame_{t}")
        frames.append(frame)

    frame_layout = {
        "scene": {
            "xaxis": {"range": x_range_global},
            "yaxis": {"range": y_range_global},
            "zaxis": {"range": z_range_global},
            "aspectmode": "manual",  # Set the desired aspect mode
            "aspectratio": {
                "x": aspect_ratio_x,
                "y": aspect_ratio_y,
                "z": aspect_ratio_z,
            },
        }
    }
    for frame in frames:
        frame.update(layout=frame_layout)

    debugpy.breakpoint()
    fig.frames = frames

    # Define animation settings
    animation_settings = go.layout.Updatemenu(
        type="buttons",
        showactive=False,
        buttons=[
            {
                "label": "Play",
                "method": "animate",
                "args": [
                    None,
                    {"frame": {"duration": 100, "redraw": True}, "fromcurrent": True},
                ],
            },
            {
                "label": "Pause",
                "method": "animate",
                "args": [
                    [None],
                    {
                        "frame": {"duration": 0, "redraw": True},
                        "mode": "immediate",
                        "transition": {"duration": 0},
                    },
                ],
            },
        ],
    )
    # Set initial frame
    fig.update(frames=frames)

    # Set layout and animation settings
    fig.update_layout(
        title="Position Data Animation",
        updatemenus=[animation_settings],
        scene=dict(
            xaxis_title="Position X",
            yaxis_title="Position Y",
            zaxis_title="Position Z",
        ),
    )

    # Show the plot
    fig.show()


def main():
    # Read the CSV file
    filename = "/Users/AndrewTon/Downloads/MSD_wt.csv"
    time_spacing = 3  # minutes
    T = pd.read_csv(filename)

    # Set plot settings
    plt.rcParams["figure.figsize"] = (10, 6)
    plt.rcParams["font.size"] = 24

    # Loop over unique samples
    unique_samples = T["sample"].unique()
    for sample_num in range(len(unique_samples)):
        # Filter data for the current sample
        sample_data = T[T["sample"] == unique_samples[sample_num]]
        unique_track_ids = sample_data["TrackID"].unique()
        plot_3d_movie(sample_data)
        print("leaving sample loop")

        # msd = calculate_msd(sample_data, unique_track_ids)
        # plt.plot(np.arange(len(msd)) * time_spacing, msd, linewidth=2)
        # plt.xlabel("Time (min)")
        # plt.ylabel("MSD/a0")

    # plot_power_law_fits(msd, time_spacing)
    # plt.xscale("log")
    # plt.yscale("log")

    debugpy.breakpoint()
    print("closing program")

    plt.show()


if __name__ == "__main__":
    main()
