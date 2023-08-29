import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from scipy.spatial import Voronoi, Delaunay
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
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


def validTrackIDs(T, startFrame, windowSize):
    # loop over unique cell track IDs
    # calculate range of time points
    # return trackIDs that persist through the requested time range
    #   (start frame to start frame + windowSize)
    unique_track_ids = T["TrackID"].unique()
    unique_track_ids = [int(k) for k in unique_track_ids]
    lifetimes = []
    trackStartStopTimes = np.zeros([len(unique_track_ids), 3])

    for i, trackID in enumerate(unique_track_ids):
        filtered_data = T[T["TrackID"] == trackID]
        track_start = filtered_data["Time"].min()
        track_end = filtered_data["Time"].max()
        trackStartStopTimes[i] = [
            trackID,
            track_start,
            track_end,
        ]
        lifetimes.append(track_end - track_start + 1)

    validIDs = []
    for i in range(1, T["Time"].max() - windowSize + 1):
        counter = 0
        for j in trackStartStopTimes:
            if j[1] <= i and j[2] >= i + windowSize:
                counter += 1
                if i == startFrame:
                    validIDs.append(j[0])
        print(
            "starting at frame i = ",
            i,
            ", detected ",
            counter,
            " tracks that persist the duration of ",
            windowSize,
            " frames!",
        )

    # plot lifetimes
    plt.hist(lifetimes, bins=20, edgecolor="black")
    plt.xlabel("cell track lifetime")
    plt.ylabel("Frequency")
    return np.array(validIDs)


def create_voronoi_plot(cells, fig, ax):
    # fig, ax = plt.subplots()
    for i, cell in enumerate(cells):
        polygon = cell["vertices"]
        ax.fill(*zip(*polygon), facecolor="none", edgecolor="k", lw=0.2)
        ax.text(cell["original"][0], cell["original"][1], str(i))
    return None


def plot_tri(ax, points, tri, distanceThreshold):
    # plot edges of the delaunay triangulation, but omit any edges with distance greater than a distance threshold
    edges = collect_edges(tri)
    x = np.array([])
    y = np.array([])
    z = np.array([])
    for i, j in edges:
        if np.sqrt(np.sum((points[i] - points[j]) ** 2)) < distanceThreshold:
            x = np.append(x, [points[i, 0], points[j, 0], np.nan])
            y = np.append(y, [points[i, 1], points[j, 1], np.nan])
            z = np.append(z, [points[i, 2], points[j, 2], np.nan])
    ax.plot3D(x, y, z, color="k", lw="1")
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], color="r", alpha=0.5)


def collect_edges(tri):
    edges = set()

    def sorted_tuple(a, b):
        return (a, b) if a < b else (b, a)

    # Add edges of tetrahedron (sorted so we don't add an edge twice, even if it comes in reverse order).
    for i0, i1, i2, i3 in tri.simplices:
        edges.add(sorted_tuple(i0, i1))
        edges.add(sorted_tuple(i0, i2))
        edges.add(sorted_tuple(i0, i3))
        edges.add(sorted_tuple(i1, i2))
        edges.add(sorted_tuple(i1, i3))
        edges.add(sorted_tuple(i2, i3))
    return edges


def calculateNeighborExchanges(T, trackIDs, startFrame, windowSize):
    # for tracks labeled by trackIDs,
    #  compute neighbor exchanges that occur between frames, starting at startFrame and ending at startFrame + windowSize
    filtered_data = T[T["TrackID"].isin(trackIDs)]
    globalXLim = min(filtered_data["posx"]), max(filtered_data["posx"])
    globalYLim = min(filtered_data["posy"]), max(filtered_data["posy"])
    globalZLim = min(filtered_data["posz"]), max(filtered_data["posz"])
    # plot_3d_movie(filtered_data)

    DelaunayEdgeDistanceThreshold = 10  # microns
    distances = []
    numberNearestNeighbors = []
    filteredNeighborsPerFrame = []

    debugpy.breakpoint()

    for frame in range(startFrame, startFrame + windowSize):
        frame_data = filtered_data[filtered_data["Time"] == frame]
        points = frame_data[["posx", "posy", "posz"]].values
        print("len points = ", len(points))
        tri = Delaunay(points)
        """fig = plt.figure()
        ax = plt.axes(projection="3d")
        ax.set_xlim3d(globalXLim[0], globalXLim[1])
        ax.set_ylim3d(globalYLim[0], globalYLim[1])
        ax.set_zlim3d(globalZLim[0], globalZLim[1])
        plot_tri(ax, points, tri, DelaunayEdgeDistanceThreshold)
        """
        indptr_nb, nbs = tri.vertex_neighbor_vertices

        nbDict = {i: nbs[indptr_nb[i] : indptr_nb[i + 1]] for i in range(len(points))}
        nbDictFiltered = {}

        # loop through the dictionary; store neighbor distances for a histogram, and count/store neighbors that are below a distance threshold
        for point, neighbors in nbDict.items():
            thresholdedNeighbors = []
            for neighbor in neighbors:
                neighborDistance = np.sqrt(
                    np.sum((tri.points[point] - tri.points[neighbor]) ** 2)
                )
                distances.append(neighborDistance)
                if neighborDistance < DelaunayEdgeDistanceThreshold:
                    thresholdedNeighbors.append(neighbor)
            numberNearestNeighbors.append(len(thresholdedNeighbors))
            nbDictFiltered[point] = thresholdedNeighbors

        filteredNeighborsPerFrame.append(nbDictFiltered)

    # use filtered neighbor dictionaries to calculate neighbor exchanges
    for frame, filteredNeighborDict in enumerate(filteredNeighborsPerFrame[:-1]):
        nextFrame = frame + 1
        for point, neighbors in filteredNeighborDict.items():
            diff = set(neighbors) ^ set(filteredNeighborsPerFrame[nextFrame][point])
            if diff:
                print("testing set diff")
                debugpy.breakpoint()

    plt.figure()
    plt.hist(distances, bins=50, range=[0, 50], color="blue", alpha=0.7)
    plt.xlabel(r"Distance between points ($\mu$m)")
    plt.ylabel("Counts")

    plt.figure()
    plt.hist(numberNearestNeighbors, bins=15, range=[0, 15], color="blue", alpha=0.7)
    plt.xlabel(r"Neighbors")
    plt.ylabel("Counts")

    debugpy.breakpoint()
    print("done calculating neighbor exchanges")


def main():
    # Read the CSV file
    filename = "/Users/AndrewTon/Downloads/MSD_wt3.csv"
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

        # plot_3d_movie(sample_data)

        # get the trackIDs of tracks that persist from frame 1 to half the movie duration
        halfMovieDuration = int(sample_data["Time"].max() / 2)
        persistingTrackIDs = validTrackIDs(sample_data, 1, halfMovieDuration)

        calculateNeighborExchanges(
            sample_data, persistingTrackIDs, 1, halfMovieDuration
        )

        print("leaving sample loop")
        debugpy.breakpoint()

        # msd = calculate_msd(sample_data, unique_track_ids)
        # plt.plot(np.arange(len(msd)) * time_spacing, msd, linewidth=2)
        # plt.xlabel("Time (min)")
        # plt.ylabel("MSD/a0")

    # plot_power_law_fits(msd, time_spacing)
    # plt.xscale("log")
    # plt.yscale("log")

    print("closing program")

    plt.show()


if __name__ == "__main__":
    main()
