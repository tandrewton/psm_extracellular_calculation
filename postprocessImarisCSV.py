import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from scipy.spatial import Voronoi, Delaunay
from scipy.stats import mannwhitneyu
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from plotly.subplots import make_subplots
from scipy.optimize import curve_fit
from collections import defaultdict
import sys, os, glob
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
    # return trackIDs that persist through the requested time range
    #   (start frame to start frame + windowSize)
    unique_track_ids = T["TrackID"].unique()
    unique_track_ids = [int(k) for k in unique_track_ids]
    lifetimes = []
    trackStartStopTimes = np.zeros([len(unique_track_ids), 3])

    # loop over unique cell track IDs, calculate the starting and ending Time of each trackID
    for i, trackID in enumerate(unique_track_ids):
        filtered_data = T[T["TrackID"] == trackID]
        track_start = filtered_data["Time"].min()
        track_end = filtered_data["Time"].max()
        trackStartStopTimes[i] = [
            int(trackID),
            int(track_start),
            int(track_end),
        ]
        lifetimes.append(track_end - track_start + 1)

    # calculate # tracks along all windows of size windowSize
    # only record the trackIDs for a window beginning at startFrame and ending at startFrame+windowSize
    validIDs = []
    for i in range(1, T["Time"].max() - windowSize + 1):
        counter = 0
        for j in trackStartStopTimes:
            # check if track starts before i and ends after i+windowSize
            if j[1] <= i and j[2] > i + windowSize:
                counter += 1
                if i == startFrame:
                    validIDs.append(j[0])

    # plot lifetimes
    # plt.hist(lifetimes, bins=20, edgecolor="black")
    # plt.xlabel("cell track lifetime")
    # plt.ylabel("Frequency")
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
    ax.plot3D(x, y, z, color="k", lw="0.3")
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], color="r", alpha=0.5, s=5)


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


def interpolateMissingTracks(T, frameNum, missingTrackIDs):
    newRows = []
    for id in missingTrackIDs:
        # if one of frameNum-1 or frameNum+1 does not exist, it should just default to the existing one, which is interpolation when both exist, and copying if only one exists
        rowsToInterpolate = T[
            (T["TrackID"] == id)
            & ((T["Time"] == frameNum - 1) | (T["Time"] == frameNum + 1))
        ].copy()
        # store average between the two rows in column_averages, then add it to the dataframe
        column_averages = {}
        for column in rowsToInterpolate.columns:
            unique_values = rowsToInterpolate[column].unique()
            # print(rowsToInterpolate[column])
            # Check if the column has more than one unique value
            if len(unique_values) > 1:
                # Calculate the average of the column
                average = np.mean(unique_values)
                column_averages[column] = average
            elif len(unique_values) == 1:
                # If only one unique value, store it directly
                column_averages[column] = unique_values[0]
        newRows.append(column_averages)
    return newRows


def calculateSpeedsPerFrame(T):
    # output: np array of speeds, each entry is the average speed of a given frame
    speedPerFrame = []
    frames = T["Time"].unique()
    for frame in frames:
        filtered_data = T[T["Time"] == frame]
        speeds = np.mean(filtered_data["speed"])
        # speeds = np.mean(np.sqrt(filtered_data["velx"]**2 + filtered_data["velz"]**2))
        speedPerFrame.append(speeds)
    return np.array(speedPerFrame)


def calculateSpeedsPerEmbryo(T):
    # output: np array of speeds, each entry is the average speed of a given frame
    speedPerEmbryo = np.array([])
    frames = T["Time"].unique()
    for frame in frames:
        filtered_data = T[T["Time"] == frame]
        # speeds = np.mean(np.sqrt(filtered_data["velx"]**2 + filtered_data["velz"]**2))
        speedPerEmbryo = np.concatenate(
            (speedPerEmbryo, filtered_data["speed"].to_numpy())
        )
    return np.array(speedPerEmbryo)


def calculateNeighborExchanges(T, trackIDs, startFrame, windowSize):
    # for tracks labeled by trackIDs, which are guaranteed to be born by startFrame and to die after startFrame + windowSize,
    #  compute neighbor exchanges that occur between frames, starting at startFrame and ending at startFrame + windowSize
    filtered_data = T[T["TrackID"].isin(trackIDs)]
    globalXLim = min(filtered_data["posx"]), max(filtered_data["posx"])
    globalYLim = min(filtered_data["posy"]), max(filtered_data["posy"])
    globalZLim = min(filtered_data["posz"]), max(filtered_data["posz"])
    # plot_3d_movie(filtered_data)

    DelaunayEdgeDistanceThreshold = 12  # microns
    distances = []
    numberNearestNeighbors = []
    filteredNeighborsPerFrame = []
    neighborExchangePerFrame = []

    for frame in range(startFrame, startFrame + windowSize):
        frame_data = filtered_data[filtered_data["Time"] == frame]
        # calculate which trackIDs are not in frame_data
        missingTrackIDs = set(filtered_data["TrackID"]) - set(frame_data["TrackID"])
        #  use interpolation to generate a row for each of those missing trackIDs in this frame
        newRows = pd.DataFrame(interpolateMissingTracks(T, frame, missingTrackIDs))
        # add the new rows, then sort by trackID values in order to have a consistent order
        frame_data = pd.concat([frame_data, newRows], ignore_index=True)
        frame_data_sorted = frame_data.sort_values(by="TrackID")
        frame_data_raw = T[T["Time"] == frame]
        """
        # to compare between filtered+interpolated vs raw data,
        #   compare frame_data_sorted and T[T["Time"] == frame]
        # First histogram
        plt.subplot(2, 3, 1)
        plt.hist(frame_data_sorted["velx"], bins=20, color="blue", alpha=0.7)
        plt.xlim(np.min(frame_data_raw["velx"]), np.max(frame_data_raw["velx"]))
        plt.title("x all")

        # Second histogram
        plt.subplot(2, 3, 4)
        plt.hist(frame_data_raw["velx"], bins=20, color="green", alpha=0.7)
        plt.xlim(np.min(frame_data_raw["velx"]), np.max(frame_data_raw["velx"]))
        plt.title("x filtered")

        # Third histogram
        plt.subplot(2, 3, 2)
        plt.hist(frame_data_sorted["vely"], bins=20, color="red", alpha=0.7)
        plt.xlim(np.min(frame_data_raw["vely"]), np.max(frame_data_raw["vely"]))
        plt.title("y all")

        # Fourth histogram
        plt.subplot(2, 3, 5)
        plt.hist(frame_data_raw["vely"], bins=20, color="purple", alpha=0.7)
        plt.xlim(np.min(frame_data_raw["vely"]), np.max(frame_data_raw["vely"]))
        plt.title("y filtered")

        # Fifth histogram
        plt.subplot(2, 3, 3)
        plt.hist(frame_data_sorted["velz"], bins=20, color="orange", alpha=0.7)
        plt.xlim(np.min(frame_data_raw["velz"]), np.max(frame_data_raw["velz"]))
        plt.title("z all")

        # Sixth histogram
        plt.subplot(2, 3, 6)
        plt.hist(frame_data_raw["velz"], bins=20, color="pink", alpha=0.7)
        plt.xlim(np.min(frame_data_raw["velz"]), np.max(frame_data_raw["velz"]))
        plt.title("z filtered")

        # Adjust layout and show the histograms
        plt.tight_layout()
        

        # label="all tracks, N = " + str(len(frame_data_raw["posx"])),
        # label="filtered tracks, N = " + str(len(frame_data_sorted["posx"])),
        mpl.rcParams["legend.fontsize"] = 10
        plt.legend()
        plt.show()
        # in red, show all tracks
        # in black, show filtered tracks
        """

        points = frame_data_sorted[["posx", "posy", "posz"]].values
        nanRowInds = frame_data_sorted.loc[frame_data_sorted["posx"].isna()].index
        tri = Delaunay(points)
        """
        fig = plt.figure(1)
        ax = plt.axes(projection="3d")
        ax.set_xlim3d(globalXLim[0], globalXLim[1])
        ax.set_ylim3d(globalYLim[0], globalYLim[1])
        ax.set_zlim3d(globalZLim[0], globalZLim[1])
        plot_tri(ax, points, tri, DelaunayEdgeDistanceThreshold)
        # Remove all ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        # Optional: Remove grey pane background
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor('white')
        ax.yaxis.pane.set_edgecolor('white')
        ax.zaxis.pane.set_edgecolor('white')
        # Remove the axes spines
        ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        # Optional: Remove grid lines
        ax.grid(False)


        debugpy.breakpoint()
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
        neighborExchanges = 0
        nextFrame = frame + 1
        for point, neighbors in filteredNeighborDict.items():
            diff = set(neighbors) ^ set(filteredNeighborsPerFrame[nextFrame][point])
            if diff:
                neighborExchanges += 1
        neighborExchangePerFrame.append(neighborExchanges)

    # plt.figure()
    # plt.hist(distances, bins=50, range=[0, 50], color="blue", alpha=0.7)
    # plt.xlabel(r"Distance between points ($\mu$m)")
    # plt.ylabel("Counts")

    # plt.figure()
    # plt.hist(numberNearestNeighbors, bins=15, range=[0, 15], color="blue", alpha=0.7)
    # plt.xlabel(r"Neighbors")
    # plt.ylabel("Counts")

    print("done calculating neighbor exchanges")

    return np.array(neighborExchangePerFrame), numberNearestNeighbors


def main():
    folder = "/Users/AndrewTon/Documents/YalePhD/projects/ZebrafishSegmentation/psm_extracellular_calculation/cellmotionproc"
    masterSpreadsheetsFolder = "/Users/AndrewTon/Documents/YalePhD/projects/ZebrafishSegmentation/psm_extracellular_calculation/masterspreadsheets_packingfrac_cellshape"

    # list of spreadsheets
    files = glob.glob(os.path.join(folder, "*.csv"))
    # list of genotypes sorted in same order as files
    # name of file does not have correct genotype info, but genotype is correct in last csv column

    genotypes = []
    # NE_filenames for Miriam, to ease sorting/plotting in R in her pipeline
    NE_filenames = []

    NE_rate = []
    exploded_NE_rate = []
    NE_std = []
    exploded_NE_std = []
    neighborCount = []
    neighborCount_std = []
    cell_speeds = []
    cell_speed_std = []
    cell_speeds_per_genotype = defaultdict(lambda: np.array([]))
    cell_NE_per_genotype = defaultdict(lambda: np.array([]))
    nSkip = 0  # skip every nSkip frames, i.e. keep Time=1, Time=1+n, Time=1+2n
    time_spacing = 3 * (nSkip + 1)  # minutes
    for filename in files:
        print("current file is: ", filename)
        T = pd.read_csv(filename)

        genotypes.append(T["genotype"][0])
        NE_filenames.append(filename.split("/")[-1])
        print("current genotype is ", T["genotype"][0])
        # filter dataset to subsample time
        firstAndSkippedFramesCondition = (T.Time % (nSkip + 1) == 0) | (T.Time == 1)
        T = T[firstAndSkippedFramesCondition]

        stack = [T["Time"].min()]
        for unique_time in np.sort(T["Time"].unique()):
            prevTime = stack.pop()
            if unique_time - prevTime > 1:
                T.loc[T["Time"] == unique_time, "Time"] = prevTime + 1
            stack.append(prevTime + 1)
        if "sample" not in T:
            T["sample"] = "itg_cdh_1"

        # Set plot settings
        plt.rcParams["figure.figsize"] = (10, 6)
        plt.rcParams["font.size"] = 24

        # Loop over unique samples
        unique_samples = T["sample"].unique()
        for sample_num in range(len(unique_samples)):
            print("current sample is :", unique_samples[sample_num])
            # Filter data for the current sample
            sample_data = T[T["sample"] == unique_samples[sample_num]]
            # plot_3d_movie(sample_data)

            # get the trackIDs of tracks that persist from frame 1 to half the movie duration
            # halfMovieDuration = int(sample_data["Time"].max() / 2)
            halfMovieDuration = int(sample_data["Time"].max() - 5)
            print(halfMovieDuration, sample_data["Time"].max())
            # halfMovieDuration = int(sample_data["Time"].max()) - 2
            persistingTrackIDs = validTrackIDs(sample_data, 1, halfMovieDuration)

            neighborExchangePerFrame, numNeighbors = calculateNeighborExchanges(
                sample_data, persistingTrackIDs, 1, halfMovieDuration
            )

            print(f"numNeighbors={np.mean(numNeighbors)}, {np.std(numNeighbors)}")

            print(f"neighborExchangePerFrame={neighborExchangePerFrame}")

            # divide by len track IDs to get per cell
            # divide by time spacing to get per minute
            neighborExchangePerFrame = (
                neighborExchangePerFrame / len(persistingTrackIDs) / time_spacing
            )

            print(
                "exchange per cell per minute : ",
                np.mean(neighborExchangePerFrame),
                " +/- ",
                np.std(neighborExchangePerFrame),
                "number of cells : ",
                len(persistingTrackIDs),
                "number of neighbors per cell : ",
                np.mean(numNeighbors),
                " +/- ", 
                np.std(numNeighbors),
                "total_time : ",
                time_spacing * halfMovieDuration,
            )

            NE_rate.append(np.mean(neighborExchangePerFrame))
            NE_std.append(np.std(neighborExchangePerFrame))
            exploded_NE_rate.append(neighborExchangePerFrame)
            exploded_NE_std.append(0)
            neighborCount.append(np.mean(numNeighbors))
            neighborCount_std.append(np.std(numNeighbors))

            # grab current list of NEs, concatenate new NEs onto the value and add back to dict
            current_NE = cell_NE_per_genotype[T["genotype"][0]]
            current_NE = current_NE[~np.isnan(current_NE)]
            cell_NE_per_genotype[T["genotype"][0]] = np.concatenate(
                (
                    current_NE,
                    neighborExchangePerFrame[~np.isnan(neighborExchangePerFrame)],
                )
            )

            speedsPerEmbryo = calculateSpeedsPerEmbryo(sample_data)
            cell_speeds.append(np.nanmean(speedsPerEmbryo))
            cell_speed_std.append(np.nanstd(speedsPerEmbryo))

            current_speeds = cell_speeds_per_genotype[T["genotype"][0]]
            current_speeds = current_speeds[~np.isnan(current_speeds)]
            cell_speeds_per_genotype[T["genotype"][0]] = np.concatenate(
                (current_speeds, speedsPerEmbryo[~np.isnan(speedsPerEmbryo)])
            )

        NE_current_file = pd.DataFrame({"NE": neighborExchangePerFrame})
        NE_current_file["genotype"] = T["genotype"][0]
        split_path = filename.split('cellmotionproc')
        NE_folder_filename = split_path[0]+'cellmotionproc/NE_data'+split_path[1]
        NE_current_file.to_csv(NE_folder_filename[:-4]+"_NE_.csv")

        nn_current_file = pd.DataFrame({"numNeighbor": neighborCount})
        nn_current_file["genotype"] = T["genotype"][0]
        split_path = filename.split('cellmotionproc')
        nn_folder_filename = split_path[0]+'cellmotionproc/nn_data'+split_path[1]
        nn_current_file.to_csv(nn_folder_filename[:-4]+"_nn_.csv")

    # exchange genotype names for shorter names for plotting simplicity
    genotype_replacements = {
        'wt': 'Wild type',
        'MZitg': r'itg$\alpha$5',
        'itga5--': r'itg$\alpha$5',
        'MZitgcdh': r'cdh2; itg$\alpha$5',
        'itga5--cdh2--': r'cdh2; itg$\alpha$5',
        'cdh2--fbn2b--': 'cdh2; Fbn2b',
        'cdh2--fn1a--fn1b--': 'cdh2; Fn1a;1b',
        'fn1a--fn1b--cdh2--': 'cdh2; Fn1a;1b',
        'cdh2MO-fn1a--fn1b--fbn2b--': 'cdh2-MO; Fn1a;1b; Fbn2b', 
        'cdh2MO-fbn2b--fn1a--fn1b--': 'cdh2-MO; Fn1a;1b; Fbn2b', 
        'cdh2MOfbn2b--fn1a--fn1b--': 'cdh2-MO; Fn1a;1b; Fbn2b', 
        'fbn2b': 'Fbn2b', 
        'fbn2b--': 'Fbn2b',
        'fbn2b--fn1a--fn1b--': 'Fn1a;1b; Fbn2b', 
        'fbn2b--_fn1a--_fn1b--': 'Fn1a;1b; Fbn2b', 
        'fn1a--fn1b--': 'Fn1a;1b'
    }

    # Function to filter genotypes based on the replacements dictionary
    def filter_genotype(genotype):
        return genotype_replacements.get(genotype, genotype)

    # Apply the filter function to each genotype in the list
    genotypes = [filter_genotype(genotype) for genotype in genotypes]
    
    NE_df = pd.DataFrame({'Genotype': genotypes, 'NE': NE_rate, 'NE_std': NE_std, 'filename': NE_filenames})
    exploded_NE_df = pd.DataFrame({'Genotype': genotypes, 'NE': exploded_NE_rate, 'NE_std': exploded_NE_std, 'filename': NE_filenames})
    nn_df = pd.DataFrame({'Genotype': genotypes, 'nn': neighborCount, 'nn_std': neighborCount_std, 'filename': NE_filenames})

    speed_df = pd.DataFrame({'Genotype': genotypes, 'speed': cell_speeds, 'speed_std': cell_speed_std})
    combined_df = pd.concat([NE_df, speed_df["speed"]], axis=1, sort=False)

    #plot_ordering = ['wt', 'MZitg', 'cdh2', 'MZitgcdh', 'cdh2--fbn2b--', 'cdh2--fn1a--fn1b--', 'cdh2MO-fn1a--fn1b--fbn2b--', 'fbn2b', 'fbn2b--fn1a--fn1b--', 'fn1a--fn1b--', 'wt_PZ']
    plot_ordering = ['Wild type', 'Fbn2b', 'Fn1a;1b', 'Fn1a;1b; Fbn2b', r'itg$\alpha$5', 'cdh2', 'cdh2; Fbn2b', 'cdh2; Fn1a;1b', 'cdh2-MO; Fn1a;1b; Fbn2b', r'cdh2; itg$\alpha$5']
    color_by_genotype = ['#B92137', '#EA8F19', '#ECC21C', '#009800','#00FF00', '#0072BD', '#00B0F6', '#DB90FF','#FF62BC', '#7E2F8E', '#969696']

    NE_df['Genotype'] = pd.Categorical(NE_df['Genotype'], categories=plot_ordering, ordered=True)
    NE_df = NE_df.sort_values('Genotype')
    NE_df = NE_df.dropna()
    exploded_NE_df['Genotype'] = pd.Categorical(exploded_NE_df['Genotype'], categories=plot_ordering, ordered=True)
    exploded_NE_df = exploded_NE_df.sort_values('Genotype')
    exploded_NE_df = exploded_NE_df.dropna()
    nn_df['Genotype'] = pd.Categorical(nn_df['Genotype'], categories=plot_ordering, ordered=True)
    nn_df = nn_df.sort_values('Genotype')
    nn_df = nn_df.dropna()

    speed_df["Genotype"] = pd.Categorical(
        speed_df["Genotype"], categories=plot_ordering, ordered=True
    )
    speed_df = speed_df.sort_values("Genotype")
    speed_df = speed_df.dropna()

    combined_df["Genotype"] = pd.Categorical(
        combined_df["Genotype"], categories=plot_ordering, ordered=True
    )
    combined_df = combined_df.sort_values("Genotype")
    combined_df = combined_df.dropna()

    plt.figure(figsize=(12,8))
    exploded_NE_df = exploded_NE_df.explode('NE')
    exploded_NE_df['NE'] = exploded_NE_df['NE'].astype(float)  # Converting 'NE' to float for plotting
    plt.scatter(exploded_NE_df['Genotype'], exploded_NE_df['NE'])
    plt.xticks(rotation=75)
    plt.xlabel("Genotype")
    plt.ylabel("NE rate (per cell per min)")
    plt.tight_layout()  # Automatically adjust subplot params

    plt.figure(figsize=(12,8))
    plt.scatter(nn_df['Genotype'], nn_df['nn'])
    plt.xticks(rotation=75)
    plt.xlabel("Genotype")
    plt.ylabel("n_neighbor/cell")
    plt.tight_layout()  # Automatically adjust subplot params

    for genotype in plot_ordering:
        print(r'Neighbor exchanges: cdh2; itg$\alpha$5,', genotype, mannwhitneyu(exploded_NE_df[exploded_NE_df["Genotype"] == r'cdh2; itg$\alpha$5']["NE"], exploded_NE_df[exploded_NE_df["Genotype"] == genotype]["NE"], alternative='greater'))
        print(r'Speeds: cdh2; itg$\alpha$5,', genotype, mannwhitneyu(speed_df[speed_df["Genotype"] == r'cdh2; itg$\alpha$5']["speed"], speed_df[speed_df["Genotype"] == genotype]["speed"], alternative='greater'))

    for key in cell_NE_per_genotype.keys():
        print(key)
        print(
            r"Neighbor exchanges: cdh2; itg$\alpha$5,",
            key,
            mannwhitneyu(
                cell_NE_per_genotype["MZitgcdh"],
                cell_NE_per_genotype[key],
                alternative="greater",
            ),
        )
        print(
            r"Speeds: cdh2; itg$\alpha$5,",
            key,
            mannwhitneyu(
                cell_speeds_per_genotype["MZitgcdh"],
                cell_speeds_per_genotype[key],
                alternative="greater",
            ),
        )

    NE_df.to_csv('NE_df.csv', index=False)
    exploded_NE_df.to_csv('exploded_NE_df.csv', index=False)
    nn_df.to_csv('nn_df.csv', index=False)

    # Convert 'Category' to a categorical type and assign numerical values
    NE_df['GenotypeNum'] = NE_df['Genotype'].astype('category').cat.codes
    exploded_NE_df['GenotypeNum'] = exploded_NE_df['Genotype'].astype('category').cat.codes
    nn_df['GenotypeNum'] = nn_df['Genotype'].astype('category').cat.codes
    offset_width = 0.15  # Adjust the spacing between points in the same category
    grouped = NE_df.groupby('GenotypeNum')
    offsets = {cat: np.linspace(-offset_width/2, offset_width/2, num=len(group)) for cat, group in grouped}
    NE_df['Offset'] = NE_df.apply(lambda row: offsets[row['GenotypeNum']][grouped.groups[row['GenotypeNum']].get_loc(row.name)], axis=1)
    nn_df['Offset'] = nn_df.apply(lambda row: offsets[row['GenotypeNum']][grouped.groups[row['GenotypeNum']].get_loc(row.name)], axis=1)
    

    plt.figure(figsize=(12, 8))
    plt.errorbar(
        NE_df["GenotypeNum"] + NE_df["Offset"],
        NE_df["NE"],
        yerr=NE_df["NE_std"],
        fmt="o",
        capsize=5,
    )
    # Set the x-ticks and their labels to the original string values
    plt.xticks(range(len(NE_df["Genotype"].unique())), NE_df["Genotype"].unique())
    plt.xticks(rotation=75)
    plt.xlabel("Genotype")
    plt.ylabel(r"NE rate (per cell per min)")
    plt.tight_layout()  # Automatically adjust subplot params

    plt.figure(figsize=(12, 8))
    plt.errorbar(
        nn_df["GenotypeNum"] + nn_df["Offset"],
        nn_df["nn"],
        yerr=nn_df["nn_std"],
        fmt="o",
        capsize=5,
    )
    # Set the x-ticks and their labels to the original string values
    plt.xticks(range(len(nn_df["Genotype"].unique())), nn_df["Genotype"].unique())
    plt.xticks(rotation=75)
    plt.xlabel("Genotype")
    plt.ylabel(r"n_neighbor/cell")
    plt.tight_layout()  # Automatically adjust subplot params

    plt.figure(figsize=(12, 8))
    plt.scatter(speed_df["Genotype"], speed_df["speed"])
    plt.xticks(rotation=75)
    plt.xlabel("Genotype")
    plt.ylabel(r"Speed $(\mu m/min)$")
    plt.tight_layout()  # Automatically adjust subplot params

    # Convert 'Category' to a categorical type and assign numerical values
    speed_df["GenotypeNum"] = speed_df["Genotype"].astype("category").cat.codes
    offset_width = 0.1  # Adjust the spacing between points in the same category
    grouped = speed_df.groupby("GenotypeNum")
    offsets = {
        cat: np.linspace(-offset_width / 2, offset_width / 2, num=len(group))
        for cat, group in grouped
    }
    speed_df["Offset"] = speed_df.apply(
        lambda row: offsets[row["GenotypeNum"]][
            grouped.groups[row["GenotypeNum"]].get_loc(row.name)
        ],
        axis=1,
    )

    plt.figure(figsize=(12, 8))
    plt.errorbar(
        speed_df["GenotypeNum"] + speed_df["Offset"],
        speed_df["speed"],
        yerr=speed_df["speed_std"],
        fmt="o",
        capsize=5,
    )
    # Set the x-ticks and their labels to the original string values
    plt.xticks(range(len(speed_df["Genotype"].unique())), speed_df["Genotype"].unique())
    plt.xticks(rotation=75)
    plt.xlabel("Genotype")
    plt.ylabel(r"Speed $(\mu m/min)$")
    plt.tight_layout()  # Automatically adjust subplot params

    plt.figure(figsize=(12, 8))
    sns.scatterplot(data=combined_df, x="speed", y="NE", hue="Genotype")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
    plt.xlabel(r"Speed $(\mu m/min)$")
    plt.ylabel(r"NE rate (per cell per min)")
    plt.tight_layout()  # Automatically adjust subplot params

    packingFractionCSV = pd.read_csv(
        masterSpreadsheetsFolder + "/packingfrac_transverse_summary.csv"
    )
    phi_genotypes = [
        filter_genotype(genotype) for genotype in packingFractionCSV["genotype"]
    ]
    phi_df = pd.DataFrame(
        {"Genotype": phi_genotypes, "mean": packingFractionCSV["mean"]}
    )
    phi_df["Genotype"] = pd.Categorical(
        phi_df["Genotype"], categories=plot_ordering, ordered=True
    )
    phi_df = phi_df.sort_values("Genotype")
    phi_df = phi_df.dropna()
    plt.figure(figsize=(12, 8))
    plt.scatter(phi_df["Genotype"], phi_df["mean"])
    plt.xticks(rotation=75)
    plt.xlabel("Genotype")
    plt.ylabel(r"$\phi$")
    plt.tight_layout()  # Automatically adjust subplot params

    shapeCSV = pd.read_csv(
        masterSpreadsheetsFolder + "/cellprops_summary_transverse.csv"
    )
    shape_genotypes = [filter_genotype(genotype) for genotype in shapeCSV["genotype"]]
    shape_df = pd.DataFrame(
        {"Genotype": shape_genotypes, "mean": shapeCSV["Circularity_mean"]}
    )
    shape_df["Genotype"] = pd.Categorical(
        shape_df["Genotype"], categories=plot_ordering, ordered=True
    )
    shape_df = shape_df.sort_values("Genotype")
    shape_df = shape_df.dropna()
    plt.figure(figsize=(12, 8))
    plt.scatter(shape_df["Genotype"], shape_df["mean"])
    plt.xticks(rotation=75)
    plt.xlabel("Genotype")
    plt.ylabel(r"C")
    plt.tight_layout()  # Automatically adjust subplot params

    debugpy.breakpoint()
    # need to test colors
    genotype_color_dict = {plot_ordering[i]: color_by_genotype[i] for i in range(len(plot_ordering))}

    fig, axs = plt.subplots(4, 1, sharex=True, figsize=(15, 15))
    plt.rcParams.update({"font.size": 30})
    plt.subplots_adjust(hspace=0)
    axs[0].scatter(phi_df["Genotype"], phi_df["mean"], c=[genotype_color_dict[i] for i in phi_df["Genotype"]])
    axs[0].set_ylabel(r"$\phi$")
    axs[0].set_ylim(0.5, 1.0)
    axs[0].set_yticks([0.6, 0.8, 1.0])
    axs[1].scatter(shape_df["Genotype"], shape_df["mean"], c=[genotype_color_dict[i] for i in shape_df["Genotype"]])
    axs[1].set_ylabel(r"C")
    axs[1].set_ylim(0.77, 0.93)
    axs[1].set_yticks([0.8, 0.84, 0.88, 0.92])
    axs[2].errorbar(speed_df['GenotypeNum']+speed_df['Offset'], speed_df['speed'], yerr=speed_df['speed_std'], fmt='o', capsize=0, zorder=1, ecolor='black')
    axs[2].scatter(speed_df['GenotypeNum']+speed_df['Offset'], speed_df['speed'], c=[genotype_color_dict[i] for i in speed_df["Genotype"]], zorder=2)
    axs[2].set_ylabel(r"v $(\mu m/min)$")
    axs[2].set_ylim(0,1)
    axs[2].set_yticks([0.25, 0.5, 0.75])
    axs[3].errorbar(NE_df['GenotypeNum']+NE_df['Offset'], NE_df['NE'], yerr=NE_df['NE_std'], fmt='o', capsize=0, zorder=1, ecolor='black')
    axs[3].scatter(NE_df['GenotypeNum']+NE_df['Offset'], NE_df['NE'], c=[genotype_color_dict[i] for i in NE_df["Genotype"]], zorder=2)
    #axs[3].scatter(exploded_NE_df['GenotypeNum'], exploded_NE_df['NE'])
    axs[3].set_ylabel(r"NE $(cell \cdot min)^{-1}$")
    axs[3].set_ylim(0, 0.22)
    # axs[3].set_yticks([0, 0.1, 0.2])
    axs[3].set_xticks(range(len(NE_df["Genotype"].unique())))
    axs[3].set_xticklabels(NE_df["Genotype"].unique())
    axs[3].tick_params(axis="x", rotation=60)
    plt.xlabel("Genotype")
    plt.tight_layout()  # Automatically adjust subplot params
    plt.savefig("stackedExperimentalOrderedPlots.eps", bbox_inches="tight")
    plt.show()

    phi_df.to_csv('phi_df_experimental.csv')
    shape_df.to_csv('shape_df_experimental.csv')
    speed_df.to_csv('speed_df_experimental.csv')
    NE_df.to_csv('NE_df_experimental.csv')

    print("done!")

    debugpy.breakpoint()


if __name__ == "__main__":
    main()
