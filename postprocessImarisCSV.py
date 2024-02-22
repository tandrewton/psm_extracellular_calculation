import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from scipy.spatial import Voronoi, Delaunay
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from plotly.subplots import make_subplots
from scipy.optimize import curve_fit
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

            """for i in unique_values:
                if i != i and column != "...10":
                    print(
                        f"detected nan with id {id}, frameNum {frameNum}, column {column}"
                    )
            """
        newRows.append(column_averages)
    return newRows

def calculateSpeedsPerFrame(T):
    # output: np array of speeds, each entry is the average speed of a given frame
    speedPerFrame = []
    frames = T["Time"].unique()
    for frame in frames:
        filtered_data = T[T["Time"] == frame]
        speeds = np.mean(filtered_data["speed"])
        #speeds = np.mean(np.sqrt(filtered_data["velx"]**2 + filtered_data["velz"]**2))
        speedPerFrame.append(speeds)
    return np.array(speedPerFrame)

def calculateSpeedsPerEmbryo(T):
    # output: np array of speeds, each entry is the average speed of a given frame
    speedPerEmbryo = np.array([])
    frames = T["Time"].unique()
    for frame in frames:
        filtered_data = T[T["Time"] == frame]
        #speeds = np.mean(np.sqrt(filtered_data["velx"]**2 + filtered_data["velz"]**2))
        speedPerEmbryo = np.concatenate((speedPerEmbryo, filtered_data["speed"].to_numpy()))
    return np.array(speedPerEmbryo)


def calculateNeighborExchanges(T, trackIDs, startFrame, windowSize):
    # for tracks labeled by trackIDs, which are guaranteed to be born by startFrame and to die after startFrame + windowSize,
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
        fig = plt.figure()
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
        neighborExchanges = 0
        nextFrame = frame + 1
        for point, neighbors in filteredNeighborDict.items():
            diff = set(neighbors) ^ set(filteredNeighborsPerFrame[nextFrame][point])
            if diff:
                neighborExchanges += 1
        neighborExchangePerFrame.append(neighborExchanges)

    plt.figure()
    plt.hist(distances, bins=50, range=[0, 50], color="blue", alpha=0.7)
    plt.xlabel(r"Distance between points ($\mu$m)")
    plt.ylabel("Counts")

    plt.figure()
    plt.hist(numberNearestNeighbors, bins=15, range=[0, 15], color="blue", alpha=0.7)
    plt.xlabel(r"Neighbors")
    plt.ylabel("Counts")

    print("done calculating neighbor exchanges")

    return np.array(neighborExchangePerFrame)


def main():
    folder = "/Users/AndrewTon/Documents/YalePhD/projects/ZebrafishSegmentation/psm_extracellular_calculation/cellmotionproc"
    #wt_filename1 = folder + "MSD_wt.csv"
    #wt_filename2 = folder + "MSD_wt3.csv"
    #cdh_filename = folder + "cdh_MSD.csv"
    #itg_cdh_filename = folder + "itgcdh_230831.csv"
    # should just glob all file names ending in .csv
    '''wt_filename1 = folder + "wt_052312b .csv"
    wt_filename2 = folder + "wt_041211b .csv"
    wt_filename3 = folder + "wt_033111 .csv"
    cdh_filename1 = folder + "cdh2_20231106_1 .csv"
    cdh_filename2 = folder + "cdh2_20231106_2 .csv"
    cdh_fn_filename1 = folder + "cdh2--fn1a--fn1b--_20231113_2 .csv"
    itga5_filename1 = folder + "itga5_20231104_2 .csv"
    itga5_filename2 = folder + "itga5_20231202_3 .csv"
    fbn_filename1 = folder + "fbn2b_231011_em2 .csv"
    MZitg_filename1 = folder + "MZitg_231030_em1 .csv"
    MZitg_cdh_filename1 = folder + "MZitgcdh_230921 .csv"
    MZitg_cdh_filename2 = folder + "MZitgcdh2_230831 .csv"
    itg_cdh_filename = folder + "itgcdh_230831.csv"
    '''

    # list of spreadsheets
    files = glob.glob(os.path.join(folder, "*.csv"))
    # list of genotypes sorted in same order as files
    # name of file does not have correct genotype info, but genotype is correct in last csv column
    genotypes = []

    NE_rate = []
    NE_std = []
    cell_speeds = []
    cell_speed_std = []
    nSkip = 0  # skip every nSkip frames, i.e. keep Time=1, Time=1+n, Time=1+2n
    time_spacing = 3 * (nSkip + 1)  # minutes
    for filename in files:
        print("current file is: ", filename)
        T = pd.read_csv(filename)

        #if (T["genotype"][0] == "wt_PZ"):
        #    continue

        genotypes.append(T["genotype"][0])
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
            halfMovieDuration = int(sample_data["Time"].max() / 2)
            print(halfMovieDuration, sample_data["Time"].max())
            # halfMovieDuration = int(sample_data["Time"].max()) - 2
            persistingTrackIDs = validTrackIDs(sample_data, 1, halfMovieDuration)

            neighborExchangePerFrame = calculateNeighborExchanges(
                sample_data, persistingTrackIDs, 1, halfMovieDuration
            )

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
                "total_time : ",
                time_spacing*halfMovieDuration
            )


            NE_rate.append(np.mean(neighborExchangePerFrame))
            NE_std.append(np.std(neighborExchangePerFrame))
            
            speedsPerEmbryo = calculateSpeedsPerEmbryo(sample_data)
            cell_speeds.append(np.nanmean(speedsPerEmbryo))
            cell_speed_std.append(np.nanstd(speedsPerEmbryo))


            """
            plt.figure()
            msd = calculate_msd(sample_data, sample_data["TrackID"].unique())
            plt.plot(
                np.arange(len(msd)) * time_spacing,
                msd,
                linewidth=2,
                color="r",
                label="all tracks, N = " + str(len(sample_data["TrackID"].unique())),
            )
            msd_filtered = calculate_msd(sample_data, persistingTrackIDs)
            plt.plot(
                np.arange(len(msd_filtered)) * time_spacing,
                msd_filtered,
                linewidth=2,
                color="k",
                label="filtered tracks, N = " + str(len(persistingTrackIDs)),
            )
            plt.xlabel("Time (min)")
            plt.ylabel("MSD/a0")
            plt.title(unique_samples[sample_num])
            plt.legend()
            """

        # plot_power_law_fits(msd, time_spacing)
        # plt.xscale("log")
        # plt.yscale("log")

    
    NE_df = pd.DataFrame({'Genotype': genotypes, 'NE': NE_rate, 'NE_std': NE_std})
    speed_df = pd.DataFrame({'Genotype': genotypes, 'speed': cell_speeds, 'speed_std': cell_speed_std})
    combined_df = pd.concat([NE_df, speed_df["speed"]], axis=1, sort=False)

    debugpy.breakpoint()

    plot_ordering = ['wt', 'MZitg', 'cdh2', 'MZitgcdh', 'cdh2--fbn2b--', 'cdh2--fn1a--fn1b--', 'cdh2MO-fn1a--fn1b--fbn2b--', 'fbn2b', 'fbn2b--fn1a--fn1b--', 'fn1a--fn1b--', 'wt_PZ']
    NE_df['Genotype'] = pd.Categorical(NE_df['Genotype'], categories=plot_ordering, ordered=True)
    NE_df = NE_df.sort_values('Genotype')

    speed_df['Genotype'] = pd.Categorical(speed_df['Genotype'], categories=plot_ordering, ordered=True)
    speed_df = speed_df.sort_values('Genotype')

    combined_df['Genotype'] = pd.Categorical(combined_df['Genotype'], categories=plot_ordering, ordered=True)
    combined_df = combined_df.sort_values('Genotype')

    plt.figure(figsize=(12,8))
    plt.scatter(NE_df['Genotype'], NE_df['NE'])
    plt.xticks(rotation=45)
    plt.xlabel('Genotype')
    plt.ylabel('NE rate (per cell per min)')
    plt.tight_layout()  # Automatically adjust subplot params

    # Convert 'Category' to a categorical type and assign numerical values
    NE_df['GenotypeNum'] = NE_df['Genotype'].astype('category').cat.codes
    offset_width = 0.1  # Adjust the spacing between points in the same category
    grouped = NE_df.groupby('GenotypeNum')
    offsets = {cat: np.linspace(-offset_width/2, offset_width/2, num=len(group)) for cat, group in grouped}
    NE_df['Offset'] = NE_df.apply(lambda row: offsets[row['GenotypeNum']][grouped.groups[row['GenotypeNum']].get_loc(row.name)], axis=1)

    plt.figure(figsize=(12,8))
    plt.errorbar(NE_df['GenotypeNum']+NE_df['Offset'], NE_df['NE'], yerr=NE_df['NE_std'], fmt='o', capsize=5)
    # Set the x-ticks and their labels to the original string values
    plt.xticks(range(len(NE_df['Genotype'].unique())), NE_df['Genotype'].unique())
    plt.xticks(rotation=45)
    plt.xlabel('Genotype')
    plt.ylabel(r"NE rate (per cell per min)")
    plt.tight_layout()  # Automatically adjust subplot params

    plt.figure(figsize=(12,8))
    plt.scatter(speed_df['Genotype'], speed_df['speed'])
    plt.xticks(rotation=45)
    plt.xlabel('Genotype')
    plt.ylabel(r"Speed $(\mu m/min)$")
    plt.tight_layout()  # Automatically adjust subplot params

    # Convert 'Category' to a categorical type and assign numerical values
    speed_df['GenotypeNum'] = speed_df['Genotype'].astype('category').cat.codes
    offset_width = 0.1  # Adjust the spacing between points in the same category
    grouped = speed_df.groupby('GenotypeNum')
    offsets = {cat: np.linspace(-offset_width/2, offset_width/2, num=len(group)) for cat, group in grouped}
    speed_df['Offset'] = speed_df.apply(lambda row: offsets[row['GenotypeNum']][grouped.groups[row['GenotypeNum']].get_loc(row.name)], axis=1)

    plt.figure(figsize=(12,8))
    plt.errorbar(speed_df['GenotypeNum']+speed_df['Offset'], speed_df['speed'], yerr=speed_df['speed_std'], fmt='o', capsize=5)
    # Set the x-ticks and their labels to the original string values
    plt.xticks(range(len(speed_df['Genotype'].unique())), speed_df['Genotype'].unique())
    plt.xticks(rotation=45)
    plt.xlabel('Genotype')
    plt.ylabel(r"Speed $(\mu m/min)$")
    plt.tight_layout()  # Automatically adjust subplot params

    plt.figure(figsize=(12,8))
    sns.scatterplot(data=combined_df, x="speed", y="NE", hue="Genotype")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.xlabel(r"Speed $(\mu m/min)$")
    plt.ylabel(r"NE rate (per cell per min)")
    plt.tight_layout()  # Automatically adjust subplot params
    plt.show()


if __name__ == "__main__":
    main()
