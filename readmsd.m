% read miriam's msd files
set(0,'DefaultFigureWindowStyle','docked')
set(groot, 'defaultAxesFontSize', 24)

filename = "/Users/AndrewTon/Downloads/MSD_wt.csv";
T = readtable(filename);
%filename = "/Users/AndrewTon/Downloads/MSD_wt3.csv";
%filename = "/Users/AndrewTon/Downloads/cdh_MSD.csv";


%% loop over trackIDs and calculate or plot msd
figure(1); clf; hold on;

uniqueSamples = unique(T.sample);

for sampleNum = 1:length(uniqueSamples)
    % get table rows where sample col matches the unique sample ID
    sampleTable = T(T.sample == string(uniqueSamples{sampleNum}),:);
    uniqueTrackIDs = unique(sampleTable{:,"TrackID"});
    msd = nan(length(uniqueTrackIDs), max(sampleTable.normTime));
    for ii=1:length(uniqueTrackIDs)
        % for each trackID in this sample, record msd and normtime
        normtime = sampleTable(sampleTable.TrackID == uniqueTrackIDs(ii),:).normTime;
        msd_cell = sampleTable(sampleTable.TrackID == uniqueTrackIDs(ii),:).MSD;
        if (length(msd(ii,:)) > length(msd_cell))
            msd(ii,:) = padarray(msd_cell', [0 length(msd(ii,:)) - length(msd_cell)] , nan, "post");
        end
    end
    msd = mean(msd, 'omitnan')/(pi*25);
    msd = msd(~isnan(msd));
    plot(msd,'linewidth',2)
    xlabel('Time')
    ylabel('MSD/a0')
    
    % plot a power law
    % Define the power law model
    powerLawModel = @(A, B, x) A * x.^B;
    
    halfInd = round(length(msd)/2);
    % Perform the fit
    initialGuess = [0, 0]; % Initial guess for coefficients A and B
    fitResult = fit((halfInd:length(msd))', msd(halfInd:end)', powerLawModel, 'StartPoint', initialGuess);
    %plot(fitResult, 'k--');% sprintf('Fit: A=%0.2f, B=%0.2f', fitResult.A, fitResult.B));
    sprintf('Fit: A=%0.2f, B=%0.2f', fitResult.A, fitResult.B)
end
