%% for DP simulations
close all; clear; 
set(0,'DefaultFigureWindowStyle','docked')
fontsize(gcf,28,"points")
folder = "DP_simulation_frames/";
rng(1)

cellDiam = 75; % in pixels
windowLengthPixels = round(2*cellDiam);

% for the original image
figure(1); 
% for the histograms
figure(2); hold on; tiledlayout(1, length(windowLengthPixels));

for ss=windowLengthPixels
    windowPackingFractions = [];
    windowN = [];
    randWindowPackingFractions = [];
    randWindowN = [];
    firstFrame = 20;
    lastFrame = 24;
    for ff=firstFrame:lastFrame
        filename = "testdata8fr"+ff;
        raw_file = filename+".tif";
        bd_file = filename+"_bd.tif";
        I = imbinarize(im2gray(imread(folder+raw_file)));
        psm_mask = ~imbinarize(im2gray(imread(folder+bd_file)));
        if (ff == firstFrame)
            figure(1); hold on;
            imshow(I); hold on;
        end
        % find indices of largest area blob in psm_mask
        psm_mask = bwareafilt(psm_mask,1);
    
%         [pct, n] = windowImageCountsFractional(I, ss, psm_mask);
%         windowPackingFractions = [windowPackingFractions (1-pct)];
%         windowN = [windowN n];
        
        [pct, n] = randWindowImageCountsFractional(I, ss, psm_mask, 1000);
        randWindowPackingFractions = [randWindowPackingFractions (1-pct)];
        randWindowN = [randWindowN n];
    end
    binedges = 0:0.02:1;

%     % transform to weighted bin counts
%     [b,~] = discretize(windowPackingFractions, binedges);
%     c = accumarray(b', windowN');
%     windowBinCounts=[c;zeros(numel(binedges)-1-numel(c),1)];

    [b,~] = discretize(randWindowPackingFractions, binedges); % get bincounts a for each bin, hist inds b for each element
    c = accumarray(b', randWindowN'); % into inds b, accumulate weights N
    randWindowBinCounts = [c;zeros(numel(binedges)-1-numel(c),1)]; % zero pad to right size

    figure(2);
    nexttile; hold on;

    %histogram('binedges', binedges, 'bincounts', windowBinCounts','normalization', 'pdf');
    histogram('binedges', binedges, 'bincounts', randWindowBinCounts','normalization', 'pdf');

    %histogram(windowPackingFractions, binedges, 'normalization', 'pdf')
    %histogram(randWindowPackingFractions, binedges, 'normalization', 'pdf')
    xlim([0 1])
    xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    xlabel("\phi")
    ylabel("P(\phi)")
    set(gca,"Fontsize", 26)
    axis square
    title(string(windowLengthPixels/cellDiam))
end
figure(2); linkaxes
ylim([0 20])
