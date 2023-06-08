% trying to use windowImageCountsFractional

%function [pct, n] = windowImageCountsFractional(image, sideLen, psm, cutoff)
% Count positive fraction in tiling windows
% Inputs: 
%   image: 2D logical matrix
%   sideLen: Integer of square window side length
%   psm: Mask of relevant tissue, only windows overlapping mask are counted
%   cutoff: Optional argument, minimum window overlap with PSM included
% Outputs:
%   pct: Fraction of positive pixels in windows
%   n: Measured window size as fraction of max. Ie if window only half overlaps psm n=.5
%%
close all; clear; 
set(0,'DefaultFigureWindowStyle','docked')
fontsize(gcf,28,"points")
folder = "mc_simulation_frames/";
rng(1)

cutoff = 1e-10;
cellDiam = 40; % in pixels
windowLengthPixels = round([0.75 1 1.25]*cellDiam);

% for the original image
figure(1); 
% for the histograms
figure(2); hold on; tiledlayout(1, length(windowLengthPixels));

for ss=windowLengthPixels
    windowPackingFractions = [];
    windowN = [];
    randWindowPackingFractions = [];
    randWindowN = [];
    firstFrame = 1;
    lastFrame = 1;
    for ff=firstFrame:lastFrame
        raw_file = "MC_cycle"+ff+".tif";
        bd_file = "MC_cycle"+ff+"_bd.tif";
        I = imbinarize(im2gray(imread(folder+raw_file)));
        psm_mask = imbinarize(im2gray(imread(folder+bd_file)));
        if (ff == firstFrame)
            figure(1); hold on;
            imshow(I); hold on;
        end
        % find indices of largest area blob in psm_mask
        psm_mask = bwareafilt(psm_mask,1);
    
        [pct, n] = windowImageCountsFractional(I, ss, psm_mask);
        windowPackingFractions = [windowPackingFractions (1-pct)];
        windowN = [windowN n];
        
        [pct, n] = randWindowImageCountsFractional(I, ss, psm_mask, 1000);
        randWindowPackingFractions = [randWindowPackingFractions (1-pct)];
        randWindowN = [randWindowN n];
    end
    binedges = 0:0.05:1;

    % transform to weighted bin counts
    [b,~] = discretize(windowPackingFractions, binedges);
    c = accumarray(b', windowN');
    windowBinCounts=[c;zeros(numel(binedges)-1-numel(c),1)];

    [b,~] = discretize(randWindowPackingFractions, binedges); % get bincounts a for each bin, hist inds b for each element
    c = accumarray(b', randWindowN'); % into inds b, accumulate weights N
    randWindowBinCounts = [c;zeros(numel(binedges)-1-numel(c),1)]; % zero pad to right size

    figure(2);
    nexttile; hold on;

    histogram('binedges', binedges, 'bincounts', windowBinCounts','normalization', 'pdf');
    histogram('binedges', binedges, 'bincounts', randWindowBinCounts','normalization', 'pdf');

    %histogram(windowPackingFractions, binedges, 'normalization', 'pdf')
    %histogram(randWindowPackingFractions, binedges, 'normalization', 'pdf')
    xlim([0 1])
    xlabel("\phi")
    ylabel("P(\phi)")
    set(gca,"Fontsize", 26)
    axis square
    title(string(ss/windowLengthPixels(2)))
end
figure(2); linkaxes
