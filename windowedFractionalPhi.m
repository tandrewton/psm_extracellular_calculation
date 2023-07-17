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
fontsize(gcf,14,"points")
folder = "images_for_andrew/";
imageSets = dir(fullfile(folder, '*itga5--cdh2--_13*'));
%imageSets = dir(fullfile(folder, '*cdh2_*'));

% for the histograms
figure(2); hold on; tiledlayout(1, length(imageSets));

for gg=1:length(imageSets)
    genotypeInfo = imageSets(gg).name;
    imageFile = genotypeInfo;
    
    randWindowPackingFractions = [];
    randWindowN = [];
    
    psm_cells_all = readstack(folder + imageFile, 2, 4);
    psm_raw_all = readstack(folder + imageFile, 3, 4);
    psm_mask_all = readstack(folder + imageFile, 4, 4);
    pixelSeparation = 0.2075665; % microns
    voxelDepthRatio = 1/pixelSeparation; % FIJI metadata
    windowLengthPixels = round(20/pixelSeparation); % 20 microns
    
    %imshow(permute(imresize(squeeze(psm_raw_all(50,:,:)),[510 67*voxelDepthRatio],'nearest'),[2 1 3]));
    handle2DSlices = true;
    if (handle2DSlices)
        firstFrame = 1;
        %lastFrame = firstFrame + 10;
        lastFrame = length(psm_raw_all(:,1,1));
        
        for ff=firstFrame:lastFrame
            ff
            [yDim xDim zDim] = size(psm_raw_all);
            % read in image and rescale for visibility
            psm_raw = permute(imresize(squeeze(psm_raw_all(ff,:,:)),[xDim zDim*voxelDepthRatio],'nearest'),[2 1 3]);
            psm_mask = permute(imresize(squeeze(psm_mask_all(ff,:,:)),[xDim zDim*voxelDepthRatio],'nearest'), [2 1 3]);
    
            I = psm_raw > 1;
            %figure(1); 
            %imshow(I)
            if (mod(ff,5) == 0 && false)
                figure(3+ff);
                imshow(permute(imresize(squeeze(psm_cells_all(ff,:,:)),[xDim zDim*voxelDepthRatio],'nearest'), [2 1 3]))
            end
            [pct, n] = randWindowImageCountsFractional_cuboid(I, windowLengthPixels, psm_mask, 1000);
            randWindowPackingFractions = [randWindowPackingFractions 1-pct'];
            randWindowN = [randWindowN n'];
        end
    else
        % handling the 3D image all at once
        I = psm_raw_all > 1;
        windowLengthPixels = [windowLengthPixels round(windowLengthPixels/voxelDepthRatio)];
        psm_mask = psm_mask_all;

        [pct, n] = randWindowImageCountsFractional_cuboid(I, windowLengthPixels, psm_mask, 100000);
        randWindowPackingFractions = [randWindowPackingFractions 1-pct];
        randWindowN = [randWindowN n];
    end

    % transform to weighted histogram
    binedges = 0:0.02:1;
    % get hist inds b for each bin
    [b,~] = discretize(randWindowPackingFractions, binedges); 
    % accumulate weights N into inds b
    c = accumarray(b', squeeze(randWindowN)); 
    % zero pad
    randWindowBinCounts = [c;zeros(numel(binedges)-1-numel(c),1)]; 
    
    %figure(); hold on;
    figure(2);
    nexttile; hold on;
    title(genotypeInfo(1:end-4),'interpreter', 'none')
    histogram('binedges', binedges, 'bincounts', randWindowBinCounts','normalization', 'pdf');
    xlim([0 1])
    xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ylim([0 10])
    xlabel("\phi")
    ylabel("P(\phi)")
    set(gca,"Fontsize", 14)
    axis square
end
