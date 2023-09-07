%% for DP simulations
close all; clear; 
set(0,'DefaultFigureWindowStyle','docked')
fontsize(gcf,28,"points")
%folder = "DP_simulation_frames/";
if ismac
    folder = "/Users/AndrewTon/Documents/YalePhD/projects/dpm/output/cells/psm/";
elseif ispc
    folder = "C:/Users/atata/projects/dpm/output/cells/psm/";
else
    disp('Platform not supported')
end

rng(1)
k_ecm = "0.005";
k_off = "1.0";
v0_arr = ["0.02" "0.04" "0.08"];
att = "0.001";

cellDiam = 50; % in pixels
windowLengthPixels = round(2*cellDiam);

% for the original image
figure(1); 
% for the histograms
figure(2); hold on; tiledlayout(1, length(windowLengthPixels));

for v0=v0_arr
    for ss=windowLengthPixels
        windowPackingFractions = [];
        windowN = [];
        randWindowPackingFractions = [];
        randWindowN = [];
        firstFrame = 10; % make sure to adjust these as necessary
        lastFrame = 25;
        for ff=firstFrame:lastFrame
            %filename = "psmtest/psmpsmtest8_seed1fr"+ff;
            filename = "psm_calA01.0_phi0.74_tm10.0_v0"+v0+"_t_abp1.0k_ecm"+k_ecm+"k_off"+k_off+"/"...
                + "psm_N40_dur1000_att"+att+"_sd1fr"+ff;
            raw_file = filename+".tif";
            bd_file = filename+"_bd.tif";
            I_color = imread(folder+raw_file);
            I = imbinarize(im2gray(I_color), 0.8);
            %I = imdilate(I, strel("disk",1)); % let edges bleed by 1 pixel
            psm_mask = ~imbinarize(im2gray(imread(folder+bd_file)));
            if (ff == lastFrame)
                figure(1); hold on;
                imshow(I_color);
                figure(2); hold on;
                imshow(I);
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
    
        figure(3);
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
        %title(string(windowLengthPixels/cellDiam))
    end
    figure(3); linkaxes
    ylim([0 21.1])
end