clear;close all;
set(0,'DefaultFigureWindowStyle','docked')

folder = "images_for_andrew/";
att = "0.0";
numSeeds = 10;
fileheader = "last_frame_PSM_images/";
windowSize = [1/2 2/2 3/2 4/2 5/2 6/2];

% for each of the synthetic datasets, specify cell diameter (roughly) here
cellDiameter = 175;

for ii=1:length(windowSize)
    % a is the multiplier for window size = a * cell diameter
    a = windowSize(ii);
    phi_arr_all_seeds = [];

    for jj=1:numSeeds
        synthetic_data = fileheader+"last_frame_PSM_sim_att"+att+"_sd"+jj+".tiff";

        % read in image
        I = im2gray(imread(synthetic_data));

        % split I into windows
        % window sideLen should be comparable to cell diameter
        cellDiameterPixels = cellDiameter;
        im_arr = windowImage(I, a*cellDiameterPixels);
        phi_arr = zeros(size(im_arr));

        for i=1:length(im_arr)
            phi_arr(i) = calcSubPhi(im_arr{i});
        end

        if (jj==Inf)
            % show some configurations, for examples
            figure();
            imshow(im_arr{1+floor(length(im_arr)/10)})
            figure()
            imshow(im_arr{1+floor(length(im_arr)/5)})
            figure()
            imshow(im_arr{1+floor(length(im_arr)/2)})
        end
        phi_arr_all_seeds = [phi_arr_all_seeds phi_arr];
    end
        
    figure();
    histogram(phi_arr_all_seeds,40,'Normalization','pdf')
    xlabel('$\phi$','interpreter','latex', 'fontsize', 24)
    ylabel('$P(\phi)$','interpreter','latex', 'fontsize', 24)
    title("Window size = "+a+"d")
    ylim([0 25])
    xlim([0 1])
    saveas(gcf, "dist_phi_att_"+att+"_window_"+a+".png")
end


function images = windowImage(image, sideLen)
    % break an image into subImages of size sideLen x sideLen
    images = {};
    %start from (0,0) and window along x until greater than image size,
    %then increase y and window along x again.
    origin = [1 1];
    while (origin(1)+sideLen < length(image(1,:)))
        while (origin(2)+sideLen < length(image(:,1)))
            images{end+1} = [image(origin(2):origin(2)+sideLen, ...
                origin(1):origin(1)+sideLen)];
            origin(2) = origin(2) + sideLen;
        end
        origin(2) = 1;
        origin(1) = origin(1) + sideLen;
    end
end

function image_arr = breakNWindows(image, N)
    % divide image into N rows and N columns. 
    % each (row,col) is going to be a subImage
    % leftover data is not formed into any row or column.
    image_arr = {};
    wid = floor(numel(image(:,1))/N);
    len = floor(numel(image(1,:))/N);
    for i=1:N
        for j=1:N
            image_arr{end+1} = image((i-1)*wid+1:i*wid, (j-1)*len+1:j*len);
        end
    end
end

function phi = calcSubPhi(subImage)
    numPixels = numel(subImage);
    %numBrightPixels = sum(subImage>50, 'all');
    %phi = numBrightPixels/numPixels;
    numBlackPixels = sum(subImage<255,'all');
    phi = numBlackPixels/numPixels;
end

% function calcSubPhi
%
% calcSubPhi takes in a window of pixels and computes the number of bright
%   pixels. It returns the number of bright pixels divided by the total
%   number of pixels. 
%
% assumption: dark pixels and bright pixels are defined by a threshold.
%       we are not worrying about pixel bleeding.
%
% input: mxn array of pixel intensities
% output: extracellular fraction phi