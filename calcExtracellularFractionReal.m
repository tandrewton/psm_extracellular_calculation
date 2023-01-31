clear;close all;
set(0,'DefaultFigureWindowStyle','docked')

folder = "images_for_andrew/";
wt_file = "wt/wt_10_2_crop-slice.tif";
cdh_file = "cdh2/cdh2_11_1-1-slice.tif";
fbn_mo_file = "fbn2b_MO/fbn2bMO_15_1-slice.tif";
fbn1_fbn2_file = "fn1a--fn1b--/fn1a--fn1b--_12_1-slice.tif";
fn_cdh_MO_file1 = "fn1a--fn1b--cdh2--fbn2b_MO/fbn2bMO_fn1a--fn1b--cdh2--_10_2-slice.tif";
fn_cdh_MO_file2 = "fn1a--fn1b--cdh2--fbn2b_MO/fbn2bMO_fn1a--fn1b--cdh2--_11_2-slice.tif";

bd_file = "wt_11_1_psmarea.tif";
raw_file = "wt_11_1_rawbinary-1.tif";

numFrames = 1;
phi_arr_all_frames = [];

%windowSize = [1/2 2/2 3/2 4/2 5/2 6/2];
windowSize = [1/2];
cm = colormap(parula(length(windowSize)));

% for each of the synthetic datasets, specify cell diameter (roughly) here
cellDiameter = 250;

for ii=1:length(windowSize)
    hist_fig_id = 10;
    % a is the multiplier for window size = a * cell diameter
    a = windowSize(ii);
    phi_arr_all_frames = [];

    for jj=1:numFrames
        psm_raw = imread(folder+raw_file,jj);
        psm_boundary = imread(folder+bd_file,jj);
        imshow(psm_raw)

        B = bwboundaries(psm_boundary,'noholes');
        B = B{1};
        
        % read in image
        I = psm_raw;
        %need to invert grayscale values (0->255, 255->0)

        % split I into windows
        % window sideLen should be comparable to cell diameter
        cellDiameterPixels = cellDiameter;
        im_arr = windowImage(I, round(a*cellDiameterPixels), B);
        %im_arr = windowImage(I, round(a*cellDiameterPixels));
        phi_arr = zeros(size(im_arr));

        for i=1:length(im_arr)
            phi_arr(i) = calcSubPhi(im_arr{i});
        end

        if (jj==1)
            % show some configuratitoons, for examples
            figure();
            imshow(im_arr{1+floor(length(im_arr)/10)})
            figure()
            imshow(im_arr{1+floor(length(im_arr)/5)})
            figure()
            imshow(im_arr{1+floor(length(im_arr)/2)})
        end
        phi_arr_all_frames = [phi_arr_all_frames phi_arr];
    end
    figure(hist_fig_id); hold on;
    [N,edges] = histcounts(phi_arr_all_frames, 'Normalization','pdf');
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    plot(edges, N,'DisplayName', "P($\phi |$ a="+a+"d)", 'linewidth', 3, 'Color', cm(ii,:));
    %histogram(phi_arr_all_seeds,40,'Normalization','pdf')
end
figure(hist_fig_id); hold on;
xlabel('$\phi$','interpreter','latex', 'fontsize', 24)
ylabel('$P(\phi)$','interpreter','latex', 'fontsize', 24)
title("all window sizes")
legend('fontsize', 12, 'interpreter','latex', 'location', 'Best')
saveas(gcf, "dist_phi_psm.png")

% function images = windowImage(image, sideLen)
%     % break an image into subImages of size sideLen x sideLen
%     images = {};
%     %start from (0,0) and window along x until greater than image size,
%     %then increase y and window along x again.
%     origin = [1 1];
%     while (origin(1)+sideLen < length(image(1,:)))
%         while (origin(2)+sideLen < length(image(:,1)))
%             images{end+1} = [image(origin(2):origin(2)+sideLen, ...
%                 origin(1):origin(1)+sideLen)];
%             origin(2) = origin(2) + sideLen;
%         end
%         origin(2) = 1;
%         origin(1) = origin(1) + sideLen;
%     end
% end

function images = windowImage(image, sideLen, boundary)
    % break an image into subImages of size sideLen x sideLen, with
    % optional boundary argument for reducing # windows based on
    % intersection with boundary
    images = {};
    %start from (0,0) and window along x until greater than image size,
    %then increase y and window along x again.
    origin = [1 1];
    while (origin(1)+sideLen < length(image(1,:)))
        while (origin(2)+sideLen < length(image(:,1)))
            if (exist('boundary', 'var'))
                %fprintf("origin = [%d %d], size(images) = %d\n", origin(1), origin(2), length(images));
                % boundary option: discard subimages overlapping w/ boundary
                % construct boundary polyshape and window polyshape
                polybd = polyshape(boundary);
                xpos = [origin(1); origin(1); origin(1)+sideLen; origin(1)+sideLen];
                ypos = [origin(2); origin(2)+sideLen; origin(2)+sideLen; origin(2)];
                polywindow = polyshape(xpos, ypos);
                overlapArea = area(intersect(polywindow, polybd));
                %if (overlapArea < area(polywindow)) % strict
                if (overlapArea < 1e-10) % lenient
                    % window overlap is too small, discard window
                    origin(2) = origin(2) + sideLen;
                    continue;
                else
                    images{end+1} = [image(origin(2):origin(2)+sideLen, origin(1):origin(1)+sideLen)];
                end
            else
                images{end+1} = [image(origin(2):origin(2)+sideLen, origin(1):origin(1)+sideLen)];
            end
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