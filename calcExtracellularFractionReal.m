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

numFrames = 10;
phi_arr_all_frames = [];

% mark what kind of boundary we're using
isNoBoundary = false;
isLenientBoundary = false;
isStrictBoundary = true;
assert(isNoBoundary+isLenientBoundary+isStrictBoundary == 1);

% adjust filename depending on boundary
if (isNoBoundary)
    filenameModifier = "no_bd_";
elseif (isLenientBoundary)
    filenameModifier = "lenient_bd_";
elseif (isStrictBoundary)
    filenameModifier = "strict_bd_";
end

windowSize = [1/2 2/2 3/2 4/2 5/2 6/2];
%windowSize = [2/2];
cm = colormap(parula(length(windowSize)));

% for each dataset, specify cell diameter in pixels (roughly) here
cellDiameter = 50;
for startFrame=[1,11,21,41]
    for ii=1:length(windowSize)
        hist_fig_id = 999+startFrame+numFrames;
        % a is the multiplier for window size = a * cell diameter
        a = windowSize(ii);
        phi_arr_all_frames = [];
    
        for jj=startFrame:startFrame+numFrames-1
            psm_raw = imread(folder+raw_file,jj);
            psm_boundary = imread(folder+bd_file,jj);
    
            if (ii == 1)
                figure()
                imshow(psm_raw)
            end
    
            B = bwboundaries(psm_boundary,'noholes');
            B = B{1};
            
            % read in image
            I = psm_raw;
            %need to invert grayscale values (0->255, 255->0)
    
            % split I into windows
            % window sideLen should be comparable to cell diameter
            cellDiameterPixels = cellDiameter;
            if (isNoBoundary)
                im_arr = windowImage(I, round(a*cellDiameterPixels));
            elseif (isLenientBoundary)
                im_arr = windowImage(I, round(a*cellDiameterPixels), B, isStrictBoundary);
            elseif (isStrictBoundary)
                im_arr = windowImage(I, round(a*cellDiameterPixels), B, isStrictBoundary);
            end
            phi_arr = zeros(size(im_arr));
    
            for i=1:length(im_arr)
                phi_arr(i) = calcSubPhi(im_arr{i});
            end
    
            if (false)
                % show some configurations, for examples
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
        [sigma,mu] = std(phi_arr_all_frames);
        plot(edges, N,'DisplayName', "P($\phi |$ a="+a+"d), $\mu$ = "+...
            sprintf('%0.3f',mu)+"$\pm$"+sprintf('%0.3f',sigma), ...
            'linewidth', 3, 'Color', cm(ii,:));
        %histogram(phi_arr_all_seeds,40,'Normalization','pdf')
        trueMean = mean(phi_arr_all_frames);
        fprintf("true mu = %f\n", trueMean);
    end
    figure(hist_fig_id); hold on;
    xlabel('$\phi$','interpreter','latex', 'fontsize', 24)
    ylabel('$P(\phi)$','interpreter','latex', 'fontsize', 24)
    title("frames "+startFrame+"-"+(startFrame+numFrames-1))
    legend('fontsize', 12, 'interpreter','latex', 'location', 'Best')
    xlim([0 1])
    saveas(gcf, "dist_phi_psm_"+filenameModifier+startFrame+"-"+(startFrame+numFrames-1)+".png")
end

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

function images = windowImage(image, sideLen, boundary, isStrict)
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
                if (isStrict)
                    cutoffArea = area(polywindow); %strict
                else
                    cutoffArea = 1e-10; %lenient
                end
                if (overlapArea < cutoffArea)
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

function phi = calcSubPhi(subImage)
    numPixels = numel(subImage);
    numBlackPixels = sum(subImage<255,'all');
    phi = numBlackPixels/numPixels;
end