clear;close all;
set(0,'DefaultFigureWindowStyle','docked')

folder = "images_for_andrew/";
%att = "0.01";
%att_arr = ["0.0" "0.01" "0.1"];
att_arr = ["0.0"];
numSeeds = 20;
firstSeed = 14;
%fileheader = "last_frame_PSM_images/";
fileheader = "mc_simulation_frames/";
%windowSize = [1/2 2/2 3/2 4/2 5/2 6/2];
windowSize = [1/2 2/2 3/2 4/2];
cm = colormap(parula(length(windowSize)));

% for each of the synthetic datasets, specify cell diameter (roughly) here
cellDiameter = 30;

% mark what kind of boundary we're using
isStrictBoundary = true;
assert(isStrictBoundary == 1);

% adjust filename depending on boundary
filenameModifier = "strict_bd_";

for aa=1:length(att_arr)
    att = att_arr(aa);
    hist_fig_id = 10 + aa;
    for ii=1:length(windowSize)
        windowSize
        % a is the multiplier for window size = a * cell diameter
        a = windowSize(ii);
        phi_arr_all_seeds = [];
    
        for jj=1:numSeeds
            %synthetic_data = fileheader+"last_frame_PSM_sim_att"+att+"_sd"+jj+".tiff";
            %synthetic_boundary = fileheader+"last_frame_PSM_sim_att"+att+"_sd"+jj+"_bd.tif";
            %synthetic_data = "diskPacking"+jj+".tif";
            %synthetic_boundary = "diskPacking"+jj+"_bd.tif";
            synthetic_data = fileheader+"MC_cycle"+(jj+firstSeed)+".tif";
            synthetic_boundary = fileheader+"MC_cycle"+(jj+firstSeed)+"_bd.tif";
            boundary = im2gray(imread(synthetic_boundary));
            B = bwboundaries(~(imbinarize(boundary)),'noholes');
            assert(length(B)==1); % only 1 boundary allowed
            
            % read in image
            I = im2gray(imread(synthetic_data));
    
            % split I into windows
            % window sideLen should be comparable to cell diameter
            cellDiameterPixels = cellDiameter;
            im_arr = windowImage(I, round(a*cellDiameterPixels), B{1});
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
            phi_arr_all_seeds = [phi_arr_all_seeds phi_arr];
        end
        figure(hist_fig_id); hold on;
        [N,edges] = histcounts(phi_arr_all_seeds, 'Normalization','pdf');
        edges = edges(2:end) - (edges(2)-edges(1))/2;
        [sigma,mu] = std(phi_arr_all_seeds);
        plot(edges, N,'DisplayName', "P($\phi |$ a="+a+"d), $\mu$ = "+...
            sprintf('%0.3f',mu)+"$\pm$"+sprintf('%0.3f',sigma), ...
            'linewidth', 3, 'Color', cm(ii,:));
    end
    figure(hist_fig_id); hold on;
    xlabel('$\phi$','interpreter','latex', 'fontsize', 24)
    ylabel('$P(\phi)$','interpreter','latex', 'fontsize', 24)
    %title("att="+att)
    legend('fontsize', 12, 'interpreter','latex', 'location', 'Best')
    xlim([0 1])
    %saveas(gcf, "dist_phi_synthetic_"+filenameModifier+"att_"+att+".png")
    saveas(gcf, "diskPackingDistribution.png")
end

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
                cutoffArea = area(polywindow); %strict
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
