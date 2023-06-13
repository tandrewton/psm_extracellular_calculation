function [pct, n] = randWindowImageCountsFractional(image, sideLen, psm, numSamples)
% Count positive fraction in tiling windows
% Inputs: 
%   image: 2D logical matrix
%   sideLen: Integer of square window side length
%   psm: Mask of relevant tissue, only windows overlapping mask are counted
%   cutoff: Optional argument, minimum window overlap with PSM included
%   numSamples: Number of window locations (with random positions)
% Outputs:
%   pct: Fraction of positive pixels in windows
%   n: Measured window size as fraction of max. Ie if window only half overlaps psm n=.5 
if ~exist('numSamples','var')
    numSamples=1000;
end
    cutoff=1e-10;
    pct=[];
    n=[];
    gap=sideLen-1;
    psm=logical(psm);
    % continue sampling windows until there are numSamples data points
    while (length(pct) < numSamples)
        % pick random windows within the image boundaries
        r=randi([sideLen, size(image,1)- sideLen]);
        c=randi([sideLen, size(image,2)- sideLen]);
        assert(r+gap < size(image,1));
        assert(c+gap < size(image,2));

        % relpsm is the section of psm (mask) according to the window with
        % origin at (r,c)
        relpsm=psm(r:r+gap, c:c+gap); 
        if sum(relpsm, "all") == sideLen^2 % box of interest entirely within tissue mask
            relim = image(r:r+gap, c:c+gap);
            pct=[pct, sum(relim, 'all')/sideLen^2];
            n=[n, 1];
        elseif sum(relpsm, "all")/(sideLen^2)>=cutoff % box of interest partially within tissue mask
            [row, col]=find(relpsm);
            row=row+r-1;
            col=col+c-1; 
            % row, col are now the indices of image that are both within the box and within the tissue mask
            ind=sub2ind(size(image), row, col);
            pct=[pct, sum(image(ind))/length(ind)];
            % pct[end] is now the total # of white pixels (within overlap region)
            % divided by the total # pixels (within overlap region)
            % so its the packing fraction of a subset of a window
            n=[n, length(ind)/(sideLen^2)];
        else
            % window does not overlap enough with the mask, skip
        end
    end
end
