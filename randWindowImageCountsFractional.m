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
        r=randi([1, size(image,1)- sideLen]);
        c=randi([1, size(image,2)- sideLen]);
        assert(r+gap < size(image,1));
        assert(c+gap < size(image,2));
        % rare edge cases where image, psm have different sizes. can either
        % ignore this frame or try to shift the window a bit, choosing to
        % shift for now.
        if (r+gap > size(psm,1))
            r = r - abs(size(image,1) - size(psm,1));
        end
        if (c+gap > size(psm,2))
            c = c - abs(size(image,2) - size(psm,2));
        end

        % relpsm is the section of psm (mask) according to the window with
        % origin at (r,c)
        relpsm=psm(r:r+gap, c:c+gap); 
        relim = image(r:r+gap, c:c+gap);
        if sum(relpsm, "all") == sideLen^2 % box of interest entirely within tissue mask
            pct=[pct, sum(relim, 'all')/sideLen^2];
            n=[n, 1];
        elseif sum(relpsm, "all")/(sideLen^2)>=cutoff % box of interest partially within tissue mask
            relpix = relim(relpsm);
            pct=[pct, sum(relpix)/length(relpix)];
            % pct[end] is now the total # of white pixels (within overlap region)
            % divided by the total # pixels (within overlap region)
            % so its the packing fraction of a subset of a window
            %n=[n, length(ind)/(sideLen^2)];
            n=[n, length(relpix)/(sideLen^2)];
        else
            % window does not overlap enough with the mask, skip
        end
    end
end
