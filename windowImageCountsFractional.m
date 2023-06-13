function [pct, n] = windowImageCountsFractional(image, sideLen, psm)
% Count positive fraction in tiling windows
% Inputs: 
%   image: 2D logical matrix
%   sideLen: Integer of square window side length
%   psm: Mask of relevant tissue, only windows overlapping mask are counted
%   cutoff: Optional argument, minimum window overlap with PSM included
% Outputs:
%   pct: Fraction of positive pixels in windows
%   n: Measured window size as fraction of max. Ie if window only half overlaps psm n=.5 
    cutoff=1e-10;
    assert(cutoff > 0) % if cutoff = 0, could get 0/0 nans
    pct=[];
    n=[];
    %start from (0,0) and window along x until greater than image size,
    %then increase y and window along x again.
    r=1;
    c=1;
    gap=sideLen-1;
    psm=logical(psm);
    while (r+gap < size(image, 1))
        while (c+gap < size(image, 2))
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
                end
            c = c + sideLen;
        end
        c = 1;
        r = r + sideLen;
    end
end
