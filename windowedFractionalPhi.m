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
folder = "images_for_andrew/";
imageFile = "cdh2_10_1.tif";

windowLengthPixels = 60;

randWindowPackingFractions = [];
firstFrame = 21;
lastFrame = 50;
for ff=firstFrame:lastFrame
    psm_raw = imread(folder+imageFile, 3 + (ff-1)*4);
    psm_mask = imread(folder+imageFile, 4 + (ff-1)*4);
    % read in image
    I = psm_raw > 1;
    figure(1); 
    imshow(I)

    [pct, n] = randWindowImageCountsFractional(I, windowLengthPixels, psm_mask, 1000);
    randWindowPackingFractions = [randWindowPackingFractions 1-pct];
end
figure(); hold on;
histogram(randWindowPackingFractions, 'normalization', 'pdf')
xlim([0 1])
ylim([0 5])
xlabel("\phi")
ylabel("P(\phi)")
set(gca,"Fontsize", 26)
axis square
