function imstack=readstack(fname, desiredchannel, totchannel)
info=imfinfo(fname);
totimages=numel(imfinfo(fname));
zslices=totimages/totchannel;
v=info.BitDepth;
rows=info.Height;
cols=info.Width;
imstack=zeros(rows, cols, zslices, strcat('uint', num2str(v)));
series=totchannel*(1:zslices-1)+desiredchannel;
imstack(:, :, 1)=imread(fname, desiredchannel);
j=1;
for k=series
    j=j+1;
    imstack(:, :, j)=imread(fname, k);
end
end