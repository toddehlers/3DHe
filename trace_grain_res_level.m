function [BW_fill_smooth,area,perimeter,max_l,min_l,orientation,center]=trace_grain_res_level(image,ps,level,smooth)

% read image with crystal parallel to c
BW = im2bw(image,level/1.2);
% invert image
BW = ~BW;
BW = 1-BW;
BW = (BW == 0);
% fill holes
BW_fill = imfill(BW, 'holes');
if smooth==1
    % smooth outline
    seD = strel('diamond',1);
    BW_fill_smooth = imerode(BW_fill,seD);
    BW_fill_smooth = imerode(BW_fill_smooth,seD);
    BW_fill_smooth = imdilate(BW_fill_smooth,seD);
    BW_fill_smooth = imdilate(BW_fill_smooth,seD);
    
else
    BW_fill_smooth=BW_fill;
end
% find wrongly selected areas
par_area=regionprops(BW_fill_smooth,'Area');
if length(par_area)>1
    [C,idx]=max(struct2array(par_area));
    par_pixels=regionprops(BW_fill_smooth,'PixelIdxList');
    for i=1:length(par_area)
        if i~=idx
            BW_fill_smooth(par_pixels(i).PixelIdxList)=0;
        end
    end
end
% get outline
BW_outline = bwperim(BW_fill_smooth);
% plot original image with estimated outline
I_out=image;
I_out(BW_outline)=255;
imshow(I_out), title('outlined original image');
% properties of image
area=regionprops(BW_fill_smooth,'Area');
area.Area1=area.Area*(ps)^2;
perimeter=regionprops(BW_fill_smooth,'Perimeter');
perimeter.Perimeter1=perimeter.Perimeter*(ps);
max_l=regionprops(BW_fill_smooth,'MajorAxisLength');
max_l.MajorAxisLength1=max_l.MajorAxisLength*(ps);
min_l=regionprops(BW_fill_smooth,'MinorAxisLength');
min_l.MinorAxisLength1=min_l.MinorAxisLength*(ps);
orientation=regionprops(BW_fill_smooth,'Orientation');
center=regionprops(BW_fill_smooth,'Centroid');