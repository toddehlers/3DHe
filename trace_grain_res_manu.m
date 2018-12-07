function [BW,area,perimeter,max_l,min_l,orientation,center]=trace_grain_res_manu(image,ps)

% draw outline
BW = roipoly(image);
% plot original image with estimated outline
imshow(I), title('outlined original image');
% properties of image
area=regionprops(BW,'Area');
area.Area1=area.Area*(ps)^2;
perimeter=regionprops(BW,'Perimeter');
perimeter.Perimeter1=perimeter.Perimeter*ps;
max_l=regionprops(BW,'MajorAxisLength');
max_l.MajorAxisLength1=max_l.MajorAxisLength*ps;
min_l=regionprops(BW,'MinorAxisLength');
min_l.MinorAxisLength1=min_l.MinorAxisLength*ps;
orientation=regionprops(BW,'Orientation');
center=regionprops(BW,'Centroid');