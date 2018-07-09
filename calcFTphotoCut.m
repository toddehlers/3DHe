function [Ft_num,mass_num,matrix_3D]=calcFTphotoCut(mineral,file1,file2,R,res,r232_238,r147_238)
% script to calculate Ft-values and mass of (intact, broken) apatite or zircon 
% grains from microscope pictures taken parallel and perpendicular to the 
% crystallographic c-axis (Jan/2017)
% please contact Christoph Glotzbach (christoph.glotzbach@uni-tuebingen.de) 
% if you have any problems, or you want to publish results

% required input
%mineral='ap'; % either ap for apatite or zr for zircon
%file1='FT01_parC.tif'; % name of file that shows the crystal parallel to C
%R1=0.4; % resolution of image in microns/pixel
%file2='FT01_perC.tif'; % name of file that shows the crystal perpendicular to C
%R2=0.4; % resolution of image in microns/pixel
%res=0.5; % reduce resolution, value between 1 and 0
%r232_238=0; % ratio of 232Th/238U atoms in grain; 
                % 0 if not yet known, in this case a mean of measured ratios is used
%r147_238=0;     % ratio of 147Sm/238U atoms in grain;
                % 0 if not yet known, in this case a mean of measured ratios is used


% do not change from here on, 

% read image with crystal parallel to c
smooth=1; % 1 for smoothing, 0 for original image
disp(' ')
disp('Find outline of grain')
[image1,area1,perimeter1,max_l1,min_l1,orientation1,level1,center]=trace_grain_res(file1,R,res,smooth);
str='1';
while str=='1' | str=='2'
    prompt = 'press (return) if OK, press (1) or (2) to repeat tracing with lower or higher treshold or do manual tracing (3): ';
    str = input(prompt,'s');
    switch str
        case '1'
            level1=level1*0.75;
            [image1,area1,perimeter1,max_l1,min_l1,orientation1,center]=trace_grain_res_level(file1,R,res,level1,smooth);
        case '2'
            level1=level1/0.75;
            [image1,area1,perimeter1,max_l1,min_l1,orientation1,center]=trace_grain_res_level(file1,R,res,level1,smooth);
        case '3'
            [image1,area1,perimeter1,max_l1,min_l1,orientation1,center]=trace_grain_res_manu(file1,R,res);
    end
end
% show automatically defined c-axis and change if neccessary, end by
% double-click
slope=tand(orientation1.Orientation);
intercept=center.Centroid(1)-slope*(size(image1,1)-center.Centroid(2));
x1=center.Centroid(1)-100;
x2=center.Centroid(1)+100;
y1=size(image1,1)-(slope*x1+intercept);
y2=size(image1,1)-(slope*x2+intercept);
if y1>size(image1,1) | y2>size(image1,1) | y1<0 | y2<0
    y1=200;
    y2=200;
end
h = imline(gca, [x1 x2],[y1 y2]);
setColor(h,[0 1 0]);
id = addNewPositionCallback(h,@(pos) title(mat2str(pos,3)));
% After observing the callback behavior, remove the callback.
% using the removeNewPositionCallback API function.
removeNewPositionCallback(h,id);
position = wait(h);

% recalculate the c-axis
disp('automatic orientation:')
orientation1.Orientation
orientation1.Orientation=-atand((position(1,2)-position(2,2))/(position(1,1)-position(2,1)));
disp('manually corrected orientation:')
orientation1.Orientation

disp(' ')



% read image with crystal perpendicular to c
disp(' ')
disp('Find outline of grain')
[image2,area2,perimeter2,max_l2,min_l2,orientation2,level2]=trace_grain_res(file2,R,res,smooth);
str='1';
while str=='1' | str=='2'
    prompt = 'press (return) if OK, press (1) or (2) to repeat tracing with lower or higher treshold or do manual tracing (3): ';
    str = input(prompt,'s');
    switch str
        case '1'
            level2=level2*0.75;
            [image2,area2,perimeter2,max_l2,min_l2,orientation2]=trace_grain_res_level(file2,R,res,level2,smooth);
        case '2'
            level2=level2/0.75;
            [image2,area2,perimeter2,max_l2,min_l2,orientation2]=trace_grain_res_level(file2,R,res,level2,smooth);
        case '3'
            [image2,area2,perimeter2,max_l2,min_l2,orientation2]=trace_grain_res_manu(file2,R,res);
    end
end
disp(' ')

% rotate image1 by orientation1 and determine starting point and length
% perpendicular to the c-axis
image1_rot=imrotate(image1,-orientation1.Orientation);
pixels=regionprops(image1_rot,'PixelList');
image1_rot_first_x=min(pixels.PixelList(:,1));
image1_rot_first_y=min(pixels.PixelList(:,2));
image1_rot_last_y=max(pixels.PixelList(:,2));
min_l1.MinorAxisLength=image1_rot_last_y-image1_rot_first_y;

% find orientation of parallel image(s) in relation to perpendicular image
prec=2;
for i=prec:prec:180
    BW_rot = imrotate(image2,i);
    pixels=regionprops(BW_rot,'PixelList');
    dist_x(i/2)=max(pixels.PixelList(:,1))-min(pixels.PixelList(:,1));
end
[dist,angle_par2b]=min(abs(dist_x-min_l1.MinorAxisLength));


% make 3D model of two pictures

% rotate image2 by angle_par2b and determine width
image2_rot = imrotate(image2,angle_par2b*prec);
image2_rot_width=max(find(sum(image2_rot)>1*(R/res)))-min(find(sum(image2_rot)>1*(R/res)));
% make 3D model
max_x=round(max_l1.MajorAxisLength)+400;
max_y=round(min_l1.MinorAxisLength)+100;
max_z=round(max_l2.MajorAxisLength)+200;
matrix_3D=zeros(max_x,max_y,max_z);
% go along the x-axis of image1_rot and add ones to 3D matrix with shape of
% image2_rot
area=0;
for i=1:size(image1_rot,2)
    if sum(image1_rot(:,i))>3
        i_first(i)=find(image1_rot(:,i),1,'first');
        i_last(i)=find(image1_rot(:,i),1,'last');
        size_prism(i)=i_last(i)-i_first(i);
        % resize image2_rot
        resize_factor=(i_last(i)-i_first(i))/image2_rot_width;
        image2_rot_resize=imresize(image2_rot,resize_factor,'nearest');
        % get area
        area=area+perimeter1.Perimeter1*resize_factor*(R/res);
        % cut image
        pixels=regionprops(image2_rot_resize,'PixelList');
        if size(pixels,1)>1
            for j=1:size(pixels,1)
                size_pixels(j)=size(pixels(j).PixelList,1);
            end
            [C,idx]=max(size_pixels);
            for j=1:size(pixels,1)
                if j~=idx
                    image2_rot_resize(pixels(j).PixelList(:,2),pixels(j).PixelList(:,1))=0;
                end
            end
            clear pixels
            pixels=regionprops(image2_rot_resize,'PixelList');
        end
        minx=min(pixels.PixelList(:,1));
        miny=min(pixels.PixelList(:,2));
        lengthx=max(pixels.PixelList(:,1))-minx;
        lengthy=max(pixels.PixelList(:,2))-miny;
        image2_rot_resize_crop=imcrop(image2_rot_resize,[minx,miny,lengthx,lengthy]);
        % add to 3D matrix
        x=i-image1_rot_first_x+200;
        y=50+i_first(i)-image1_rot_first_y:size(image2_rot_resize_crop,2)+50+i_first(i)-image1_rot_first_y-1;
        z=round(size(matrix_3D,3)/2)-round(size(image2_rot_resize_crop,1)/2):size(image2_rot_resize_crop,1)+round(size(matrix_3D,3)/2)-round(size(image2_rot_resize_crop,1)/2)-1;
        %y=50+i_first(i)-image1_rot_first_y:size(image2_rot_resize_crop,1)+50+i_first(i)-image1_rot_first_y-1;
        %z=round(size(matrix_3D,3)/2)-round(size(image2_rot_resize_crop,2)/2):size(image2_rot_resize_crop,2)+round(size(matrix_3D,3)/2)-round(size(image2_rot_resize_crop,2)/2)-1;
        matrix_3D(x,y,z)=image2_rot_resize_crop';
    end
end

% add missing pieces in the c-axis perpendicular picture
disp(' ')
disp('Add pieces grinded or polished away?')
str='0';
prompt = 'Press (1) to manually add missing pieces, press (0) if not applicable: ';
str = input(prompt,'s');
if str=='1'
    matrix_3D_broken=zeros(max_x,max_y,max_z);
    BW=roipoly(image2_rot);
    image2_rot_broken=image2_rot;
    image2_rot_broken(BW==1)=1;
    pixels=regionprops(image2_rot_broken,'PixelList');
    image2_rot_broken_width=max(pixels.PixelList(:,1))-min(pixels.PixelList(:,1));
    for i=1:size(image1_rot,2)
        if sum(image1_rot(:,i))>3
            i_first(i)=find(image1_rot(:,i),1,'first');
            i_last(i)=find(image1_rot(:,i),1,'last');
            size_prism(i)=i_last(i)-i_first(i);
            % resize image2_rot
            resize_factor=(i_last(i)-i_first(i))/image2_rot_broken_width;
            image2_rot_broken_resize=imresize(image2_rot_broken,resize_factor);
            % get area
            area=area+perimeter1.Perimeter1*resize_factor*(R/res);
            % cut image
            pixels=regionprops(image2_rot_broken_resize,'PixelList');
            if size(pixels,1)>1
                for j=1:size(pixels,1)
                    size_pixels(j)=size(pixels(j).PixelList,1);
                end
                [C,idx]=max(size_pixels);
                for j=1:size(pixels,1)
                    if j~=idx
                        image2_rot_resize(pixels(j).PixelList(:,2),pixels(j).PixelList(:,1))=0;
                    end
                end
                clear pixels
                pixels=regionprops(image2_rot_broken_resize,'PixelList');
            end
            minx=min(pixels.PixelList(:,1));
            miny=min(pixels.PixelList(:,2));
            lengthx=max(pixels.PixelList(:,1))-minx;
            lengthy=max(pixels.PixelList(:,2))-miny;
            image2_rot_broken_resize_crop=imcrop(image2_rot_broken_resize,[minx,miny,lengthx,lengthy]);
            % add to 3D matrix
            x=i-image1_rot_first_x+200;
            y=50+i_first(i)-image1_rot_first_y:size(image2_rot_broken_resize_crop,2)+50+i_first(i)-image1_rot_first_y-1;
            z=round(size(matrix_3D,3)/2)-round(size(image2_rot_broken_resize_crop,1)/2):size(image2_rot_broken_resize_crop,1)+round(size(matrix_3D,3)/2)-round(size(image2_rot_broken_resize_crop,1)/2)-1;
            %y=50+i_first(i)-image1_rot_first_y:size(image2_rot_resize_crop,1)+50+i_first(i)-image1_rot_first_y-1;
            %z=round(size(matrix_3D,3)/2)-round(size(image2_rot_resize_crop,2)/2):size(image2_rot_resize_crop,2)+round(size(matrix_3D,3)/2)-round(size(image2_rot_resize_crop,2)/2)-1;
            matrix_3D_broken(x,y,z)=image2_rot_broken_resize_crop';
        end
    end
end
disp(' ')

% smooth matrix 
matrix_3D_smooth=imdilate(matrix_3D,strel3d(6));
matrix_3D_smooth=imerode(matrix_3D_smooth,strel3d(6));

if exist('matrix_3D_broken')
    matrix_3D_broken_smooth=imdilate(matrix_3D_broken,strel3d(6));
    matrix_3D_smooth=imerode(matrix_3D_broken_smooth,strel3d(6));
end

% get surface area and volume
volume=sum(matrix_3D_smooth(:))*(R/res)^3;
% a = fv.vertices(fv.faces(:, 2), :) - fv.vertices(fv.faces(:, 1), :);
% b = fv.vertices(fv.faces(:, 3), :) - fv.vertices(fv.faces(:, 1), :);
% c = cross(a, b, 2);
% area = 1/2 * sum(sqrt(sum(c.^2, 2)));


% calculate the Ft based on mean stopping distances from Ketcham et al.
% 2011

% give ratio of U238/Th232 and U238/Sm147 of grain, if not available
% average values are assumed (can be used to estimate FT for grain
% selection)

% first take average 
switch mineral
    case 'ap'
        mass_num=volume*1E-4^3*3.2*1E6
        sd_238U=18.81;
        sd_235U=21.8;
        sd_232Th=22.25;
        sd_147Sm=5.93;
        if r232_238==0 
            r232_238=0.4444;
            r147_238=0.0615;
        end
    case 'zr'
        mass_num=volume*1E-4^3*4.65*1E6
        sd_238U=15.55;
        sd_235U=18.05;
        sd_232Th=18.43;
        sd_147Sm=4.76;
        if r232_238==0 
            r232_238=0.271;
            r147_238=0;
        end
end

% get relative contribution of U,Th,Sm on alpha-decay
lambda238=1.55125E-10;
lambda235=9.84500E-10;
lambda232=4.94750E-11;
lambda147=6.53900E-12;

U238=1;
U235=U238/136.1384*7*lambda235;
Th232=U238*r232_238*6*lambda232;
Sm147=U238*r147_238*1*lambda147;
U238=U238*8*lambda238;
radio_nukl_sum=sum([U238 U235 Th232 Sm147]);
U238=U238/radio_nukl_sum;
U235=U235/radio_nukl_sum;
Th232=Th232/radio_nukl_sum;
Sm147=Sm147/radio_nukl_sum;

% select randomly ne pixel and distort He particles according to the mean 
% stopping distance and random phi and theta values
ne=200000; % number of random events
ne_nuklide=[round(U238*ne) round(U235*ne) round(Th232*ne) round(Sm147*ne)];
sd_nuklide=[sd_238U sd_235U sd_232Th sd_147Sm];
if str=='0'
    pixels=find(matrix_3D_smooth==1);
    pixels_n=randsample(pixels,ne,'true');
    [xi,yi,zi]=ind2sub(size(matrix_3D_smooth),pixels_n);
    matrix_He=zeros(size(matrix_3D_smooth));
else
    pixels=find(matrix_3D_broken_smooth==1);
    pixels_n=randsample(pixels,ne,'true');
    [xi,yi,zi]=ind2sub(size(matrix_3D_broken_smooth),pixels_n);
    matrix_He=zeros(size(matrix_3D_broken_smooth));
end
pixels_n=randsample(pixels,ne,'true');
theta=randsample(360,ne,'true');
phi=randsample(180,ne,'true');
n=0;
hold on
for j=1:4
    for i=1:ne_nuklide(j)
        x = xi(i) + round((sd_nuklide(j)/(R/res) * sind(phi(i)) * cosd(theta(i))));
        y = yi(i) + round((sd_nuklide(j)/(R/res) * sind(phi(i)) * sind(theta(i))));
        z = zi(i) + round((sd_nuklide(j)/(R/res) * cosd(phi(i))));
        if x>0 & y>0 & z>0 & x<=size(matrix_He,1) & y<=size(matrix_He,2) & z<=size(matrix_He,3)
            matrix_He(x,y,z)=matrix_He(x,y,z)+1;
        end
%         % visualize decay
%         quiver3(yi(i),xi(i),zi(i),y-yi(i),x-xi(i),z-zi(i))
%         M(i)=getframe(gcf);
    end
end

% % rotate figure and capture
% for i=201:242
%     camorbit(0.9,-0.9);
%     M(i)=getframe(gcf);
% end
% for i=243:300
%     camorbit(0,-1);
%     M(i)=getframe(gcf);
% end
% 
% movie2avi(M,'AlphaDecayMovie.avi')


% make slice through 3D matrix perpendicular to c-axis
% figure
% slice_c_per=150;
% xtick=(R/res):(R/res):size(matrix_3D,1)*(R/res);
% ytick=(R/res):(R/res):size(matrix_3D,2)*(R/res);
% imagesc(xtick,ytick,squeeze(matrix_3D(1:size(matrix_3D,1),1:size(matrix_3D,2),slice_c_per))+squeeze(matrix_He(1:size(matrix_3D,1),1:size(matrix_3D,2),slice_c_per)))
% figure
% contour(ytick',xtick',squeeze(matrix_3D(1:size(matrix_3D,1),1:size(matrix_3D,2),slice_c_per)),[.5,.5],'k')
% % make slice through 3D matrix perpendicular to c-axis
% slice_c_par=100;
% ytick=(R1/res):(R1/res):size(matrix_3D,2)*(R1/res);
% ztick=(R1/res):(R1/res):size(matrix_3D,3)*(R1/res);
% imagesc(ytick,ztick,squeeze(matrix_3D(slice_c_par,1:size(matrix_3D,2),1:size(matrix_3D,3)))+squeeze(matrix_He(slice_c_par,1:size(matrix_3D,2),1:size(matrix_3D,3))))
% he_map=squeeze(matrix_He(slice_c_par,1:size(matrix_3D,2),1:size(matrix_3D,3)));
% grain_map=squeeze(matrix_3D(slice_c_par,1:size(matrix_3D,2),1:size(matrix_3D,3)));

% determine Ft, how much He-particle retained in the grain
if str=='0'
    Ft_num=sum(matrix_He(matrix_3D_smooth==1))/ne
else
    Ft_num=sum(matrix_He(matrix_3D_smooth==1))/sum(matrix_3D_smooth(pixels_n)==1)
end

