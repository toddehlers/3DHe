% script to load files from folder to calculate FTs from all possible
% combinations of images taken parallel and perpendicular to the
% crystallographic c-axis
% written by C. Glotzbach, 04/2017

%clear
%clc

% read input parameter table (IPT) which should have 14 columns and X number
% of rows depending on the number of pictures
% Table Columns are...
% 1. Path: Path to the file of picture taken from the crystal
% 2. Pixel_size: original resolution of image in microns per pixel
% 3. Mode: 'normal - 0' mode for grains or broken grains (example: Test-1),
%           'cut - 1' mode for grains grained and polished (example: Test-2)
% 4. Mineral: 'ap - 0' apatite grain
%              'zr - 1' zircon grain
% 5. Ratio_232_238: ratio of 232Th and 238U, measured in mol; if 0, mean values are used
% 6. Ratio_147_238: ratio of 147Sm and 238U, measured in mol; if 0 mean values are used
% 7. Orientation2C: defines the orientation of the crystal c-axis to the field of view
%         parallel - 0, perpendicular - 1
% 8. Shape: hexagonal - 0, ellipsoid - 1, cylinder - 2, block - 3 (analytical equations from Ketcham et al. 2011)
% 9. Nr_pyramids: number of pyramids, required for analytical calculation
% 10. Broken_tips: true or false
% 11. Zoned: are radio nuclides inhomogenously distributed, true or false
% 12. Rim_width: width of zoned rim
% 13. Ratio_Rim_Core: radio nuclide ratio between rim and core (1:
%     homogenous distribution, <1: rim depleted, >1: rim concentrated)

input_file = fullfile(grain_folder, input_file);
output_file = fullfile(grain_folder, output_file);
%
%
% display(input_file)
% display(output_file)
% display(grain_folder)
%
% exit

IPT = readtable(input_file);

% open image files and resize to 2 microns (ps=pixel size)
per = 0; par = 0;
par_I = []; per_I = [];
ps = 2;


for i = 1:size(IPT, 1)
    % create image from txt-file
    outline=readtable(char(IPT.coordinateFile(i)));
    X=[outline.Var1(1:2:end); flipud(outline.Var1(2:2:end))];
    Y=[outline.Var2(1:2:end); flipud(outline.Var2(2:2:end))];
    BW=poly2mask(X,Y,max(Y)+50,max(X)+50);
    switch char(IPT.orientation(i))
        case 'parallel'
            par = par + 1;
            image_name = ['image_par' num2str(par)];
            images.(image_name)=BW;
            par_I = [par_I i];
        case 'perpendicular'
            per = per + 1;
            image_name = ['image_per' num2str(per)];
            images.(image_name)=BW;
            per_I = [per_I i];
    end
end

o = 0;
done_par = [];
done_per = [];
for j = 1:par
    for k = 1:per
        o = o + 1;
        [path1, name1, ext1] = fileparts(char(IPT.coordinateFile(par_I(j))));
        [path2, name2, ext2] = fileparts(char(IPT.coordinateFile(per_I(k))));
        picture_combination(o,:) = string(strcat(name1, ext1,'<>', name2, ext2));
        if any(done_par == j)
            image_par = images.(['image_par' num2str(j)]);
            image_par_broken = images.(['image_par_broken' num2str(j)]);
        else
            image_par = [];
            image_par_broken = [];
        end
        if any(done_per == k)
            image_per = images.(['image_per' num2str(k)]);
        else
            image_per = [];
        end
        switch char(IPT.mode(1))
            case 'normal'
            % calculate Fts of intact or broken grains
            [Results(o), vol_photo, images.(['image_par' num2str(j)]), images.(['image_par_broken' num2str(j)]), images.(['image_per' num2str(k)])] = calcFTphoto(IPT, images.(['image_par' num2str(j)]), images.(['image_per' num2str(k)]), ps, image_par, image_par_broken, image_per);
            case 'cut'
            % calculate Fts of grains cut by polishing parallel to the c-axis
            [Results(o), vol_photo] = calcFTphotoCut(IPT, images.(['image_par' num2str(j)]), images.(['image_per' num2str(k)]));
        end
        done_par = [done_par j];
        done_per = [done_per k];
        %% calculate sphere-equivalent radius (ser)

        % calculate the volume
        V_num(o) = sum(vol_photo(:)) * (ps)^3;

        % calculate surface area
        fv = isosurface(vol_photo,0.5);
        a = fv.vertices(fv.faces(:, 2), :) - fv.vertices(fv.faces(:, 1), :);
        b = fv.vertices(fv.faces(:, 3), :) - fv.vertices(fv.faces(:, 1), :);
        c = cross(a, b, 2);
        A_num(o) = 1/2 * sum(sqrt(sum(c.^2, 2))) * (ps)^2;


    end
end

% calculate the mass
Mass_num=V_num*1E-4^3*3.2*1E6;

% sphere-equivalent radius
SER_num=3*V_num./A_num;

%% plot grain model
%plot_grain(vol_photo,0.5)

%% write output table
output_table=table(picture_combination,[Results.Ft_num]',A_num',V_num',SER_num',[Results.Ft_ana]',[Results.A_ana]',[Results.V_ana]',[Results.SER_ana]',...
   'VariableNames',{'Picture_combination' 'Ft_num' 'Area_num' 'Volume_num' 'SER_num' 'Ft_ana' 'Area_ana' 'Volume_ana' 'SER_ana'});
% add mean values
mean_values={'Mean',mean([Results.Ft_num]),mean(A_num),mean(V_num),mean(SER_num),mean([Results.Ft_ana]),mean([Results.A_ana]),mean([Results.V_ana]),mean([Results.SER_ana])};
output_table=[output_table;mean_values];
% add standard devation of mean values
std_values={'Standard_deviation',std([Results.Ft_num]),std(A_num),std(V_num),std(SER_num),std([Results.Ft_ana]),std([Results.A_ana]),std([Results.V_ana]),std([Results.SER_ana])};
output_table=[output_table;std_values];
% write table to csv file
writetable(output_table, output_file)
