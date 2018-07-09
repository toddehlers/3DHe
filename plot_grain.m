function [ output_args ] = plot_grain(volume,isoval)

cut_surf_x=size(volume,1);
cut_surf_y=size(volume,2);
cut_surf_z=size(volume,3);

figure
hold on
Ds = smooth3(volume(1:cut_surf_x,1:cut_surf_y,1:cut_surf_z));
hiso = patch(isosurface(Ds,isoval),...
   'FaceColor',[1,.75,.65],...
   'EdgeColor','none');
isonormals(Ds,hiso)
hcap = patch(isocaps(volume(1:cut_surf_x,1:cut_surf_y,1:cut_surf_z),isoval),...
   'FaceColor','interp',...
   'EdgeColor','none');
colormap hsv
view(35,30), axis tight; axis equal
daspect([1,1,1])
lightangle(45,30);
lighting gouraud
hcap.AmbientStrength = 0.6;
hiso.SpecularColorReflectance = 0;
hiso.SpecularExponent = 50;
xlabel('microns')
ylabel('microns')
zlabel('microns')


end

