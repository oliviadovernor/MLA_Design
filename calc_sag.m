function [sag] = calc_sag(R,N_pixels_obj,N_pixels_mla,interval,lenscentre)
% calculating lens sagitta
% R = 36.63e-3; % radius of curvature [mm]
% interval - can be the manufacturing resolution or camera pixel size [m]
% lenscentre - xy coordinate of MLA centre [pixels]

d = (-N_pixels_obj/2):(N_pixels_obj/2-1);
d = d*interval;
lenscentre = lenscentre*interval;

s = R-sqrt(R^2-(((2*N_pixels_mla/sqrt(3)))/2*interval).^2); %take diagonal to make sag extend over all hexagon area

% tilt in X direction
[x,y] = meshgrid(d,d);
rho = sqrt(((x-lenscentre(1)).^2)+((y-lenscentre(2)).^2));
sag = sqrt(R^2-rho.^2)+s-R;

