function[thickness] = get_physical_mask_7hex(D_mla,D_cam,N_pixels_obj,N_pixels_mla,x_mla,y_mla,curvature)

% Calculating sagitta of a single u-lens
R = curvature*1e-3; % radius of curvature [mm]
interval = D_cam;

% Generating MLA+mask
N_views = numel(x_mla);
[px,py] = meshgrid((-N_pixels_obj/2):(N_pixels_obj/2-1),(-N_pixels_obj/2):(N_pixels_obj/2-1));
u = zeros(size(px,1),size(px,2));
v = u;
mla_pattern = zeros(N_pixels_obj,N_pixels_obj);

total_inpolygon = zeros(size(px,1),size(px,2));
[x,y] = meshgrid((-N_pixels_obj/2):(N_pixels_obj/2-1),(-N_pixels_obj/2):(N_pixels_obj/2-1));


for i =1:N_views
    hexgon_shape = nsidedpoly(6,'Center',[x_mla(i) y_mla(i)],'SideLength',N_pixels_mla/sqrt(3));
    points_on_hexagon_X(i,:) = hexgon_shape.Vertices(:,1);
    points_on_hexagon_Y(i,:) = hexgon_shape.Vertices(:,2);
    in_polygon = inpolygon(px,py,points_on_hexagon_X(i,:),points_on_hexagon_Y(i,:));
    u = u+in_polygon.*x_mla(i)/N_pixels_mla;
    v = v+in_polygon.*y_mla(i)/N_pixels_mla;
    lens_shape = calc_sag(R,N_pixels_obj,N_pixels_mla,interval,[x_mla(i) y_mla(i)]);
    lens_shape = lens_shape.*in_polygon;
    mla_pattern = mla_pattern+lens_shape;  
    total_inpolygon = total_inpolygon+in_polygon;
end

thickness = mla_pattern;
% thickness = (max(thickness,[],'all')-thickness).*total_inpolygon;

% thickness = thickness - not(total_inpolygon)*max(mask_pattern+sag,[],'all');