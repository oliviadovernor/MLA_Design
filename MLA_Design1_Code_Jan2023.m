% Generate 15x15 grid for Powerphotonics
% Boya
clear; clc; 

lambda = 640e-9;
k=2*pi/lambda; % wavenumber
n_imm=1.518; % refractive index of immersion medium
n_s=1.33; % refractive index of sample
%NA= 1.27; %water objective numerical aperture
NA=1.49; % oil objective numerical aperture
f_obj=2e-3; % focal length of objective lens (oil) 
%f_obj=2e-3; % focal length of objective lens (water) 
D_cam=1e-6; % pixel size of xy grid powerphotonics

D=2*NA*f_obj; % diameter of objective lens

%% 1. 3 Hex, 2.0 pitch, no tilt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_tube=200e-3; % focal length of tube lens
f_fl=175e-3; % focal length of fourier lens (f4)
f_mla=100e-3; % focal length of MLA
D_mla=2.0e-3; % diameter of u-lens, between flat edges (pitch)

M_relay=(f_fl/f_tube); % magnification of the relay system
Mag=((f_tube/f_obj)*(f_mla/f_fl));
size_bfp=D*M_relay*(n_s/n_imm); % diameter of relayed BFP
N_pixels_obj=round(size_bfp/D_cam); % number of pixels in objective plane (one dimension)
rho_max=size_bfp/2;
N_mla=size_bfp/D_mla;
N_pixels_mla=round(D_mla/D_cam); % number of pixels per microlens 
[x_mla,y_mla]=get_MLAcentres(N_mla,N_pixels_mla); % xy-coordinate of u-lens
curvature =45.67; % unit in mm depends on wavelengh 
thickness1 = get_physical_mask(D_mla,D_cam,N_pixels_obj,N_pixels_mla,x_mla,y_mla,curvature);
thickness1 = thickness1*1e6;
figure(1); imagesc(thickness1); title('3 Hex, 2.0 mm pitch')
cbar = colorbar; 
cbar.Label.String = 'Z height (um)'; 
min(thickness1(:));
T1Correction = abs(min(thickness1(:)));
T1Max = abs(max(thickness1(:)));

%% 2. 3 Hex, 1.9 pitch, slight adaptation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_tube=200e-3; % focal length of tube lens
f_fl=175e-3; % focal length of fourier lens (f4)
f_mla=100e-3; % focal length of MLA
D_mla=1.9e-3; % diameter of u-lens, between flat edges (pitch)

M_relay=(f_fl/f_tube); % magnification of the relay system
Mag=((f_tube/f_obj)*(f_mla/f_fl));
size_bfp=D*M_relay*(n_s/n_imm); % diameter of relayed BFP
N_pixels_obj=round(size_bfp/D_cam); % number of pixels in objective plane (one dimension)
rho_max=size_bfp/2;
N_mla=size_bfp/D_mla;
N_pixels_mla=round(D_mla/D_cam); % number of pixels per microlens 
[x_mla,y_mla]=get_MLAcentres(N_mla,N_pixels_mla); % xy-coordinate of u-lens
curvature =45.67; % unit in mm
thickness3 = get_physical_mask(D_mla,D_cam,N_pixels_obj,N_pixels_mla,x_mla,y_mla,curvature);
thickness3 = thickness3*1e6;
figure(2); imagesc(thickness3); title('3 Hex, 1.9 mm pitch')
cbar = colorbar; 
cbar.Label.String = 'Z height (um)'; 
min(thickness3(:));
T3Correction = abs(min(thickness3(:)));
T3Max = abs(max(thickness3(:)));


%% 1. 7 Hex, 1.5 pitch, no tilt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_tube=200e-3; % focal length of tube lens
f_fl=175e-3; % focal length of fourier lens
f_mla=100e-3; % focal length of MLA
D_mla=1.5e-3; % diameter of u-lens, between flat edges (pitch)

M_relay=(f_fl/f_tube); % magnification of the relay system
size_bfp=D*M_relay*(n_s/n_imm); % diameter of relayed BFP
N_pixels_obj=round(size_bfp/D_cam); % number of pixels in objective plane (one dimension)
rho_max=size_bfp/2;
N_mla=size_bfp/D_mla;
N_pixels_mla=round(D_mla/D_cam); % number of pixels per microlens 
[x_mla,y_mla]=get_MLAcentres(N_mla,N_pixels_mla); % xy-coordinate of u-lens
curvature = 45.67;
pixel_size= (D_cam/Mag)*1000000000; % nm 
thickness7 = get_physical_mask_7hex(D_mla,D_cam,N_pixels_obj,N_pixels_mla,x_mla,y_mla,curvature);
thickness7 = thickness7*1e6;
figure(3); imagesc(thickness7); title('7 Hex, 1.5 mm pitch')
cbar = colorbar; 
cbar.Label.String = 'Z height (um)';
min(thickness7(:));
T7Correction = abs(min(thickness7(:)));
T7Max = abs(max(thickness7(:)));

%% 2. 7 Hex, 1.4 pitch,slight adaptation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_tube=200e-3; % focal length of tube lens
f_fl=175e-3; % focal length of fourier lens
f_mla=100e-3; % focal length of MLA
D_mla=1.4e-3; % diameter of u-lens, between flat edges (pitch)

M_relay=(f_fl/f_tube); % magnification of the relay system
size_bfp=D*M_relay*(n_s/n_imm); % diameter of relayed BFP
N_pixels_obj=round(size_bfp/D_cam); % number of pixels in objective plane (one dimension)
rho_max=size_bfp/2;
N_mla=size_bfp/D_mla;
N_pixels_mla=round(D_mla/D_cam); % number of pixels per microlens 
[x_mla,y_mla]=get_MLAcentres(N_mla,N_pixels_mla); % xy-coordinate of u-lens
curvature = 45.67;
pixel_size= (D_cam/Mag)*1000000000; % nm 
thickness4 = get_physical_mask_7hex(D_mla,D_cam,N_pixels_obj,N_pixels_mla,x_mla,y_mla,curvature);
thickness4 = thickness4*1e6;
figure(4); imagesc(thickness4); title('7 Hex, 1.4 mm pitch')
cbar = colorbar; 
cbar.Label.String = 'Z height (um)';
min(thickness4(:));
T4Correction = abs(min(thickness4(:)));
T4Max = abs(max(thickness4(:)));

%% Downsample grid from 1um to 10 um to fit on PowerPhotonics grid 

hex3_2mm_10um = thickness1(10:10:end, 6:10:end); 
figure(5); imagesc(hex3_2mm_10um); title('downsampled grid 10um 2mm 3hex')
cbar = colorbar; 
cbar.Label.String = 'Z height (um)'; 
min(hex3_2mm_10um(:));
DownsampledMin1 = abs(min(hex3_2mm_10um(:))); 
DownsampledMax1 = abs(max(hex3_2mm_10um(:))); 

hex3_19mm_10um = thickness3(6:10:end, 6:10:end); 
figure(6); imagesc(hex3_19mm_10um); title('downsampled grid 10um 1.9mm 3hex')
cbar = colorbar; 
cbar.Label.String = 'Z height (um)';
min(hex3_19mm_10um(:)); 
DownsampledMin3 = abs(min(hex3_19mm_10um(:)));
DownsampledMax3 = abs(max(hex3_19mm_10um(:))); 

hex7_15mm_10um = thickness7(6:10:end, 6:10:end); 
figure(7); imagesc(hex7_15mm_10um); title('Downsampled grid 10um 1.5mm 7 hex')
cbar = colorbar; 
cbar.Label.String = 'Z height (um)'; 
min(hex7_15mm_10um(:)); 
DownsampledMin7 = abs(min(hex7_15mm_10um(:))); 
DownsampledMax7 = abs(max(hex7_15mm_10um(:))); 

hex7_14mm_10um = thickness4(6:10:end, 6:10:end); 
figure(8); imagesc(hex7_14mm_10um); title('Downsampled grid 10um 1.4mm 7 hex')
cbar = colorbar; 
cbar.Label.String = 'Z height (um)'; 
min(hex7_14mm_10um(:)); 
DownsampledMin4 = abs(min(hex7_14mm_10um(:))); 
DownsampledMax4 = abs(max(hex7_14mm_10um(:)));


%%%Figure for thesis (manufacture imposed resolution)
%data = thickness4(:); %1um sampling 
%data = hex7_14mm_10um(:); downsampled data 10um
%data = round(hex7_14mm_10um*1000)/1000;% Rounds z values to 0.001 um (z = 1nm resolution) 
%y_value = 229;
% Extract the row from the data matrix
%row_data = data(y_value,:);
% Plot x against z
%figure(11);
%plot(row_data);
%xlim([0, 456]);
%xlabel('x (um)');
%ylabel('Z height (um)');
%title('Plot of Z height for y = 229');

%% Zemax modelling needs to be in mm not um 
ExtractedLensT4 = hex7_14mm_10um(:); % units in um
ExtractedLensT4 = ExtractedLensT4/1000; %units in mm
MLAT4_Zemax_mm = reshape(ExtractedLensT4',457^2,1); 
writematrix(MLAT4_Zemax_mm, '7MLA_14mm_pitch');
max_zemaxT4 = abs(max(ExtractedLensT4(:)));
min_zemaxT4 = abs(min(ExtractedLensT4(:)));

ExtractedLensT3 = hex3_19mm_10um(1:350,1:350); % units in um
ExtractedLensT3 = ExtractedLensT3/1000; %units in mm
MLAT3_Zemax_mm = reshape(ExtractedLensT3',350^2,1); 
writematrix(MLAT3_Zemax_mm, 'Lens_3MLA_19mm_pitch');
max_zemaxT3 = abs(max(ExtractedLensT3(:)));
min_zemaxT3 = abs(min(ExtractedLensT3(:)));

ExtractedLensT7 = hex7_15mm_10um(:); % units in um
ExtractedLensT7 = ExtractedLensT7/1000; %units in mm
MLAT7_Zemax_mm = reshape(ExtractedLensT7',457^2,1); 
writematrix(MLAT7_Zemax_mm, '7MLA_15mm_pitch');
max_zemaxT7 = abs(max(ExtractedLensT7(:)));
min_zemaxT7 = abs(min(ExtractedLensT7(:)));

%% Full 15x15 layout
full_map = NaN(1500,1500);
spacing1 = floor((1500-457*2)/4);
full_map((1:456)+(spacing1),(1:457)+(spacing1)) = hex3_2mm_10um;
full_map((1:457)+(spacing1),(1:457)+(spacing1*3)+457) = hex7_15mm_10um;
full_map((1:457)+(spacing1*3)+457,(1:457)+(spacing1)) = hex7_14mm_10um;
full_map((1:457)+(spacing1*3)+457,(1:457)+(spacing1*3)+457) = hex3_19mm_10um;

figure(9)
imagesc(full_map)
title('Full 15x15 map')
cbar = colorbar; 
cbar.Label.String = 'Z height (um)'; 

%% Generating LightForge format - xyz in um, first line counting in 10s
size_map = size(full_map,1);
spacing = 10; 
coordinate_map = [];

XYZ_grid = zeros(size_map+1,size_map+1); %add zeros to first row and first column 
XYZ_grid(1,:) = 0:10:(size_map*10); %adds array of 1500 in x 
XYZ_grid(:,1) = 0:10:(size_map*10); %and y to XYZ_grid
XYZ_grid(2:(size_map+1),2:(size_map+1)) = full_map; % data populated from the second row and second column of XYZ_grid
XYZ_grid = round(XYZ_grid*1000)/1000;% Rounds z values to 0.001 um (z = 1nm resolution) 
DesignMLA = XYZ_grid;

figure(10)
visible_data = DesignMLA(2:end, 2:end);
imagesc(visible_data)
title('Powerphotonics - MLA Design')
% Adjust colormap and color axis range
cbar = colorbar; 
cbar.Label.String = 'Z height (um)'; 
caxis([min((visible_data(:))), max((visible_data(:)))]); %#ok<CAXIS> % Set the color axis range to the minimum and maximum values of the data

writematrix(DesignMLA,'v2_MLADesign_f88mm.csv')
writematrix(DesignMLA,'v2_MLADesign_f88mm.dat')
