clear;clc;close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation/Creation of fake restart file from 
% reanalysis product to a regional model with hybrid z-sigma coordinates
% Anna Katavouta, NOC, Liverpool 09/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%addpath '/vkamino/scratch/accord/SEAsia_R36/RESTART_SEAsiaR36/Inpaint_nans'
addpath '/Volumes/Elements/SEAsia_R36/RESTART_INT/Inpaint_nans'
% the code requires the inpaint_nan functions that are available for matlab
% (need to be downloaded: https://uk.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans
%% read coordinates and mask of your regional model
file='domain_cfg.nc';
lat_regional=ncread(file,'nav_lat');
lon_regional=ncread(file,'nav_lon');
e3t_regional=ncread(file,'e3t_0');
e3u_regional=ncread(file,'e3u_0');
e3v_regional=ncread(file,'e3v_0');
lev=ncread(file,'nav_lev');

file='mesh_mask.nc';
maskt_regional=double(ncread(file,'tmask'));
masku_regional=double(ncread(file,'umask'));
maskv_regional=double(ncread(file,'vmask'));
maskt_regional(maskt_regional==0)=nan;masku_regional(masku_regional==0)=nan;maskv_regional(maskv_regional==0)=nan;

%% read coordinates and create the mask for the reanlysis product (if 
% the product has its own mask even better)
file_data='global-reanalysis-phy-001-030-daily_1637057175950.nc';
lat_reanal=ncread(file_data,'latitude');
lon_reanal=ncread(file_data,'longitude');
Depth_reanal=ncread(file_data,'depth');
mask_reanal=ncread(file_data,'so');mask_reanal(~isnan(mask_reanal))=1;

%% estimate depths from e3 level thickness
Depth_regional(:,:,1)=(e3t_regional(:,:,1)./2).*maskt_regional(:,:,1);
for zz=2:size(e3t_regional,3)
    Depth_regional(:,:,zz)=nansum((e3t_regional(:,:,1:zz-1).*maskt_regional(:,:,1:zz-1)),3)+(e3t_regional(:,:,zz)./2).*maskt_regional(:,:,zz);
end

%% Interpolation/create of restart file
%fields
field_2D=string( {'sshn'} );
field_3D=string( {'un';'vn';'tn';'sn'} );

%interpolate the 2D fields
for ii=1:length(field_2D)
   Temp_out_2D=Int_RESTART_2D_reanalysis(lat_regional,lon_regional,maskt_regional,lat_reanal,lon_reanal,mask_reanal,field_2D(ii),e3t_regional,lev,file_data);
end

%interpolate the 3D fields
for ii=1:length(field_3D)    
    if strcmp(field_3D(ii),'sn') || strcmp(field_3D(ii),'tn')
       mask_in=maskt_regional;
    end
    if strcmp(field_3D(ii),'un')
       mask_in=masku_regional;
    end
    if strcmp(field_3D(ii),'vn')
       mask_in=maskv_regional;
    end
    Temp_out_3D=Int_RESTART_3D_reanalysis(lat_regional,lon_regional,mask_in,Depth_regional,lat_reanal,lon_reanal,mask_reanal,Depth_reanal,field_3D(ii),e3t_regional,lev,file_data);
end
