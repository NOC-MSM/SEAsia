clear;clc;close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Anna Katavouta NOC 07/2021
% create river forcing for NEMO from 
% JRA55 rivers dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% choose year
year=2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read SEAsia domain (or your domain)
file='domain_cfg.nc'
lat_region=ncread(file,'nav_lat');
lon_region=ncread(file,'nav_lon');
e1t_region=ncread(file,'e1t');
e2t_region=ncread(file,'e2t');
% define your land mask, alternative use the actual mask file from nemo
Land_region=single(ncread(file,'top_level'));Land_region(Land_region==0)=nan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read JRA55 rivers
file=['/projectsa/NEMO/slwa/global_rivers/JRA55-do-1-5-0/friver_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_',num2str(year),'0101-',num2str(year),'1231.nc'];
% for simplicity/rapid extract only your rough subdomain from the global JRA55 domain, you do not have to.
ilat=200; nlat=300; ilon=250; nlon=400;
x1=ncread(file,'lon',[ilon],[nlon]);
lat=ncread(file,'lat',[ilat],[nlat])';
lon=x1;lon(x1>180)=lon(x1>180)-360; %wrap to -180 180 your longitude
friver=ncread(file,'friver',[ilon ilat 1],[nlon nlat Inf]); % in kg m-2 s-1
friver(friver==0)=NaN;

% cell areas
file='/projectsa/NEMO/slwa/global_rivers/JRA55-do-1-5-0/areacello_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr.nc';
a1=ncread(file,'areacello',[ilon ilat],[nlon nlat]);
friver=friver.*a1; % kg/s
rho0=1000; %you can use exact water density but 1/1000 for freshwater will do
river_outlow=friver./rho0; % m3/s 
mask=friver;mask(~isnan(mask))=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% find coast in our domain
Re=6367456.0*pi/180;
dr=mean([e1t_region(:) ; e2t_region(:)]);
dx=mean(e1t_region(:))./(Re.*cos(mean(lon_region(:)).*180/pi));
dy=mean(e2t_region(:))./Re;

%get mask coast_region
mask_coast= isnan(Land_region);
mask_sea= Land_region==1;

mask_coast = [mask_coast(1,:); mask_coast; mask_coast(end,:)] ;
mask_coast = [mask_coast(:,1), mask_coast, mask_coast(:,end)] ;

mask_coast = mask_coast(1:end-2,2:end-1) + mask_coast(3:end,2:end-1)      ...
          + mask_coast(2:end-1,1:end-2) + mask_coast(2:end-1,3:end) ;
mask_region_C = (mask_coast>0 & +mask_sea==1) ;

% No rivers on boundaries
mask_region_C(1,:)=0; mask_region_C(end,:)=0;mask_region_C(:,1)=0;mask_region_C(:,end)=0;
mask_region_C(2,:)=0; mask_region_C(end-1,:)=0;mask_region_C(:,2)=0;mask_region_C(:,end-1)=0;

mask_region=double(mask_region_C);
mask_region(isnan(Land_region))=nan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% loop over time as different mask and sometimes rivers in some datasets
[LLAT LLON]=meshgrid(lat,lon);

for tt=1:size(river_outlow,3)
    River_remap=remap_river(mask(:,:,tt),LLON,LLAT,mask_region,lon_region,lat_region,river_outlow(:,:,tt));
    River(:,:,tt)=River_remap;
    tt
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% spread large rivers to adjacent sea points
% define a threshold for maximum river volume flux in one grid cell
THRESOLD= 1500 %in m^3/s so it depends from your model gridcell/resolution how large it should be

River_uns=River;
mask_spread=mask_region; mask_spread(~isnan(mask_spread))=1;
for tt=1:size(River,3)
    [Indi,Indj]=find(squeeze(River(:,:,tt))>THRESOLD);
    Rad_s=2;%choose how many points to spread intitially, with each iteration this 
    %radius will grow until all large rivers have been spread
    
    while Indi>0
    for ii=1:length(Indi)
         Dummy2=mask_spread(Indi(ii)-Rad_s:Indi(ii)+Rad_s,Indj(ii)-Rad_s:Indj(ii)+Rad_s);
         np=nansum(nansum(Dummy2,2),1);
         [Di, Dj]=find(~isnan(Dummy2));
         Ri_SUM=nansum(nansum(River(Indi(ii)-Rad_s:Indi(ii)+Rad_s,Indj(ii)-Rad_s:Indj(ii)+Rad_s,tt),2),1);
         for jj=1:length(Di)
             River(Indi(ii)-(Rad_s+1)+Di(jj),Indj(ii)-(Rad_s+1)+Dj(jj),tt)=Ri_SUM/np;
         end
    end
    % repeat until all the large rivers have been spread
    Rad_s=Rad_s+1;
    display([num2str(Rad_s),' with time ',num2str(tt)])
    [Indi,Indj]=find(squeeze(River(:,:,tt))>THRESOLD);
    end
    
    % check that you included all river output and everythign consistent 
    if (nansum(squeeze(nansum(River_uns(:,:,tt),1)),2)-nansum(squeeze(nansum(River(:,:,tt),1)),2))>0.000001
        tt
        display (['Error total runoff in t=',num2str(tt)])
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% rivers from m^3/s to to Kg m^(-2) s^(-1) that nemo wants
for tt=1:size(River,3)
    rorunoff(:,:,tt)=River(:,:,tt).*rho0./(e1t_region.*e2t_region);
end    
rorunoff(isnan(rorunoff))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% write data	
output=(['rivers_JRA55_y',num2str(year),'.nc']);
nx=size(lon_region,1);
ny=size(lon_region,2);
tt=1:size(rorunoff,3);

nccreate(output,'time_counter','Dimensions',{'time_counter',Inf},'Format','netcdf4_classic');
nccreate(output,'nav_lat','Dimensions',{'x',nx,'y',ny},'Format','netcdf4_classic');
nccreate(output,'nav_lon','Dimensions',{'x',nx,'y',ny},'Format','netcdf4_classic');
nccreate(output,'rorunoff','Dimensions',{'x',nx,'y',ny,'time_counter',Inf},'Format','netcdf4_classic');

ncwrite(output,'time_counter',tt);
ncwrite(output,'nav_lon',lon_region);
ncwrite(output,'nav_lat',lat_region);
ncwrite(output,'rorunoff',rorunoff);
ncwriteatt(output, 'rorunoff', 'long_name','river runoff' );
ncwriteatt(output, 'rorunoff','units','Kg m^-2 s^-1');
