function Temp_out=Int_RESTART_2D_reanalysis(lat_h,lon_h,mask_h,lat_c1,lon_c1,mask_c,field_2D,e3t,lev,file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation of 2D field in hybrid s-z coordinates 
% Anna Katavouta, NOC, Liverpool 09/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read the field and interpolate
[lat_c lon_c]=meshgrid(lat_c1,lon_c1);

Temp_in=ncread(file,'zos');
Temp_out=griddata(double(lon_c),double(lat_c),double(mask_c(:,:,1).*Temp_in),double(lon_h),double(lat_h),'linear');

%flood it to deal with any problems with mismatched mask/land between your
%model and the data
Temp_out=inpaint_nans(Temp_out,2);

%if you want do not want to mask your flooded field with the nemo mask comment the
%line below (the fields will look/saved flooded)
Temp_out=Temp_out.*mask_h(:,:,1);Temp_out(isnan(Temp_out))=0;

%% read extra / set up the file to write / write the file
x=size(lon_h,1);y=size(lon_h,2);z=length(lev);
filename='MYRESTART.nc'

if ~isfile(filename) 
    time=1;
    ntime=0;
    nccreate(filename,'nav_lat', 'Dimensions',{'x',x,'y',y})
    ncwrite(filename,'nav_lat',lat_h);
    nccreate(filename,'nav_lon', 'Dimensions',{'x',x,'y',y})
    ncwrite(filename,'nav_lon',lon_h);
    nccreate(filename,'nav_lev', 'Dimensions',{'nav_lev',z})
    ncwrite(filename,'nav_lev',lev);
    nccreate(filename,'time_counter', 'Dimensions',{'time_counter',Inf},'Datatype','double')
    ncwrite(filename,'time_counter',time);
    nccreate(filename,'ntime', 'Dimensions',{'time_counter',Inf},'Datatype','double')
    ncwrite(filename,'ntime',ntime);
    nccreate(filename,'e3t_b', 'Dimensions',{'x',x,'y',y,'nav_lev',z,'time_counter',Inf},'Datatype','double')
    ncwrite(filename,'e3t_b',e3t);
    nccreate(filename,'e3t_n', 'Dimensions',{'x',x,'y',y,'nav_lev',z,'time_counter',Inf},'Datatype','double')
    ncwrite(filename,'e3t_n',e3t);
end

nccreate(filename,field_2D, 'Dimensions',{'x',x,'y',y,'time_counter',Inf},'Datatype','double')
ncwrite(filename,field_2D,Temp_out);

end

