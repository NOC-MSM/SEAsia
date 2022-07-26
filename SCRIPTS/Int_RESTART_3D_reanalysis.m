function Temp_out=Int_RESTART_3D_reanalysis(lat_h,lon_h,mask_h,Depth_h,lat_c1,lon_c1,mask_c,Depth_c,field_3D,e3t,lev,file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation of 3D field in hybrid s-z coordinates 
% Anna Katavouta, NOC, Liverpool 09/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read the field and interpolate
%those names are specific to copernicus data 
%they will need to be changed based of the data you are using
if strcmp(field_3D,'un')
    name_read=string( {'uo'} );
end
if strcmp(field_3D,'vn')
    name_read=string( {'vo'} );
end
if strcmp(field_3D,'sn')
    name_read=string( {'so'} );
end
if strcmp(field_3D,'tn')
    name_read=string( {'thetao'} );
end
Temp_in=ncread(file,[name_read]);
[lat_c lon_c]=meshgrid(lat_c1,lon_c1);

    % 1. interpolate horizontaly
    Temp2=double(squeeze(Temp_in.*mask_c));
    for kk=1:size(Temp2,3)
        Temp_in_21(:,:,kk)=griddata(double(lon_c),double(lat_c),double(Temp2(:,:,kk)),double(lon_h),double(lat_h),'linear');
        %flood only for T and S not for velocities
        if strcmp(field_3D,'tn') || strcmp(field_3D,'sn')
           Temp_in_2(:,:,kk)=inpaint_nans(Temp_in_21(:,:,kk),2);
        else
           Temp_in_2(:,:,kk)=Temp_in_21(:,:,kk);
        end
    end
    
    %2. interpolate vertically 
    for k1=1:size(Temp_in_2,1)
        for k2=1:size(Temp_in_2,2)
            Temp_d=double(squeeze(Depth_h(k1,k2,:)).*squeeze(mask_h(k1,k2,:)));Temp_d(isnan(Temp_d))=[];
            if length(Temp_d)>0
                TT1=interp1(Depth_c,double(squeeze(Temp_in_2(k1,k2,:))),Temp_d);
                if length(TT1)<length(lev)
                 TT1(length(TT1)+1:length(lev))=nan;
                end
                 k_start=find(Temp_d<Depth_c(1));%extrapolate value at the surface though
                 if length(k_start)>0
                 TT1(k_start)=TT1(nanmax(k_start+1));
                 end
             else
                TT1(1:size(mask_h,3))=nan;
            end
            Temp_out(k1,k2,:)=TT1;
        end
    end
    
   %flood again !not necessary but just in case
   %uncomment the lines below if you want to flood again
   %if strcmp(field_3D,'tn') || strcmp(field_3D,'sn')
   %     for kv=1:size(Temp_out,3)
   %         Temp_out(:,:,kv)=inpaint_nans(Temp_out(:,:,kv),2);
   %     end
   % end

    %if you want do not waht to mask your flooded field with the nemo mask comment the
    %line below (your fields will look / be saved flooded)
    %Temp_out=mask_h.*Temp_out;
    
    Temp_out(isnan(Temp_out))=0;

%% set up  and write the netcdfi RESTART file
x=size(lon_h,1);y=size(lon_h,2);z=size(Temp_out,3);
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
nccreate(filename,field_3D, 'Dimensions',{'x',x,'y',y,'nav_lev',z,'time_counter',Inf},'Datatype','double')
ncwrite(filename,field_3D,Temp_out);

end

