function River=remap_river(mask,LLON,LLAT,mask_region,lon_region,lat_region,river_in)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remap rivers
% Anna Katavouta, NOC, Liverpool 07/2021
% mask : the land mask of the river dataset from which you extract (e.g. JRA-55)
% LLON : longitude of the river dataset (e.g., JRA-55)
% LLAT : latitude of the river dataset
% mask_region : the land mask of your regional model domain
% lon_region : the longitude of your regional model
% lat_region : the latitude of your regional model
% river_in : the river ourflow of your dataset in m^3/s (e.g., JRA-55)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    River=zeros(size(lon_region,1),size(lon_region,2));
    %% move everything to the regional model coast
    mask_log=(~isnan(mask));
    lon_log=LLON(mask_log);
    lat_log=LLAT(mask_log);

    %keep only what is in your domain;
    lat_log(lat_log>nanmax(lat_region(:)))=nan;lat_log(lat_log<nanmin(lat_region(:)))=nan;
    lon_log(isnan(lat_log))=nan;
    lon_log(lon_log<nanmin(lon_region(:)))=nan;lon_log(lon_log>nanmax(lon_region(:)))=nan;
    lat_log(isnan(lon_log))=nan;
    lon_save=lon_log;
    lat_log(isnan(lat_log))=[];
    lon_log(isnan(lon_log))=[];

    %find nearest neighbour on the regional coast grid
    lon_coast=lon_region;lat_coast=lat_region;
    lon_coast(find(mask_region~=1))=nan;lat_coast(find(mask_region~=1))=nan;
    lon_coast(isnan(lon_coast))=[];lat_coast(isnan(lat_coast))=[];
    for ii=1:length(lon_log)
       dxf = [lon_log(ii), lat_log(ii)];
       G = [lon_coast(:),lat_coast(:)];
       d =  sum( (G-dxf).^2, 2);
       [minDist, idxMinDist] = min(d);
       solution = G(idxMinDist,:);
       CC(ii,:)=solution;
       [Ic Jc]=find(lon_region==solution(1));
       ind_ci(ii)=Ic(1);
       [Ic Jc]=find(lat_region==solution(2));
       ind_cj(ii)=Jc(1); 
    end    
    
    %% restruct your matrix and
    %ensure that you add catchmens that go to the same point   
    R_log1=river_in; 
    R_log=R_log1(mask_log);
    R_log(isnan(lon_save))=[];
    for ii=1:length(ind_ci)
       if ~isnan(ind_ci(ii)) 
           River(ind_ci(ii),ind_cj(ii))=River(ind_ci(ii),ind_cj(ii))+(R_log(ii));
        end      
    end

    % check that you included all river output
    R_log(isnan(ind_ci))=nan;R_log(isnan(ind_cj))=nan;
    if (nansum(R_log)-nansum(squeeze(nansum(River,1)),2))>0.000001
       display (['Error total runoff in t=',num2str(tt)])
    end

end
