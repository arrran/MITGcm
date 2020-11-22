clear all

%% Directory:
 cd ('/mnt/data/malyaral/MITgcm-checkpoint66j/exp/ch_1/run');
 %%
 addpath (genpath('/mnt/data/malyaral/OneDrive/MATLAB/MITgcm'))


 %% Basic data load:
 dt=1800;
startdate=datenum(1992,01,01);
 
 
XC=rdmds('XC');% xaxis
YC=rdmds('YC');% yaxis
RC=rdmds('RC');% zaxis
% mask_surf=rdmds('maskInC');% 
% mask_surf(mask_surf==0)=NaN; 

hFacC=rdmds('hFacC');%
mask_ocean_surf(:,:)=hFacC(:,:,1);
mask_ocean_surf(mask_ocean_surf==0)=NaN;

shelf=rdmds('shelficemassinit');
shelf_mask=shelf;
shelf_mask(shelf>0)=1;
shelf_mask(shelf==0)=NaN;

%
hFacC=rdmds('hFacC');%
mask_ocean_surf(:,:)=hFacC(:,:,1);
mask_ocean_surf(mask_ocean_surf==0)=NaN;
mask_ocean_surf(mask_ocean_surf<1)=1;

% figure
% pcolor (mask_ocean_surf')
% shading flat
 %
    
         colorOrder = ...
[  0 0 1 % 1 BLUE dark
   0 1 0 % 2 GREEN (pale)
   1 0 0 % 3 RED
   0 1 1 % 4 CYAN
   1 0 1 % 5 MAGENTA (pale)
   1 1 0 % 6 YELLOW (pale)
   0 0 0 % 7 BLACK
   0 0.75 0.75 % 8 TURQUOISE
   0 0.5 0 % 9 GREEN (dark)
   0.75 0.75 0 % 10 YELLOW (dark)
   1 0.50 0.25 % 11 ORANGE
   0.75 0 0.75 % 12 MAGENTA (dark)
   0.7 0.7 0.7 % 13 GREY
   0.8 0.7 0.6 % 14 BROWN (pale)
   0.6 0.5 0.4  % 15 BROWN (dark)
   0.3 0.5 0.8 % 16
   0.3010 0.7450 0.9330 % 17light blue
   0.4 0 0.9 % 18 purple
   245/245 66/245 105/245]; % 19 


%
nx=50;
ny=50;
nz=30;

for i=1:nx
    for j=1:ny
        for k=1:nz
            lat_plot(i,j,k)=YC(i,j);
            lon_plot(i,j,k)=XC(i,j);
            depth_plot(i,j,k)=RC(k);
        end
    end
end
 
%%
clear time_char
for i=1:10   
   time_options(i,1)=216*(i-1); 
   
   time_char(i,:)=[year(startdate+time_options(i,1)*dt/(3600*24)) , month(startdate+time_options(i,1)*dt/(3600*24)), day(startdate+time_options(i,1)*dt/(3600*24)) ];
end

%%
 NF=5;
 time=time_options(NF);


%%
  
 % Plots:
 % Horizontal plots of 2d variables
      % Pcolor in nx ny
     plot_depth          =[1];
     plot_eta            =[1]; 
     plot_shelf_melt     =[1];
     plot_si_area        =[1];
     plot_si_heff   =[1];
      
     plot_EXFhs          =[1];
     plot_EXFhl          =[1];
     plot_EXFpreci       =[1];
     plot_EXFswnet       =[1];
     plot_EXFlwnet       =[1];
     plot_EXFatemp       =[1];
     plot_EXFqnet        =[1];     
     plot_EXF_hflux      =[1];
     
     plot_SHIuStar       =[1];
     plot_SHIgammT       =[1];
     plot_SHIgammS       =[1];
     plot_SHIfwFlx       =[1];
     plot_SHIhtFlx       =[1];
     plot_oceFreez       =[1];
     
     plot_KPPdiffKzS     =[1 70];       % true level
     plot_KPPdiffKzT     =[1 70];       % true level
     plot_KPPhbl         =[1 1];
     
     plot_SIarea         =[1];
     plot_SIuice         =[1];
     plot_SIvice         =[1];
     plot_SIqnet         =[1];
     plot_SIempmr        =[1];
     plot_SIatmQnt       =[1];
     
 % Horizontal plots of 3d variables
     % Pcolor in nx ny 
     plot_t_ic           =[0 1];        % true level
     plot_s_ic           =[1 25];        % true level
     plot_t              =[1 25];        % true level
     plot_s              =[1 25];        % true level
     
     plot_v              =[1 3];        % true level
     plot_u              =[1 3];       % true level
     plot_vel            =[1 22 1 3];   % true level step scale (do not change step)
     plot_SI_velocity    =[1 1 3];   % true level step scale (do not change step)
     
     % Pcolor in Lat Lon 
     plot_t_geo          =[1 25];         % true level
     plot_s_geo          =[1 1];         % true level
     plot_vel_geo        =[1 22 2 3];    % true level step scale
     
     plot_depth_geo      =[1];
     
 % Vertical Sections - Meridional
    % Depth with meridional section (define I)
    plot_depth_sec_m     =[1 23];        % true section
    
    plot_s_sec_m_geo     =[1 25];
    plot_t_sec_m_geo     =[1 25];
    plot_u_sec_m_geo     =[1 70];

    plot_t_sec_m         =[1 45];
    plot_s_sec_m         =[1 25];
    plot_u_sec_m         =[1 25];
 
  % Vertical Sections - Zonal
    % Depth with Zonal section (define J)
    plot_depth_sec_z     =[1 63];        % true section
    
    plot_s_sec_z_geo     =[1 22];
    plot_t_sec_z_geo     =[1 22];
    plot_v_sec_z_geo     =[1 71];

    plot_t_sec_z         =[1 60];
    plot_s_sec_z         =[1 71];
    plot_v_sec_z         =[1 71];
    
    % Along front section
    plot_t_sec_front     =[1];
    plot_s_sec_front     =[1];
    plot_mixing_sec_front     =[1];
   %TS diagram
   plot_ts               =[1];
    
  % Time Series
  plot_time_t              =[1 1 NF; 17 77 01; 63 73 01;  63 73 25; 43 79 01; ];        % true; time1; time_end;  point XYZ; 
  plot_time_s              =[1 1 NF; 17 100 01; 17 100 33; 100 68 01; 100 68 33; 44 81 01; 44 81 25];        % true; time1; time_end;  point XYZ; 
  
  % Time Series with ECCO
  plot_time_t_ecco         =[1 1 NF; 17 77 1; 63 73 1; 43 79 1];        % true; time1; time_end;  point XYZ; 
  plot_time_s_ecco         =[1 1 NF; 17 77 1; 25 77 1; 17 77 5];        % true; time1; time_end;  point XYZ; 
  plot_time_t_ecco_oisst   =[1 1 NF];
  plot_time_si_ecco_nsidc   =[1 1 NF];
  
  % Vertical profiles in time at 1 locations
  plot_tz_t              =[1 1 NF; 63 73 0];    
  plot_tz_t_offset       =[1 1 NF; 44 81 0];    
  % Vertical profiles at 1 time 
  plot_tz_loc_t              =[1 NF 0; 63 73 0; 25 77 0];    
  

 % Horizontal plots with ECCO 
 plot_t_geo_with_ecco        =[1 3];  
 plot_SI_geo_with_ecco       =[1 ];  
 
 plot_SI_geo_with_nsidc       =[1 ];  


plot_t_itp    =[1];

plot_s_itp    =[1];


%% Temp IC Horizontal 
if (plot_t_ic(1)==1)
    clear var plotter
    var= rdmds('T',0);
    figure
    plotter(:,:)=var(:,:,plot_t_ic(2));
    pcolor(plotter'); 
    caxis([-2 1])
    colorbarset(gca,'YDir', 'Reverse');
    shading flat
    title (['Temperature IC Level' num2str(plot_s_ic(2)) ])
end

%% Salt IC Horizontal 
if (plot_s_ic(1)==1)
    clear var plotter 
    var= rdmds('S',0);
    figure
    plotter(:,:)=var(:,:,plot_s_ic(2));
    pcolor(plotter'); 
    caxis([33 35])
    colorbar
    shading flat 
    title (['Salinity IC Level' num2str(plot_s_ic(2)) ])
    
end

%% Depth
if (plot_depth(1)==1)
    clear var plotter
    var=rdmds('Depth');
    figure
    plotter(:,:)=var(:,:);
    pcolor(plotter'); 
    colorbar
    shading flat    
    title ('Depth')
     caxis ([ 0 1000])
end

%% Eta
if (plot_eta(1)==1)
    clear var plotter
    var=rdmds('Eta',time);
    figure
    plotter(:,:)=var(:,:);
    pcolor(plotter'); 
    colorbar
    shading flat
    title (['Eta ' ' Date=' datestr(startdate+ (time*dt/(3600*24)))])
end

%% Ice Shelf meltrate  
if (plot_shelf_melt(1)==1)
    clear var plotter
    melt_coeff=3600*24*365/1028;
    var= rdmds('SHICE_fwFlux',time);
%     load('/mnt/data/malyaral/OneDrive/MATLAB/MITgcm/5_ross_1/IC_iceshelf_part1.mat', 'base')
%     var(isnan(base))=NaN;
    figure
    plotter(:,:)=var(:,:).*melt_coeff;
    pcolor(plotter'); 
       caxis([-1 1])
    colorbar
%      colormap (map)
    shading flat
    title (['Meltrate m/year' ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Sea ice Area
if (plot_si_area(1)==1)
 clear var plotter lat_trg depth_trg
 var=rdmds('AREA',time);
    figure
    plotter(:,:)=var(:,:);
    pcolor(plotter'); 
%     caxis([0 1])
    colorbar
    shading flat 
    title (['Sea Ice Area' ' Date=' datestr(startdate+ (time*dt/(3600*24)))])
end

%% Sea ice Area from diagnostic
if (plot_si_heff(1)==1)
 clear var plotter lat_trg depth_trg
 var=rdmds('HEFF',time);
    figure
    plotter(:,:)=var(:,:);
    pcolor(plotter'); 
     caxis([0 2])
    colorbar
    shading flat 
    title (['Sea Ice Thickness' ' Date=' datestr(startdate+ (time*dt/(3600*24)))])
end

%%
% KPP diagnostic
%
%% KPP Diagnostic - KPP hsbl 
if (plot_KPPhbl(1)==1)    
    clear var plotter    
    var= rdmds('kppDiag',time);    
    figure
%     plotter(:,:)=var(:,:,3).*mask_surf;
     plotter(:,:)=var(:,:,3);
    pcolor(plotter'); 
      caxis([0 30])
    colorbar
    colormap (RED)
    shading flat
    title (['KPP hbl' ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%%


%%
% EXF diagnostic
%
%% Surface Diagnostic - EXF hs 
if (plot_EXFhs(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
    plotter(:,:)=var(:,:,16).*mask_surf;
%     plotter(:,:)=var(:,:,18);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
%     colormap (RED)
    shading flat
    title (['Sensible heat, >0 increase T.' ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Surface Diagnostic - EXF hl 
if (plot_EXFhl(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
    plotter(:,:)=var(:,:,17).*mask_surf;
%     plotter(:,:)=var(:,:,19);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
%     colormap (RED)
    shading flat
    title (['Latent heat, >0 increase T.  Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Surface Diagnostic - EXF precipitiation
if (plot_EXFpreci(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
     figure
    plotter(:,:)=var(:,:,18).*mask_surf;
     plotter(:,:)=var(:,:,20);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
%     colormap (RED)
    shading flat
    title (['Precipitation . Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Surface Diagnostic - EXF sw net
if (plot_EXFswnet(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
    plotter(:,:)=var(:,:,19).*mask_surf;
%     plotter(:,:)=var(:,:,21);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['SW net upward (-350 0). <0 increases T . Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Surface Diagnostic - EXF lw net
if (plot_EXFlwnet(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
    plotter(:,:)=var(:,:,20).*mask_surf;
%     plotter(:,:)=var(:,:,22);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['LW net upward (-20 170). >0 decreases T . Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Surface Diagnostic - EXF atemp
if (plot_EXFatemp(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
%     plotter(:,:)=var(:,:,21).*mask_surf-273.15;
     plotter(:,:)=var(:,:,21);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['Temperature 2m C. Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Surface Diagnostic - EXF qnet
if (plot_EXFqnet(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
    plotter(:,:)=var(:,:,22).*mask_surf;
%     plotter(:,:)=var(:,:,24);
    pcolor(plotter'); 
%      caxis([-200 0])
    colorbar
    shading flat
    title (['Net upward heat flux. >0 decreases T.  Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Surface Diagnostic - EXF hflux
if (plot_EXF_hflux(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
    plotter(:,:)=(var(:,:,16)+var(:,:,17)+var(:,:,18)+var(:,:,19)).*mask_surf;
%     plotter(:,:)=var(:,:,18)+var(:,:,19)+var(:,:,20)+var(:,:,21);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['Net upward heat flux (-250 to 600). >0 decreases T.  Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%%
% Shelfice diagnostic
%
%% Surface Diagnostic - SHI ustar
if (plot_SHIuStar(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
    plotter(:,:)=var(:,:,4).*shelf_mask;
%      plotter(:,:)=var(:,:,4);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['Ice shelf - U star. Date=' datestr(startdate+ (time*dt/(3600*24)))])
end

%% Surface Diagnostic - SHI gamma T
if (plot_SHIgammT(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
    plotter(:,:)=var(:,:,5).*shelf_mask;
%     plotter(:,:)=var(:,:,5);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['Ice shelf - gamma T. Date=' datestr(startdate+ (time*dt/(3600*24)))])
end

%% Surface Diagnostic - SHI gamma S
if (plot_SHIgammS(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
    plotter(:,:)=var(:,:,6).*shelf_mask;
%     plotter(:,:)=var(:,:,6);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['Ice shelf - gamma S. Date=' datestr(startdate+ (time*dt/(3600*24)))])
end

%% Surface Diagnostic - SHI fw flux
if (plot_SHIfwFlx(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
    plotter(:,:)=var(:,:,7);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['Ice shelf - FW flux. Date=' datestr(startdate+ (time*dt/(3600*24)))])
end

%% Surface Diagnostic - SHI heat flux
if (plot_SHIhtFlx(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
    plotter(:,:)=var(:,:,8).*shelf_mask;
%     plotter(:,:)=var(:,:,8);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['Ice shelf - Heat flux. Date=' datestr(startdate+ (time*dt/(3600*24)))])
end

%% Surface Diagnostic - oceFreez
if (plot_oceFreez(1)==1)    
    clear var plotter    
    var= rdmds('surfDiag',time);    
    figure
%     plotter(:,:)=var(:,:,17).*shelf_mask;
    plotter(:,:)=var(:,:,17);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['oce Freez. Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end
%%
% Kpp diagnostic
%
%% KPP diffKzT Meridional Section Map
if (plot_KPPdiffKzT(1)==1)
    clear var plotter lat_trg depth_trg mask_ocean_l
    var= rdmds('KPPdiffKzT',time);  
    figure
    mask_ocean_l=hFacC(plot_KPPdiffKzT(2),:,:); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;  
    plotter(:,:)= var(plot_KPPdiffKzT(2),:,:).*mask_ocean_l;
    
    lat_trg(:,:)=    lat_plot(plot_KPPdiffKzT(2),:,:);
    depth_trg(:,:)=depth_plot(plot_KPPdiffKzT(2),:,:);
    pcolor(lat_trg,depth_trg,plotter); 
%      caxis([0.00001 0.001])
    colorbar
    shading flat
        hold on
%    v=[-2, -1.9, -1.8 -1.5,0,1,2,3];
%    [C,h]=contour(lat_trg,depth_trg,plotter,'Color','w');
%    clabel(C,h,v);
    
    title (['Kpp diffKzT meridional section I=' num2str(plot_KPPdiffKzT(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% KPP diffKzS Meridional Section Map
if ( plot_KPPdiffKzS(1)==1)
    clear var plotter lat_trg depth_trg mask_ocean_l
    var= rdmds('KPPdiffKzS',time);  
    figure
    mask_ocean_l=hFacC(plot_KPPdiffKzS(2),:,:); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;  
    plotter(:,:)= var(plot_KPPdiffKzS(2),:,:).*mask_ocean_l;
    
    lat_trg(:,:)=    lat_plot(plot_KPPdiffKzS(2),:,:);
    depth_trg(:,:)=depth_plot(plot_KPPdiffKzS(2),:,:);
    pcolor(lat_trg,depth_trg,plotter); 
%      caxis([0.00001 0.001])
    colorbar
    shading flat
        hold on
%    v=[-2, -1.9, -1.8 -1.5,0,1,2,3];
%    [C,h]=contour(lat_trg,depth_trg,plotter,'Color','w');
%    clabel(C,h,v);
    
    title (['Kpp diffKzs meridional section I=' num2str(plot_KPPdiffKzS(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% KPP Diagnostic - hbl 
if ( plot_KPPhbl(1)==1)    
    clear var plotter    
    var= rdmds('KPPhbl',time);  
    
        mask_ocean_l=hFacC(:,:,plot_KPPhbl(2)); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;
    
    
    figure 
       plotter(:,:)=var(:,:,plot_KPPhbl(2)).*mask_ocean_l;
    
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['KPP hbl. Date=' datestr(startdate+ (time*dt/(3600*24)))])
end


%%
% Seaice diagnostic
%
%% SeaIce Diagnostic - SI area 
if (plot_SIarea(1)==1)    
    clear var plotter    
    var= rdmds('siDiag',time);    
    figure 
plotter(:,:)=var(:,:,1).*mask_surf;
%     plotter(:,:)=var(:,:,1);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['Sea Ice Area. Date=' datestr(startdate+ (time*dt/(3600*24)))])
end

%% SeaIce Diagnostic - SI uice 
if (plot_SIuice(1)==1)    
    clear var plotter    
    var= rdmds('siDiag',time);    
    figure 
    plotter(:,:)=var(:,:,2).*mask_surf;
%     plotter(:,:)=var(:,:,2);
    pcolor(plotter'); 
%     caxis([-0.5 0.5])
%     colormap (RED)
    colorbar
    shading flat
    title (['Sea Ice Uice. Date=' datestr(startdate+ (time*dt/(3600*24)))])
end

%% SeaIce Diagnostic - SI vice 
if (plot_SIvice(1)==1)    
    clear var plotter    
    var= rdmds('siDiag',time);    
    figure 
    plotter(:,:)=var(:,:,3).*mask_surf;
%     plotter(:,:)=var(:,:,3);
    pcolor(plotter'); 
%      caxis([-0.5 0.5])
    colorbar
    shading flat
    title (['Sea Ice Vice. Date=' datestr(startdate+ (time*dt/(3600*24)))])
end

%% SeaIce Diagnostic - SI qnet 
if (plot_SIqnet(1)==1)    
    clear var plotter    
    var= rdmds('siDiag',time);    
    figure 
    plotter(:,:)=var(:,:,4).*mask_surf;
%     plotter(:,:)=var(:,:,4);
    pcolor(plotter'); 
%      caxis([-500 500])
    colorbar
    shading flat
    title (['Sea Ice: Ocean surface heat flux. >0 decreases T (freezing). Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% SeaIce Diagnostic - SI empmr
if (plot_SIempmr(1)==1)    
    clear var plotter    
    var= rdmds('siDiag',time);    
    figure 
    plotter(:,:)=var(:,:,5).*mask_surf;
%     plotter(:,:)=var(:,:,5);
    pcolor(plotter'); 
%      caxis([-5 1])
    colorbar
    shading flat
    title (['Sea Ice: Ocean surface fw flux. >0 increases Salt (freezing). Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% SeaIce Diagnostic - SI empmr
if (plot_SIatmQnt(1)==1)    
    clear var plotter    
    var= rdmds('siDiag',time);    
    figure 
    plotter(:,:)=var(:,:,6).*mask_surf;
%     plotter(:,:)=var(:,:,6);
    pcolor(plotter'); 
%     caxis([-5 1])
    colorbar
    shading flat
    title (['Sea Ice: Net atmospheric heat flux. >0 decreases T (in ice free region). Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end


%%
% Horizontal Slices
%
%% Temp Horizontal Slice
if (plot_t(1)==1)
    clear var plotter mask_ocean_l
    var= rdmds('T',time);
    figure
     
    mask_ocean_l=hFacC(:,:,plot_t(2)); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;
    plotter(:,:)=var(:,:,plot_t(2)).*mask_ocean_l;
%     plotter(:,:)=var(:,:,plot_t(2));
    pcolor(plotter'); 
%              caxis([-2 +4])
%       caxis([-2.5 2])
    colorbar
    shading flat

    title (['Temperature Depth ' num2str(-RC(1,1,plot_t(2)))  ' Date=' datestr(startdate+ (time*dt/(3600*24)))  ])
    

     
     
     
end

%% Salt Horizontal Slice
if (plot_s(1)==1)
    clear var plotter mask_ocean_l
    var= rdmds('S',time);
    figure
    mask_ocean_l=hFacC(:,:,plot_s(2)); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
mask_ocean_l(mask_ocean_l>0)=1;
    plotter(:,:)=var(:,:,plot_s(2)).*mask_ocean_l;
%      plotter(:,:)=var(:,:,plot_s(2));
    pcolor(plotter'); 
       caxis([33 35])
    colorbar
    shading flat
    title (['Salinity Level ' num2str(plot_s(2))  ' Date=' datestr(startdate+ (time*dt/(3600*24)))  ])
end

%% V Horizontal 
if (plot_v(1)==1)
    clear var plotter mask_ocean_l
    var= rdmds('V',time);
    figure
    mask_ocean_l=hFacC(:,:,plot_v(2)); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;
    plotter(:,:)=var(:,:,plot_v(2)).*mask_ocean_l;
%     plotter(:,:)=var(:,:,plot_v(2));
    pcolor(plotter'); 
%     caxis([33 35])
    colorbar
    shading flat
    title (['V (meridional) Level ' num2str(plot_v(2))  ' Date=' datestr(startdate+ (time*dt/(3600*24)))  ])
end

%% U Horizontal 
if (plot_u(1)==1)
    clear var plotter mask_ocean_l
    var= rdmds('U',time);
    figure
    mask_ocean_l=hFacC(:,:,plot_u(2)); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;
    plotter(:,:)=var(:,:,plot_u(2)).*mask_ocean_l;
%     plotter(:,:)=var(:,:,plot_u(2));
    pcolor(plotter'); 
%     caxis([33 35])
    colorbar
    shading flat
    title (['U (zonal) Level ' num2str(plot_u(2))  ' Date=' datestr(startdate+ (time*dt/(3600*24)))  ])
end

%% Velocity Horizontal 
if (plot_vel(1)==1)
    clear var1 var2 plotter1 plotter 2
    var1= rdmds('U',time);
    var2= rdmds('V',time);
    figure
    plotter1(:,:)=var1(:,:,plot_vel(2));
    plotter2(:,:)=var2(:,:,plot_vel(2));
     step=plot_vel(3); scale=plot_vel(4);
    quiver(plotter1(1:step:end,1:step:end)',plotter2(1:step:end,1:step:end)',scale); 
    xlim ([0 210]);
    ylim ([0 200]);
%     caxis([33 35])
%     colorbar
    shading flat
    title (['Velocity Level ' num2str(plot_vel(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24)))  ])
end

%%
% Maps
%
%% Temp Horizontal - Map
if (plot_t_geo(1)==1)
    clear var plotter lat_trg lon_trg mask_ocean_l lon_p lat_p
    var= rdmds('T',time);
    figure
    mask_ocean_l=hFacC(:,:,plot_t_geo(2)); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;
    plotter(:,:)=var(:,:,plot_t_geo(2)).*mask_ocean_l;
    
%     plotter(:,:)=var(:,:,plot_t_geo(2));
    lat_trg(:,:)=lat_plot(:,:,1);
    lon_trg(:,:)=lon_plot(:,:,1);
    pcolor(lon_trg,lat_trg,plotter);  
%      caxis([-2 2])
    colorbar
     shading flat
        hold on
       v=[-2, -1.9, -1.8 -1.5, -1, 0, 1,2,3,4,5,6,7,8,9];
   [C,h]=contour(lon_trg,lat_trg,plotter,v,'Color','w');
   clabel(C,h,v);
          hold on
        front=[ 32 77; 33 77; 34 77; 35 77; 36 77; 37 77; 38 77; 39 77; 40,77; ...
        41 76; 42 76; 43 76; 44 76; 45 76; 46 76; 47 76;48 76; 49 76; 50 76 ; ...
        51 75; 52 74; 53 74;...
        54 73; 55 73; 56 73; 57 73; 58 73; 59 73; 60 73; 61 73; 62 73; 63 73; 64 73; 65 73; 66 73; 67 73;...
        68 72; 69 72; 70 72; 72 70; 73 70; 73 69; ...
        74 68; 75 68; 76 68; 77 68; 78 68; 79 68; 80 68; 81 68;...
        82 67; 83 67; 84 67; 85 67; 86 67; 87 67; 88 67; 89 67; 90 67; 91 67; 92 67; 93 67; 94 67; 95 67; 96 67; ...
        97 66; 98 66; 99 66; 100 66; 101 66; 102  66; 103  66; 104  66; 105 66; 106 66; 107 66; 108 66; 109 66; 110 66; 111 66; 112 66; 113 66; ...
        114 67; 115 67; 116 67; 117 67; 118 67; 119 68; 120 68; 121 69; 122 70; 123 71; 124 72; 125 73; 126 74 ];
    points=size (front,1 );    
    for p=1:points
        lon_p(p,: )= lon_plot(front(p,1),front(p,2)-2,: );
        lat_p(p,: )= lat_plot(front(p,1),front(p,2)-2,: );
    end
    

%     plot (lon_p, lat_p,'color','k','LineWidth',2 )
    title (['Temperature Depth ' num2str(RC(1,1, plot_t_geo(2))) ' Date=' datestr(startdate+ (time*dt/(3600*24)))  ])
end

%% Salt Horizontal - Map
if (plot_s_geo(1)==1)
    clear var plotter lat_trg lon_trg mask_ocean_l
    var= rdmds('S',time);
    figure
    mask_ocean_l=hFacC(:,:,plot_s_geo(2)); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;
    plotter(:,:)=var(:,:,plot_s_geo(2)).*mask_ocean_l;
    
%     plotter(:,:)=var(:,:,plot_s_geo(2));
    lat_trg(:,:)=lat_plot(:,:,1);
    lon_trg(:,:)=lon_plot(:,:,1);
    pcolor(lon_trg,lat_trg,plotter);  
    caxis([33 35])
    colorbar
    shading flat
    
        hold on
       v=[32:0.5:34 34.2:0.2:35];
   [C,h]=contour(lon_trg,lat_trg,plotter,v,'Color','w');
   clabel(C,h,v);
    
    
    title (['Salinity Depth ' num2str(RC(1,1, plot_s_geo(2))) ' Date=' datestr(startdate+ (time*dt/(3600*24)))  ])
end

%% Velocity Horizontal Map
if (plot_vel_geo(1)==1)
    clear var1 var2 plotter1 plotter 2 mask_ocean_l
    var1= rdmds('U',time);
    var2= rdmds('V',time);
    figure
    mask_ocean_l=hFacC(:,:,plot_vel_geo(2)); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;
    plotter1(:,:)=var1(:,:,plot_vel_geo(2)).*mask_ocean_l;
    plotter2(:,:)=var2(:,:,plot_vel_geo(2)).*mask_ocean_l;    
    
%     plotter1(:,:)=var1(:,:,plot_vel_geo(2));
%     plotter2(:,:)=var2(:,:,plot_vel_geo(2));
    lat_trg(:,:)=lat_plot(:,:,1);
    lon_trg(:,:)=lon_plot(:,:,1);
    step=plot_vel_geo(3); scale=plot_vel_geo(4);
    quiver(lon_trg(1:step:end,1:step:end),lat_trg(1:step:end,1:step:end), plotter1(1:step:end,1:step:end) ,plotter2(1:step:end,1:step:end),scale); 
    xlim ([160 230]);
    ylim ([-85 -65]);
    shading flat
    title (['Velocity Level ' num2str(plot_vel_geo(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24)))  ])
end


%% Depth - Map
if (plot_depth_geo(1)==1)
    clear var plotter lat_trg lon_trg
    var=rdmds('Depth'); 
    figure
    plotter(:,:)=var(:,:);
    lat_trg(:,:)=lat_plot(:,:,1);
    lon_trg(:,:)=lon_plot(:,:,1);
    pcolor(lon_trg,lat_trg,plotter);  
    colorbar
    shading flat
    title ('Depth')
end
%%
% Sections
%
%%
% Meridional Sections
%
%% Depth with Meridional Section
if (plot_depth_sec_m(1)==1)
    clear var plotter
    var=rdmds('Depth');
    figure
    plotter(:,:)=var(:,:);
    pcolor(plotter'); 
    colorbar
%     caxis ([0 50])
    shading flat 
    hold on
    plot ([plot_depth_sec_m(2) plot_depth_sec_m(2)], [1 ny],'color','w' )
    title (['Position of meridional section I=' num2str(plot_depth_sec_m(2))  ])
end
%% Temp Meridional Section Map
if (plot_t_sec_m_geo(1)==1)
    clear var plotter lat_trg depth_trg mask_ocean_l
    var= rdmds('T',time);
    figure
    mask_ocean_l=hFacC(plot_t_sec_m_geo(2),:,:); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;  
    plotter(:,:)= var(plot_t_sec_m_geo(2),:,:).*mask_ocean_l;
    
%     plotter(:,:)=         var(plot_t_sec_m_geo(2),:,:);
    lat_trg(:,:)=    lat_plot(plot_t_sec_m_geo(2),:,:);
    depth_trg(:,:)=depth_plot(plot_t_sec_m_geo(2),:,:);
    pcolor(lat_trg,depth_trg,plotter); 
%     caxis([-2 1])
    colorbar
    shading flat
        hold on
   v=[-2, -1.9, -1.8 -1.5,0,1,2,3];
   [C,h]=contour(lat_trg,depth_trg,plotter,v,'Color','w');
   clabel(C,h,v);
    
    title (['Temperature meridional section I=' num2str(plot_t_sec_m_geo(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Salt Meridional Section Map
if (plot_s_sec_m_geo(1)==1)
    clear var plotter lat_trg depth_trg mask_ocean_l
    var= rdmds('S',time);
    figure
    mask_ocean_l=hFacC(plot_s_sec_m_geo(2),:,:); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;  
    plotter(:,:)=var(plot_s_sec_m_geo(2),:,:).*mask_ocean_l;
    
%     plotter(:,:)=         var(plot_s_sec_m_geo(2),:,:);
    lat_trg(:,:)=    lat_plot(plot_s_sec_m_geo(2),:,:);
    depth_trg(:,:)=depth_plot(plot_s_sec_m_geo(2),:,:);
    pcolor(lat_trg,depth_trg,plotter); 
%     caxis([33 35])
    colorbar
    shading flat
    title (['Salinity meridional section I=' num2str(plot_s_sec_m_geo(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% U Meridional Section Map
if (plot_u_sec_m_geo(1)==1)
    clear var plotter lat_trg depth_trg mask_ocean_l
    var= rdmds('U',time);
    figure
    mask_ocean_l=hFacC(plot_u_sec_m_geo(2),:,:); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;  
    plotter(:,:)= var(plot_u_sec_m_geo(2),:,:).*mask_ocean_l;    
    
%     plotter(:,:)=         var(plot_u_sec_m_geo(2),:,:);
    lat_trg(:,:)=    lat_plot(plot_u_sec_m_geo(2),:,:);
    depth_trg(:,:)=depth_plot(plot_u_sec_m_geo(2),:,:);
    pcolor(lat_trg,depth_trg,plotter); 
%     caxis([33 35])
    colorbar
    shading flat
    title (['U meridional section I=' num2str(plot_u_sec_geo(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Temp Meridional Section
if (plot_t_sec_m(1)==1)
    clear var plotter 
    var= rdmds('T',time);
    figure
    plotter(:,:)=var(plot_t_sec_m(2),:,:);
    pcolor(plotter'); 
%     caxis([-2.5 1])
%   ylim ([-500 0]);
    colorbar
    shading flat
    set(gca,'YDir', 'Reverse');
  
    
    title (['Temperature meridional section I=' num2str(plot_t_sec_m(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Salt Meridional Section
if (plot_s_sec_m(1)==1)
    clear var plotter
    var= rdmds('S',time);
    figure 
    plotter(:,:)=var(plot_s_sec_m(2),:,:);
    pcolor(plotter'); 
    caxis([33 35])
    colorbar
    shading flat
    set(gca,'YDir', 'Reverse');
    title (['Salinity meridional section I=' num2str(plot_s_sec_m_geo(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% U Meridional Section
if (plot_u_sec_m(1)==1)
    clear var plotter
    var= rdmds('U',time);
    figure 
    plotter(:,:)=var(plot_u_sec_m(2),:,:);
    pcolor(plotter'); 
%     caxis([33 35])
    colorbar
    shading flat
    set(gca,'YDir', 'Reverse');
    title (['U meridional section I=' num2str(plot_u_sec_geo(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end



%%
% Zonal Sections
%
%% Depth with Zonal Section
if (plot_depth_sec_z(1)==1)
    clear var plotter
    var=rdmds('Depth');
    figure
    plotter(:,:)=var(:,:);
    pcolor(plotter'); 
    colorbar
%     caxis ([0 50])
    shading flat 
    hold on
    plot ([1 nx],[plot_depth_sec_z(2) plot_depth_sec_z(2)],'color','w' )
    title (['Position of zonal section J=' num2str(plot_depth_sec_z(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% T Zonal Section Map
if (plot_t_sec_z_geo(1)==1)
    clear var plotter lon_trg depth_trg
    var= rdmds('T',time);
    figure
    plotter(:,:)=         var(:,plot_t_sec_z_geo(2),:);
    lon_trg(:,:)=    lon_plot(:,plot_t_sec_z_geo(2),:);
    depth_trg(:,:)=depth_plot(:,plot_t_sec_z_geo(2),:);
    pcolor(lon_trg,depth_trg,plotter); 
      caxis([-2.3 -1])
 ylim ([-500 0])
    colorbar
    shading flat
            hold on
   v=[-2, -1.9, -1.8 -1.5];
   [C,h]=contour(lon_trg,depth_trg,plotter,v,'Color','w');
   clabel(C,h,v);
    title (['Temperature zonal section J=' num2str(plot_t_sec_z_geo(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Salt Zonal Section Map
if (plot_s_sec_z_geo(1)==1)
    clear var plotter lon_trg depth_trg
    var= rdmds('S',time);
    figure
    plotter(:,:)=         var(:,plot_s_sec_z_geo(2),:);
    lon_trg(:,:)=    lon_plot(:,plot_s_sec_z_geo(2),:);
    depth_trg(:,:)=depth_plot(:,plot_s_sec_z_geo(2),:);
    pcolor(lon_trg,depth_trg,plotter); 
    caxis([33 35])
    colorbar
    shading flat
    title (['Salinity zonal section J=' num2str(plot_s_sec_z_geo(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% V Zonal Section Map

if (plot_v_sec_z_geo(1)==1)
    clear var plotter lon_trg depth_trg
    var= rdmds('V',time);
    figure
    plotter(:,:)=         var(:,plot_v_sec_z_geo(2),:);
    lon_trg(:,:)=    lon_plot(:,plot_v_sec_z_geo(2),:);
    depth_trg(:,:)=depth_plot(:,plot_v_sec_z_geo(2),:);
    pcolor(lon_trg,depth_trg,plotter); 
%     caxis([33 35])
    colorbar
    shading flat
    title (['V zonal section J=' num2str(plot_v_sec_z_geo(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ])
end

%% Temp Zonal Section
if (plot_t_sec_z(1)==1)
    clear var plotter
    var= rdmds('T',time);
    figure
    plotter(:,:)=var(:,plot_t_sec_z(2),:);
    pcolor(plotter'); 
      caxis([-2.5 -1])
%       xlim([ 63 93])
    colorbar
    shading flat
    set(gca,'YDir', 'Reverse');
    title (['Temperature section at J =' num2str(plot_t_sec_z(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24))) ] )
end

%% Salt Zonal Section
if (plot_s_sec_z(1)==1)
    clear var plotter
    var= rdmds('S',time);
    figure 
    plotter(:,:)=var(:,plot_s_sec_z(2),:);
    pcolor(plotter'); 
    caxis([33 35])
    colorbar
    shading flat
    set(gca,'YDir', 'Reverse');
     title (['Salinity section at J ' num2str(plot_s_sec_z(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24)))] )
end


%% V Zonal Section
if (plot_v_sec_z(1)==1)
    clear var plotter
    var= rdmds('V',time);
    figure 
    plotter(:,:)=var(:,plot_v_sec_z(2),:);
    pcolor(plotter'); 
%     caxis([33 35])
    colorbar
    shading flat
    set(gca,'YDir', 'Reverse');
     title (['V velocity section at J ' num2str(plot_v_sec_z(2)) ' Date=' datestr(startdate+ (time*dt/(3600*24)))] )
end

%%
% Time series
%
%% Temp Time Series
if (plot_time_t(1,1)==1)
    clear var plotter name z time_stamp
    %Points loop
    points=size (plot_time_t,1 )-1;
    for p=1:points
       % Time loop 
        for t=plot_time_t(1,2):plot_time_t(1,3)
            time1=time_options(t);
                var= rdmds('T',time1);
                plotter(p,t)=var(plot_time_t(1+p,1),plot_time_t(1+p,2),plot_time_t(1+p,3));
                time_stamp(t,:)=(startdate+ (time1*dt/(3600*24)));
        end
    end
    %
    figure
    for p=1:points
        hold on
    z(p)=plot (time_stamp,plotter(p,:));
    set(z(p),'Color',colorOrder(p,:));
    name(p,1:(length([num2str(plot_time_t(1+p,1)) ', '  num2str(plot_time_t(1+p,2)) ', ' num2str(plot_time_t(1+p,3))])))= [num2str(plot_time_t(1+p,1)) ', '  num2str(plot_time_t(1+p,2)) ', ' num2str(plot_time_t(1+p,3))];
    end
    set(gca, 'Xlim', [time_stamp(plot_time_t(1,2),:) time_stamp(plot_time_t(1,3),:)])  
    legend ( [z(1:p)], [name(1:p,:)])
    set(gca,'XTick', [ datenum(1992,1,1) datenum(1992,2,1) datenum(1992,3,1) datenum(1992,4,1) datenum(1992,5,1) ...
    datenum(1992,6,1) datenum(1992,7,1) datenum(1992,8,1) datenum(1992,9,1) datenum(1992,10,1) datenum(1992,11,1) datenum(1992,12,1)...
    datenum(1993,1,1) datenum(1993,2,1) datenum(1993,3,1) datenum(1993,4,1) datenum(1993,5,1) datenum(1993,6,1) ...
    datenum(1993,7,1) datenum(1993,8,1) datenum(1993,9,1) datenum(1993,10,1) datenum(1993,11,1) datenum(1993,12,1) ...
    datenum(1994,1,1) datenum(1994,2,1) datenum(1994,3,1) datenum(1994,4,1) datenum(1994,5,1) datenum(1994,6,1) ...
    datenum(1994,7,1) datenum(1994,8,1) datenum(1994,9,1) datenum(1994,10,1) datenum(1994,11,1) datenum(1994,12,1) ...
    datenum(1995,1,1) datenum(1995,2,1) datenum(1995,3,1) datenum(1995,4,1) datenum(1995,5,1) datenum(1995,6,1) ...
    datenum(1995,7,1) datenum(1995,8,1) datenum(1995,9,1) datenum(1995,10,1) datenum(1995,11,1) datenum(1995,12,1) ...
    ] );
    datetick('x', 'mm/yy', 'keeplimits', 'keepticks')  %dd/mm
    title (['Temperature Time Series'   ])
    ylim([-650 0])
end
%% Salt Time Series
if (plot_time_s(1,1)==1)
    clear var plotter name z time_stamp
    %Points loop
    points=size (plot_time_s,1 )-1;
    for p=1:points
       % Time loop 
        for t=plot_time_s(1,2):plot_time_s(1,3)
            time1=time_options(t);
                var= rdmds('S',time1);
                plotter(p,t)=var(plot_time_s(1+p,1),plot_time_s(1+p,2),plot_time_s(1+p,3));
                time_stamp(t,:)=(startdate+ (time1*dt/(3600*24)));
        end
    end
    %
    figure
    for p=1:points
        hold on
    z(p)=plot (time_stamp,plotter(p,:));
    set(z(p),'Color',colorOrder(p,:));
     name(p,1:length([num2str(plot_time_s(1+p,1)) ', '  num2str(plot_time_s(1+p,2)) ', ' num2str(plot_time_s(1+p,3))]))= [num2str(plot_time_s(1+p,1)) ', '  num2str(plot_time_s(1+p,2)) ', ' num2str(plot_time_s(1+p,3))];
    end
     legend ( [z(1:p)], [name(1:p,:)])
      set(gca, 'Xlim', [time_stamp(plot_time_s(1,2),:) time_stamp(plot_time_s(1,3),:)]) 
    set(gca,'XTick', [ datenum(1992,1,1)  datenum(1992,4,1) datenum(1992,7,1)  datenum(1992,10,1) ...
    datenum(1993,1,1)  datenum(1993,4,1) datenum(1993,7,1) datenum(1993,10,1) ...
    datenum(1994,1,1) datenum(1994,4,1)  datenum(1994,7,1)  datenum(1994,10,1) ...
    datenum(1995,1,1)  datenum(1995,4,1)  datenum(1995,7,1) datenum(1995,10,1) ...
    datenum(1996,1,1)  datenum(1996,4,1)  datenum(1996,7,1) datenum(1996,10,1) ...
    datenum(1997,1,1)  datenum(1997,4,1)  datenum(1997,7,1) datenum(1997,10,1) ...
    datenum(1998,1,1)  datenum(1998,4,1)  datenum(1998,7,1) datenum(1998,10,1) ...
    datenum(1999,1,1)  datenum(1999,4,1)  datenum(1999,7,1) datenum(1999,10,1) ...
    datenum(2000,1,1)  datenum(2000,4,1)  datenum(2000,7,1) datenum(2000,10,1) ...
    ] );
    datetick('x', 'mm/yy', 'keeplimits', 'keepticks')  %dd/mm
    title (['Salinity Time Series'   ])
end
%%

%% Temp Time Series of Vert profiles
if (plot_tz_t(1,1)==1)
    clear var plotter name z depth 2 profile2
       % Time loop 
        for t=plot_tz_t(1,2):plot_tz_t(1,3)
            time1=time_options(t);
                var= rdmds('T',time1);
                plotter(t,:)=var(plot_tz_t(2,1),plot_tz_t(2,2),:);
                time_stamp(t,:)=(startdate+ (time1*dt/(3600*24)));
        end
   %
    mask_ocean_z(:)=hFacC(plot_tz_t(2,1),plot_tz_t(2,2),:); 
    mask_ocean_z(mask_ocean_z==0)=NaN;
    mask_ocean_z(mask_ocean_z>0)=1;
    %
    figure
    for t=plot_tz_t(1,2):plot_tz_t(1,3)
    hold on
    depth2(:)=depth_plot (plot_tz_t(2,1),plot_tz_t(2,2),:);
    profile2(:)=plotter(t,:).*mask_ocean_z(1,:);
    z(t)=plot (profile2,depth2);
    name(t,:)= [datestr(time_stamp(t,:) )  ];
    end
    legend ( [z(1:t)], [name(1:t,:)])
    title (['Temperature Profiles. Point ' num2str(plot_tz_t(2,1) ) ', ' num2str(plot_tz_t(2,2) ) ' in time'   ])
end  

 %% Temp Time Series of Vert profiles (with offset)
if (plot_tz_t_offset(1,1)==1)
    clear var plotter name z depth 2 profile2 mask_ocean_z
       % Time loop 
        for t=plot_tz_t_offset(1,2):plot_tz_t_offset(1,3)
            time1=time_options(t);
                var= rdmds('T',time1);
                plotter(t,:)=var(plot_tz_t_offset(2,1),plot_tz_t_offset(2,2),:);
                time_stamp(t,:)=(startdate+ (time1*dt/(3600*24)));
        end
   %
    mask_ocean_z(:)=hFacC(plot_tz_t_offset(2,1),plot_tz_t_offset(2,2),:); 
    mask_ocean_z(mask_ocean_z==0)=NaN;
    mask_ocean_z(mask_ocean_z>0)=1;
    %
    figure
    for t=plot_tz_t_offset(1,2):plot_tz_t_offset(1,3)
    hold on
    depth2(:)=depth_plot (plot_tz_t_offset(2,1),plot_tz_t_offset(2,2),:);
    profile2(:)=plotter(t,:).*mask_ocean_z(1,:)+0.1*t;
    z(t)=plot (profile2,depth2);
    name(t,:)= [datestr(time_stamp(t,:) )  ];
    end
    legend ( [z(1:t)], [name(1:t,:)])
    title (['Temperature Profiles. Point ' num2str(plot_tz_t_offset(2,1) ) ', ' num2str(plot_tz_t_offset(2,2) ) ' in time'   ])
end  
    
    
    %% Temp Time Series of Vert profiles (at 1 time)
if (plot_tz_loc_t(1,1)==1)
    clear var plotter name z depth 2 profile2
    %Points loop
    points=size (plot_tz_loc_t,1 )-1;
    time1=time_options(plot_tz_loc_t(1,2)   );
    for p=1:points
                var= rdmds('T',time1);
                plotter(p,:)=var(plot_tz_loc_t(1+p,1),plot_tz_loc_t(1+p,2),:);
    end
    %

    figure
    for p=1:points
        %
        mask_ocean_z(:)=hFacC(plot_tz_loc_t(1+p,1),plot_tz_loc_t(1+p,2),:); 
        mask_ocean_z(mask_ocean_z==0)=NaN;
        mask_ocean_z(mask_ocean_z>0)=1;
        %    
        hold on

    depth2(:)=depth_plot (plot_tz_loc_t(1+p,1),plot_tz_loc_t(1+p,2),:);
    profile2(:)=plotter(p,:).*mask_ocean_z;
    z(p)=plot (profile2,depth2);
    set(z(p),'Color',colorOrder(p,:));
    name(p,:)= [num2str(plot_tz_loc_t(1+p,1)) ', '  num2str(plot_tz_loc_t(1+p,2)) ];
    end
   
    legend ( [z(1:p)], [name(1:p,:)])
   
    title (['Temperature Profiles at time ' datestr(startdate+ (time1*dt/(3600*24)))  ])
end


%% Temp Frontal Section (inside cavity)
if (plot_t_sec_front(1)==1)
    clear var plotter points
    var= rdmds('T',time);
    front=[13 80; 14 80;15 80;16 80;17 80;18 80;19 80;20 80;21 80;22 80; ...
        23 79; 24 79; 25 79; 26 79; ...
        27 78; 28 78; 29 78; ...
        30 77; 31 77; 32 77; 33 77; 34 77; 35 77; 36 77; 37 77; 38 77; 39 77; 40,77; ...
        41 76; 42 76; 43 76; 44 76; 45 76; 46 76; 47 76;48 76; 49 76; 50 76 ; ...
        51 75; 52 74; 53 74;...
        54 73; 55 73; 56 73; 57 73; 58 73; 59 73; 60 73; 61 73; 62 73; 63 73; 64 73; 65 73; 66 73; 67 73;...
        68 72; 69 72; 70 72; 72 70; 73 70; 73 69; ...
        74 68; 75 68; 76 68; 77 68; 78 68; 79 68; 80 68; 81 68;...
        82 67; 83 67; 84 67; 85 67; 86 67; 87 67; 88 67; 89 67; 90 67; 91 67; 92 67; 93 67; 94 67; 95 67; 96 67; ...
        97 66; 98 66; 99 66; 100 66; 101 66; 102  66; 103  66; 104  66; 105 66; 106 66; 107 66; 108 66; 109 66; 110 66; 111 66; 112 66; 113 66; ...
        114 67; 115 67; 116 67; 117 67; 118 67; 119 68; 120 68; 121 69; 122 70; 123 71; 124 72; 125 73; 126 74 ];
    points=size (front,1 );    
    %
    for p=1:points
        mask_ocean_l=hFacC(front(p,1),front(p,2)-3,:,:); 
        mask_ocean_l(mask_ocean_l==0)=NaN;
        mask_ocean_l(mask_ocean_l>0)=1;  
        lon_p(p,: )= lon_plot(front(p,1),front(p,2)-3,: );
        depth_p(p,:)= depth_plot(front(p,1),front(p,2)-3,: );
        temp_p(p,:)= var(front(p,1),front(p,2)-3,: ).*mask_ocean_l;
    end
   %
   figure
   pcolor(lon_p,depth_p,temp_p); 
   hold on
   v=[-2, -1.9, -1.8 -1.5];
   [C,h]=contour(lon_p,depth_p,temp_p,v,'Color','w');
   clabel(C,h,v);
%        caxis([-2.5 -1])
   ylim([-1000 0])
   colorbar
   shading flat
   title (['Temperature section along ice shelf front' ' Date=' datestr(startdate+ (time*dt/(3600*24))) ] )
end

%%

%% Temp Frontal Section (outside cavity)
if (plot_t_sec_front(1)==1)
    clear var plotter points
    var= rdmds('T',time);
    front=[13 80; 14 80;15 80;16 80;17 80;18 80;19 80;20 80;21 80;22 80; ...
        23 79; 24 79; 25 79; 26 79; ...
        27 78; 28 78; 29 78; ...
        30 77; 31 77; 32 77; 33 77; 34 77; 35 77; 36 77; 37 77; 38 77; 39 77; 40,77; ...
        41 76; 42 76; 43 76; 44 76; 45 76; 46 76; 47 76;48 76; 49 76; 50 76 ; ...
        51 75; 52 74; 53 74;...
        54 73; 55 73; 56 73; 57 73; 58 73; 59 73; 60 73; 61 73; 62 73; 63 73; 64 73; 65 73; 66 73; 67 73;...
        68 72; 69 72; 70 72; 72 70; 73 70; 73 69; ...
        74 68; 75 68; 76 68; 77 68; 78 68; 79 68; 80 68; 81 68;...
        82 67; 83 67; 84 67; 85 67; 86 67; 87 67; 88 67; 89 67; 90 67; 91 67; 92 67; 93 67; 94 67; 95 67; 96 67; ...
        97 66; 98 66; 99 66; 100 66; 101 66; 102  66; 103  66; 104  66; 105 66; 106 66; 107 66; 108 66; 109 66; 110 66; 111 66; 112 66; 113 66; ...
        114 67; 115 67; 116 67; 117 67; 118 67; 119 68; 120 68; 121 69; 122 70; 123 71; 124 72; 125 73; 126 74 ];
    points=size (front,1 );    
    %
    for p=1:points
        mask_ocean_l=hFacC(front(p,1),front(p,2),:,:); 
        mask_ocean_l(mask_ocean_l==0)=NaN;
        mask_ocean_l(mask_ocean_l>0)=1;  
        lon_p(p,: )= lon_plot(front(p,1),front(p,2),: );
        depth_p(p,:)= depth_plot(front(p,1),front(p,2),: );
        temp_p(p,:)= var(front(p,1),front(p,2),: ).*mask_ocean_l;
    end
   %
   figure
   pcolor(lon_p,depth_p,temp_p); 
%     caxis ([-2.2 0])
   
   colorbar
   hold on
   v=[-2, -1.9, -1.8 -1.5:0.5:1];
   [C,h]=contour(lon_p,depth_p,temp_p,v,'Color','w');
   clabel(C,h,v);
% %       caxis([-2.5 0])
    ylim([-1000 0])
  

   shading flat
   title (['Temperature section along ice shelf front' ' Date=' datestr(startdate+ (time*dt/(3600*24))) ] )
end

%% Salt Frontal Section
if (plot_s_sec_front(1)==1)
    clear var plotter points
    var= rdmds('S',time);
    front=[13 80; 14 80;15 80;16 80;17 80;18 80;19 80;20 80;21 80;22 80; ...
        23 79; 24 79; 25 79; 26 79; ...
        27 78; 28 78; 29 78; ...
        30 77; 31 77; 32 77; 33 77; 34 77; 35 77; 36 77; 37 77; 38 77; 39 77; 40,77; ...
        41 76; 42 76; 43 76; 44 76; 45 76; 46 76; 47 76;48 76; 49 76; 50 76 ; ...
        51 75; 52 74; 53 74;...
        54 73; 55 73; 56 73; 57 73; 58 73; 59 73; 60 73; 61 73; 62 73; 63 73; 64 73; 65 73; 66 73; 67 73;...
        68 72; 69 72; 70 72; 72 70; 73 70; 73 69; ...
        74 68; 75 68; 76 68; 77 68; 78 68; 79 68; 80 68; 81 68;...
        82 67; 83 67; 84 67; 85 67; 86 67; 87 67; 88 67; 89 67; 90 67; 91 67; 92 67; 93 67; 94 67; 95 67; 96 67; ...
        97 66; 98 66; 99 66; 100 66; 101 66; 102  66; 103  66; 104  66; 105 66; 106 66; 107 66; 108 66; 109 66; 110 66; 111 66; 112 66; 113 66; ...
        114 67; 115 67; 116 67; 117 67; 118 67; 119 68; 120 68; 121 69; 122 70; 123 71; 124 72; 125 73; 126 74 ];
    points=size (front,1 );    
    %
    for p=1:points
        mask_ocean_l=hFacC(front(p,1),front(p,2),:,:); 
        mask_ocean_l(mask_ocean_l==0)=NaN;
        mask_ocean_l(mask_ocean_l>0)=1;  
        lon_p(p,: )= lon_plot(front(p,1),front(p,2),: );
        depth_p(p,:)= depth_plot(front(p,1),front(p,2),: );
        temp_p(p,:)= var(front(p,1),front(p,2),: ).*mask_ocean_l;
    end
   %
   figure
   pcolor(lon_p,depth_p,temp_p); 
   hold on
    v=[33.8, 34 ,34.1,34.2, 34.3,34.4,34.5, 34.6,34.7,34.8,34.9,35];
   [C,h]=contour(lon_p,depth_p,temp_p,v,'Color','w');
   clabel(C,h,v);
%       caxis([-2.5 0])
   ylim([-1000 0])
   colorbar
   shading flat
   title (['Salinity section along ice shelf front' ' Date=' datestr(startdate+ (time*dt/(3600*24))) ] )
end

%% Depth with Frontal Section
if (plot_depth_sec_z(1)==1)
    clear var plotter lat_trg lon_trg
    var=rdmds('Depth');
    figure
    plotter(:,:)=var(:,:).*mask_surf;
    lat_trg(:,:)=lat_plot(:,:,1);
    lon_trg(:,:)=lon_plot(:,:,1);
    pcolor(lon_trg,lat_trg,plotter);  
    caxis([ 0 1000])
    colorbar
    shading flat 
    hold on
        front=[13 80; 14 80;15 80;16 80;17 80;18 80;19 80;20 80;21 80;22 80; ...
        23 79; 24 79; 25 79; 26 79; ...
        27 78; 28 78; 29 78; ...
        30 77; 31 77; 32 77; 33 77; 34 77; 35 77; 36 77; 37 77; 38 77; 39 77; 40,77; ...
        41 76; 42 76; 43 76; 44 76; 45 76; 46 76; 47 76;48 76; 49 76; 50 76 ; ...
        51 75; 52 74; 53 74;...
        54 73; 55 73; 56 73; 57 73; 58 73; 59 73; 60 73; 61 73; 62 73; 63 73; 64 73; 65 73; 66 73; 67 73;...
        68 72; 69 72; 70 72; 72 70; 73 70; 73 69; ...
        74 68; 75 68; 76 68; 77 68; 78 68; 79 68; 80 68; 81 68;...
        82 67; 83 67; 84 67; 85 67; 86 67; 87 67; 88 67; 89 67; 90 67; 91 67; 92 67; 93 67; 94 67; 95 67; 96 67; ...
        97 66; 98 66; 99 66; 100 66; 101 66; 102  66; 103  66; 104  66; 105 66; 106 66; 107 66; 108 66; 109 66; 110 66; 111 66; 112 66; 113 66; ...
        114 67; 115 67; 116 67; 117 67; 118 67; 119 68; 120 68; 121 69; 122 70; 123 71; 124 72; 125 73; 126 74 ];
    points=size (front,1 );    
    for p=1:points
        lon_p(p,: )= lon_plot(front(p,1),front(p,2),: );
        lat_p(p,: )= lat_plot(front(p,1),front(p,2),: );
    end
    

    plot (lon_p, lat_p,'color','w' )
    title (['Position of frontal section ' ])
end

%% Depth with Frontal Section 2
if (plot_depth_sec_z(1)==1)
    clear var plotter lat_trg lon_trg
    var=rdmds('Depth');
    figure
    plotter(:,:)=var(:,:).*mask_surf;
    lat_trg(:,:)=lat_plot(:,:,1);
    lon_trg(:,:)=lon_plot(:,:,1);
    pcolor(lon_trg,lat_trg,plotter);  
    caxis([ 0 1000])
    colorbar
    shading flat 
    hold on
        front=[13 80; 14 80;15 80;16 80;17 80;18 80;19 80;20 80;21 80;22 80; ...
        23 79; 24 79; 25 79; 26 79; ...
        27 78; 28 78; 29 78; ...
        30 77; 31 77; 32 77; 33 77; 34 77; 35 77; 36 77; 37 77; 38 77; 39 77; 40,77; ...
        41 76; 42 76; 43 76; 44 76; 45 76; 46 76; 47 76;48 76; 49 76; 50 76 ; ...
        51 75; 52 74; 53 74;...
        54 73; 55 73; 56 73; 57 73; 58 73; 59 73; 60 73; 61 73; 62 73; 63 73; 64 73; 65 73; 66 73; 67 73;...
        68 72; 69 72; 70 72; 72 70; 73 70; 73 69; ...
        74 68; 75 68; 76 68; 77 68; 78 68; 79 68; 80 68; 81 68;...
        82 67; 83 67; 84 67; 85 67; 86 67; 87 67; 88 67; 89 67; 90 67; 91 67; 92 67; 93 67; 94 67; 95 67; 96 67; ...
        97 66; 98 66; 99 66; 100 66; 101 66; 102  66; 103  66; 104  66; 105 66; 106 66; 107 66; 108 66; 109 66; 110 66; 111 66; 112 66; 113 66; ...
        114 67; 115 67; 116 67; 117 67; 118 67; 119 68; 120 68; 121 69; 122 70; 123 71; 124 72; 125 73; 126 74 ];
    points=size (front,1 );    
    for p=1:points
        lon_p(p,: )= lon_plot(front(p,1),front(p,2)-3,: );
        lat_p(p,: )= lat_plot(front(p,1),front(p,2)-3,: );
    end
    

    plot (lon_p, lat_p,'color','w' )
    title (['Position of frontal section ' ])
end

%% Temp Time Series
if (plot_time_t_ecco(1,1)==1)
    clear var plotter name z time_stamp points
    %
    file_name=['THETA.0010.nc'];
    file_name2=['SIarea.0010.nc'];
    temp=ncread (file_name,'THETA',[1 1 1 1], [inf inf inf 312] );
    a=ncread (file_name2,'SIarea',[1 1  1], [inf inf 312] );
    p_ecco= [141 63 1; 145 68 1; 146 66 1; 144 99 1; ...
             121 92 1; 152 135 1; 99 89 1; 136 63 1];
    points=size(p_ecco,1);
    for p=1:points
        for t=1:312
            plotter_ecco(p,t)=temp(p_ecco(p,1),p_ecco(p,2),p_ecco(p,3),t);
            ice_ecco(p,t)=a(p_ecco(p,1),p_ecco(p,2),t);
            time_ecco(t,:)=(datenum(1992,01,01)+t*2629800/(3600*24));
        end
    end
    
    %Points loop
    p_mit= [14 79 1; 19 74 1; 14 76 1; 52 77 1; ...
                42 107 1; 123 77 1; 38 137 1; 13 85 1];
    for p=1:points
       % Time loop 
        for t=plot_time_t_ecco(1,2):plot_time_t_ecco(1,3)
            time1=time_options(t);
                var= rdmds('T',time1);
                plotter(p,t)=var(p_mit(p,1),p_mit(p,2),p_mit(p,3));
                var1=rdmds('AREA',time1);
                plotter_a(p,t)=var1(p_mit(p,1),p_mit(p,2));
                time_stamp(t,:)=(startdate+ (time1*dt/(3600*24)));
        end
    end
%    
    figure
   subplot('Position',[0.1,0.6,0.8,0.3])
    for p=1:points
    hold on
    z(p)=plot (time_stamp,plotter(p,:));
    set(z(p),'Color',colorOrder(p,:),'LineWidth',1);
    
    end
    
    for p=1:points
    hold on
    z1(p)=plot (time_ecco,plotter_ecco(p,:));
    set(z1(p),'Color',colorOrder(p,:),'LineStyle','--');
    end
    
    
    set(gca, 'Xlim', [time_stamp(plot_time_t_ecco(1,2),:) time_stamp(plot_time_t_ecco(1,3),:)])  
    legend ( [z(1:p)], 'McMurdo center', 'McMurdo East', 'McMurdo West', '176E','Shelf','East Ross','Cape Adare','Vic Current')
    set(gca,'XTick', [ datenum(1992,1,1) datenum(1992,2,1) datenum(1992,3,1) datenum(1992,4,1) datenum(1992,5,1) ...
    datenum(1992,6,1) datenum(1992,7,1) datenum(1992,8,1) datenum(1992,9,1) datenum(1992,10,1) datenum(1992,11,1) datenum(1992,12,1)...
    datenum(1993,1,1) datenum(1993,2,1) datenum(1993,3,1) datenum(1993,4,1) datenum(1993,5,1) datenum(1993,6,1) ...
    datenum(1993,7,1) datenum(1993,8,1) datenum(1993,9,1) datenum(1993,10,1) datenum(1993,11,1) datenum(1993,12,1) ...
    datenum(1994,1,1) datenum(1994,2,1) datenum(1994,3,1) datenum(1994,4,1) datenum(1994,5,1) datenum(1994,6,1) ...
    datenum(1994,7,1) datenum(1994,8,1) datenum(1994,9,1) datenum(1994,10,1) datenum(1994,11,1) datenum(1994,12,1) ...
    datenum(1995,1,1) datenum(1995,2,1) datenum(1995,3,1) datenum(1995,4,1) datenum(1995,5,1) datenum(1995,6,1) ...
    datenum(1995,7,1) datenum(1995,8,1) datenum(1995,9,1) datenum(1995,10,1) datenum(1995,11,1) datenum(1995,12,1) ...
    ] );
    datetick('x', 'mm/yy', 'keeplimits', 'keepticks')  %dd/mm
    title (['Temperature Time Series with ECCO (dashed)'   ])
    
    
    subplot('Position',[0.1,0.2,0.8,0.3])
    for p=1:points
    hold on
    z(p)=plot (time_stamp,plotter_a(p,:));
    set(z(p),'Color',colorOrder(p,:),'LineWidth',1);
    end
    
    for p=1:points
    hold on
    z1(p)=plot (time_ecco,ice_ecco(p,:));
    set(z1(p),'Color',colorOrder(p,:),'LineStyle','--');
    end
    
    set(gca, 'Xlim', [time_stamp(plot_time_t_ecco(1,2),:) time_stamp(plot_time_t_ecco(1,3),:)])  
    legend ( [z(1:p)], 'McMurdo center', 'McMurdo East', 'McMurdo West', '176E','Shelf','East Ross','Cape Adare','Vic Current')
    set(gca,'XTick', [ datenum(1992,1,1) datenum(1992,2,1) datenum(1992,3,1) datenum(1992,4,1) datenum(1992,5,1) ...
    datenum(1992,6,1) datenum(1992,7,1) datenum(1992,8,1) datenum(1992,9,1) datenum(1992,10,1) datenum(1992,11,1) datenum(1992,12,1)...
    datenum(1993,1,1) datenum(1993,2,1) datenum(1993,3,1) datenum(1993,4,1) datenum(1993,5,1) datenum(1993,6,1) ...
    datenum(1993,7,1) datenum(1993,8,1) datenum(1993,9,1) datenum(1993,10,1) datenum(1993,11,1) datenum(1993,12,1) ...
    datenum(1994,1,1) datenum(1994,2,1) datenum(1994,3,1) datenum(1994,4,1) datenum(1994,5,1) datenum(1994,6,1) ...
    datenum(1994,7,1) datenum(1994,8,1) datenum(1994,9,1) datenum(1994,10,1) datenum(1994,11,1) datenum(1994,12,1) ...
    datenum(1995,1,1) datenum(1995,2,1) datenum(1995,3,1) datenum(1995,4,1) datenum(1995,5,1) datenum(1995,6,1) ...
    datenum(1995,7,1) datenum(1995,8,1) datenum(1995,9,1) datenum(1995,10,1) datenum(1995,11,1) datenum(1995,12,1) ...
    ] );
    datetick('x', 'mm/yy', 'keeplimits', 'keepticks')  %dd/mm
    title (['Sea Ice Area Time Series with ECCO (dashed)'   ])
    
    
end

%% TS diagram - in and out of cavity
if (plot_ts(1)==1)
    clear var1 var2 plotter mask
    var1= rdmds('S',time);
    var2= rdmds('T',time);
    var=rdmds('Depth');
    p1=0;p2=0;
    for i=1:nx
        for j=1:ny
            if min(var(i,j,:))<1500
                for k=1:nz
                    mask=hFacC(i,j,k); 
                    mask(mask==0)=NaN;
                    mask(mask>0)=1;

                    if (lat_plot(i,j,k)>-78) & (lon_plot(i,j,k)<202) 
                        p1=p1+1;
                        salt_1(p1,1)=var1(i,j,k)*mask;
                        temp_1(p1,1)=var2(i,j,k)*mask;
                    end
                    if (lat_plot(i,j,k)<-78)
                        p2=p2+1;
                        salt_2(p2,1)=var1(i,j,k)*mask;
                        temp_2(p2,1)=var2(i,j,k)*mask;                
                    end
                end
            end
        end
    end
    %  
    figure 
    p(1)=plot(salt_1,temp_1,'LineStyle','none'); 
    set(p(1),'Marker','.','MarkerSize',1,'MarkerEdgeColor','r')
    hold on
    p(2)=plot(salt_2,temp_2,'LineStyle','none'); 
    set(p(2),'Marker','.','MarkerSize',1,'MarkerEdgeColor','b')
    xlim ([30 35])
    title (['TS Salinity Date=' datestr(startdate+ (time*dt/(3600*24)))] )
end

%%
%
%  Maps- Compared with ECCO and Observations
%
%% Temp Map - Compared with ECCO
if (plot_t_geo_with_ecco(1)==1)
    clear var plotter lat_trg lon_trg mask_ocean_l lon_p lat_p
    var= rdmds('T',time);
    figure
    subplot('Position',[0.05, 0.2, 0.4, 0.6])
    mask_ocean_l=hFacC(:,:,plot_t_geo_with_ecco(2)); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;
    plotter(:,:)=var(:,:,plot_t_geo_with_ecco(2)).*mask_ocean_l;
    lat_trg(:,:)=lat_plot(:,:,1);
    lon_trg(:,:)=lon_plot(:,:,1);
    pcolor(lon_trg,lat_trg,plotter);  
    caxis([-2 2])
%     colorbar
    shading flat
    hold on
    v=[-2, -1.9, -1.8 -1.5, -1, 0, 1];
    [C,h]=contour(lon_trg,lat_trg,plotter,v,'Color','w');
    clabel(C,h,v);
    title (['Temperature Depth ' num2str(-RC(1,1, plot_t_geo_with_ecco(2))) ' Date=' datestr(startdate+ (time*dt/(3600*24)))  ])
    %
    subplot('Position',[0.5, 0.2, 0.4, 0.6])
    clear plooter lat lon dep t temp
    addpath ('/mnt/data/malyaral/MITgcm-master/ECCO_LLC270_V5_nctiles_monthly')
    file_name=['THETA.0010.nc'];
    lat=ncread(file_name,'lat');
    lon=ncread(file_name,'lon');
    dep=ncread(file_name,'dep');
    for i=1:270
        for j=1:270
    if (lon(i,j)<0)
        lon(i,j)=lon(i,j)+360;
    end
        end
    end
    %
    t=(floor(NF*7/365))*12+ floor ( (NF*7-(floor(NF*7/365))*365.25)/30.4) ; 
    target=-RC(1,1, plot_t_geo_with_ecco(2));
    [val,idx]=min(abs(dep-target));
    k=(idx);
    %        
    plotter=ncread (file_name,'THETA',[1 1 k t], [inf inf 1 1] );
    pcolor (lon,lat,plotter)
    caxis([-2 2])
    xlim([159 230])
    ylim([-85 -65])
    colorbar
    shading flat
    hold on
    v=[-2, -1.9, -1.8 -1.5, -1, 0, 1];
    [C,h]=contour(lon,lat,plotter,v,'Color','w');
    clabel(C,h,v);
    startdate2=datenum(1992,01,01);
    title (['Temperature ECCO Depth ' num2str(dep(k,1)) ' Date=' datestr(startdate2+ (2629800/(3600*24)*t))  ])
end
%%


%% SIarea map - compared with ECCO
if (plot_SI_geo_with_ecco(1)==1)
    clear var plotter lat_trg lon_trg mask_ocean_l lon_p lat_p
    var= rdmds('AREA',time);
    figure
    subplot('Position',[0.05, 0.2, 0.4, 0.6])
       mask_ocean_l=hFacC(:,:,1); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;
    plotter(:,:)=var(:,:).*mask_ocean_l;
    lat_trg(:,:)=lat_plot(:,:,1);
    lon_trg(:,:)=lon_plot(:,:,1);
    pcolor(lon_trg,lat_trg,plotter);  
%     caxis([-2 2])
     colorbar
    shading flat
    
    title (['Sea Ice Area Date=' datestr(startdate+ (time*dt/(3600*24)))  ])
    %
    subplot('Position',[0.5, 0.2, 0.4, 0.6])
    clear plooter lat lon dep t temp
    addpath ('/mnt/data/malyaral/MITgcm-master/ECCO_LLC270_V5_nctiles_monthly')
    file_name=['SIarea.0010.nc'];
    lat=ncread(file_name,'lat');
    lon=ncread(file_name,'lon');
    for i=1:270
        for j=1:270
    if (lon(i,j)<0)
        lon(i,j)=lon(i,j)+360;
    end
        end
    end
    %
    t=(floor(NF*7/365))*12+ floor ( (NF*7-(floor(NF*7/365))*365.25)/30.4) ; 
    %        
    plotter=ncread (file_name,'SIarea',[1 1 t], [inf inf 1] );
    pcolor (lon,lat,plotter)
    xlim([159 230])
    ylim([-85 -65])
    colorbar
    shading flat
    startdate2=datenum(1992,01,01);
    title (['Sea Ice Area ECCO  Date=' datestr(startdate2+ (2629800/(3600*24)*t))  ])
end

%% SeaIce Diagnostic - SI velocity - Compared with ECCO
if (plot_SI_velocity(1)==1)    
    clear var plotter1 plotter2  mask_ocean_l  
    var= rdmds('siDiag',time); 
    
    figure 
    subplot('Position',[0.05, 0.2, 0.4, 0.6])
    mask_ocean_l=hFacC(:,:,1); 
    mask_ocean_l(mask_ocean_l==0)=NaN;
    mask_ocean_l(mask_ocean_l>0)=1;
    plotter1(:,:)=var(:,:,2).*mask_ocean_l;
    plotter2(:,:)=var(:,:,3).*mask_ocean_l;
    step=plot_SI_velocity(2); scale=plot_SI_velocity(3);
    lat_trg(:,:)=lat_plot(:,:,1);
    lon_trg(:,:)=lon_plot(:,:,1);  
    quiver(lon_trg(1:step:end,1:step:end),lat_trg(1:step:end,1:step:end), plotter1(1:step:end,1:step:end) ,plotter2(1:step:end,1:step:end),scale);  
    xlim ([160 230]);
    ylim ([-85 -65]);
    caxis([-0.5 0.5])
    title (['Sea Ice Velocity. Date=' datestr(startdate+ (time*dt/(3600*24)))])
    
    subplot('Position',[0.5, 0.2, 0.4, 0.6])
    clear plooter lat lon dep t temp
    addpath ('/mnt/data/malyaral/MITgcm-master/ECCO_LLC270_V5_nctiles_monthly')
    t=(floor(NF*7/365))*12+ floor ( (NF*7-(floor(NF*7/365))*365.25)/30.4) ; 
    file_name=['SIuice.0010.nc'];
    vice=-ncread (file_name,'SIuice',[1 1  t], [inf inf  1] );
    file_name=['SIvice.0010.nc'];
    uice=ncread (file_name,'SIvice',[1 1  t], [inf inf 1] );
    lat=ncread(file_name,'lat');
    lon=ncread(file_name,'lon');
    for i=1:270
        for j=1:270
    if (lon(i,j)<0)
        lon(i,j)=lon(i,j)+360;
    end
        end
    end
    
    quiver(lon(1:step:end,1:step:end), lat(1:step:end,1:step:end),uice(1:step:end,1:step:end),vice(1:step:end,1:step:end),scale); 
    xlim ([160 230]);
    ylim ([-85 -65]);
    startdate2=datenum(1992,01,01);
    title (['Sea Ice Area ECCO  Date=' datestr(startdate2+ (2629800/(3600*24)*t))  ])
end


%% Temp Time Series - Compared to ECCO and OISST
if (plot_time_t_ecco_oisst(1,1)==1)
    clear var plotter name z time_stamp points p_ecco p_mit
    %
    file_name=['THETA.0010.nc'];
    file_name2=['SIarea.0010.nc'];
    temp=ncread (file_name,'THETA',[1 1 1 1], [inf inf inf 312] );
    a=ncread (file_name2,'SIarea',[1 1  1], [inf inf 312] );
    p_ecco= [139 88 1; 145 102 21; 143 67 1; 120 72 1; 148 116 21; 135, 123 21; 135, 123,1];
    points=size(p_ecco,1);
    for p=1:points
        for t=1:312
            plotter_ecco(p,t)=temp(p_ecco(p,1),p_ecco(p,2),p_ecco(p,3),t);
            time_ecco(t,:)=(datenum(1992,01,01)+t*2629800/(3600*24));
        end
    end
    %Points loop
    p_mit= [44 81 1; 63 73 25; 19 75 1; 17 100 1; 94 71 25; 87 90 25; 87,90,1];
    for p=1:points
       % Time loop 
        for t=plot_time_t_ecco_oisst(1,2):plot_time_t_ecco_oisst(1,3)
            time1=time_options(t);
                var= rdmds('T',time1);
                plotter(p,t)=var(p_mit(p,1),p_mit(p,2),p_mit(p,3));
                time_stamp(t,:)=(startdate+ (time1*dt/(3600*24)));
        end
    end
    
    load('/mnt/data/malyaral/OneDrive/MATLAB/Satellite_SST_oisst_v1/SST_9_MIT/temp_oisst_area.mat')
   
    point_oisst=[10 193;  21 158]; %RSP TNBP
    %
    figure
    subplot('Position',[0.1,0.6,0.8,0.3])
    p=1;
    hold on
    z(p)=plot (time_stamp,plotter(p,:));
    set(z(p),'Color',colorOrder(p,:),'LineWidth',1);
    hold on
    z1(p)=plot (time_ecco,plotter_ecco(p,:));
    set(z1(p),'Color',colorOrder(p,:),'LineStyle','--');
    hold on
    plotter_oisst(:)=temp_oisst(point_oisst(p,1),point_oisst(p,2),:);
    z2(p)=plot (time_oisst,plotter_oisst);
    set(z2(p),'Color','k','LineWidth',2);
    %
      set(gca, 'Xlim', [time_stamp(plot_time_t_ecco_oisst(1,2),:) time_stamp(plot_time_t_ecco_oisst(1,3),:)])  
    set(gca,'XTick', [ datenum(1992,1,1)  datenum(1992,4,1) datenum(1992,7,1)  datenum(1992,10,1) ...
    datenum(1993,1,1)  datenum(1993,4,1) datenum(1993,7,1) datenum(1993,10,1) ...
    datenum(1994,1,1) datenum(1994,4,1)  datenum(1994,7,1)  datenum(1994,10,1) ...
    datenum(1995,1,1)  datenum(1995,4,1)  datenum(1995,7,1) datenum(1995,10,1) ...
    datenum(1996,1,1)  datenum(1996,4,1)  datenum(1996,7,1) datenum(1996,10,1) ...
    datenum(1997,1,1)  datenum(1997,4,1)  datenum(1997,7,1) datenum(1997,10,1) ...
    datenum(1998,1,1)  datenum(1998,4,1)  datenum(1998,7,1) datenum(1998,10,1) ...
    datenum(1999,1,1)  datenum(1999,4,1)  datenum(1999,7,1) datenum(1999,10,1) ...
    datenum(2000,1,1)  datenum(2000,4,1)  datenum(2000,7,1) datenum(2000,10,1) ...
    ] );
    datetick('x', 'mm/yy', 'keeplimits', 'keepticks')  %dd/mm
    title (['Ross Sea Polynya - Temperature Time Series with ECCO (dashed) and OISST (bold)'   ])
    %
    subplot('Position',[0.1,0.2,0.8,0.3])
    p=4;
    hold on
    z(p)=plot (time_stamp,plotter(p,:));
    set(z(p),'Color',colorOrder(3,:),'LineWidth',1);
    hold on
    z1(p)=plot (time_ecco,plotter_ecco(p,:));
    set(z1(p),'Color',colorOrder(3,:),'LineStyle','--');
     
    hold on
    p=2;
    plotter_oisst(:)=temp_oisst(point_oisst(p,1),point_oisst(p,2),:);
    z2(p)=plot (time_oisst,plotter_oisst);
    set(z2(p),'Color','k','LineWidth',2);
    %
    set(gca, 'Xlim', [time_stamp(plot_time_t_ecco_oisst(1,2),:) time_stamp(plot_time_t_ecco_oisst(1,3),:)])  
    set(gca,'XTick', [ datenum(1992,1,1)  datenum(1992,4,1) datenum(1992,7,1)  datenum(1992,10,1) ...
    datenum(1993,1,1)  datenum(1993,4,1) datenum(1993,7,1) datenum(1993,10,1) ...
    datenum(1994,1,1) datenum(1994,4,1)  datenum(1994,7,1)  datenum(1994,10,1) ...
    datenum(1995,1,1)  datenum(1995,4,1)  datenum(1995,7,1) datenum(1995,10,1) ...
    datenum(1996,1,1)  datenum(1996,4,1)  datenum(1996,7,1) datenum(1996,10,1) ...
    datenum(1997,1,1)  datenum(1997,4,1)  datenum(1997,7,1) datenum(1997,10,1) ...
    datenum(1998,1,1)  datenum(1998,4,1)  datenum(1998,7,1) datenum(1998,10,1) ...
    datenum(1999,1,1)  datenum(1999,4,1)  datenum(1999,7,1) datenum(1999,10,1) ...
    datenum(2000,1,1)  datenum(2000,4,1)  datenum(2000,7,1) datenum(2000,10,1) ...
    ] );
    datetick('x', 'mm/yy', 'keeplimits', 'keepticks')  %dd/mm
    title (['Terra Nova Bay - Temperature Time Series with ECCO (dashed) and OISST (bold)'   ])
end



%% Sie Ice Time Series - Compared to ECCO and NSIDC
if (plot_time_si_ecco_nsidc (1,1)==1)
    clear var plotter name z time_stamp points p_ecco p_mit plotter_a ice_ecco time_ecco
    %
    file_name2=['SIarea.0010.nc'];
    a=ncread (file_name2,'SIarea',[1 1  1], [inf inf 312] );
    p_ecco= [139 88 1; 145 102 21; 143 67 1; 120 72 1; 148 116 21; 135, 123 21;];
    points=size(p_ecco,1);
    for p=1:points
        for t=1:312
            ice_ecco(p,t)=a(p_ecco(p,1),p_ecco(p,2),t);
            time_ecco(t,:)=(datenum(1992,01,01)+t*2629800/(3600*24));
        end
    end
    %Points loop
    p_mit= [44 81 1; 63 73 25; 19 75 1; 17 100 1; 94 71 25; 87 90 25];
    for p=1:points
       % Time loop 
        for t=plot_time_si_ecco_nsidc(1,2):plot_time_si_ecco_nsidc(1,3)
            time1=time_options(t);
                var1=rdmds('AREA',time1);
                plotter_a(p,t)=var1(p_mit(p,1),p_mit(p,2));
                time_stamp(t,:)=(startdate+ (time1*dt/(3600*24)));
        end
    end    
    load('/mnt/data/malyaral/OneDrive/MATLAB/Satellite_SIC_all_datasets/Nimbus7/MIT_1/mit_nsidc.mat')
    %
    figure
     subplot('Position',[0.1,0.7,0.8,0.25])
    p=1;
    z2(p)=plot (date_nsidc(1,:),ci_nsidc(:,17,91)/100);
    set(z2(p),'Color','k','LineWidth',1);
    hold on
    z(p)=plot (time_stamp,plotter_a(p,:));
    set(z(p),'Color',colorOrder(p,:),'LineWidth',1);
    hold on
    z1(p)=plot (time_ecco,ice_ecco(p,:));
    set(z1(p),'Color',colorOrder(p,:),'LineStyle','--');
    %
    set(gca, 'Xlim', [time_stamp(plot_time_si_ecco_nsidc(1,2),:) time_stamp(plot_time_si_ecco_nsidc(1,3),:)])  
    set(gca,'XTick', [ datenum(1992,1,1)  datenum(1992,4,1) datenum(1992,7,1)  datenum(1992,10,1) ...
    datenum(1993,1,1)  datenum(1993,4,1) datenum(1993,7,1) datenum(1993,10,1) ...
    datenum(1994,1,1) datenum(1994,4,1)  datenum(1994,7,1)  datenum(1994,10,1) ...
    datenum(1995,1,1)  datenum(1995,4,1)  datenum(1995,7,1) datenum(1995,10,1) ...
    datenum(1996,1,1)  datenum(1996,4,1)  datenum(1996,7,1) datenum(1996,10,1) ...
    datenum(1997,1,1)  datenum(1997,4,1)  datenum(1997,7,1) datenum(1997,10,1) ...
    datenum(1998,1,1)  datenum(1998,4,1)  datenum(1998,7,1) datenum(1998,10,1) ...
    datenum(1999,1,1)  datenum(1999,4,1)  datenum(1999,7,1) datenum(1999,10,1) ...
    datenum(2000,1,1)  datenum(2000,4,1)  datenum(2000,7,1) datenum(2000,10,1) ...
    ] );
    datetick('x', 'mm/yy', 'keeplimits', 'keepticks')  %dd/mm
    title (['Ross Sea Polynya - SIarea Time Series with ECCO (dashed) and NSIDC (black)'   ])
    %
    subplot('Position',[0.1,0.35,0.8,0.25])
    p=4;
    hold on
    z2(p)=plot (date_nsidc(1,:),ci_nsidc(:,25,101)/100);
    set(z2(p),'Color','k','LineWidth',1);
    hold on
    z(p)=plot (time_stamp,plotter_a(p,:));
    set(z(p),'Color',colorOrder(3,:),'LineWidth',1);
    hold on
    z1(p)=plot (time_ecco,ice_ecco(p,:));
    set(z1(p),'Color',colorOrder(3,:),'LineStyle','--');
    %
    set(gca, 'Xlim', [time_stamp(plot_time_si_ecco_nsidc(1,2),:) time_stamp(plot_time_si_ecco_nsidc(1,3),:)])  
    set(gca,'XTick', [ datenum(1992,1,1)  datenum(1992,4,1) datenum(1992,7,1)  datenum(1992,10,1) ...
    datenum(1993,1,1)  datenum(1993,4,1) datenum(1993,7,1) datenum(1993,10,1) ...
    datenum(1994,1,1) datenum(1994,4,1)  datenum(1994,7,1)  datenum(1994,10,1) ...
    datenum(1995,1,1)  datenum(1995,4,1)  datenum(1995,7,1) datenum(1995,10,1) ...
    datenum(1996,1,1)  datenum(1996,4,1)  datenum(1996,7,1) datenum(1996,10,1) ...
    datenum(1997,1,1)  datenum(1997,4,1)  datenum(1997,7,1) datenum(1997,10,1) ...
    datenum(1998,1,1)  datenum(1998,4,1)  datenum(1998,7,1) datenum(1998,10,1) ...
    datenum(1999,1,1)  datenum(1999,4,1)  datenum(1999,7,1) datenum(1999,10,1) ...
    datenum(2000,1,1)  datenum(2000,4,1)  datenum(2000,7,1) datenum(2000,10,1) ...
    ] );
    datetick('x', 'mm/yy', 'keeplimits', 'keepticks')  %dd/mm
    title (['Terra Nova Bay - SIarea Time Series with ECCO (dashed) and NSIDC (black)'   ])
    %
    subplot('Position',[0.1,0.05,0.8,0.25])
    p=6;
    hold on
    z2(p)=plot (date_nsidc(1,:),ci_nsidc(:,31,74)/100);
    set(z2(p),'Color','k','LineWidth',1);
    hold on
    z(p)=plot (time_stamp,plotter_a(p,:));
    set(z(p),'Color',colorOrder(3,:),'LineWidth',1);
    hold on
    z1(p)=plot (time_ecco,ice_ecco(p,:));
    set(z1(p),'Color',colorOrder(3,:),'LineStyle','--');
    %
    set(gca, 'Xlim', [time_stamp(plot_time_si_ecco_nsidc(1,2),:) time_stamp(plot_time_si_ecco_nsidc(1,3),:)])  
    set(gca,'XTick', [ datenum(1992,1,1)  datenum(1992,4,1) datenum(1992,7,1)  datenum(1992,10,1) ...
    datenum(1993,1,1)  datenum(1993,4,1) datenum(1993,7,1) datenum(1993,10,1) ...
    datenum(1994,1,1) datenum(1994,4,1)  datenum(1994,7,1)  datenum(1994,10,1) ...
    datenum(1995,1,1)  datenum(1995,4,1)  datenum(1995,7,1) datenum(1995,10,1) ...
    datenum(1996,1,1)  datenum(1996,4,1)  datenum(1996,7,1) datenum(1996,10,1) ...
    datenum(1997,1,1)  datenum(1997,4,1)  datenum(1997,7,1) datenum(1997,10,1) ...
    datenum(1998,1,1)  datenum(1998,4,1)  datenum(1998,7,1) datenum(1998,10,1) ...
    datenum(1999,1,1)  datenum(1999,4,1)  datenum(1999,7,1) datenum(1999,10,1) ...
    datenum(2000,1,1)  datenum(2000,4,1)  datenum(2000,7,1) datenum(2000,10,1) ...
    ] );
    datetick('x', 'mm/yy', 'keeplimits', 'keepticks')  %dd/mm
    title (['Contonental Slope Central - SIarea Time Series with ECCO (dashed) and NSIDC (black)'   ])
    
end

%%
%% Temp Frontal Section (outside cavity)
if (plot_mixing_sec_front(1)==1)
    clear var plotter points
    var= rdmds('kppDiag3',time);
    front=[13 80; 14 80;15 80;16 80;17 80;18 80;19 80;20 80;21 80;22 80; ...
        23 79; 24 79; 25 79; 26 79; ...
        27 78; 28 78; 29 78; ...
        30 77; 31 77; 32 77; 33 77; 34 77; 35 77; 36 77; 37 77; 38 77; 39 77; 40,77; ...
        41 76; 42 76; 43 76; 44 76; 45 76; 46 76; 47 76;48 76; 49 76; 50 76 ; ...
        51 75; 52 74; 53 74;...
        54 73; 55 73; 56 73; 57 73; 58 73; 59 73; 60 73; 61 73; 62 73; 63 73; 64 73; 65 73; 66 73; 67 73;...
        68 72; 69 72; 70 72; 72 70; 73 70; 73 69; ...
        74 68; 75 68; 76 68; 77 68; 78 68; 79 68; 80 68; 81 68;...
        82 67; 83 67; 84 67; 85 67; 86 67; 87 67; 88 67; 89 67; 90 67; 91 67; 92 67; 93 67; 94 67; 95 67; 96 67; ...
        97 66; 98 66; 99 66; 100 66; 101 66; 102  66; 103  66; 104  66; 105 66; 106 66; 107 66; 108 66; 109 66; 110 66; 111 66; 112 66; 113 66; ...
        114 67; 115 67; 116 67; 117 67; 118 67; 119 68; 120 68; 121 69; 122 70; 123 71; 124 72; 125 73; 126 74 ];
    points=size (front,1 );    
    %
    for p=1:points
        mask_ocean_l=hFacC(front(p,1),front(p,2),:,:); 
        mask_ocean_l(mask_ocean_l==0)=NaN;
        mask_ocean_l(mask_ocean_l>0)=1;  
        lon_p(p,: )= lon_plot(front(p,1),front(p,2),: );
        depth_p(p,:)= depth_plot(front(p,1),front(p,2),: );
        temp_p(p,:)= var(front(p,1),front(p,2),:,1 ).*mask_ocean_l;
    end
   %
   figure
   pcolor(lon_p,depth_p,temp_p); 
     caxis ([0 0.01])
   
   colorbar
%    hold on
%    v=[-2, -1.9, -1.8 -1.5];
%    [C,h]=contour(lon_p,depth_p,temp_p,v,'Color','w');
%    clabel(C,h,v);
% %       caxis([-2.5 0])
    ylim([-1000 0])
  

   shading flat
   title (['Temperature section along ice shelf front' ' Date=' datestr(startdate+ (time*dt/(3600*24))) ] )
end



%% Temp like itp in time
if (plot_t_itp(1)==1)
    clear var plotter points
    var= rdmds('T',time);
    point=[73,78];
    p=1;
    %
    for t=1:NF
        time1=time_options(t);
        var= rdmds('T',time1);
        time_p(t,1:70)= startdate+ (time1*dt/(3600*24));
        depth_p(t,:)= depth_plot(point(p,1),point(p,2),: );
        temp_p(t,:)= var(point(p,1),point(p,2),: );
    end
   %
   figure
   pcolor(time_p,depth_p,temp_p); 
   hold on
   v=[-2, -1.9, -1.8 -1.5];
   [C,h]=contour(time_p,depth_p,temp_p,v,'Color','w');
   clabel(C,h,v);
%        caxis([-2.5 -1])
    ylim([-600 0])
   colorbar
   shading flat
   title (['Temperature Central ' ] )
          set(gca,'XTick', [ datenum(1992,1,1)  datenum(1992,4,1) datenum(1992,7,1)  datenum(1992,10,1) ...
    datenum(1993,1,1)  datenum(1993,4,1) datenum(1993,7,1) datenum(1993,10,1) ...
    datenum(1994,1,1) datenum(1994,4,1)  datenum(1994,7,1)  datenum(1994,10,1) ...
    datenum(1995,1,1)  datenum(1995,4,1)  datenum(1995,7,1) datenum(1995,10,1) ...
    datenum(1996,1,1)  datenum(1996,4,1)  datenum(1996,7,1) datenum(1996,10,1) ...
    datenum(1997,1,1)  datenum(1997,4,1)  datenum(1997,7,1) datenum(1997,10,1) ...
    datenum(1998,1,1)  datenum(1998,4,1)  datenum(1998,7,1) datenum(1998,10,1) ...
    datenum(1999,1,1)  datenum(1999,4,1)  datenum(1999,7,1) datenum(1999,10,1) ...
    datenum(2000,1,1)  datenum(2000,4,1)  datenum(2000,7,1) datenum(2000,10,1) ...
    ] );
    datetick('x', 'mm/yy', 'keeplimits', 'keepticks')  %dd/mm
end



%% Temp like itp in time
if (plot_s_itp(1)==1)
    clear var plotter points
%     var= rdmds('S',time);
    point=[73,78];
    p=1;
    %
    for t=1:NF
        time1=time_options(t);
        var= rdmds('S',time1);
        time_p(t,1:70)= startdate+ (time1*dt/(3600*24));
        depth_p(t,:)= depth_plot(point(p,1),point(p,2),: );
        temp_p(t,:)= var(point(p,1),point(p,2),: );
    end
   %
   figure
   pcolor(time_p,depth_p,temp_p); 
   caxis([33 35])
   hold on
   v=[33:0.2:35];
   [C,h]=contour(time_p,depth_p,temp_p,v,'Color','w');
   clabel(C,h,v);
%        caxis([-2.5 -1])
   ylim([-600 0])
   colorbar
   shading flat
   title (['Salinity ' ] )
       set(gca,'XTick', [ datenum(1992,1,1)  datenum(1992,4,1) datenum(1992,7,1)  datenum(1992,10,1) ...
    datenum(1993,1,1)  datenum(1993,4,1) datenum(1993,7,1) datenum(1993,10,1) ...
    datenum(1994,1,1) datenum(1994,4,1)  datenum(1994,7,1)  datenum(1994,10,1) ...
    datenum(1995,1,1)  datenum(1995,4,1)  datenum(1995,7,1) datenum(1995,10,1) ...
    datenum(1996,1,1)  datenum(1996,4,1)  datenum(1996,7,1) datenum(1996,10,1) ...
    datenum(1997,1,1)  datenum(1997,4,1)  datenum(1997,7,1) datenum(1997,10,1) ...
    datenum(1998,1,1)  datenum(1998,4,1)  datenum(1998,7,1) datenum(1998,10,1) ...
    datenum(1999,1,1)  datenum(1999,4,1)  datenum(1999,7,1) datenum(1999,10,1) ...
    datenum(2000,1,1)  datenum(2000,4,1)  datenum(2000,7,1) datenum(2000,10,1) ...
    ] );
    datetick('x', 'mm/yy', 'keeplimits', 'keepticks')  %dd/mm
end


