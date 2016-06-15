% % 
% See Get_EDGES_Beams_1_degree_steps.m
% %

clear;
close all;
saveforpaper=false;
longitude=-45;
v_i=150; % used later for beam - sky convolution
orientation=true;

time='2000/03/20 12:06:40.59856';
[d_year, d_month, d_day, d_hour, d_minute, d_second] = datevec(datenum(time,'yyyy/mm/dd HH:MM:SS'));
JD = juliandate(d_year,d_month,d_day,d_hour,d_minute,d_second);
T_UT1 = (JD-2451545)./36525;
ThetaGMST = 67310.54841 + (876600*3600 + 8640184.812866).*T_UT1 ...
+ .093104.*(T_UT1.^2) - (6.2*10^-6).*(T_UT1.^3);
ThetaGMST = mod((mod(ThetaGMST,86400*(ThetaGMST./abs(ThetaGMST)))/240),360);
% ThetaGMST should be near 0 or 360 degrees.

% plot the haslam sky map

load('sky_square_40to200.mat')
sky150(:,:)=sky_square(150,:,:);
[myt, myra, mydec] = getsky('C:\Users\Tom\OneDrive\ASU_Thesis\Beams\Beam_on_Sky\radio408.RaDec.fits');

myt=(myt./10);
myt150=myt*(408/150)^2.5;
RA_sky=ones(541,1)*myra;
dec_out=ones(1081,1)*mydec;
Dec_sky=dec_out';
latitudes = -90:(1/3):90;
% Get the sky data - only the dec_out and ra_out are used.
% the sky data was preprocessed, scaled and stored by frequencyto save time

figure(100);
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none')
set(figure(100),'OuterPosition',[.12 .6 .3 .4])

xmin=1;    xmax=1081;    ymin=1;    ymax=541;  xstep=135; ystep=90;
obj100=pcolor(1:1081,1:541,sky150);
title('Haslam Sky Map at 150 MHz','FontSize',16);
set(obj100,'linestyle','none');
% caxis
colorbar;
caxis([0 3000]);
ylabel(colorbar,'Kelvin', 'Rotation',90, 'FontSize' , 14);
% X axis graph points
xlim([ xmin xmax]);
set(gca,'Xtick', (xmin:xstep:xmax));
set(gca, 'XTickLabel', ( xmin:xstep:xmax));
xlabel('RA data point index','FontSize',16);
% Y axis Dec
ylim([ ymin ymax]);
set(gca,'Ytick', (ymin:ystep:ymax));
ylabel('Dec data point index','FontSize',16);

set(gca,'FontSize',14);

figure(101);
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.02 .6 .3 .4])

xmin=-180;    xmax=180;    ymin=-90;    ymax=90;
obj101=pcolor(myra,latitudes,sky150);
title('Haslam Sky Map at 150 MHz','FontSize',16);
set(obj101,'linestyle','none');

xlabel('RA (Degrees)','FontSize',16);
ylabel('Dec (degrees)','FontSize',16);

colorbar;
caxis([0 3000]);
ylabel(colorbar,'Kelvin', 'Rotation',90, 'FontSize' , 14);

% X axis RA 
xlim([ xmin xmax]);
set(gca,'Xtick', (xmin:30:xmax));
set(gca, 'XTickLabel', ([ xmin:30:xmax]));

% Y axis Dec
ylim([ ymin ymax]);
set(gca,'Ytick', ([ymin:30:ymax]));

set(gca,'FontSize',14);

% Skymap with LST axis
figure(102);
set(figure(102),'Units','normalized','MenuBar','none','ToolBar','none')
set(figure(102),'OuterPosition',[.02 .5 .3 .4])

obj102=pcolor(myra,latitudes,sky150);
title('Haslam Sky Map at 150 MHz','FontSize',16);
set(obj102,'linestyle','none');
colorbar;
caxis([0 3000]);
ylabel(colorbar,'Kelvin', 'Rotation',90, 'FontSize' , 14);

ylabel('Dec (degrees)','FontSize',16);
xlim([ xmin xmax]);
set(gca,'Xtick', (xmin:30:xmax));

% LST axis
set(gca, 'XTickLabel', ([ 12 14 16 18 20 22 0 2 4 6 8 10 12]));
xlabel('LST (Hours))','FontSize',16);

ylim([ ymin ymax]);
set(gca,'Ytick', ([ymin:30:ymax]));

set(gca,'FontSize',15);

%Skymap Spectrum at Lat=-26.7
figure(777)
set(figure(777),'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.02 .1 .3 .4])

skym26p7(:,:)=sky_square(100:200,65*3,:); % full frequency 100 to 200 MHz at lat -26 all LST values
longit=180:-1/3:-180;
datapts=1:1081;
obj777=pcolor(100:200,datapts,skym26p7(1:101,:)');set(obj777,'linestyle','none');
title('Haslam Sky Spectrum at Dec -27','FontSize',16);
colorbar;
ylabel(colorbar,'Kelvin', 'Rotation',90, 'FontSize' , 14);
caxis([100 3000]);

xlim([ 100 200]);
set(gca,'Xtick', ([100:10:200]));
ylim([ 1 1081]);
set(gca,'Ytick', ([1:135:1081]));
%set(gca, 'YTickLabel', ([ 12 14 16 18 20 22 0 2 4 6 8 10 12]));
ylabel('LST (Hours))','FontSize',16);

xlabel('Frequency (MHz)','FontSize',16);

set(gca,'FontSize',15);

%%
% Lat -26.7 temperature
figure(888)
clf(888);
set(figure(888),'Units','normalized','MenuBar','none','ToolBar','none')
set(figure(888),'OuterPosition',[.12 .15 .3 .4])
skytemp_m26p7=zeros(1,1081);
for lst=1:1081
    for freq_index=1:101
        skytemp_m26p7(lst)=skytemp_m26p7(lst) + skym26p7(freq_index,lst);
    end % frequency_index
end % lst

xmin=1;    xmax=1081;    ymin=1;    ymax=541;  xstep=90; ystep=90;
datapts=1:1081;
obj888=plot(datapts,skytemp_m26p7/101);
title('Sky Temp Averaged over Freq Range vs. LST (Dec -27)','FontSize',16);
xlim([ xmin xmax]);
set(gca,'Xtick', (xmin:xstep:xmax));
%set(gca, 'XTickLabel', ([ 12 14 16 18 20 22 0 2 4 6 8 10 12]));
set(gca, 'XTickLabel', ([ 12 10 8 6 4 2 24 22 20 18 16 14 12]));
xlabel('LST (Hours)','FontSize',16);
%set(gca,'XDir','Reverse');

ylabel('Sky Temp Sum / 101','FontSize',16);
ylim([ 0 7000]);
set(gca,'Ytick', ([0:1000:7000]));

set(gca,'FontSize',15);
grid on;
set(gca,'GridLineStyle','-');



%%

% 'EDGES_beam_data_100-200_MHz_1_degree_steps' file contains:
% 101 frequency blocks from 100 MHz to 200 MHz
% with 181 Elevation steps from  theta = 0 to 180 degrees in steps of 1
% with 360 Azimuth steps from phi = 0 to 360 degrees in steps of 1
% Note increasing values of phi are in a counterclockwise direction - as
% per spherical coordinate conventions: +z axis pointing up, +x axis coming
% at you, and + y axis pointing to the right. phi starts at +x=0 and moves
% towards +y. This is opposite to azimuthal values

load(['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/Blade_highband/beam_data_Normed_linear_Highband_Blade_90MHz-200MHz.mat']);
%load('C:\Users\Tom\OneDrive\ASU_Thesis\Beams\Beam_on_Sky\EDGES\beam_data_EDGES_linear_normed_1deg_step_without_LNA_box_June');
%load('C:\Users\Tom\OneDrive\ASU_Thesis\Beams\Beam_on_Sky\Gaussian\beam_data_square_Gaussian_4p25_sigma_linear_normed_1deg');
%load('C:\Users\Tom\OneDrive\ASU_Thesis\Beams\Beam_on_Sky\Gaussian\beam_data_square_Gaussian_30_sigma_linear_normed_1deg');
%load('C:\Users\Tom\OneDrive\ASU_Thesis\Beams\Beam_on_Sky\Dipole\beam_data_square_Dipole_linear_normalized_176p7MHz_1deg');

% Az block looks like:                  el block looks like:               Phi block looks like:
% 360  359  358  ... 0      line 1       90   90   90   90  ...  90        0     1    2    3  ...   360
% 360  359  358  ... 0      line 2       89   89   89   89  ...  89        0     1    2    3  ...   360
% ...  ...  ...  ... ...     ...        ...  ...  ...  ...  ...  ...       0     1    2    3  ...   360
% 360  359  358  ... 0      line 181    -90  -90  -90  -90  ...  -90       0     1    2    3  ...   360
% col1 col2 col3 ... col361             col1 col2 col3 col4 ... col361     col1 col2 col3 col4 ... col361

Az_sec=zeros(181,361);
El_sec=zeros(181,361);
for Az=360:-1:0  % same as Az=360:-1:0; repmat(Az,181,1);
    Az_sec(:,360-Az+1)=Az;
end
for El=90:-1:-90 % same as El=[90:-1:-90]'; repmat(tom_el,1,361);
    El_sec(90-El+1,:)=El;
end

%%
Polarization=-45; % positive numbers mean the beam is rotated clockwise, which means the matrix gets shifted to the left.
display(Polarization);
Polarization_step=Polarization; % rotation columns = Polarization number

% Initialize some variables
v_res_mod=10;
h_res_mod=15;

latitudes = -90:10/v_res_mod:90;

figure_counter=0;

%%
for lat_index=66:25:66
    
    %lat_index=1; % lat_index=1 --> lat(181) = +90, lat(1) = -90
    latitude = latitudes(lat_index);
    
    %Convert the sky map coordinates to Az and El. (% ThetaGMST should be near 0 or 360 degrees.)
    % Sticking with the coordinate system of the sky, just in Az/El coordinates
    [Az_sky, El_sky]=RaDec2AzEl_fresh_orig_mod1(RA_sky,Dec_sky, latitude,longitude,ThetaGMST);
    
    
    %%
    %  ELEVATION output from RaDec2AzEl
    %
    
    figure(1+figure_counter);
    set(gcf,'Units','normalized','MenuBar','none','ToolBar','none')
    set(figure(gcf),'OuterPosition',[.32 .2 .3 .4])
    obj1=pcolor(180:-1/3:-180,-90:1/3:90,El_sky);
    set(obj1,'linestyle','none');colorbar;
    ylabel(colorbar,'Elevation (Degrees)', 'Rotation',90, 'FontSize' , 14);
    title(['ELEVATION, Output of RaDec2AzEl, Lat= ' num2str(latitudes(lat_index))],'FontSize',16);
    
    hcb=colorbar;
    caxis([-90 90]);
    set(hcb,'YTick',[-90  -60 -30  0 30 60 90], 'FontSize', 10)
    
    xlabel('RA (degrees)');
    ylabel('Dec (Degrees)');
    
    h_xlabel=get(gca,'xlabel');
    set(h_xlabel,'FontSize',16);
    xlim([ -180 180]);
    set(gca,'Xtick', ([-180:30:180]));
    %set(gca,'XDir','reverse');
    
    h_ylabel=get(gca,'ylabel');
    set(h_ylabel,'FontSize',16);
    ylim([ -90 90]);
    set(gca,'Ytick', ([-90:30:90]))
    
    set(gca,'FontSize',13); % this sets the size of the numbers on the axis
    
    %%
    %  AZIMUTH output from RaDec2AzEl
    %
    figure(2+figure_counter);
    set(gcf,'Units','normalized','MenuBar','none','ToolBar','none')
    set(figure(gcf),'OuterPosition',[.32 .1 .3 .4])
    obj2=pcolor(myra,-90:1/3:90,Az_sky);
    set(obj2,'linestyle','none');colorbar;  ylabel(colorbar,'Azimuth (Degrees)', 'Rotation',90, 'FontSize' , 14);
    title(['AZIMUTH, Output of RaDec2AzEl,  Lat= ' num2str(latitudes(lat_index))],'FontSize',16);
    
    hcb=colorbar;
    caxis([0 360]);
    set(hcb,'YTick',[0  60 120  180 240 300 360], 'FontSize', 10)
    
    xlabel('RA (degrees)');
    ylabel('Dec (Degrees)');
    
    h_xlabel=get(gca,'xlabel');
    set(h_xlabel,'FontSize',16);
    xlim([ -180 180]);
    set(gca,'Xtick', ([-180:30:180]));
    %set(gca,'XDir','reverse');
    
    h_ylabel=get(gca,'ylabel');
    set(h_ylabel,'FontSize',16);
    ylim([ -90 90]);
    set(gca,'Ytick', ([-90:30:90]));
    
    set(gca,'FontSize',13); % this sets the size of the numbers on the axis
    
    %%
    %  Plot Beam at v=100 MHz
    %
    
    %Get the beam block of data
    %beam_sec=reshape(beam(i:i+65340), [181 361]);
    beam_sec(: , :) = beam_data_square(v_i , : , :);
    
    if orientation==true
        % Imprint a pattern on the beam
        beam_sec(1:25,2:7)=3000.0/2000;     % North or  0 / 360 degrees a tad west of north - this is skinny width
        beam_sec(1:25,82:98)=3000.0/2000;   % 270 degrees - West - (CST phi = 90)           - this is medium width
        beam_sec(1:25,165:195)=3000.0/2000; % South or 180 degrees                          - this is a bit wider
        beam_sec(1:25,250:290)=3000.0/2000; % 90 degrees - East - (CST phi = 270)           - this is the fat guy
        beam_sec(1:45,354:359)=3000.0/2000; % north again - just a little east - this is the tall skinny spike
    else
        ind=find(beam_sec>0.5);
        beam_sec(ind)=1;
    end
    
    % Polarization rotates the beam to the right position
    beam_sec(:,1:360)=circshift(beam_sec(:,1:360),[0,-Polarization]);
    beam_sec(:,361)=beam_sec(:,1);
    
    % BEAM in CST Coordinates
    figure(3+figure_counter);
    set(gcf,'Units','normalized','MenuBar','none','ToolBar','none')
    set(figure(gcf),'OuterPosition',[.62 .1 .3 .4])
    obj3=pcolor(0:360,0:180,beam_sec);
    set(obj3,'linestyle','none');colorbar; ylabel(colorbar,'Beam Magnitude (Normalized)', 'Rotation',90, 'FontSize' , 14);
    title(['EDGES Beam CST Coordinates, Lat= ' num2str(latitudes(lat_index))],'FontSize',16);
    xlabel('\phi in CST CCW Coordinates (degrees)');
    ylabel('\theta in CST Coordinates (Degrees)');
    
    h_xlabel=get(gca,'xlabel');
    set(h_xlabel,'FontSize',16);
    h_ylabel=get(gca,'ylabel');
    set(h_ylabel,'FontSize',16);
    
    xmin=0;
    xmax=360;
    xlim([ xmin xmax]);
    set(gca,'Xtick', ([xmin:30:xmax]));
    %set(gca,'XDir','reverse');
    
    ymin=0;
    ymax=180;
    ylim([ ymin ymax]);
    set(gca,'Ytick', ([ymin:30:ymax]));
    %set(gca,'YDir','reverse');
    
    set(gca,'FontSize',13); % this sets the size of the numbers on the axis
    caxis([0 1]);
    % BEAM in Az / El Coordinates
    figure(4+figure_counter);
    set(gcf,'Units','normalized','MenuBar','none','ToolBar','none')
    set(figure(gcf),'OuterPosition',[.62 .2 .3 .4])
    xmin=0;    xmax=360;    ymin=-90;    ymax=90;
    
    obj4=pcolor(360:-1:0,ymax:-1:ymin,beam_sec);
    set(obj4,'linestyle','none');colorbar; ylabel(colorbar,'Beam Magnitude (Normalized)', 'Rotation',90, 'FontSize' , 14);
    title(['EDGES Beam Az / El Coordinates, Lat= ' num2str(latitudes(lat_index))],'FontSize',16);
    xlabel('Azimuth (Degrees)');
    ylabel('Elevation (Degrees)');
    
    h_xlabel=get(gca,'xlabel');
    set(h_xlabel,'FontSize',16);
    h_ylabel=get(gca,'ylabel');
    set(h_ylabel,'FontSize',16);
    
    xlim([ xmin xmax]);
    set(gca,'Xtick', ([xmin:30:xmax]));
    set(gca,'XDir','reverse');
    
    ylim([ ymin ymax]);
    set(gca,'Ytick', ([ymin:30:ymax]));
    %set(gca,'YDir','reverse');
    
    set(gca,'FontSize',13); % this sets the size of the numbers on the axis
    caxis([0 1]);
    
    %%
    %   Plot the sky-beam projection
    %
    
    % Get the beam pattern to have the same corrdinate base shape as the sky
    
    figure(5+figure_counter);
    set(gcf,'Units','normalized','MenuBar','none','ToolBar','none')
    set(figure(gcf),'OuterPosition',[.62 .6 .3 .4])
    P_beam_sec=interp2(Az_sec,El_sec, beam_sec, Az_sky,El_sky);
    obj5=pcolor(180:-1/3:-180,-90:1/3:90,P_beam_sec);
    set(obj5,'linestyle','none');colorbar; ylabel(colorbar,'Beam Magnitude (Normalized to 1.0 Max)', 'Rotation',90, 'FontSize' , 14);
    title(['Projected Beam, Lat= ' num2str(latitudes(lat_index))],'FontSize',16);
    xlabel('RA (degrees)');
    ylabel('Dec (Degrees)');
    
    h_xlabel=get(gca,'xlabel');
    set(h_xlabel,'FontSize',16);
    xlim([ -180 180]);
    set(gca,'Xtick', ([-180:30:180]));
    %set(gca,'XDir','reverse');
    
    h_ylabel=get(gca,'ylabel');
    set(h_ylabel,'FontSize',16);
    ylim([ -90 90]);
    set(gca,'Ytick', ([-90:30:90]));
    
    set(gca,'FontSize',13); % this sets the size of the numbers on the axis
    caxis([0 1]);
    
    
    %%
    %  Project Beam * Skymap
    load('C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/sky_square_with_sinefactor_40to200.mat');
    % This loads, RA_sky, Dec_sky, and sky_square(frequency,Dec,RA)
    sinefactor_dec=cosd(Dec_sky);
    
    Normalization=sum(sum(sinefactor_dec.*P_beam_sec));
    P_beam_sec = P_beam_sec/Normalization;
    
    sky_with_sinefactor_150(:,:)=sky_square(150,:,:);
    image=P_beam_sec.*sky_with_sinefactor_150;
    
    figure(6+figure_counter);
    set(gcf,'Units','normalized','MenuBar','none','ToolBar','none')
    set(figure(gcf),'OuterPosition',[.32 .6 .3 .4])
    xmin=-180;    xmax=180;    ymin=-90;    ymax=90;
    
    % Fixed the image so the caxis looked more "normal"
    obj6=pcolor(myra,-90:(1/3):90,Normalization*image);
    title(['Sky Map * Beam Values / (Integral of Projected Beam over Sky) , Lat= ' num2str(latitudes(lat_index))],'FontSize',14);
    %set(gca,'XDir','reverse');
    set(obj6,'linestyle','none');
    colorbar;
    %caxis([0 3000/Normalization]); Units too small - makes graph hard to
    %interpret
    caxis([0 3000]);
    ylabel(colorbar,'Kelvin', 'Rotation',90, 'FontSize' , 14);
    
    xlabel('LST (h)');
    ylabel('Dec (degrees)');
    h_xlabel=get(gca,'xlabel');
    set(h_xlabel,'FontSize',16);
    h_ylabel=get(gca,'ylabel');
    set(h_ylabel,'FontSize',16);
    set(gca,'FontSize',15);
    
    xlim([ xmin xmax]);
    set(gca,'Xtick', ([xmin:30:xmax]));
    
    ylim([ ymin ymax]);
    set(gca,'Ytick', ([ymin:30:ymax]));
    
    xmin=-180;    xmax=180; xstep=30;
    xlim([ xmin xmax]);
    set(gca,'Xtick', (xmin:xstep:xmax));
    set(gca, 'XTickLabel', ([ 12 14 16 18 20 22 0 2 4 6 8 10 12]));
    %set(gca, 'XTickLabel', ([ 12 10 8 6 4 2 24 22 20 18 16 14 12]));
   
    % prep to save file
    mymap=colormap('jet');
    colormap(mymap);
    c=colorbar;
    axpos = get(gca,'Position');
    cpos = get(c,'Position');
    cpos(3) = 0.5*cpos(3);
    set(c,'Position',cpos);
    set(gca,'Position',axpos);
    
    ylabel(c,'Kelvin', 'Rotation',90, 'FontSize' , 7);
    set(gca,'FontSize',6);
    ylabel('Dec (degrees)', 'FontSize' , 6);
    xlabel('LST (h)', 'FontSize' , 6);
    title(' ');
      
    set(gcf, 'PaperSize', [3.15 2.30]); %Set the paper to have width 5 and height 5.
    set(gcf, 'PaperPosition', [0 0 3.15 2.30]); %Position plot at left hand corner with width 5 and height 5.
    
    if saveforpaper
        filename_mnras = ['C:\Users\Tom\OneDrive\ASU_Thesis\Written_Papers\Spectral_Index\pdf_images\' ];
        print( [filename_mnras 'beam_on_sky_projection'], '-dpdf','-r800') %Save figure
    end
    
    %%
    latitude
    
    figure_counter=figure_counter+10;
end

% P_beam_sec_sinefactor=sinefactor .* P_beam_sec;
% P_beam_sec = P_beam_sec./(sum(P_beam_sec_sinefactor(:)));
%
% %Change the frequency of the sky map to mach the frequency of the beam.
% sky(: , :) = sky_square(v_i, : , : );
% image=P_beam_sec.*sky.*sinefactor;
% Final_spectrum_v(lat_index,hr_index,v_i) = sum(image(:));







%