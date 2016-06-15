% This creates the Beam Adjustment Factor
%stop
%% 
clear;
% Blade Section


Polarization=0; %Polarization=185; Polarization=0;
use_cmb=false;
save_png_enable=true;
save_mnras_pdf_enable=true;
save_beam_factor=true;

if use_cmb
    CMB=2.725;
    cmb_factor='with_cmb_adjust';
else
    CMB=0;
    cmb_factor='without_cmb_adjust';
end

Antenna='Blade';
band='high_band';

v_initial=90;
v_final=200;
v_step=v_final-v_initial+1;

cd 'C:\Users\Tom\OneDrive\ASU_Thesis\Beams\Beam_on_Sky\Results\Blade_highband';
load(['Blade_highband_Spectrum_Pm' num2str(Polarization) '_90to200_L-26_90MHz-200MHz.mat'])
BladeLat26S(:,:)=Final_spectrum_v(65,:,v_initial:v_final); % data available for v_init to v_final MHz
frequency_150=150;
BladeLat26S_flat = zeros(361,v_step); % simulation data for the beam
BladeLat26S_flat_LST_binned=zeros(72,v_step);      % LST from 1/6 hr to 23 5/6 hrs
BladeLat26S_flat_LST_freq_binned=zeros(72,250); % (90.2 to 189.8 in 0.4 MHz steps)
adjustment_hb=zeros(72,250); % 72 LST values, 250 frequency values 90.2 to 189.8
adjustment_lb=zeros(72,125); % 72 LST values, 250 frequency values 90.2 to 189.8

% Scale back down to the original lower Haslam Temperatures of 408 so the CMB can be
% subtracted before scaling to 150 MHz and then added back in at 150 MHz
% Step 1 - Original Haslam Map has the CMB. So scale back the 90 to 190 MHz Spectra to 408 MHz and then subtract the CMB

for frequency = 90:190 % undo the frequency scaling
    frequency_factor=(frequency/408)^2.5;
    for LST=1:361
        BladeLat26S_flat(LST,frequency-89) =  frequency_factor * BladeLat26S(LST,frequency-89)-CMB;
    end
end
% Step 2 Scale everything back to 150 MHz and then add back the CMB
BladeLat26S_flat=BladeLat26S_flat*(408/frequency_150)^2.5 + CMB;

% Step 3 Normalize at 150 MHz

Norm=BladeLat26S_flat(:,frequency_150-89); % index 1 --> 90 MHz thus index 161 --> 150 MHz
for frequency = 90:190 % normalize
    for LST=1:361
        BladeLat26S_flat(LST,frequency-89) =  BladeLat26S_flat(LST,frequency-89)/Norm(LST);
    end
end


adjustment_hb_lst_0to24=zeros(361,v_step);
adjustment_hb_lst_0to24(1:180,:) = BladeLat26S_flat(180:-1:1,:); % LST goes from 12 to 0 then from 24 to 12
adjustment_hb_lst_0to24(181:361,:) = BladeLat26S_flat(361:-1:181,:);
%
% Put 360 LST bins of 4 mins per bin into 72 bins of 20 mins per bin

for frequency =90:190
    for LST=1:72
        BladeLat26S_flat_LST_binned(LST,frequency-89) = sum(BladeLat26S_flat(5*(LST-1)+1:5*(LST-1)+1+5,frequency-89))/6.0;
    end
end


for LST=1:72 % interpolate into 0.4MHz frequency slot
    vq90_190 = interp1(90:190,BladeLat26S_flat_LST_binned(LST,1:101),90.2:0.4:189.8);
    adjustment_hb(LST,:) = vq90_190;
end

% OK, adjustment is listed from LST going from 12 hrs to 24 hrs and then
% from 0 hrs to 12 hrs. Need to fix that
% Old LST went from 180 degrees lat (GMST = 0), downward to 0 degrees lat,
% then to -180 degrees lat --> LST 12 6 0/24 18 12

adjustment_shift(1:36,:) = adjustment_hb(36:-1:1,:); % LST goes from 12 to 0 then from 24 to 12
adjustment_shift(37:72,:) = adjustment_hb(72:-1:37,:);
adjustment_freq=adjustment_shift;
adjustment_hb_final=adjustment_freq; % save variable in a unique name before we clobber it in the low band

if save_beam_factor
    cd 'C:\Users\Tom\OneDrive\ASU_Thesis\Beams\EDGES\high_band\BeamFactor'
    save(['BladeBeamAdjustmentHighBand_P' num2str(Polarization) '_' cmb_factor '.mat'],'adjustment_freq');
end

% Plot the basic beam correction map
figure(1);
clf(1);
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.01 .05 .30 .35])
imagesc(180:-1:-180, 1:101, BladeLat26S_flat(:,1:101)');
mymap=colormap('jet');
colormap(mymap);
c=colorbar;

axpos = get(gca,'Position');
cpos = get(c,'Position');
cpos(3) = 0.5*cpos(3);
set(c,'Position',cpos);
set(gca,'Position',axpos);

set(gca, 'XTickLabel', (''));
set(gca,'FontSize',14);
set(gca,'Xtick', (5:6:71));
set(gca, 'XTickLabel', (1.5:2:23.5));

% LST axis
xmin=-180; xmax=180;
xlim([ xmin xmax]);

caxis([ 0.90 1.10]);
set(c,'YTick',0.90:0.05:1.15)
set(c,'YTickLabel',{ '0.90' '0.95' '1.00' '1.05' '1.10' '1.15'})
    
set(gca,'Xtick', (xmin:30:xmax));
set(gca, 'XTickLabel', ([ 12 14 16 18 20 22 0 2 4 6 8 10 12]));
set(gca,'FontSize',15);

ylabel(c,'Beam Factor', 'Rotation',90, 'FontSize' , 14);

ylabel('Frequency', 'FontSize' , 14);
set(gca,'Ytick', (1:10:101));
set(gca, 'YTickLabel', ([ 90:10:190]));

xlabel('LST (hr)', 'FontSize' , 14);
title(['Beam Factor Relative to 150 MHz'], 'FontSize' , 16);

if save_mnras_pdf_enable
    set(gcf, 'PaperSize', [3.15 2.90]); %Set the paper to have width 3.15 and height 2.9.
    ylabel(c,'Beam Factor', 'Rotation',90, 'FontSize' , 7);
    ylabel('Frequency (MHz)', 'FontSize' , 7);
    caxis([ 0.90 1.10]);
    set(c,'YTick',0.90:0.05:1.10)
    set(c,'YTickLabel',{ '0.90' '0.95' '1.00' '1.05' '1.10'})
    xlabel('LST (hr)', 'FontSize' , 7);
    set(gca,'FontSize',7);
    title(' ');
    set(gcf, 'PaperPosition', [0 0 3.15 2.90]); %Position plot at left hand corner with width 5 and height 5.
    filename_mnras = ['C:\Users\Tom\OneDrive\ASU_Thesis\Written_Papers\Spectral_Index\pdf_images\' Antenna ];
    print( [filename_mnras '_Pm' num2str(Polarization) '_BeamFactor_' cmb_factor ], '-dpdf','-r800') %Save figure
end

if save_png_enable
    set(gcf, 'PaperSize', [7.5 6]); %Set the paper to have width 7.5 and height 6.
    set(gca,'FontSize',14);
    ylabel(c,'Beam Factor', 'Rotation',90, 'FontSize' , 14);
    ylabel('Frequency (MHz)', 'FontSize' , 14);
    xlabel('LST (hr)', 'FontSize' , 14);
    title([Antenna ' Beam Factor Relative to 150 MHz'], 'FontSize' , 16);
    set(gcf, 'PaperPosition', [0 0 7.5 6]); %Position plot at left hand corner with width 7.5 and height 6.
    savefilename=['C:\Users\Tom\OneDrive\ASU_Thesis\Beams\EDGES\' band '\BeamFactor\' Antenna];
    saveas(gcf, [savefilename '_Pm' num2str(Polarization) '_BeamFactor_' cmb_factor '.png'], 'png');
end

figure(11);
clf(11);
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.01 .40 .30 .35])

imagesc(1:361, 1:101, adjustment_hb_lst_0to24(:,1:101)');
mymap=colormap('jet');
colormap(mymap);
c=colorbar;

axpos = get(gca,'Position');
cpos = get(c,'Position');
cpos(3) = 0.5*cpos(3);
set(c,'Position',cpos);
set(gca,'Position',axpos);

set(gca, 'XTickLabel', (''));
set(gca,'FontSize',14);

% LST axis
xmin=1; xmax=361;
xlim([ xmin xmax]);

caxis([ 0.90 1.10]);
set(c,'YTick',0.90:0.05:1.15)
set(c,'YTickLabel',{ '0.90' '0.95' '1.00' '1.05' '1.10' '1.15'})
    
set(gca,'Xtick', (xmin:30:xmax));
set(gca, 'XTickLabel', ([  0 2 4 6 8 10 12 14 16 18 20 22 24]));
set(gca,'FontSize',15);

ylabel(c,'Beam Factor', 'Rotation',90, 'FontSize' , 14);

ylabel('Frequency', 'FontSize' , 14);
set(gca,'Ytick', (1:10:101));
set(gca, 'YTickLabel', ([ 90:10:190]));

xlabel('LST (hr)', 'FontSize' , 14);
title(['Beam Factor Relative to 150 MHz'], 'FontSize' , 16);

if save_mnras_pdf_enable
    set(gcf, 'PaperSize', [3.15 2.90]); %Set the paper to have width 3.15 and height 2.9.
    ylabel(c,'Beam Factor', 'Rotation',90, 'FontSize' , 7);
    ylabel('Frequency (MHz)', 'FontSize' , 7);
    caxis([ 0.90 1.10]);
    set(c,'YTick',0.90:0.05:1.10)
    set(c,'YTickLabel',{ '0.90' '0.95' '1.00' '1.05' '1.10'})
    xlabel('LST (hr)', 'FontSize' , 7);
    set(gca,'FontSize',7);
    title(' ');
    set(gcf, 'PaperPosition', [0 0 3.15 2.90]); %Position plot at left hand corner with width 5 and height 5.
    filename_mnras = ['C:\Users\Tom\OneDrive\ASU_Thesis\Written_Papers\Spectral_Index\pdf_images\' Antenna ];
    print( [filename_mnras '_Pm' num2str(Polarization) '_BeamFactor_LST_0-24_' cmb_factor], '-dpdf','-r800') %Save figure
end

if save_png_enable
    set(gcf, 'PaperSize', [7.5 6]); %Set the paper to have width 7.5 and height 6.
    set(gca,'FontSize',14);
    ylabel(c,'Beam Factor', 'Rotation',90, 'FontSize' , 14);
    ylabel('Frequency (MHz)', 'FontSize' , 14);
    xlabel('LST (hr)', 'FontSize' , 14);
    title([Antenna ' Beam Factor Relative to 150 MHz'], 'FontSize' , 16);
    set(gcf, 'PaperPosition', [0 0 7.5 6]); %Position plot at left hand corner with width 7.5 and height 6.
    savefilename=['C:\Users\Tom\OneDrive\ASU_Thesis\Beams\EDGES\' band '\BeamFactor\' Antenna];
    saveas(gcf, [savefilename '_Pm' num2str(Polarization) '_BeamFactor_LST_0-24_' cmb_factor '.png'], 'png');
end

stop
%%
% Fourpoint Section
clear;
Antenna='Fourpoint';
band='high_band';
save_png_enable=true;
save_mnras_pdf_enable=false;

cd 'C:\Users\Tom\OneDrive\ASU_Thesis\Beams\Beam_on_Sky\Results\EDGES\November 2014 Mesh Study';
load('EDGES_Spectrum_Pm0_100to200_L-26_Run07_pass4.mat')
FourpointLat26S(:,:)=Final_spectrum_v(65,:,100:200);
frequency_150=150;
FourpointLat26S_flat = zeros(361,101);
FourpointLat26S_flat_LST_binned=zeros(72,101);      % LST from 1/6 hr to 23 5/6 hrs
adjustment=zeros(72,250); % 72 LST values, 250 frequency values 90.2 to 189.8 25 values per 10 MHz

for frequency = 100:200
    frequency_factor=(frequency/frequency_150)^2.5;
    for LST=1:361
        FourpointLat26S_flat(LST,frequency-99) =  frequency_factor * FourpointLat26S(LST,frequency-99)/FourpointLat26S(LST,frequency_150-99);
    end
end

% 
% Put 360 LST bins of 4 mins per bin into 72 bins of 20 mins per bin

for frequency = 100:200
    for LST=1:72
        FourpointLat26S_flat_LST_binned(LST,frequency-99) = sum(FourpointLat26S_flat(5*(LST-1)+1:5*(LST-1)+1+5,frequency-99))/6.0;
    end
end


for LST=1:72 % interpolate into 0.4MHz frequency slot
    vq = interp1(100:200,FourpointLat26S_flat_LST_binned(LST,1:101),100.2:0.4:199.8);
    adjustment(LST,:)=vq;
end

% OK, adjustment is listed from LST going from 12 hrs to 24 hrs and then
% from 0 hrs to 12 hrs. Need to fix that
% Old LST went from 180 degrees lat (GMST = 0), downward to 0 degrees lat,
% then to -180 degrees lat --> LST 12 6 0/24 18 12

adjustment_shift(1:36,:) = adjustment(36:-1:1,:);
adjustment_shift(37:72,:) = adjustment(72:-1:37,:);
% 
adjustment_freq(1:72,1:25)=1; % 90 to 100 MHz is set to 1.0 as no data is available
adjustment_freq(1:72,26:250)=adjustment_shift(1:72,1:225); % 100 MHz to 190 MHz data goes here

cd 'C:\Users\Tom\OneDrive\ASU_Thesis\Beams\EDGES\v1\Fourpoint_2015'
save('FourpointBeamAdjustmentHighBand.mat','adjustment_freq');

% Plot the basic beam correction map
figure(2);
clf(2);
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.01 .05 .30 .35])
imagesc(180:-1:-180, 1:length(FourpointLat26S_flat(1,:)), FourpointLat26S_flat');
colormap('default');
c=colorbar;

axpos = get(gca,'Position');
cpos = get(c,'Position');
cpos(3) = 0.5*cpos(3);
set(c,'Position',cpos);
set(gca,'Position',axpos);

set(gca, 'XTickLabel', (''));
set(gca,'FontSize',14);
set(gca,'Xtick', (5:6:71));
set(gca, 'XTickLabel', (1.5:2:23.5));

% LST axis
xmin=-180; xmax=180;
xlim([ xmin xmax]);
set(gca,'Xtick', (xmin:30:xmax));
set(gca, 'XTickLabel', ([ 12 14 16 18 20 22 0 2 4 6 8 10 12]));
set(gca,'FontSize',15);
caxis([ 0.90 1.15]);

ylabel(c,'Beam Factor', 'Rotation',90, 'FontSize' , 14);

ylabel('Frequency', 'FontSize' , 14);
set(gca,'Ytick', (1:10:101));
set(gca, 'YTickLabel', ([ 100:10:200]));

xlabel('LST (hr)', 'FontSize' , 14);
title([Antenna ' Beam Factor Relative to 150 MHz'], 'FontSize' , 16);

if save_mnras_pdf_enable
    set(gcf, 'PaperSize', [3.15 2.90]); %Set the paper to have width 3.15 and height 2.9.
    ylabel(c,'Beam Factor', 'Rotation',90, 'FontSize' , 8);
    ylabel('Frequency)', 'FontSize' , 8);
    xlabel('LST (hr)', 'FontSize' , 8);
    title(' ');
    set(gcf, 'PaperPosition', [0 0 3.15 2.90]); %Position plot at left hand corner with width 5 and height 5.
    filename_mnras = ['C:\Users\Tom\OneDrive\ASU_Thesis\Written_Papers\Spectral_Index\mnras_numbered_images\' Antenna ];
    print( [filename_mnras '_BeamFactor'], '-dpdf','-r800') %Save figure
end

if save_png_enable
    set(gcf, 'PaperSize', [7.5 6]); %Set the paper to have width 7.5 and height 6.
    set(gca,'FontSize',14);
    ylabel(c,'Beam Factor', 'Rotation',90, 'FontSize' , 14);
    ylabel('Frequency (MHz)', 'FontSize' , 14);
    xlabel('LST (hr)', 'FontSize' , 14);
    title([Antenna ' Beam Factor Relative to 150 MHz'], 'FontSize' , 16);
    set(gcf, 'PaperPosition', [0 0 7.5 6]); %Position plot at left hand corner with width 7.5 and height 6.
    savefilename=['C:\Users\Tom\OneDrive\ASU_Thesis\Beams\EDGES\' band '\v1\' Antenna];
    saveas(gcf, [savefilename '_BeamFactor.png'], 'png');
end


%%