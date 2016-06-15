% %
% % Make the spectrum for an Antenna beam at a single latitude
% %
clear; close all;
datestr(now)

%Antenna ='Blade_Finite_Ground'  ;  variation={'5p35x5p35m_with_wings'};
%Antenna ='Blade_Finite_Ground'  ;  variation={'20x20m_no_wings'};
%Antenna = 'EDGES'  ;  variation = {'Run07_pass4'};
Antenna = 'Blade_highband'; variation = {'90MHz-200MHz'};

% If[antennaType == "Blade_highband", 
%   filename = 
%    "Blade_highband\Blade_highband_Spectrum_Pm" <> polarization <> 
%     "_100to200_" <> variation <> ".mat"];
% If[antennaType == "EDGES",      
%   filename = 
%    "EDGES\November 2014 Mesh Study\EDGES_Spectrum_Pm" <> 
%     polarization <> "_100to200_L"  <> userlatitude <> "_" <> 
%     variation <> ".mat"];


% Antenna = 'DARE'; variation={'2015_clipped'};
% Antenna = 'EDGES'; variation={ 'Run07_pass4'};
% Antenna = 'CST_LongDipole';
% Antenna = 'CST_Finite_Dipole';
% Antenna='EDGESPanels';
% Antenna = 'Gaussian';
% Antenna = 'Dipole'
% Antenna = 'Flat'
% Antenna='HFSS_Panels'
% Antenna='EDGES_newcap_run1';
% Antenna='EDGES_newcap_run2';
% Antenna='EDGES_Nov_cap_run1';
% Antenna = 'DipoleLong' % Theoretical Finite Dipole
% variation = {'4_05' '4_10', '8_05' '8_10' '16_05' '16_10' '32_05' '32_10' }
% variation = { '32_05_05'  '32_05_10'  '32_05_25'  '32_05_75' }
% variation={ '32_01_01'  '32_025_025'  '32_05_025' };
% variation={ '32_01_50'  '32_025_50' };
% variation={ 'Run2_L33p4_D0p05_G0p05' 'Run3_L33p4_D0p05_G0p5' };
% variation={'Nov_2014_Run48-56'};


%variation={'2delta'};

% Initialize some variables
        if strcmp(Antenna,'Dare')
            v_initial=40;
            v_final=119;
        else
            v_initial=90;
            v_final=200;
        end
        latitudes = -90:1:90;
        longitudes=180:-1:-180;
        frequency_steps = v_final-v_initial+1;
        Final_spectrum_v=zeros(length(latitudes),length(longitudes),frequency_steps);
        
        v_i=v_initial;

%% For the finite dipole variations
if strcmp(Antenna,'CST_Finite_Dipole')
    cd 'C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/CST_Finite_Dipole';
    files=dir('*.mat');
    temp=struct2cell(files);
    num_files=length(temp(1,:));
    clear variation;
    counter=0;
    for i=1:num_files
        if strfind(files(i).name,'L33')
            if isempty(strfind(files(i).name, 'oct'))
                counter=counter+1;
                variation{counter}=(files(i).name);
            end
        end
    end
    for i=1:num_files
        
        if strfind(files(i).name, 'oct')
            counter=counter+1;
            variation{counter}=(files(i).name);
            
        end
    end
    for i=1:length(variation)
        variation=strrep(variation, 'beam_data_CST_Finite_Dipole_lin_norm_1deg_','');
        variation=strrep(variation, '.mat','');
    end
end
num_runs=length(variation);
%%
% this may take 30 seconds
load('C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/sky_square_with_sinefactor_40to200.mat');

for index=1:length(variation)
    display(variation(index));
        
    tic
    % the date chosen sets ThetaGMST to 360 degrees.
    time='2000/03/20 12:06:40.59856';
    [d_year, d_month, d_day, d_hour, d_minute, d_second] = datevec(datenum(time,'yyyy/mm/dd HH:MM:SS'));
    JD = juliandate(d_year,d_month,d_day,d_hour,d_minute,d_second);
    T_UT1 = (JD-2451545)./36525;
    ThetaGMST = 67310.54841 + (876600*3600 + 8640184.812866).*T_UT1 ...
        + .093104.*(T_UT1.^2) - (6.2*10^-6).*(T_UT1.^3);
    ThetaGMST = mod((mod(ThetaGMST,86400*(ThetaGMST./abs(ThetaGMST)))/240),360);
    % ThetaGMST = 360 degrees
    
    % Get the sky data - only the dec_out and ra_out are used.
    % the sky data was preprocessed, scaled and stored by frequencyto save time
    tic
    
    % This loads, RA_sky, Dec_sky, and sky_square(frequency,Dec,RA)
    sinefactor_dec=cosd(Dec_sky);
    toc
    % step size of 1/3 degree - [541 x 1081] (Dec x RA)
    % RA block looks like:                                 Dec block looks like:
    % 180.00  179.67 ... 0  ... -179.67 -180.00  line 1     -90.0 -90.0  -90.0  ... -90.0
    % 180.00  179.67 ... 0  ... -179.67 -180.00  line 2     -89.7 -89.7  -89.7  ... -89.7
    % ...  ...  ...  ...    ...  ...      ...     ...         ...  ...    ...   ...  ...
    % 180.00  179.67 ... 0  ... -179.67 -180.00  line 271      0    0      0    ...   0
    % ...  ...  ...  ...    ...  ...      ...     ...         ...  ...    ...   ...  ...
    % 180.00  179.67 ... 0  ... -179.67 -180.00  line 540     90    90    90    ...   90
    % col1    col2   col541 ... col1080  col1081 line 541    col1  col2  col3   ... col1081
    
    % SkyMap Data
    % index(Dec,RA)                            Corresponding Values of index
    % (1,1)     (1,2)    ...   (1,1081)   |    (-90.00,180)     (-90.00,179.67)  ... (-90.00,-180)
    % (2,1)     (2,2)    ...   (2,1081)   |    (-89.67,180)     (-89.67,179.67)  ... (-89.67,-180)
    %  ...       ...     ...     ...      |     ...                 ...          ...
    % (541,1)   (541,2)  ...  (541,1081)  |    (+90.00,180)     (+90.00,179.67)  ... (+90.00,-180)
    
    %load('C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/EDGES/beam_data_EDGES_linear_normed_1deg_step_without_LNA_box_June.mat');
    %load('C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/EDGES/beam_data_square_EDGES_linear_normalized');
    %load('C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/EDGES/beam_data_EDGES_linear_normed_1deg_step_with_one_LNA_box_June.mat');
    tic
    if (strcmp(Antenna,'EDGES') || strcmp(Antenna,'EDGESPanels'))
        %load(['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/' Antenna  '/beam_data_' Antenna '_linear_normed_1deg_step_without_LNA_box_June.mat']);
        %load(['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/' Antenna  '/beam_data_' Antenna '_linear_normed_1deg_step_' variation{index} '.mat']);
        load(['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/' Antenna  '/November 2014 Mesh Study/beam_data_' Antenna '_linear_normed_1deg_step_Nov_2014_Mesh_Study_' variation{index} '.mat']);
        % Use with new top cap and delta studies
         % load(['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/' Antenna  '/November 2014 Mesh Study/beam_data_' Antenna '_linear_normed_1deg_step_Nov_2014_TopCap_' variation{index} '.mat']);

    end
    if strcmp(Antenna,'Blade_Finite_Ground')
        load(['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/' Antenna '/beam_data_Normed_linear_Highband_Blade_Ground_' variation{index} '.mat']);
    end
    if strcmp(Antenna,'Blade_highband')
        load(['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/' Antenna '/beam_data_Normed_linear_Highband_Blade_' variation{index} '.mat']);
    end
    if strcmp(Antenna,'DARE')
        load(['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/Dare/DAREbeam2015/beam_data_DARE_linear_normed_1deg_step_Sep_' variation{index} '.mat']);
    end
    if strcmp(Antenna,'CST_LongDipole')
        load(['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/' Antenna  '/beam_data_' Antenna '_linear_normed_1deg_step_Sep' variation{index} '.mat']);
    end
    
    if strcmp(Antenna,'CST_Finite_Dipole')
        cd (['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/', Antenna ]);
        
        load( ['beam_data_CST_Finite_Dipole_lin_norm_1deg_' variation{index} '.mat']);
    end
    
    if strcmp(Antenna , 'Gaussian')
        load('C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/Gaussian/beam_data_square_Gaussian_30_sigma_linear_normed_1deg.mat');
    end
    if strcmp(Antenna, 'Dipole')
        load('C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/Dipole/beam_data_square_Dipole_linear_normalized_176p7MHz_1deg.mat');
    end
    if strcmp(Antenna, 'Flat')
        Beam_Width=160;
        load(['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/Flat/beam_data_square_Flat_linear_normalized_' num2str(Beam_Width) 'FWHM_1deg.mat']);
    end
    if strcmp(Antenna, 'HFSS_Panels')
        load('C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/HFSS_Panels/beam_data_EDGES_linear_normed_1deg_step_HFSS_Panels_July.mat');
    end
    if strcmp(Antenna, 'EDGES_Nov_cap_run1')
        load('C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/EDGES/beam_data_EDGES_linear_normed_1deg_step_Nov_Run1.mat');
    end
    if strcmp(Antenna, 'DipoleLong')
        load('C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_sky/DipoleLong/testbeam3_with_ground_normalized.mat');
    end
    toc
    tic
    Az_sec=zeros(181,361);
    El_sec=zeros(181,361);
    for Az=360:-1.0:0
        Az_sec(:,360-Az+1)=Az;
    end
    for El=90:-1:-90
        El_sec(90-El+1,:)=El;
    end
    toc
    % step size of 1 - [181 x 361]
    % Az block looks like:                  el block looks like:
    % 360  359  358  ... 0      line 1       90   90   90   90  ...  90
    % 360  359  358  ... 0      line 2       89   89   89   89  ...  89
    % ...  ...  ...  ... ...     ...        ...  ...  ...  ...  ...  ...
    % 360  359  358  ... 0      line 181    -90  -90  -90  -90  ...  -90
    % col1 col2 col3 ... col361             col1 col2 col3 col4 ... col361
    
    % Beam Data
    % index(El,Az)                            Corresponding Values of index
    % (1,1)     (1,2)    ...   (1,361)   |    (+90,360)     (+90,359)  ... (+90,0)
    % (2,1)     (2,2)    ...   (2,361)   |    (+89,360)     (+89,359)  ... (+89,0)
    %  ...       ...     ...     ...     |     ...              ...          ...
    % (181,1)   (181,2)  ...  (181,361)  |    (-90,360)     (-90,359)  ... (-90,0)
    toc
    
    for Polarization=[0 5 180 185 ] % as defined as a counterclockwise rotation
        
        display(Polarization);
        Polarization_step=Polarization; % rotation columns = Polarization number
        
        %%
        tic
        % The goal of the for loop is to make a pcolor plot at different frequencies
        % in the range between 100MHz and 190MHz. The variable v is the frequency.
        for latitude_index=65:65
            latitude = latitudes(latitude_index);
            
            for longitude_index=1:361
                longitude = longitudes(longitude_index);
                %Convert the sky map coordinates to Az and El.
                [Az_sky, El_sky]=RaDec2AzEl_fresh_orig_mod1(RA_sky, Dec_sky, latitude,longitude,ThetaGMST);
                
                for i=1:1:frequency_steps
                    
                    %Get the beam block of data
                    beam_sec(:,:) = beam_data_square(v_i,:,:);
                    
                    %  Polarization rotates the beam to the right position -
                    % UNCOMMENT TO ACTIVATE
                                     beam_sec(:,1:360)=circshift(beam_sec(:,1:360), [0,+Polarization_step]);
                                     beam_sec(:,361)=beam_sec(:,1);
                    
                    % zero out anything below the horizon
                    % ind = find(el_sec<0);
                    % beam_sec(ind) = 0;
                    
                    %Find the beam value for the (RA,Dec) pairs that have been mapped to Az/El earlier
                    P_beam_sec=interp2(Az_sec,El_sec, beam_sec, Az_sky,El_sky);
                    Normalization=sum(sum(sinefactor_dec.*P_beam_sec));
                    P_beam_sec = P_beam_sec/Normalization;
                    %Change the frequency of the sky map to mach the frequency of the beam.
                    sky_with_sinefactor(:,:)=sky_square(v_i,:,:);
                    image=P_beam_sec.*sky_with_sinefactor;
                    Final_spectrum_v(latitude_index,longitude_index,v_i) = sum(image(:));
                    v_i=v_i+1;
                    
                end % i loop over specified frequencies
                v_i=v_initial;
            end % hour index
            
            % Save the Final Spectrum data
            if strcmp(Antenna, 'Flat')
                savefilename=['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/Results/' Antenna '/' Antenna '_' num2str(Beam_Width) '_Spectrum_Pm' ...
                    num2str(Polarization) '_' num2str(v_initial) 'to' num2str(v_final) '_L' num2str(latitude) '_noLNA_June_1deg.mat'];
            elseif strcmp(Antenna, 'HFSS_Panels')
                savefilename=['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/Results/' Antenna '/' Antenna '_Spectrum_Pm' ...
                    num2str(Polarization) '_' num2str(v_initial) 'to' num2str(v_final) '_L' num2str(latitude) '_noLNA_July_1deg.mat'];
            elseif strcmp(Antenna,'CST_LongDipole')
                savefilename=['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/Results/' Antenna '/' Antenna '_Spectrum_Pm' ...
                    num2str(Polarization) '_' num2str(v_initial) 'to' num2str(v_final) '_L' num2str(latitude) '_noLNA_Sep' variation{index}  '_1deg.mat'];
            elseif strcmp(Antenna,'CST_Finite_Dipole')
                savefilename=['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/Results/' Antenna '/' Antenna '_Spectrum_Pm' ...
                    num2str(Polarization) '_' num2str(v_initial) 'to' num2str(v_final) '_L' num2str(latitude) '_' variation{index}  '_1deg.mat'];
            elseif strcmp(Antenna, 'EDGES_Nov_cap_run1')
                savefilename=['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/Results/EDGES_newcap/EDGES_Nov_cap_run1_Spectrum_Pm' ...
                    num2str(Polarization) '_' num2str(v_initial) 'to' num2str(v_final) '_L' num2str(latitude) '_noLNA_Aug_1deg.mat'];
            elseif strcmp(Antenna, 'DARE')
                savefilename=['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/Results/Dare/DARE_Spectrum_Pm' ...
                    num2str(Polarization) '_' num2str(v_initial) 'to' num2str(v_final) '_L' num2str(latitude) '_' variation{index} '.mat'];
            elseif strcmp(Antenna, 'Blade_Finite_Ground')
                savefilename=['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/Results/' Antenna '/' Antenna '_Spectrum_Pm' ...
                    num2str(Polarization) '_' num2str(v_initial) 'to' num2str(v_final) '_L' num2str(latitude) '_' variation{index} '.mat'];
            elseif strcmp(Antenna, 'Blade_highband')
                savefilename=['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/Results/' Antenna '/' Antenna '_Spectrum_Pm' ...
                    num2str(Polarization) '_' num2str(v_initial) 'to' num2str(v_final) '_L' num2str(latitude) '_' variation{index} '.mat'];
           else
                savefilename=['C:/Users/Tom/OneDrive/ASU_Thesis/Beams/Beam_on_Sky/Results/' Antenna '/November 2014 Mesh Study/' Antenna '_Spectrum_Pm' ...
                    num2str(Polarization) '_' num2str(v_initial) 'to' num2str(v_final) '_L' num2str(latitude) '_' variation{index} '.mat'];
            end
            
            save(savefilename,  'Final_spectrum_v');
            
            
            display(sprintf('Polarization of %s degrees at latitude %s degrees has finished.\n',  num2str(Polarization), num2str(latitude)));
        end % latitude index
        toc
        
        %
        
    end % Polarization loop
end % variation loop
datestr(now)



