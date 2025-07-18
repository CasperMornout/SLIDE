%% Script to perform GBS/GBO identification, using the SLIDE method, on SEM-DIC/EBSD/Profilometry data of a Zn coating.

% This script uses the same data as in the publication "SLIDE: Automated 
% identification and quantification of grain boundary sliding and opening 
% in 3D" (https://doi.org/10.1016/j.scriptamat.2025.116861)

% Author: C.J.A. Mornout, c.j.a.mornout@tue.nl
% // Eindhoven University of Technology, Hoefnagels Group

% Date: 17-07-2025
% the latest version of this code can be found on
% https://github.com/CasperMornout/SLIDE

% MTEX is required to use this code
% This code has been tested to run for MTEX Version 5.11.2

%% Initialization

close all; clearvars; clc;

% setting MTEX preferences
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoplane');

% add paths in the working directory
addpath(genpath(pwd))

% load struct with all aligned datasets
% The file has been split up to cope with Github file size limitations
File1 = load('Data/ZnCoating_Data_1.mat');
File2 = load('Data/ZnCoating_Data_2.mat');

ZnCoating_DataStruct = struct;
ZnCoating_DataStruct.ebsd = File1.ZnCoating_DataStruct.ebsd;
ZnCoating_DataStruct.grains = File1.ZnCoating_DataStruct.grains;
ZnCoating_DataStruct.CONFOCAL = File2.ZnCoating_DataStruct.CONFOCAL;
ZnCoating_DataStruct.IBSE = File2.ZnCoating_DataStruct.IBSE;

clear File1 File2

%% User Settings

userSettings = struct;

%%%% Settings about grain boundary smoothing %%%%%
userSettings.smoothfactor = 50; % the number of times the grains will be smoothed

%%%% Settings about Integration Path Length (IPL) %%%%
userSettings.IPL_min = 900; % distance of IPL_min parameter in nanometer
userSettings.IPL_min_IBSE = 1800; % distance of normal path length for IBSE data in nanometer
userSettings.IPL_min_CONF = 1800; % distance of normal path length for Confocal data in nanometer

%%%% Experimental Settings %%%%
userSettings.inc_chosen = [1 2 3 4 5]; % increment of deformation. Only for increment 5, Confocal Data is available
userSettings.residualTH = 0.4; % Residual threshold above which data is removed
userSettings.strainTH = 0.1; % strain threshold below which data is removed
userSettings.distPureSliding = 0.1; % length in micrometer along which sliding needs to be identified (if length is lower, the datapoint is removed)

%%%% Opening Identification Settings %%%%
userSettings.openingThreshold = 0.2745; % brightness below which opening is identified
userSettings.numPixUnder = 3; % number of consecutive pixels the IBSE needs to be below the opening threshold to be considered opening
userSettings.nrSegments = 11; % number of line segments used to get average IBSE response per GB segment (line segment spacing = IBSE pixelsize)

%%%% Plotting Settings %%%%
userSettings.plottingIntermediateResults = 1; % select 1 if you want to plot intermediate results

%% PART 0 - Unpacking the aligned data

% load the ebsd and grains variables
ebsd    = ZnCoating_DataStruct.ebsd;
grains  = ZnCoating_DataStruct.grains;

% if selected, plot the ebsd and the grains
if userSettings.plottingIntermediateResults

    close all; h = newMtexFigure;
    plot(ebsd,ebsd.orientations)
    hold on
    plot(grains.boundary)
    mtexTitle('Orientations + GB')

    nextAxis

    plot(ebsd,ebsd.prop.inc_5_Eequi)
    hold on
    plot(grains.boundary,'lineColor','red')
    clim([0.02 1])
    set(gca,'colorscale','log')
    colormap viridis
    mtexColorbar
    mtexTitle('Effective Strain + GB')

    saveas(gcf,'Figures/Data_overview.png')

end

%% PART 1 - Assign GB Normals/Tangentials to Neighbouring EBSD Datapoints

% define distance from the GBs in which we need to assign datapoints a GB
% tangential and normal (corresponds to IPL_min parameter set by user)
userSettings.padding = userSettings.IPL_min + 100;

% Smooth Grains
grains_smooth = grains;
for i = 1:userSettings.smoothfactor; grains_smooth = smooth(grains_smooth); end

% assign the GB normals and tangentials to neighbouring EBSD datapoints
% that are close enough to the GB
boundary_info = assignGB_n_t(ebsd,grains_smooth,userSettings.padding);

if userSettings.plottingIntermediateResults
    % indiviudal plots showing GB normals for each grain, as assigned to
    % the neighbouring pixels
    SLIDE_Plot_GB_Tangentials_and_Normals(boundary_info,grains,ebsd)

    % plot of smoothened GB
    close all
    figure('WindowState','maximized')
    plot(grains_smooth.boundary,'lineWidth',8,'lineColor','red','DisplayName','Smooth GB')
    hold on
    plot(grains.boundary,'lineWidth',2,'lineColor','k','DisplayName','Raw GB')
    saveas(gcf,'Figures/Grains_Smooth.png')
end

%% PART 2.1 - In-Plane Datapointwise SLIDE Execution

% SLIDE_Struct will contain all information and will gradually be filled
SLIDE_Struct = struct;

% Loop over the desired strain increments
for strainInc = userSettings.inc_chosen

    % Adjust the settings to point to the correct strain increment
    userSettings.currentInc = strainInc;

    % Assign strain and displacement fields to feed into SLIDE
    % For this, we use the 'ebsd' map from MTEX
    ebsd.prop.dispUfield = ebsd.prop.(['inc_' num2str(userSettings.currentInc) '_U']);
    ebsd.prop.dispVfield = ebsd.prop.(['inc_' num2str(userSettings.currentInc) '_V']);
    ebsd.prop.strainfield = ebsd.prop.(['inc_' num2str(userSettings.currentInc) '_Eequi']);
    ebsd.prop.U = vector3d(ebsd.prop.dispUfield,ebsd.prop.dispVfield,zeros(size(ebsd.prop.dispUfield)));
    
    % run SLIDE and save information in "SLIDE_struct"
    % This struct will gradually be filled up with the following data:
    % 2.1) Datapoint-wise SLIDE execution, sorted by grain
    % 2.2) GB-Integration of SLIDE results, sorted by grain
    % 3) Opening/Sliding Identification from IBSE data, sorted by grain
    % 4) 3D Sliding displacements from optical profilometry, sorted by grain
    
    % Execute Datapointwise SLIDE for all points close to a GB
    SLIDE_Struct.(['inc_' num2str(userSettings.currentInc)]) = SLIDE_Datapointwise(ebsd,grains,boundary_info);

end

%% PART 2.2 - GB-Integration to Obtain Sliding Displacements for GB

for strainInc = userSettings.inc_chosen

    % Adjust the settings to point to the correct strain increment
    userSettings.currentInc = strainInc;

    GB_integration_input = struct;
    
    GB_integration_input.ebsd = ebsd;
    GB_integration_input.grains = grains_smooth;
    GB_integration_input.SLIDE_Struct = SLIDE_Struct.(['inc_' num2str(userSettings.currentInc)]);
    GB_integration_input.boundary_info = boundary_info;
    GB_integration_input.userSettings = userSettings;
    
    % Excecute GB-integration to obtain sliding displacements for the GB
    SLIDE_Struct.(['inc_' num2str(userSettings.currentInc)]) = SLIDE_GB_integration(GB_integration_input);
    
    clear GB_integration_input

end

%% PART 3 - Identify Opening and Sliding with In-Beam SE Data

for strainInc = userSettings.inc_chosen

    userSettings.currentInc = strainInc;
    IBSE_input = struct;
    
    % Read the X/Y maps for IBSE Forward-Deformed Alignment
    % These fields show, for each pixel in the undeformed EBSD dataset, the corresponding X and Y
    % location in the IBSE dataset, taking into account forward deformed alignment
    IBSE_input.xIBSE = ZnCoating_DataStruct.IBSE.(['inc_' num2str(userSettings.currentInc)]).xIBSE;
    IBSE_input.yIBSE = ZnCoating_DataStruct.IBSE.(['inc_' num2str(userSettings.currentInc)]).yIBSE;

    % add additional input
    IBSE_input.filepath = ['ZnCoating_IBSE_Increment' num2str(userSettings.currentInc) '.tif'];
    IBSE_input.psize =  ZnCoating_DataStruct.IBSE.(['inc_' num2str(userSettings.currentInc)]).psize;
    IBSE_input.userSettings = userSettings;
    IBSE_input.ebsd = ebsd;
    IBSE_input.grains = grains_smooth;
    IBSE_input.SLIDE_Struct = SLIDE_Struct.(['inc_' num2str(userSettings.currentInc)]);
    IBSE_input.boundary_info = boundary_info;
    
    % Identify Opening and sliding with in-beam SE data
    SLIDE_Struct.(['inc_' num2str(userSettings.currentInc)]) = SLIDE_Identify_opening(IBSE_input);
    
    clear IBSE_input

end

%% PART 4 - Compute dZ (3D Sliding Displacement) Using Optical Profilometry
% Only works for the final increment of deformation (increment 5), because
% only for this increment, confocal data is available

userSettings.currentInc = 5;

CONF_input = struct;

% Read the Confocal data
CONF_input.Confocal_Z = ZnCoating_DataStruct.CONFOCAL.confocal_z;
CONF_input.Confocal_X = ZnCoating_DataStruct.CONFOCAL.confocalX;
CONF_input.Confocal_Y = ZnCoating_DataStruct.CONFOCAL.confocalY;

% Read the X/Y maps for Confocal Forward-Deformed Alignment
% These fields show, for each EBSD pixel, the corresponding X and Y
% location in the Confocal dataset
CONF_input.xCONF = ZnCoating_DataStruct.CONFOCAL.xCONF;
CONF_input.yCONF = ZnCoating_DataStruct.CONFOCAL.yCONF;

% add additional inputs
CONF_input.userSettings = userSettings;
CONF_input.ebsd = ebsd;
CONF_input.grains = grains_smooth;
CONF_input.SLIDE_Struct = SLIDE_Struct.(['inc_' num2str(userSettings.currentInc)]);
CONF_input.boundary_info = boundary_info;

% Identify Out-of-plane displacements with profilometry data
SLIDE_Struct.(['inc_' num2str(userSettings.currentInc)]) = SLIDE_Profilometry(CONF_input);

clear CONF_input

%% PART 5 - Post Processing Data
% Combine all the grain-based data for the individual grains into one
% grains2d variable

grains_SLIDE = struct;

for strainInc = userSettings.inc_chosen

    userSettings.currentInc = strainInc;

    input = struct;
    
    input.grains = grains_smooth;
    input.userSettings = userSettings;
    input.boundary_info = boundary_info;
    input.SLIDE_Struct = SLIDE_Struct.(['inc_' num2str(userSettings.currentInc)]);

    if userSettings.currentInc == 5; input.CONFdata = 1; 
    else; input.CONFdata = 0; end
    
    % combine SLIDE results
    grains_SLIDE.(['inc_' num2str(userSettings.currentInc)]) = combineSLIDEresults(input);
    
    clear input

end

%% Compile and Save Results

% Generate a struct to save only the necessary data
SLIDE_Results = struct;

% append the SLIDE results and smooth GBs to the struct
SLIDE_Results.grains_SLIDE = grains_SLIDE;
SLIDE_Results.grains = grains_smooth;

% append the original data to the struct
SLIDE_Results.grains_raw = grains;
SLIDE_Results.IBSE = ZnCoating_DataStruct.IBSE;
SLIDE_Results.CONFOCAL = ZnCoating_DataStruct.CONFOCAL;

% append the user settinsg to the struct
SLIDE_Results.userSettings = userSettings;

% Save Results
save('Data/ZnCoating_Data_Results.mat','SLIDE_Results')