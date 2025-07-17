%% Script to plot results from the SLIDE Method

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
load('Data/ZnCoating_Data_Results.mat');

%% Unpack stuff

userSettings = SLIDE_Results.userSettings;
ebsd_SLIDE = SLIDE_Results.ebsd_SLIDE;
grains_SLIDE = SLIDE_Results.grains_SLIDE;

%% 1 - Overview Plot for SLIDE (3D)

% similar to Figure 5 from the publication
% Only works for increment 5 because profilometry data is not available for
% the other strain increments.

userSettings.currentInc = 5;

input = struct;
input.grains = grains_SLIDE.(['inc_' num2str(userSettings.currentInc)]);
input.ebsd = ebsd_SLIDE;
input.userSettings = userSettings;

SLIDE_plotResults3D(input)

% For saving the figures properly, the export fig toolbox is used
% https://nl.mathworks.com/matlabcentral/fileexchange/23629-export_fig
export_fig('Figures/Results/Overview_3D', '-m2', '-p0.01', '-png')

clear input

%% 2 - Overview Plot for SLIDE (2D)

% example plot for a when 3D (profilometry) data is unavailable
% works for any strain increment

userSettings.currentInc = 1;

input = struct;
input.grains = grains_SLIDE.(['inc_' num2str(userSettings.currentInc)]);
input.ebsd = ebsd_SLIDE;
input.userSettings = userSettings;

SLIDE_plotResults2D(input)
export_fig(['Figures/Results/Overview_2D_Increment_' num2str(userSettings.currentInc)], '-m2', '-p0.01', '-png')

clear input

close all

%% 3 - Evolution Plot for SLIDE

% plots evolution of strain and SLIDE results in 2D

input = struct;
input.grains = grains_SLIDE;
input.ebsd = ebsd_SLIDE;
input.userSettings = userSettings;
input.incList = [1 2 3 4 5]; % do not change

SLIDE_plotEvolution2D(input)
export_fig('Figures/Results/Evolution_SLIDE_2D', '-m2', '-p0.01', '-png')

clear input

%% 4 - Evolution Plot for GBS/GBO

% plots evolution of GBS/GBO along with the IBSE images

input = struct;
input.grains = grains_SLIDE;
input.ebsd = ebsd_SLIDE;
input.userSettings = userSettings;
input.incList = [1 2 3 4 5]; % do not change

SLIDE_plotGBS_GBO_2D(input)
export_fig('Figures/Results/Evolution_GBS_GBO_2D', '-m2', '-p0.01', '-png')

clear input

