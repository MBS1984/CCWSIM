%% This script performs 2D multiple-point simulation for binary variables using CCWSIM algorithm.
% 
% Manuscript: "CCWSIM: An efficient and fast wavelet-based CCSIM for binary 
% characterizing in large-scale geological domains"
%--------------------------------------------------------------------------
% Author: Mojtaba Bavandsavadkoohi
% E-mail: Mojtaba.Bavand@inrs.ca
% Institut National de la Recherche Scientifique
% Centre Eau Terre Environnement, Québec, Québec, Canada
%--------------------------------------------------------------------------
% SIMULATION INPUT PARAMETERS:
%
%   ti             : A 2D training image to be imported. Three TIs with different sizes are available:
%                     ti_channel, ti_stoneWall, and ti_gang.
%   hd             : For unconditional simulation, a NaN matrix with the size of 'ti' should be defined, for example, NaN(1000).
%                     For conditional simulation, a simulation grid with conditional data should be imported.
%                     Diffrent conditional data are available according to
%                     available TIs.
%   CT             : Determination of the Co-Template area size. Normally, [2 2] is a suitable size
%                     for finding hard data in the Co-Template area.       
%   OL             : Size of the overlap area between the desired pattern and the previously simulated 
%                     pattern. This parameter should be determined based on wavelet decomposition level. For detailed information,
%                     users are invited to see the README file.                                 
%   cand           : Number of candidate patterns selected based on minimum similarity distance.
%   fc             : 0 means no need for facies matching. Otherwise, the facies proportion should be provided  
%                     (e.g., for a 2-facies TI, this matrix should be like [0.28 0.72]).
%   mrp            : Multiple random path flag (1:on, 0:off). 
%   T_vibration    : Multiple template size flag (1:on, 0:off).
%   prop           : The proportion of sorted candidates considered for finding the best-matched pattern based on 
%                     hard data inside the template and in the Co-template area.
%   real_numb      : Number of realizations.
%   wavelet_level  : Number of wavelet decompositions. This script can handle up to three levels of the DWT.
%--------------------------------------------------------------------------------------------------------------------------------------------                                 
% Output Parameters:
% - C: Simulated grid.
%---------------------------------------------------------------------------------------------------------------------------------------------
% This program requires the following subroutines:

% - hd_resize_2D: It changes the simulation grid size based on wavelet decomposition level and relocates the hard data.
% - mincut_func: It finds the best minimum error boundary in simulation.
% - mincut: It finds the best minimum error boundary in simulation.
% - combine_2D: It merges the desired pattern with the previous pattern based on the best minimum error boundary. 
% - hist_cat: It matches the facies proportion for categorical variables.
% - CCWSIM_main: It simulates the categorical binary variables in the simulation grid.
% - CCWSIM_2D: It shows the simulation results and runtime.
%%
% Clear workspace, command window, and close all figures
clear; clc; close all;

% Delete existing progress bars
progressBars = findall(0, 'Type', 'figure', 'Tag', 'TMWWaitbar');
delete(progressBars);

% Load training image
load('ti_channel.mat');

% Load conditional hard data
load('hd_channel_400.mat');

% Create a nan simulation grid for unconditional simulation
% hd = nan(512);   

% Define simulation input parameters
CT = [2, 2];           
OL = 8;                 
T = 32;                
cand = 20;             
fc = 0;                
mrp = 1;                
T_vibration = 0;        
prop = 0.1;             
real_numb = 1;          
wavelet_level = 1;      

% Call the simulation function
[C] = CCWSIM_2D(ti, hd, T, OL, CT, prop, cand, fc, mrp, T_vibration, wavelet_level, real_numb);

