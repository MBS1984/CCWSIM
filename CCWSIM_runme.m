% ---------------------------------------------------------------------------------------------------------
% This script can perform 2D multiple-point simulation for 
% Categorical variables through CCWSIM algorithm.
 
% Manuscript: "CCWSIM: An Efficient and Fast Wavelet-Based CCSIM for Categorical 
% Characterization of Large-Scale Geological Domains"
%--------------------------------------------------------------------------
% Author: Mojtaba Bavandsavadkoohi
% E-mail: Mojtaba.Bavand@inrs.ca
% Institut National de la Recherche Scientifique
% centre Eau Terre Environnement, Québec, Québec, Canada
%--------------------------------------------------------------------------
% INPUT PARAMETERS:

%   ti             : A 2D training image should be imported. Three TIs with size are available:
%                     ti_Channel, ti_StoneWall and ti_Delta
%   hd             : for unconditional simulation according to the size of ti, a NaN matrix should be defined for example NaN(2000),               
%                     for conditional simulation, a simulation grid with conditional data should be imported.                                                               
%                     Three HDs are available according to available TIs,respectively: cnd_Channel, cnd_StoneWall and cnd_Delta.
%                     Also, to test the performance of the algorithm in the case of dense hard data,
%                     a file called cnd_Channel_Dense is prepared. 
%   CT             : determination of the Co-Template area size. normally, [2 2] is a suitable size
%                     for finding hard data in Co-Template area.       
%   OL             : size of the overlap area between the desired pattern and previous simulated 
%                     pattern. This parameter should be determined based on wavelet decomposition level. For detail, information
%                     users are invited to see the README file.                                 
%   cand           : number of candidate patterns that are selected based on minimum similarity distance.
%   fc             : 0 equal to no need for facies matching. Otherwise, the facies proportion should be provided  
%                     for example: for a 2 facies TI, this matrix should be like [0.28 0.72])
%   mrp            : Multiple random path flag (1:on, 0:off). 
%   T_vibration    : Multiple template size flag (1:on, 0:off)
%   prop           : the proportion of sorted candidates that are considered for finding best-matched pattern based on 
%                     hard data inside the template and in the co-template area.
%   real_numb      : number of realizations
%   wavelet_level  : number of wavelet decomposition. This script can handle up to three levels of the DWT.
%--------------------------------------------------------------------------------------------------------------------------------------------                                 
% Output Parameters
% - C: output Simulation grid 
% - error: The number mismatch HD ( zero means no mismatch)
%---------------------------------------------------------------------------------------------------------------------------------------------
% This program needs the following subroutines:
% - hd_resize_2D: It changes the simulation grid size based on wavelet decomposition level and relocates the hard data.
% - mincut_func: It finds the best minimum error boundary in simulation.
% - mincut: It finds the best minimum error boundary in simulation.
% - combine_2D: It merges desire pattern with previous pattern based on best minimum error boundary. 
% - hist_cat: It matches the facies proportion for categorical variable.
% - CCWSIM_main: It simulates the categorical binary variables in simulation grid.
% - CCWSIM_2D: It shows the simulation results and runtime. 
%--------------------------------------------------------------------------
%%
% Clear workspace, command window, and close all figures
clear; clc; close all;

% Delete existing progress bars
progressBars = findall(0, 'Type', 'figure', 'Tag', 'TMWWaitbar');
delete(progressBars);

% Load training image
% load('ti_channel.mat');
load('ti_Delta_sund.mat');

% Load conditional hard data
% load('hd_channel_400.mat');
load('hd_Delta_sund_500.mat');
% Create a nan simulation grid for unconditional simulation
% hd = nan(512);   

% Define simulation input parameters
CT = [2, 2];            % Co-Template size
OL = 8;                 % Overlap region
T = 32;                 % Template size
cand = 20;              % Number of Candidate pattern
fc = 0;                 % 0 means no facies matching otherwise proportion should be provided. 
mrp = 1;                % Markov random path
T_vibration = 0;        % Multiple template size flag (1:on, 0:off)
prop = 0.1;             % The proportion of TI which will be scanned to find the matched hard data.
real_numb = 1;          % Nuber of realizations
wavelet_level = 1;      % Wavelet level ( can be 1, 2 or 3)

% Call the simulation function
[C] = CCWSIM_2D(ti, hd, T, OL, CT, prop, cand, fc, mrp, T_vibration, wavelet_level, real_numb);

