clc;
clear;
close all;

% Define constellation parameters
P = 6;                % Number of orbit planes
S = 11;               % Satellites per orbit
semiMajorAxis = 6371000 + 780000; % Orbital radius in meters (R_E + altitude)
inclination = 86.4;   % degrees

% Define Ground Station parameters
gs_lat = 45.07;
gs_long = 7.69;
gs_alt = 0;

% Call the WP2 function
[closest_sat_coords, closest_sat_indices, closest_sat_elevations_discrete, closest_sat_dists, num_times] = WP2_function(P, S, semiMajorAxis, inclination, gs_lat, gs_long, gs_alt);

%% Call the WP3 function
Sample_Rate = 400; % Hz
seed = 60;
[Sat_Ang_time, full_fadedWave, full_inputWave, full_channel, full_state, full_time] = WP3_function(Sample_Rate, num_times, closest_sat_elevations_discrete, seed);
