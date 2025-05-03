clc;
clear;
close all;

% Define constellation parameters
P = 6;                % Number of orbit planes
S = 11;               % Satellites per orbit
semiMajorAxis = 6371000 + 780000; % Orbital radius in meters (R_E + altitude)
inclination = 86.4;   % degrees

% Call the WP2 function
[closest_sat_coords, closest_sat_indices, closest_sat_elevations_discrete, closest_sat_dists] = ...
    WP2(P, S, semiMajorAxis, inclination, gs_lat, gs_long);

% Call the WP3 function

