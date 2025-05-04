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

%% HANDOVER FUNCTION

Handover_Type = 0;

[closest_sat_elevations_discrete_2] = Handover_2(closest_sat_elevations_discrete,closest_sat_indices);

% if ho 1 in wp3 closest_sat_1 else closest_sat_2

%% Call the WP3 function

CarrierFrequency = 3.8e9;
Sample_Rate = 400; % Hz
seed = 65;

if Handover_Type == 0 % old HO
    [Sat_Ang_time, full_fadedWave, full_inputWave, full_channel, full_state, full_time, all_channel_gains] = WP3_function(Sample_Rate, num_times, closest_sat_elevations_discrete, seed, CarrierFrequency);

else
    [Sat_Ang_time, full_fadedWave, full_inputWave, full_channel, full_state, full_time, all_channel_gains] = WP3_function(Sample_Rate, num_times, closest_sat_elevations_discrete_2, seed, CarrierFrequency);

end

%% Call WP4 function
[Bit_rate,channel_gain_avg,L] = WP4_function(all_channel_gains, closest_sat_dists,CarrierFrequency);
