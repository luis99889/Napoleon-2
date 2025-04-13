clc;
clear;
close all;

% Initializing the satellite scenario
startTime = datetime(2025,3,02,8,23,0);
stopTime = startTime + hours(4);
sampleTime = 60;
sc = satelliteScenario(startTime,stopTime,sampleTime);

% Define constellation parameters
P = 6;                % Number of orbit planes
S = 11;                % Satellites per orbit
semiMajorAxis = 7511000; % meters
inclination = 86.4;     % degrees

% Add satellites using the constellation function
sat = createConstellation(sc, P, S, semiMajorAxis, inclination);

% Turin's coordinates
gu_lat = 45.07;
gu_long = 7.69;
gs = groundStation(sc, gu_lat, gu_long);

% Extract Positions and Analyze Connections
times = startTime : seconds(sampleTime) : stopTime;

% Define dimensions
num_times = numel(times);
num_sats = numel(sat);

% Initialize array to store positions (3D: 3 coords × times × satellites)
all_sat_pos = zeros(3, num_times, num_sats);

all_sat_pos = states(sat, "CoordinateFrame","ecef");


% This loop is not necessary since states(sat) works as well

% % Loop over each time step
% for j = 1:num_times
%     sat_states = states(sat, times(j));        % Get states for all satellites at times(j)
%     all_sat_pos(:, j, :) = sat_states(1:3, :); % Store positions (x, y, z)
% end
% 

gu_alt = 0;

% Ground station position: ECEF 
ecef_pos = lla2ecef([gu_lat, gu_long, gu_alt]);
gs_ecef_pos = repmat(ecef_pos', 1, numel(times));


% utc_times = datevec(times); % Converts datetime to [Y, M, D, H, M, S]
% 
% gs_eci = zeros(3, num_times);
% for j = 1:num_times
%     gs_eci(:, j) = ecef2eci(gs_ecef_pos(:, j), utc_times(j, :));
% end
% gs_eci = ecef2eci(gs_ecef, utc_times);

% % % % gs_eci = ecef2eci(gs_ecef, times);

% Access status
ac = access(sat, gs);
intvls = accessIntervals(ac);

% status = accessStatus(ac, times);
% status = 
% % 
% for j = 1:num_times
%     status_temp = accessStatus(ac, times(j));        % Get states for all satellites at times(j)
%     status(:, j, :) = status_temp(1:3, :); % Store positions (x, y, z)
% end

% Initialize status array to store logical access status at each time
status = false(num_times, num_sats);  % Preallocate with logical false

% Extract status for each time step
for j = 1:num_times
    status(j, :) = accessStatus(ac, times(j));  % Store logical status for each satellite
end


tau=1;

% Find best satellite at each time
for j = 1:numel(times)
    active_sats = find(status(j, :));
    if ~isempty(active_sats)
        distances = zeros(1, length(active_sats));
        for k = 1:length(active_sats)
            sat_idx = active_sats(k);
            sat_pos = all_sat_pos(:, j, sat_idx);
            gs_pos = gs_ecef_pos(:, j);
            distances(k) = norm(gs_pos - sat_pos);
        end
        [min_dist, idx] = min(distances);
        best_sat = active_sats(idx);
        fprintf('At time %s, best satellite is %d with distance %.2f km\n', ...
                string(times(j)), best_sat, min_dist/1000);
    else
        fprintf('At time %s, no satellite has access\n', string(times(j)));
    end
end

% Visualize (optional)
% show(sat);
% groundTrack(sat, LeadTime=7200);
play(sc);


function sat = createConstellation(sc, P, S, semiMajorAxis, inclination)
    % Total number of satellites
    total_sats = P * S;
    
    % Define constant orbital parameters for all satellites
    eccentricity = zeros(1, total_sats);           % Circular orbits with eccentricitt set to zero
    argumentOfPerigee = zeros(1, total_sats);      % Arbitrary for circular orbits, set to 0
    semiMajorAxis = repmat(semiMajorAxis, 1, total_sats); % Same for all
    inclination = repmat(inclination, 1, total_sats);     % Same for all
    
    % Compute RAAN: repeat each plane's RAAN for S satellites
    RAAN = repelem((0:P-1) * 180 / P, S);
    
    % Compute true anomaly: repeat the sequence [0, 360/S, ..., (S-1)*360/S] for P planes
    trueAnomaly = repmat((0:S-1) * 360 / S, 1, P);
    
    % Create the satellite array
    sat = satellite(sc, semiMajorAxis, eccentricity, inclination, RAAN, ...
                    argumentOfPerigee, trueAnomaly);
end
