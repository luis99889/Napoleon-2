clc;
clear;
close all;

% Example MATLAB script to compute satellite orbit parameters

% %% Define Constants and Inputs
% % Standard gravitational parameter (Earth): mu [km^3/s^2]
% mu = 3.986e5;  
% % Earth radius [km]
% RE = 6371;  
% % Satellite orbit altitude [km]
% h0 = 550;    
% % Orbital inclination [radians] (example: 45°)
% iK = deg2rad(45);  
% % Earth rotation angular speed [rad/s]
% omega_E = 7.2921159e-5;  
% 
% % Time instance: e.g., compute parameters 1 hour into the pass, t in seconds
% t = 3600;  
% % Initial angular position (set by convention)
% alpha0 = 0;    
% 
% %% Compute Orbit Parameters
% [omega, alpha, latitude, longitude] = computeOrbitParameters(mu, RE, h0, iK, omega_E, t, alpha0);
% 
% fprintf('Angular Speed (ω) [rad/s]: %f\n', omega);
% fprintf('Angular Position (α) [rad]: %f\n', alpha);
% fprintf('Sub-satellite Latitude [deg]: %f\n', rad2deg(latitude));
% fprintf('Sub-satellite Longitude [deg]: %f\n', rad2deg(longitude));
% 
% %% Function Definitions
% 
% function [omega, alpha, latitude, longitude] = computeOrbitParameters(mu, RE, h0, iK, omega_E, t, alpha0)
%     % Computes the key orbit parameters based on the geometrical model.
%     %
%     % Inputs:
%     %   mu      - Standard gravitational parameter (km^3/s^2)
%     %   RE      - Earth radius (km)
%     %   h0      - Satellite orbit altitude (km)
%     %   iK      - Orbital inclination (rad)
%     %   omega_E - Earth rotation angular speed (rad/s)
%     %   t       - Time elapsed (s)
%     %   alpha0  - Initial angular position (rad)
%     %
%     % Outputs:
%     %   omega     - Angular speed [rad/s]
%     %   alpha     - Angular position on the orbit at time t [rad]
%     %   latitude  - Sub-satellite point latitude [rad]
%     %   longitude - Sub-satellite point longitude [rad]
% 
%     % 1. Angular speed (using Kepler's third law):
%     %    ω = √(µ / (RE + h0)³)
%     omega = sqrt(mu / (RE + h0)^3);
% 
%     % 2. Angular position of the satellite:
%     %    α(t) = ω * t + α₀
%     alpha = omega * t + alpha0;
% 
%     % 3. Sub-satellite latitude:
%     %    λ = arcsin( sin(iK) * sin(α) )
%     latitude = asin(sin(iK) * sin(alpha));
% 
%     % 4. Compute intermediate phase for longitude:
%     %    eϕ = atan2( cos(iK)* sin(α), cos(α) )
%     ephi = atan2(cos(iK) * sin(alpha), cos(alpha));
% 
%     % 5. Sub-satellite longitude considering Earth's rotation:
%     %    ϕ(t) = mod( eϕ - ω_E * t + π, 2π ) - π
%     longitude = mod(ephi - omega_E * t + pi, 2*pi) - pi;
% end
% 
% % %% Create a tle file
% % 
% % % Define computed or assumed orbital parameters:
% % satName      = 'Satellite 01';
% % satelliteID  = '25544';       % Example NORAD catalog number
% % classification = 'U';         % Unclassified
% % 
% % % For epoch, suppose you computed:
% % epochStr = '25061.34958333'; % Format: YYDDD.FFFFFF (year 2025, 61st day, fraction of day)
% % 
% % % Example numerical parameters (all values must be formatted to match TLE field widths)
% % firstDerivative   = '.00001234';  % Options for drag term, etc.
% % secondDerivative  = '00000-0';
% % BSTAR             = '32123-4';
% % ephemerisType     = '0';
% % elementSetNumber  = '9991';
% % 
% % % Create line 1 (using proper spacing and fixed widths)
% % line1 = sprintf('1 %s%s %s   %s  %s  %s %s %s', ...
% %     satelliteID, classification, '98067A', epochStr, firstDerivative, secondDerivative, BSTAR, elementSetNumber);
% % 
% % % For line 2, define:
% % inclination    = 45.0000;    % degrees
% % RAAN           = 21.5993;    % degrees, can be computed from your model
% % eccentricity   = '0000000';  % circular orbit
% % argPerigee     = 0.0000;     % degrees (for circular orbit, arbitrary)
% % meanAnomaly    = 11.5090;    % degrees, computed orbital position
% % meanMotion     = 15.48157;   % revolutions per day, computed from angular speed
% % 
% % % Create line 2. Note that eccentricity is inserted as a formatted string.
% % line2 = sprintf('2 %s %8.4f %8.4f %7s %8.4f %8.4f %11.8f', ...
% %     satelliteID, inclination, RAAN, eccentricity, argPerigee, meanAnomaly, meanMotion);
% % 
% % % Optional: you could add a renormalization step to compute checksum for each line.
% % 
% % % Assemble full TLE content, including title line:
% % tleContent = sprintf('%s\n%s\n%s\n', satName, line1, line2);
% % 
% % % Write the TLE content to a file with a .tle extension
% % tleFileName = 'satellite_01.tle';
% % fid = fopen(tleFileName, 'wt');
% % if fid == -1
% %     error('Unable to open file for writing: %s', tleFileName);
% % end
% % fprintf(fid, '%s', tleContent);
% % fclose(fid);
% % 
% % fprintf('TLE file "%s" created successfully.\n', tleFileName);
% 
% %% Initializing the satellite
% startTime = datetime(2025,3,02,8,23,0);
% stopTime = startTime + hours(4);
% sampleTime = 60;
% sc = satelliteScenario(startTime,stopTime,sampleTime);
% 
% %% Define the Orbital (Keplerian) Parameters
% % These parameters are defined in the Geocentric Celestial Reference Frame (GCRF).
% % Example values (adjust as needed) 
% semiMajorAxis = [7171000, 7171000, 7171000, 7171000];             % km, approx. Earth radius (6371 km) + altitude (400 km)
% eccentricity = [0, 0, 0, 0];                 % Circular orbit: 0
% inclination = [90, 90, 90, 90];                 % Inclination in degrees
% rightAscensionOfAscendingNode = [295, 295, 297, 297]; % RAAN in degrees, example value
% argumentOfPerigee = [235, 45, 235, 45];           % Argument of perigee in degrees, example value
% trueAnomaly = [0, 0, 0, 0];                  % True anomaly in degrees
% 
% %% Add the Satellite Using the Geometrical Parameters
% % The syntax below uses the six orbital parameters:
% sat = satellite(sc, semiMajorAxis, eccentricity, inclination, ...
%                 rightAscensionOfAscendingNode, argumentOfPerigee, trueAnomaly);
% 
% gu_lat = 90;
% gu_long = 90;
% 
% groundStation(sc, gu_lat, gu_long);
% 
% % ac = access(sat,gs);
% % intvls = accessIntervals(ac);
% % 
% % %% Visualize the Scenario
% % v = satelliteScenarioViewer(sc);
% % play(sc);
% 
% show(sat)
% groundTrack(sat,LeadTime=7200)
% 
% % tleFile = "C:\Users\Nuovo Utente\Uni_Project\04_WPs\WP2\satellite_01.tle";
% % sat = satellite(sc, tleFile);
% 
% % Calculate the total duration of the scenario in seconds
% duration = stopTime - startTime;
% totalSeconds = seconds(duration);
% 
% % Create the GroundTrack object for the first satellite
% gt = groundTrack(sat(1), 'LeadTime', totalSeconds);
% 
% % Extract arrays of latitude and longitude of the SSP
% lat = gt.Latitude;    % Latitude array (degrees)
% lon = gt.Longitude;   % Longitude array (degrees)
% 
% % (Optional) Get the corresponding time array
% time = startTime + seconds(gt.Time);
% 
% % Example really important about LEOs
% % openExample('satcom/AnalyzeNTNCoverageAndCapacityExample')
% 
% 
% % investigate the geoTrajectory function to obtain the satellite
% % coordinates through time


%%  Grok generated

clc;
clear;
close all;

% Initializing the satellite
startTime = datetime(2025,3,02,8,23,0);
stopTime = startTime + hours(4);
sampleTime = 60;
sc = satelliteScenario(startTime,stopTime,sampleTime);

% % Define the Orbital (Keplerian) Parameters
% semiMajorAxis = [7171000, 7171000, 7171000, 7171000]; % meters (not km)
% eccentricity = [0, 0, 0, 0];
% inclination = [90, 90, 90, 90];
% rightAscensionOfAscendingNode = [290, 295, 297, 300];
% argumentOfPerigee = [0, 90, 180, 270];
% trueAnomaly = [0, 0, 0, 0];
% 
% 
% % Add the Satellite
% sat = satellite(sc, semiMajorAxis, eccentricity, inclination, ...
%                 rightAscensionOfAscendingNode, argumentOfPerigee, trueAnomaly);


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
    eccentricity = zeros(1, total_sats);           % Circular orbits
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
