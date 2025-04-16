clc;
clear;
close all;

% Scenario
startTime = datetime(2020,5,1,11,36,0);
stopTime = startTime + days(1);
sampleTime = 60; % in seconds
sc = satelliteScenario(startTime, stopTime, sampleTime);

% Ground station (Torino)
lat = 45.07; 
lon = 7.69;
gs = groundStation(sc, lat, lon);

% Orbital parameters
a = 7171e3; % semi-major axis [m] (Earth radius + 800 km altitude)
e = 0;      % eccentricity
i = 45;     % inclination [deg]
argPeriapsis = 0; % [deg]
trueAnomalies = linspace(0, 360, 8); % 7 sat => 7 step => 8 points, exclude last
trueAnomalies(end) = [];

% Costellation parameters
numPlanes = 3;
numSatsPerPlane = 7;

% RAANs equally spaced
RAANs = linspace(0, 360, numPlanes+1); 
RAANs(end) = [];

% Add satellites
satArray = [];
for p = 1:numPlanes
    RAAN = RAANs(p);
    for s = 1:numSatsPerPlane
        ta = trueAnomalies(s);
        sat = satellite(sc, a, e, i, RAAN, argPeriapsis, ta);
        satArray = [satArray; sat];
        access(sat, gs); % Access analysis
    end
end

% Visualizzazione
play(sc);
