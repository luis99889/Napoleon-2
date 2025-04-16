clc;
clear;
close all;

% Scenario
startTime = datetime(2020,5,1,9,36,0);
stopTime = startTime + hours(5);
sampleTime = 60; % in seconds
sc = satelliteScenario(startTime, stopTime, sampleTime);

% Ground station (Torino)
lat = 45.07; 
lon = 7.69;
gs = groundStation(sc, lat, lon);

% Parametri orbitali
a = 7171e3; % semi-major axis [m] (Earth radius + 800 km altitude)
e = 0;      % eccentricity
i = 90;     % inclinazione [deg] - orbite polari (passano per i poli)
argPeriapsis = 0; % [deg]
trueAnomalies = linspace(0, 360, 8); % 7 sat => 7 step => 8 points, exclude last
trueAnomalies(end) = [];

% Parametri della costellazione
numPlanes = 3; % Numero di orbite parallele
numSatsPerPlane = 7; % Numero di satelliti per orbita

% RAANs distribuiti uniformemente
RAANs = linspace(0, 90, numPlanes+1); 
RAANs(end) = [];

% Aggiungi i satelliti
satArray = [];
for p = 1:numPlanes
    RAAN = RAANs(p); % Ogni piano avrà un RAAN diverso
    for s = 1:numSatsPerPlane
        ta = trueAnomalies(s); % Anomalia vera per i satelliti
        sat = satellite(sc, a, e, i, RAAN, argPeriapsis, ta); % Aggiungi satelliti con inclinazione di 90°
        satArray = [satArray; sat];
        access(sat, gs); % Analisi di accesso al Ground User
    end
end

% Visualizzazione
play(sc);
