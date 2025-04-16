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
a = 7171e3; % semi-major axis [m]
e = 0;      
i = 90;     
argPeriapsis = 0; 
numPlanes = 3;
numSatsPerPlane = 5;

% Anomalie vere per i satelliti
trueAnomalies = linspace(0, 360, numSatsPerPlane + 1);
trueAnomalies(end) = [];

% RAANs distribuiti uniformemente
RAANs = linspace(0, 90, numPlanes+1); 
RAANs(end) = [];

% Offset di fase per sfalsare i satelliti tra le orbite
phaseOffset = 360 / (numSatsPerPlane * numPlanes); 

% Aggiungi i satelliti
satArray = [];
for p = 1:numPlanes
    RAAN = RAANs(p);
    for s = 1:numSatsPerPlane
        ta = mod(trueAnomalies(s) + (p-1)*phaseOffset, 360);
        sat = satellite(sc, a, e, i, RAAN, argPeriapsis, ta);
        satArray = [satArray; sat];
        access(sat, gs);
    end
end

% Visualizzazione
play(sc);

