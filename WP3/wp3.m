clear; clc; close all;

% Create an ITU-R P.681-11 channel
chan = p681LMSChannel;


% Environment type
chan.Environment = "Urban";
% Carrier frequency (in Hz)
chan.CarrierFrequency = 3.8e9;
% Elevation angle with respect to ground plane (in degrees)
chan.ElevationAngle = 45;
% Speed of movement of ground terminal (in m/s)
chan.MobileSpeed = 3;
% Direction of movement of ground terminal (in degrees)
chan.AzimuthOrientation = 0;

% Sampling rate (in Hz)
chan.SampleRate = 400;

chan.InitialState = "Good";

chan.FadingTechnique = "Filtered Gaussian noise"; 

seed = 72;
chan.RandomStream = "mt19937ar with seed";
chan.Seed = seed;


% Set random number generator with seed
rng(seed);
% Channel duration (in seconds)
chanDur = 500;
% Random input waveform
numSamples = floor(chan.SampleRate*chanDur)+1;
in = complex(randn(numSamples,1),randn(numSamples,1));
% Pass the input signal through channel
[fadedWave,channelCoefficients,sampleTimes,stateSeries] = step(chan,in);

distanceVector = chan.MobileSpeed * sampleTimes;

figure(1)
plot(distanceVector, 20*log10(abs(in)), distanceVector, 20*log10(abs(fadedWave)))
title(['Power Profile of Waveform over Distance (Speed = ' num2str(chan.MobileSpeed) ' m/s)'])
legend('Input Waveform', 'Faded Waveform')
xlabel('Distance (in meters)')
ylabel('Power (in dB)')


figure(2)
plot(distanceVector, 20*log10(abs(channelCoefficients)))
title('Channel Gain over Distance')
xlabel('Distance (in meters)')
ylabel('Path Gain (in dB)')


figure(3)
plot(distanceVector, stateSeries)
title('State Series of Channel over Distance')
axis([0 distanceVector(end) -0.5 1.5])
xlabel('Distance (in meters)')
ylabel('State')


% generate 5 different objects, quantize