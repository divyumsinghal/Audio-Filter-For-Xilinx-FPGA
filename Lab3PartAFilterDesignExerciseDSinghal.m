close all;
clear;
clc;

load('Part1.mat', 'Part1');             % loads the filter object back into workspace

Fs = 20000;
b = Part1.Numerator;
Order = length(b) - 1;

% Display the filter order
disp(['Filter Order: ', num2str(Order)]);
fvtool(b, 1, 'Fs', Fs)

% Frequency response
freqz(Part1.Numerator);
title('Frequency Response of Part1 Filter');

% Impulse response
impz(Part1.Numerator);
title('Impulse Response of Part1 Filter');

% Open Filter Analyzer for Part1
filterAnalyzer(Part1);

% Get filter information
filterInfo = info(Part1);
disp(filterInfo);
