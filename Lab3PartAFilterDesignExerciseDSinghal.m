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
