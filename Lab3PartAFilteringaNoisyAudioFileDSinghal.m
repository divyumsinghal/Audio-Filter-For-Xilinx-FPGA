close all;
clear;
clc;

load('BandstopFilterHd_Standard.mat',    'Hd_Standard');              % loads the filter object back into workspace
load('BandstopFilterHd_Standard_8.mat',  'Hd_Standard_8');            % loads the filter object back into workspace
load('BandstopFilterHd_Standard_16.mat', 'Hd_Standard_16');           % loads the filter object back into workspace
load('BandstopFilterHd_Standard_32.mat', 'Hd_Standard_32');           % loads the filter object back into workspace
load('BandstopFilterHd_Standard_64.mat', 'Hd_Standard_64');           % loads the filter object back into workspace

[x, Fs] = audioread('speech_38.wav');   % Read audio file

%% Plot to find out which frequency to remove from signal 
nfft    = 2^10;
X       = fft(x, nfft);
fstep   = Fs / nfft;
fvec    = fstep * (0:nfft/2 - 1);
fresp   = 2 * abs(X(1:nfft/2));

figure;
plot(fvec, fresp);
title('Single-Sided Amplitude Spectrum of x(t)');
xlabel('Frequency (Hz)');
ylabel('|X(f)|');
grid on;

%% Go over the vector to find the limits for the BandStop Filter
max_val = 0;
max_2 = 0;
max_fs = 0;

for i = 2:length(fresp)-1
    if fresp(i) > fresp(i-1) && fresp(i) > fresp(i+1)
        if fresp(i) > max_val
            max_2   = max_val;
            max_val = fresp(i);
            max_fs  = i;
        elseif fresp(i) > max_2
            max_2   = fresp(i);
        end
    end
end

fprintf('Sampled at %.2f Hz\n', Fs);
fprintf('Peak at %.2f Hz, amplitude = %.2f\n', fvec(max_fs), fresp(max_fs));
fprintf('Second largest local max amplitude = %.2f\n', max_2);

% thresholds
th_stop = max_2;        % The stop is at the max_Signal / 2 spot (not the max_introduced / 2)
th_pass = max_2 / 3;

% Search left (towards lower freqs)
Fstop1 = max_fs;
for i = max_fs:-1:1
    if fresp(i) < th_stop
        Fstop1 = i;
        break;
    end
end

Fpass1 = Fstop1;
for i = Fstop1:-1:1
    if fresp(i) < th_pass
        Fpass1 = i;
        break;
    end
end

% Search right (towards higher freqs)
Fstop2 = max_fs;
for i = max_fs:length(fresp)
    if fresp(i) < th_stop
        Fstop2 = i;
        break;
    end
end

Fpass2 = Fstop2;
for i = Fstop2:length(fresp)
    if fresp(i) < th_pass
        Fpass2 = i;
        break;
    end
end

% Convert indices to frequencies
f_pass1 = fvec(Fpass1);
f_stop1 = fvec(Fstop1);
f_stop2 = fvec(Fstop2);
f_pass2 = fvec(Fpass2);

fprintf('\nDetected band edges:\n');
fprintf('Fpass1 = %.2f Hz\n', f_pass1);
fprintf('Fstop1 = %.2f Hz\n', f_stop1);
fprintf('Fstop2 = %.2f Hz\n', f_stop2);
fprintf('Fpass2 = %.2f Hz\n', f_pass2);

hold on;
xline(f_pass1,      'g--',  'Fpass1');
xline(f_stop1,      'r--',  'Fstop1');
xline(fvec(max_fs), 'k-',   'Peak');
xline(f_stop2,      'r--',  'Fstop2');
xline(f_pass2,      'g--',  'Fpass2');
legend('Spectrum','Fpass1','Fstop1','Peak','Fstop2','Fpass2');


%% Filter the noisy audio
y       = filter(Hd_Standard, x);       % Apply FIR filter as desgined in fD
y_8     = filter(Hd_Standard_8, x);     % Apply FIR filter as desgined in fD
y_8     = double(y_8);                  % Convert
y_16    = filter(Hd_Standard_16, x);    % Apply FIR filter as desgined in fD
y_16    = double(y_16);                 % Convert
y_32    = filter(Hd_Standard_32, x);    % Apply FIR filter as desgined in fD
y_32    = double(y_32);                 % Convert
y_64    = Hd_Standard_64(x);    % Apply FIR filter as desgined in fD
y_64    = double(y_64);                 % Convert


% Analyze filtered signal
Y       = fft(y, nfft);
Y_8     = fft(y_8, nfft);
Y_16    = fft(y_16, nfft);
Y_32    = fft(y_32, nfft);
Y_64    = fft(y_64, nfft);

fresp_filtered      = 2*abs(Y(1:nfft/2));
fresp_filtered_8    = 2*abs(Y_8(1:nfft/2));
fresp_filtered_16   = 2*abs(Y_16(1:nfft/2));
fresp_filtered_32   = 2*abs(Y_32(1:nfft/2));
fresp_filtered_64   = 2*abs(Y_64(1:nfft/2));

figure;
plot(fvec, fresp, 'r-', 'DisplayName', 'Original', 'LineWidth', 1);
hold on;

plot(fvec, fresp_filtered,      'b-',   'DisplayName', 'Filtered',      'LineWidth', 1);
plot(fvec, fresp_filtered_8,    'g--',  'DisplayName', 'Filtered_8',    'LineWidth', 1);
plot(fvec, fresp_filtered_16,   'm--',  'DisplayName', 'Filtered_16',   'LineWidth', 1);
plot(fvec, fresp_filtered_32,   'c--',  'DisplayName', 'Filtered_32',   'LineWidth', 1);
plot(fvec, fresp_filtered_64,   'c-*',  'DisplayName', 'Filtered_64',   'LineWidth', 1);

title('Frequency Spectrum Before and After Filtering');
xlabel('Frequency (Hz)');
ylabel('|X(f)|');
legend;
grid on;

%% Save and listen to the filtered audio
audiowrite('speech_filtered.wav',       y,      Fs);
audiowrite('speech_filtered_8.wav',     y_8,    Fs);
audiowrite('speech_filtered_16.wav',    y_16,   Fs);
audiowrite('speech_filtered_32.wav',    y_32,   Fs);
audiowrite('speech_filtered_64.wav',    y_64,   Fs);

% Listen to both
pause_dur = length(x)/Fs + 1;

sound(x, Fs);       % Original
pause(pause_dur);   % Pause

sound(y, Fs);       % Filtered
pause(pause_dur);   % Pause

sound(y_8, Fs);     % Filtered
pause(pause_dur);   % Pause

sound(y_16, Fs);    % Filtered
pause(pause_dur);   % Pause

sound(y_32, Fs);    % Filtered

%% Analyse the Filters

b = Hd_Standard.Numerator;
Order = length(b) - 1;
disp(['Filter Order: ', num2str(Order)]);
fvtool(b, 1, 'Fs', Fs)

b = Hd_Standard_8.Numerator;
Order = length(b) - 1;
disp(['Filter Order: ', num2str(Order)]);
fvtool(b, 1, 'Fs', Fs)

b = Hd_Standard_16.Numerator;
Order = length(b) - 1;
disp(['Filter Order: ', num2str(Order)]);
fvtool(b, 1, 'Fs', Fs)

b = Hd_Standard_32.Numerator;
Order = length(b) - 1;
disp(['Filter Order: ', num2str(Order)]);
fvtool(b, 1, 'Fs', Fs)

b = Hd_Standard_64.Numerator;
Order = length(b) - 1;
disp(['Filter Order: ', num2str(Order)]);
fvtool(b, 1, 'Fs', Fs)