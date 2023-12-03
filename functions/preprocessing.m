function [] = preprocessing(RegisteredDataFolder, Params)
% This function will use bandstop filters and a low-pass filter to remove
% unwanted frequencies from raw phase traces.

% INPUTS
%    RegisteredDataFolder       path to the registered data folder
%    Params                     a structure that stores all the parameters

% Copyright (C) 2023  Huakun Li
%   This library is free software; you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation; either
%   version 2.1 of the License, or (at your option) any later version.
% 
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
% 
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301 USA

%% Set up parameters and paths
pre_stimulus_frames = Params.pre_stimulus_frames; % the number of pre-stimulus frames
Fs = Params.Fs; % frame rate

AnalysisDataFolder = [RegisteredDataFolder, '/analysis/'];

%% load data
disp('Loading extracted ORG signals...');
load([AnalysisDataFolder,'/ORG_signals.mat']); % load raw phase traces

%% process signals extracted from the OS & RPE band
[b1,a1] = butter(3,[2.5 5.5] / (Fs/2), 'stop'); % ~3--5Hz, heartbeat
[b2,a2] = butter(3,[5 11] / (Fs/2), 'stop'); % ~6--10Hz, secondary harmonics of heartbeat
[b3,a3] = butter(6, 10 / (Fs/2), 'low'); % ~10Hz, low-pass filter 

OSRPE_ISOS_signal = OSRPE_ISOS_signal.'; % transpose for filtering along the correct dimention
bandstop_filtered_phase = filtfilt(b1,a1,OSRPE_ISOS_signal);
bandstop_filtered_phase = filtfilt(b2,a2,bandstop_filtered_phase);
lowpass_filtered_phase = filtfilt(b3,a3,bandstop_filtered_phase);

% transpose to the original shape
OSRPE_ISOS_signal = OSRPE_ISOS_signal.';
bandstop_filtered_phase = bandstop_filtered_phase.';
lowpass_filtered_phase = lowpass_filtered_phase.';

std_baseline = std(lowpass_filtered_phase(:, 51 : pre_stimulus_frames), [], 2); % standard deviation in pre-stimulus period
ROI_index = find(std_baseline < 0.06); % index of phase traces that are stable in the pre-stimulus period

% retain stable phase traces
OSRPE_index = OSRPE_index(ROI_index, :); 
OSRPE_ISOS_raw = OSRPE_ISOS_signal(ROI_index, :);
OSRPE_ISOS_bandstop = bandstop_filtered_phase(ROI_index, :);
OSRPE_ISOS_lowpass = lowpass_filtered_phase(ROI_index, :);
%% process phase traces extracted from the ISOS layer
mean_ISOS_ELM_raw = -mean_ELM_ISOS_raw; % revert the target layer and the reference layer

% find the breathing frequency
range_breath = round([0.8,2.4] / (Fs/size(mean_ISOS_ELM_raw,2))) + 1; % 0.8~2.2Hz, breath
range_breath = range_breath(1) : range_breath(2);
fft_ISOS_ELM = abs(fft(mean_ISOS_ELM_raw, [], 2));
peak_index_breath = find((fft_ISOS_ELM(range_breath) > 1.1 * fft_ISOS_ELM(range_breath - 1)) & ...
                 (fft_ISOS_ELM(range_breath) > 1.1 * fft_ISOS_ELM(range_breath + 1))) + range_breath(1) - 1;

% remove residul artifacts induced by breathing
mean_ISOS_ELM_bandstop = mean_ISOS_ELM_raw;
for filter_i = 1:length(peak_index_breath)
    peak_i = peak_index_breath(filter_i);
    [b1,a1] = butter(3,[(peak_i-1)/(size(OSRPE_ISOS_signal,2)/2) - 0.2/(Fs/2), (peak_i-1)/(size(OSRPE_ISOS_signal,2)/2) + 0.2/(Fs/2)], 'stop'); % breath
    mean_ISOS_ELM_bandstop = filtfilt(b1,a1,mean_ISOS_ELM_bandstop);
end

[b1,a1] = butter(3,[2.2 5.5] / (Fs/2), 'stop'); % ~3--5Hz, heartbeat
mean_ISOS_ELM_bandstop = filtfilt(b1, a1, mean_ISOS_ELM_bandstop);
mean_ISOS_ELM_bandstop = filtfilt(b2, a2, mean_ISOS_ELM_bandstop);
mean_ISOS_ELM_lowpass = filtfilt(b3, a3, mean_ISOS_ELM_bandstop);

% save filtered phase traces
disp('Saving filtered ORG signals...');
save([AnalysisDataFolder,'/Filtered_ORG_signals.mat'], ...
    'OSRPE_ISOS_raw', 'OSRPE_ISOS_bandstop', 'OSRPE_ISOS_lowpass', 'OSRPE_index', ...
    'mean_ISOS_ELM_raw', 'mean_ISOS_ELM_bandstop', 'mean_ISOS_ELM_lowpass', ...
    '-v7.3');

end
