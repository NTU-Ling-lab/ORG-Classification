function classification_prolonged_recording(RegisteredDataFolder, Params)
% This function is used for processing phase traces extracted from
% prolonged recordings at a sampling rate of 25 B-scans per second.

% The script will extract Type-I and Type-II signals using a pre-trained
% support vector machine (SVM) model and further extract multilayerd ORG
% signals.

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


%% set up parameters and paths
phase_to_OPL = 850 / (4*pi);
pre_stimulus_frames = Params.pre_stimulus_frames; % the number of pre-stimulus frames
Fs = Params.Fs; % frame rate

AnalysisDataFolder = [RegisteredDataFolder, '/analysis/'];
FigureFolder = [RegisteredDataFolder, '/figures/'];
%% load data
load([AnalysisDataFolder,'/ORG_signals.mat']);
load([RegisteredDataFolder,'/seg_mask_v2.mat']);

%% project phase traces into the spatiotemporal feature space
% extract distances to the upper boundary of the BrM
dis_BrM = zeros(size(OSRPE_index,1),1);
for i = 1:length(dis_BrM)
    line_num = OSRPE_index(i,2);
    BrM_loc = find(label(:,line_num)==3, 1, 'last');
    dis_BrM(i) = BrM_loc - OSRPE_index(i,1);
end

% bandstop filters to remove residual artifacts induced by heartbeat and breathing
OSRPE_ISOS_signal_bandstop = Wrap_bandstop_filters(OSRPE_ISOS_signal); 

std_baseline = std(OSRPE_ISOS_signal_bandstop(:, 1: pre_stimulus_frames), [], 2); % standard deviation in pre-stimulus period
ROI_index = find(std_baseline < 0.3); % index of phase traces that are stable in the pre-stimulus period

% retain stable phase traces
OSRPE_ISOS_signal = OSRPE_ISOS_signal(ROI_index,:);
OSRPE_ISOS_signal_bandstop = OSRPE_ISOS_signal_bandstop(ROI_index,:);
dis_BrM = dis_BrM(ROI_index,:);

% interpolate to 200 Hz frame rate
downsample_time = 1:8:993;
interpolated_phase = zeros(size(OSRPE_ISOS_signal_bandstop,1), 1000);
for i = 1:size(interpolated_phase, 1)
    interpolated_phase(i,:)= interp1(downsample_time, ...
    OSRPE_ISOS_signal_bandstop(i,101:225), 1:1000, 'linear');
end

% smoothing intepolated phase traces using a gaussian filter
interpolated_phase = smoothdata(interpolated_phase, 2, "gaussian", 140);
cropped_phase = interpolated_phase(:,51:930); % crop both ends exhibiting high variance due to the edge effects
norm_phase = (cropped_phase - repmat(mean(cropped_phase, 2), [1, size(cropped_phase,2)])) ./ ...
             repmat(std(cropped_phase,[],2),[1,size(cropped_phase,2)]); % normalize individual phase trace     

load('./SVM/PCA_coeff.mat');
score = norm_phase * coeff; % project onto the principal component space
X = [score(:,1:2), dis_BrM]; % concaract temporal and spatial features
norm_X = X;

% normalize each feature
Transform_matrix(3,2) = 28;
for i = 1:size(norm_X,2)
    norm_X(:,i) = (norm_X(:,i) - Transform_matrix(i,1)) / Transform_matrix(i,2);
end

%% classification using the pre-trained SVM
load('./SVM/SVM_Mdl.mat');
recflag =  predict(Mdl, norm_X);

% get index of the Type-I (OS) and Type-II signal
OS_index = (recflag == 2);
RPE_index = (recflag == 1);

OS_ISOS_signal = OSRPE_ISOS_signal(OS_index, :);
RPE_ISOS_signal = OSRPE_ISOS_signal(RPE_index, :);

%% extract multilayed dynamics
samp_time = ([1: size(OSRPE_ISOS_signal, 2)] - pre_stimulus_frames) / Fs;
mean_OS_ISOS_raw = (angle(mean(exp(1i * OS_ISOS_signal), 1))); % average across different pixels and then get the phase component
mean_RPE_ISOS_raw = unwrap((angle(mean(exp(1i * RPE_ISOS_signal), 1))), [], 2);

% bandstop filtering
mean_ELM_BrM_bandstop = Wrap_bandstop_filters(unwrap(mean_ELM_BrM_raw, [], 2));
mean_ISOS_BrM_bandstop = Wrap_bandstop_filters(mean_ISOS_BrM_raw);
mean_OS_ISOS_bandstop = Wrap_bandstop_filters(mean_OS_ISOS_raw);
mean_RPE_ISOS_bandstop = Wrap_bandstop_filters(mean_RPE_ISOS_raw);

% get multilayerd ORG siganls
mean_ISOS_ELM_bandstop = unwrap(mean_ISOS_BrM_bandstop - mean_ELM_BrM_bandstop, [], 2);
mean_RPE_ELM_bandstop = unwrap(mean_RPE_ISOS_bandstop + mean_ISOS_ELM_bandstop, [], 2);
mean_BrM_RPE_bandstop = (-(mean_ISOS_BrM_bandstop + mean_RPE_ISOS_bandstop));

% plot multilayered ORG signals
figure; hold on;
plot(samp_time, mean_RPE_ELM_bandstop * phase_to_OPL, 'Color', [0.6,0.6,0.6]);
plot(samp_time, mean_BrM_RPE_bandstop * phase_to_OPL, 'Color', [0.51,0.65,0.98]);
plot(samp_time, mean_ISOS_ELM_bandstop * phase_to_OPL, 'Color', [0.35,0.75,0.35]);
plot(samp_time, mean_OS_ISOS_bandstop * phase_to_OPL, 'Color', [0.75,0.3,0.32]);
xlim([-2.5, 20]);
ylim([-50, 280]);
xlabel('Time (s)');
ylabel('OPL change (nm)');
legend('SRS', 'RPE', 'IS', 'OS', 'Location', 'northwest');
hold off;
drawnow;
pause(1);

saveas(gcf, [FigureFolder,'/extracted_signals.tif']);

save([AnalysisDataFolder, '/Clustered_signals.mat'], ...
     'samp_time', 'mean_OS_ISOS_bandstop', ...
     'mean_RPE_ELM_bandstop', 'mean_ISOS_ELM_bandstop', 'mean_BrM_RPE_bandstop');
end


function filtered_traces = Wrap_bandstop_filters(raw_traces)
% This function use bandstop filters to remove residual artifacts
% induced by heartbeat and breathing

% set up bandpass filters
[b1, a1] = butter(6, [0.23, 0.245], 'stop');
[b2, a2] = butter(3, [0.465, 0.49], 'stop');
[b3, a3] = butter(3, [0.71, 0.725], 'stop');

raw_traces = raw_traces'; % transpose for filtering along the correct dimention

filtered_traces = filtfilt(b1, a1, raw_traces);
filtered_traces = filtfilt(b2, a2, filtered_traces);
filtered_traces = filtfilt(b3, a3, filtered_traces);
filtered_traces = filtered_traces'; % transpose to the original shape
end