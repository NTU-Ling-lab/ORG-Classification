function [] = SVM_classification(RegisteredDataFolder, Params)
% This function will extract Type-I and Type-II signals using a pre-trained SVM.

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
load([AnalysisDataFolder,'/Filtered_ORG_signals.mat']);
load([RegisteredDataFolder,'/seg_mask_v2.mat']);

%% project phase traces into the spatiotemporal feature space
% extract distances to the upper boundary of the BrM
dis_BrM = zeros(size(OSRPE_index,1),1);
for i = 1:length(dis_BrM)
    line_num = OSRPE_index(i,2);
    BrM_loc = find(label(:,line_num)==3, 1, 'last');
    dis_BrM(i) = BrM_loc - OSRPE_index(i,1);
end

cropped_OSRPE_ISOS_lowpass = OSRPE_ISOS_lowpass(:,51:930); % crop both ends exhibiting high variance due to the edge effects
norm_phase = (cropped_OSRPE_ISOS_lowpass - repmat(mean(cropped_OSRPE_ISOS_lowpass, 2), [1, size(cropped_OSRPE_ISOS_lowpass,2)])) ./ ...
             repmat(std(cropped_OSRPE_ISOS_lowpass,[],2),[1,size(cropped_OSRPE_ISOS_lowpass,2)]); % normalize individual phase trace

load('./SVM/PCA_coeff.mat');
score = norm_phase * coeff; % project onto the principal component space
X = [score(:,1:2), dis_BrM]; % concaract temporal and spatial features

% normalize each feature
norm_X = X;
for i = 1:size(norm_X,2)
    norm_X(:,i) = (norm_X(:,i) - Transform_matrix(i,1)) / Transform_matrix(i,2);
end

%% classification using the pre-trained SVM
load('./SVM/SVM_Mdl.mat');
recflag = predict(Mdl, norm_X);

% show different clusters in the spatiotemporal feature space
figure1 = figure;
axes1 = axes('Parent',figure1);
hold on;
scatter3(norm_X(recflag == 2,1), norm_X(recflag == 2,2), norm_X(recflag == 2,3), ...
    10, 'filled', 'MarkerFaceColor', [0.75,0.3,0.32], 'MarkerFaceAlpha', .6);
scatter3(norm_X(recflag == 1,1), norm_X(recflag == 1,2), norm_X(recflag == 1,3), ...
    10, 'filled', 'MarkerFaceColor', [0.3,0.45,0.7], 'MarkerFaceAlpha', .6);
scatter3(norm_X(recflag == 3,1), norm_X(recflag == 3,2), norm_X(recflag == 3,3), ...
    10, 'filled', 'MarkerFaceColor', [0.92,0.82,0.77], 'MarkerFaceAlpha', .6);
scatter3(norm_X(recflag == -1,1), norm_X(recflag == -1,2), norm_X(recflag == -1,3), ...
    10, 'filled', 'MarkerFaceColor', [0.65,0.65,0.65], 'MarkerFaceAlpha', .2);
zlim([0,1]);
view(axes1,[200 35]);
% legend('Type-I', 'Type-II');
hold off;
drawnow;
pause(1);
saveas(gcf, [FigureFolder,'/spatiotemporal space.tif']);

% get index of the Type-I (OS) and Type-II signal
OS_index = (recflag == 2);
RPE_index = (recflag == 1);

%% plot Type-I (OS) and Type-II signal and the SRS dynamics
samp_time = ([1: size(OSRPE_ISOS_lowpass, 2)] - pre_stimulus_frames) / Fs; % sampling time
OS_ISOS_bandstop = OSRPE_ISOS_bandstop(OS_index, :);
RPE_ISOS_bandstop= OSRPE_ISOS_bandstop(RPE_index, :);

mean_OS_ISOS_bandstop= angle(mean(exp(1i * OS_ISOS_bandstop), 1)); % average across different pixels and then get the phase component
mean_RPE_ISOS_bandstop = angle(mean(exp(1i * RPE_ISOS_bandstop), 1));

% incorporate the Type-II signal with the dynamics between ISOS and ELM to
% get the dynamics of SRS (RPE to ISOS)
mean_RPE_ELM_bandstop = mean_RPE_ISOS_bandstop + mean_ISOS_ELM_bandstop;

figure; hold on;
plot(samp_time, phase_to_OPL * mean_OS_ISOS_bandstop, 'Color', [0.75,0.3,0.32]);
plot(samp_time, phase_to_OPL * mean_RPE_ISOS_bandstop, 'Color', [0.3,0.45,0.7]);
plot(samp_time, phase_to_OPL * mean_RPE_ELM_bandstop, 'Color', [0.5,0.5,0.5]);
xlim([-1,4]);
xlabel('Time (s)');
ylabel('\DeltaOPL (nm)');
legend('Type-I (OS)', 'Type-II', 'SRS', 'Location', 'northwest');
drawnow;
pause(1);
saveas(gcf, [FigureFolder,'/extracted signals.tif']);

%% plot the spatial distribution of Type-I and Type-II signals in the structural image
col_OS = OSRPE_index(OS_index,2); row_OS = OSRPE_index(OS_index,1);
col_RPE = OSRPE_index(RPE_index,2); row_RPE = OSRPE_index(RPE_index,1);

figure;
imshow(mean_bscan_log,[35 70],'border','tight'); hold on;
scatter(col_OS, row_OS,'r*');
scatter(col_RPE, row_RPE,'b+');
hold off;
drawnow;
pause(1);
saveas(gcf, [FigureFolder,'/spatial distribution.tif']);

% save the extracted dynamics
save([AnalysisDataFolder,'/Clustered_signals.mat'], 'samp_time', 'mean_OS_ISOS_bandstop', 'mean_RPE_ELM_bandstop');
end
