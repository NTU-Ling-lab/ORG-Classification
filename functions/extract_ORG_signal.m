function extract_ORG_signal(RegisteredDataFolder, Params)
% This function will extract phase-based ORG signals from registered
% complex-valued OCT data

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
pre_stimulus_frames = Params.pre_stimulus_frames;   % the number of pre-stimulus frames

AnalysisDataFolder = [RegisteredDataFolder, '/analysis/'];
FigureFolder = [RegisteredDataFolder, '/figures/'];

if ~exist(AnalysisDataFolder) || ~exist(FigureFolder)
    mkdir(AnalysisDataFolder);
    mkdir(FigureFolder);
end

%% Load data
disp('Loading registered data...');
load([RegisteredDataFolder, 'I_reg.mat']);
load([RegisteredDataFolder, 'seg_mask_v2.mat']);

% the variable 'label' in seg_mask_v2.mat stores segmentation results
%        label == 1: ELM
%        label == 2: IS/OS
%        label == 3; OS & RPE
%        label == 4: BrM
%        label == 0; others

% show the OCT structural image overlaid with segmentation results 
figure;
imshowpair(mean_bscan_log, label);
title('OCT structural image overlaid with segmentation');
drawnow;
pause(1);

%% Extract signals from different paired layers
disp('Extracting ORG signals...');
for refer_label_index = [2,4] % reference layer: IS/OS and BrM
    
    % get indices of all pixels in the reference layer
    refer_mask = (label == refer_label_index);
    refer_mask_indice = find(refer_mask(:));
    refer_mask_pixel_num = length(refer_mask_indice);
    refer_row = zeros(refer_mask_pixel_num,1);
    refer_col = zeros(refer_mask_pixel_num,1);
    for i = 1:refer_mask_pixel_num
        [refer_row(i), refer_col(i)] = ind2sub(size(label),refer_mask_indice(i));
    end

    switch refer_label_index
        case 2 % take the IS/OS as the refernce layer
            target_label_indice = [1,3]; % signals from the ELM layer and the OS & RPE layer will be extracted
        case 4 % take the BrM as the reference layer
            target_label_indice = [1,2]; % signals from the ELM layer and the BrM layer will be extracted
    end
    
    for target_label_index = target_label_indice

        % get indices of all pixels in the target layer
        target_mask = (label==target_label_index);
        target_mask_indice = find(target_mask(:));
        target_mask_pixel_num = length(target_mask_indice);
        target_row = zeros(target_mask_pixel_num,1);
        target_col = zeros(target_mask_pixel_num,1);
        for i = 1:target_mask_pixel_num
            [target_row(i), target_col(i)] = ind2sub(size(label),target_mask_indice(i));
        end
    
        % Extract ORG signals from the target layer with respect to the reference layer
        paired_phase = zeros(target_mask_pixel_num, size(I,3));
        for col = 1:size(I, 2)
            ROI_index_target = find(target_col == col); % find pixels in the target layer located at the colume 'col'
            ROI_index_refer = find((refer_col >= col - 2) & (refer_col <= col + 2)); % find corresponding referecne region
            if ~isempty(ROI_index_target) && ~isempty(ROI_index_refer)

                % get complex-valued signals from the reference region
                I_refer = zeros(length(ROI_index_refer), size(I,3));
                for i = 1:length(ROI_index_refer)
                    I_refer(i,:) = squeeze(I(refer_row(ROI_index_refer(i)), refer_col(ROI_index_refer(i)), :)).';
                end
        
                for i = ROI_index_target(1):ROI_index_target(end)
                    I_target = squeeze(I(target_row(i),target_col(i),:)).'; % get complex-valued signals from target pixels
                    I_target = repmat(I_target,[length(ROI_index_refer),1]);
                    I_target_refer = I_target .* conj(I_refer); % canceling phase drift by self-referencing
                    I_target_refer = I_target_refer .* conj(repmat(mean(I_target_refer(:,1:pre_stimulus_frames), 2), [1,size(I,3)])); % reference to pre-stimulus frames
                    paired_phase(i,:) = mean(I_target_refer,1); % average across the reference region
                end
            end
        end

        % Assign extraced phase traces to different variables named
        % according to the target layer and the reference layer.
        % Nomenclature rule: "a_b_signal" represents that the target
        % layer is "a" and the reference layer is "b"
        switch refer_label_index
            case 2
                switch target_label_index
                    case 1
                        ELM_ISOS_signal = paired_phase;
                    case 3
                        OSRPE_ISOS_signal = angle(paired_phase);
                        OSRPE_index = [target_row, target_col];
                end
            case 4
                switch target_label_index
                    case 1
                        ELM_BrM_signal = paired_phase;
                    case 2
                        ISOS_BrM_signal = paired_phase;
                end
        end
    end
end

% average across different pixels and then take the phase component
mean_ELM_BrM_raw = angle(mean(ELM_BrM_signal, 1));
mean_ISOS_BrM_raw = angle(mean(ISOS_BrM_signal, 1));
mean_ELM_ISOS_raw = angle(mean(ELM_ISOS_signal, 1));

% save data to the analysis folder
save([AnalysisDataFolder,'ORG_signals.mat'], ...
'OSRPE_ISOS_signal', 'OSRPE_index', ...
'mean_ELM_BrM_raw', 'mean_ISOS_BrM_raw', 'mean_ELM_ISOS_raw', '-v7.3');
end