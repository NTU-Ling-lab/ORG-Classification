% Copyright (C) 2023  Huakun Li, Tong Ling
%
%   This script provides an example for classifying ORG signals acquired in
%   prolonged recordings at a sampling rate of 25 B-scans per second using
%   a pretrained SVM model.
%
%   The script consists of two functions,
%   1) extract_ORG_signal: extracts phase-based ORG signals from registered
%   complex-valued OCT signals.
%
%   2) classification_prolonged_recording: extract Type-I and Type-II
%   signals using a pre-trained support vector machine (SVM) model and
%   further extract multilayerd ORG signals.
%
%
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


% initialization
clear; clc;
close all;
addpath('./functions/');

RegisteredDataFileName = 'example2-data';

% download example dataset
if ~(exist(RegisteredDataFileName,'dir') || exist([RegisteredDataFileName '.zip'],'file'))
    disp('Downloading example dataset (data size is 4.8 GB, may take ~5 minutes to complete)...');
    websave([RegisteredDataFileName '.zip'],'https://f2n7.c18.e2-1.dev/tonglinglab-share/ORG-Classification/example2-data.zip');
end

% unzip example dataset
if ~exist(RegisteredDataFileName,'dir')
    disp('Unzipping...');
    unzip([RegisteredDataFileName '.zip']);
end

% setup path
RegisteredDataFolder = ['./' RegisteredDataFileName '/']; % path to the registered data folder

% setup parameters
Params.pre_stimulus_frames = 125; % the number of pre-stimulus frames
Params.Fs = 25; % frame rate

% extract phase signals from registered OCT data
extract_ORG_signal(RegisteredDataFolder, Params); 

% extract multilayered ORG signals
classification_prolonged_recording(RegisteredDataFolder, Params);
