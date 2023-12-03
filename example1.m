% Copyright (C) 2023  Huakun Li, Tong Ling
%
%   This script provides an example for classifying ORG signals acquired in 
%   5 s recordings at a sampling rate of 200 B-scans per second using a
%   pretrained SVM model.
%
%   The script consists of three functions,
%   1) extract_ORG_signal: extracts phase-based ORG signals from registered
%   complex-valued OCT signals.
%
%   2) preprocessing: uses bandstop filters and a low-pass filter to remove
%   unwanted frequencies from raw phase traces.
%
%   3) SVM_classification: extracts Type-I and Type-II signals using a
%   pre-trained support vector machine (SVM) model.
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

RegisteredDataFileName = 'example1-data';

% download example dataset
if ~(exist(RegisteredDataFileName,'dir') || exist([RegisteredDataFileName '.zip'],'file'))
    disp('Downloading example dataset (data size is ~7 GB, may take 5-10 minutes to complete)...');
    websave([RegisteredDataFileName '.zip'],'https://f2n7.c18.e2-1.dev/tonglinglab-share/ORG-Classification/example1-data.zip');
end

% unzip example dataset
if ~exist(RegisteredDataFileName,'dir')
    disp('Unzipping...');
    unzip([RegisteredDataFileName '.zip']);
end

% setup path
RegisteredDataFolder = ['./' RegisteredDataFileName '/']; % path to the registered data folder

% setup parameters
Params.pre_stimulus_frames = 197; % the number of pre-stimulus frames
Params.Fs = 200; % frame rate

% extract phase signals from registered OCT data
extract_ORG_signal(RegisteredDataFolder, Params); 

% use temporal filters to process raw phase traces
preprocessing(RegisteredDataFolder, Params);

% classification using the pre-trained support vector machine (SVM)
SVM_classification(RegisteredDataFolder, Params);
