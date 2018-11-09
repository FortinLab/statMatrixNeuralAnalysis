function [Events_TS,OnTS,OffTS,TS_Raw] = EventsTS(fileDir,processorNum,recStart,rigRoom,varargin)
%% Updated 7/9/2018: Added conditional for rig used (either 212 or 206)
% This function extracts event presentation information from the ADC input
% files found in the OpenEphys data session folders 
% Example code for FPGA node (aka processor = 100):
% [Events_TS,OnTS,OffTS,TS_Raw] = EventsTS(fileDir,100,recStart,206,'downsample')

% Detect and determine # of ADC input files in set folder directory
origDir = cd; % saves current file directory
cd(fileDir) % changes to target directory
Filename = sprintf('%d_ADC',processorNum);
numfiles = length(dir([Filename,'*']));
cd(origDir) % reverts to original directory

% Adds OpenEphysTools to the current pathway since it contains supporting
% files
addpath(genpath('OpenEhysTools'))

% Cycles and extracts necessary information from detected ADC input files
for filenum = 1:numfiles
    % Extracts respective ADC input into Matlab
    myfilename = sprintf('%d_ADC%d.continuous', processorNum,filenum);
    [contData,~,info] = load_open_ephys_data_faster([fileDir '\' myfilename]);
    sampleRate = info.header.sampleRate;
    
    % Sets timestamp vector to 0 based on the starting recording
    % information
    TS_Raw = recStart:1/sampleRate:recStart+((length(contData)/sampleRate)-(1/sampleRate));
    
    % Checks sample rate to see if downsampling is necessary
    if sum(cell2mat(cellfun(@(a)(strcmp(a,'Downsample') | strcmp(a, 'downsample')), varargin, 'uniformoutput', 0)))>=1
        if sampleRate > 1000
            n_samples = sampleRate/1000;
            contData = downsample(contData,n_samples);
            TS_Raw = downsample(TS_Raw,n_samples);
        end
    end
    
    % Detects location of peaks in the digital signal aka event activation
    % If the data is recorded in the 206 behavior rig the system works by items activating on low voltage, not high.
    if rigRoom == 206
        [~, offEvents] = findpeaks([0; diff(contData)],1:length(contData), 'MinPeakProminence', 0.25);
        [~, onEvents] = findpeaks([0; diff(contData)*-1], 1:length(contData), 'MinPeakProminence', 0.25);
    end
    if rigRoom == 212
        [~, onEvents] = findpeaks([0; diff(contData)],1:length(contData), 'MinPeakProminence', 0.25);
        [~, offEvents] = findpeaks([0; diff(contData)*-1], 1:length(contData), 'MinPeakProminence', 0.25);
    end
    
    % Interpolates square pulse peak location to pinpoint timestamps
    On_field = strcat('OnPeakTS_',num2str(filenum));
    OnTS.(On_field)= interp1(TS_Raw,onEvents);
    
    Off_field = strcat('OffPeakTS_',num2str(filenum));
    OffTS.(Off_field) = interp1(TS_Raw,offEvents);
    
    % Initializes final output structure
    Events_TS.(On_field) = zeros(size(TS_Raw));
    
    % Combines timestamp indices with the raw file timestamp to create
    % index matrices
    [~,TS_index] = intersect(TS_Raw,OnTS.(On_field));
    
    % Marks timestamp locations on the final output structure for each ADC
    % file
    for i = 1:length(TS_index)
        Events_TS.(On_field)(TS_index) = 1;
    end
    
end

end




