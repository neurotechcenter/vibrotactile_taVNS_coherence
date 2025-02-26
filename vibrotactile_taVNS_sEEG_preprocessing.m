% Kara Donovan & Phil Demarest:  5/20/2022
% preprocessing script for vibrotactile data
             
clear;
clc;

%% specify subject
subjNum = 1;            

% specify re-referencing type
ref = 'car';

%% load data
load(strcat('D:\sEEG_vibrotactileProjects\vibrotactile_taVNS_data\Subject00',num2str(subjNum),'_data.mat'));
signal = rawData.signal;
sr = rawData.sr;
groups = rawData.groups;
stimulusCode = rawData.stimulusCode;
clear rawData

%% high-pass filter at 0.5 Hz
hp_cutoff = 0.5;
order = 4;
type = 'high';
[b0_hp, a0_hp] = butter(order, 2*hp_cutoff/sr, type);

hp_signal = filtfilt(b0_hp,a0_hp,double(signal'))'; % note: have to be in samples x chans for filtfilt function

%% measure line noise power and denote excessively noisy channels
line_noise = butter_bandpass_filtering(hp_signal', sr, [55 65], 4); % note: have to be in samples x chans for butter_bandpass_filtering
line_noise_power = double(abs(hilbert(line_noise)).^2);
clear line_noise

line_noise_power_ch = mean(line_noise_power);   % avg line noise power for each channel
m_power = median(line_noise_power_ch);          % median line noise power across channels
std_power = mad(line_noise_power_ch);

% select good channels
range = 3;
good_channels = find(line_noise_power_ch > (m_power-(range*std_power)) & line_noise_power_ch < (m_power+(range*std_power)));
num_chans = length(line_noise_power_ch);
bad_channels = setdiff(1:num_chans,good_channels);

% plot
figure();
stem(line_noise_power_ch); hold on
hline(m_power-(range*std_power),'g')
hline(m_power,'k')
hline(m_power+(range*std_power),'g')
vline(bad_channels,'r');   % highlight excessively noisy channels in red
title("Subject " + num2str(subjNum),'FontSize',16)
xlabel('Channel Number','FontSize',14)
ylabel('Line Noise Power','FontSize',14)

%% exclude bad channels (after visually confirming they are bad)
if subjNum == 1
    bad_channels = [66, 199];
elseif subjNum == 2
    bad_channels = 66;
elseif subjNum == 3
    bad_channels = [105, 113];
elseif subjNum == 4
    bad_channels = [64, 99, 100, 108];
elseif subjNum == 5
    bad_channels = 131;
elseif subjNum == 6
    bad_channels = [];
elseif subjNum == 7
    bad_channels = [1, 3, 34, 38, 47, 48, 61, 62];
end

hp_signal(bad_channels,:) = [];
groups(bad_channels) = [];

%% re-reference using CAR (excluding noisy channels)
if contains(ref,'car')
    final_car = mean(hp_signal);
    reref_data = hp_signal;   % chans x samples
    for ch = 1:length(groups)
        reref_data(ch,:) = hp_signal(ch,:) - final_car;
    end
elseif contains(ref,'none')
    reref_data = hp_signal;
end

%% notch filter at 60 Hz (and subsequent harmonics)
signal_preprocessed = double(multi_iirnotch_filtering(reref_data',sr,[60 120 180]))';  % multi_iirnotch_filtering requires samples x channels

%% Create Data Struct
preProcessedReRef.signal = signal_preprocessed;
preProcessedReRef.samplingRate = sr;
preProcessedReRef.channelsRemoved = bad_channels;
preProcessedReRef.stimulusCode = stimulus_code;
preProcessedReRef.groups = groups;
preProcessedReRef.refType = ref;

%% save preprocessed dat
save(['D:\sEEG_vibrotactileProjects\preprocessed_data\Subject00' num2str(subjNum) '_data.mat'],...
    'preProcessedReRef','-v7.3','-nocompression')
disp('Done saving!');

%% functions
function [ y ] = butter_bandpass_filtering( X, Fs, fBand, n)
% Author: Hohyun Cho, PhD, hohyun@wustl.edu
% X: time x channal
% Fs: sampling rate
% fBand: frequency band
% n: order of butter worth filtering
Wn = fBand;
Fn = Fs/2;
ftype = 'bandpass';
[b, a] = butter(n,Wn/Fn,ftype);
y = filtfilt(b,a,X);
end

function [filtered_signal] = multi_iirnotch_filtering(signal,srate,notch_freq)
% Author: Hohyun Cho, PhD, hohyun@wustl.edu
% define the line noise frequency and bandwidth
peak.fcenter = 60;
peak.bw      = 0.001;

% calculate the IIR-peak filter coefficients in a,b format 
peak.wo = peak.fcenter/(srate/2);  
peak.bw = peak.bw;
[peak.b,peak.a] = iirpeak(peak.wo,peak.bw);  

% define the harmonics of line noise frequency
param.filter.notch.fcenter = notch_freq;

param.filter.notch.bw      = ones(1,length(param.filter.notch.fcenter)).*0.001;

% calculate the IIR-peak filter coefficients in a,b format 
for idx = 1:length(param.filter.notch.fcenter),
    notch{idx}.wo = param.filter.notch.fcenter(idx)/(srate/2);  
    notch{idx}.bw = param.filter.notch.bw(idx);
    [notch{idx}.b,notch{idx}.a] = iirnotch(notch{idx}.wo,notch{idx}.bw);  
end

fprintf(1, '> Notch filtering signal \n');
fprintf(1,'[');
% for each channel
for idx_channel=1:size(signal,2),
    
    % get the signal for this channel
    signal_preliminary = double(signal(:,idx_channel));
    
    % remove all harmonics of line-noise
    for idx = 1:length(param.filter.notch.fcenter), 
        signal_preliminary = filtfilt(notch{idx}.b,notch{idx}.a,signal_preliminary); 
    end 
    
    % return the signal
    filtered_signal(:,idx_channel) = single(signal_preliminary);
    
    fprintf(1,'.');
end
fprintf(1,'] done\n');


end
