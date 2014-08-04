
function [] = powerFeatures_v3(data_file,bin_size,varargin)
% Created by Isaac Shamie, Jul. 2014

% Edited by Xi Jiang, Jul.30th, 2014

% data_file: .mat file directory; the file contains only a structure variable (e.g. output of ft_preprocessing)
% bin_size: duration (in seconds) of time bins in which features are calculated
% varargin: [chans, butFilterOrder, save_dir]
%   chans: index vector of desirable channels
%   butFilterOrder: Butterworth filter order for fieldtrip
%   save_dir: save directory

% some basic analysis of neural data
% such as calculating z-scores

% Extract multiple features from the dataset    

%% Set up prerequisites

data = load(data_file);

if nargin < 2
    chans = 1:length(data.label);
    butFilterOrder = 4;
    save_dir = '/home/ishamie/Desktop/projects/24_7/NY442/442_day4_bpFeatures.mat';
elseif nargin < 3
    chans = varargin{1};
    butFilterOrder = 4;
    save_dir = '/home/ishamie/Desktop/projects/24_7/NY442/442_day4_bpFeatures.mat';
elseif nargin < 4
    chans = varargin{1};
    butFilterOrder = varargin{2};
    save_dir = '/home/ishamie/Desktop/projects/24_7/NY442/442_day4_bpFeatures.mat';
else
    chans = varargin{1};
    butFilterOrder = varargin{2};
    save_dir = varargin{3};
end

values = data.trial{1}(chans,:);        % this is all of the good channels

totalTimeSamples = length(values(1,:)); % this is needed for later, when averaging over every 30 seconds
fsample = data.fsample;
% butFilterOrder = 4;                     % For Butterworth filter

clear data
disp('Hello')

%% HGP and Hilbert

highGammaRange = [180 80];
hgp = ft_preproc_bandpassfilter(values,fsample, highGammaRange, butFilterOrder, 'but');

hgp = hilbert(hgp);
HGP_envelope = abs(hgp);
HGP_phase = angle(hgp);


%HGP_bpHilbert = ft_preproc_bandpassfilter(HGP_envelope,fsample, [0.5,30], butFilterOrder, 'but'); % this is bp after hilbert to remove potential artifcacts
%zscore for HGP envelope
% [HGP_Zscore, HGP_mean, HGP_stDev] = zscore(HGP_envelope,0,2);
[HGP_Zscore, ~, ~] = zscore(HGP_envelope,0,2);

disp('Extracted HGP zscore')


save(save_dir, '*Zscore','*phase', '*envelope', '-v7.3');
clearvars -except fsample values HGP_Zscore totalTimeSamples butFilterOrder save_dir bin_size; 

%% theta and Hilbert

%same as hgp just different range
thetaRange = [7 4];
thetaWave = ft_preproc_bandpassfilter(values,fsample, thetaRange, butFilterOrder, 'but');
hil_after_theta = hilbert(thetaWave);
theta_envelope = abs(hil_after_hgp);
theta_phase = angle(hil_after_hgp);

%theta_bpHilbert = ft_preproc_bandpassfilter(theta_envelope,fsample, [0.5,30], butFilterOrder, 'but'); % this is bp after hilbert to remove potential artifcacts
[theta_Zscore, ~, ~] = zscore(theta_envelope,0,2);

save(save_dir, '*Zscore','*phase', '*envelope', '-v7.3', '-append');
clearvars -except fsample values HGP_Zscore theta_Zscore totalTimeSamples butFilterOrder save_dir bin_size; 

%% alpha and Hilbert
%same as hgp just different range
alphaRange = [12 7];
alphaWave = ft_preproc_bandpassfilter(values,fsample, alphaRange, butFilterOrder, 'but');
hil_after_alpha = hilbert(thetaWave);
alpha_envelope = abs(hil_after_alpha);
alpha_phase = angle(hil_after_alpha);

%alpha_bpHilbert = ft_preproc_bandpassfilter(alpha_envelope,fsample, [0.5,30], butFilterOrder, 'but'); % this is bp after hilbert to remove potential artifcacts
[alpha_Zscore, ~, ~] = zscore(alpha_envelope,0,2);

save(save_dir, '*Zscore','*phase', '*envelope', '-v7.3', '-append');
clearvars -except fsample values HGP_Zscore theta_Zscore alpha_Zscore totalTimeSamples save_dir bin_size; 

%% just raw LFP, no bp, no hilbert
[LFP_Zscore, ~, ~] = zscore(values,0,2);

save(save_dir, '*Zscore','*phase', '*envelope', '-v7.3', '-append');
clearvars -except fsample HGP_Zscore theta_Zscore alpha_Zscore LFP_Zscore totalTimeSamples save_dir bin_size; 

%% Score tally
% go through data every 30 seconds and get calculations on mean, variance,
% skewness, and kurtosis. create a struct with three fields and each is a
% vector of all the channels avg, var, and skewness over 30 second
% timeframes

step = fsample*bin_size; % samples/sec * 30 seconds
preassign_length = totalTimeSamples/step;
if round(preassign_length) ~= preassign_length  % when the file's duration is not divisible by the length of each "step"
    preassign_length = zeros(size(HGP_Zscore,1),1 + round(preassign_length));
else
    preassign_length = zeros(size(HGP_Zscore,1),preassign_length);
end


HGPScores = struct('mean', preassign_length, 'variance', preassign_length, 'skewness', preassign_length); %%zeros(length(channelIndices), floor(totalTimeSamples/30));
LFPScores = struct('mean', preassign_length, 'variance', preassign_length, 'skewness', preassign_length);
alphaScores = struct('mean', preassign_length, 'variance', preassign_length, 'skewness', preassign_length);
thetaScores = struct('mean', preassign_length, 'variance', preassign_length, 'skewness', preassign_length);

% HGPScores = struct('mean', {}, 'variance', {}, 'skewness', {}); %%zeros(length(channelIndices), floor(totalTimeSamples/30));
% LFPScores = struct('mean', {}, 'variance', {}, 'skewness', {});
% alphaScores = struct('mean', {}, 'variance', {}, 'skewness', {});
% thetaScores = struct('mean', {}, 'variance', {}, 'skewness', {});


for i = 1:step:totalTimeSamples
    if i+step-1 > totalTimeSamples 
        break
    end 
    
    % HGP
    meanHGP = mean(HGP_Zscore(:,i:i+step-1),2);
    varHGP = var(HGP_Zscore([i i+step-1]),0,2);
    skewHGP = skewness(HGP_Zscore([i i+step-1]),0,2);
    HGPScores.mean(:,((i-1)/step)+1) = meanHGP; HGPScores.variance(:,((i-1)/step)+1) = varHGP; HGPScores.skewness(:,((i-1)/step)+1) = skewHGP;  % Use (i-1)/step + 1 to make updates fit the structure
%     HGPScores(end+1) = temp;
    
    % LFP
    meanLFP = mean(LFP_Zscore([i i+step-1]),2);
    varLFP = var(LFP_Zscore([i i+step-1]),0,2);
    skewLFP = skewness(LFP_Zscore([i i+step-1]),0,2);
    LFPScores.mean(:,((i-1)/step)+1) = meanLFP; LFPScores.variance(:,((i-1)/step)+1) = varLFP; LFPScores.skewness(:,((i-1)/step)+1) = skewLFP;
%     LFPScores(end+1) = tempLFP;
    
    % alpha
    mean_alpha = mean(alpha_Zscore([i i+step-1]),2);
    var_alpha = var(alpha_Zscore([i i+step-1]),0,2);
    skew_alpha = skewness(alpha_Zscore([i i+step-1]),0,2);
    alphaScores.mean(:,((i-1)/step)+1) = mean_alpha; alphaScores.variance(:,((i-1)/step)+1) = var_alpha; alphaScores.skewness(:,((i-1)/step)+1) = skew_alpha;
%     alphaScores(end+1) = temp_alpha;
    
    % theta
    mean_theta = mean(theta_Zscore([i i+step-1]),2);
    var_theta = var(theta_Zscore([i i+step-1]),0,2);
    skew_theta = skewness(theta_Zscore([i i+step-1]),0,2);
    thetaScores.mean(:,((i-1)/step)+1) = mean_theta; thetaScores.variance(:,((i-1)/step)+1) = var_theta; thetaScores.skewness(:,((i-1)/step)+1) = skew_theta;
%     thetaScores(end+1) = temp_theta;

end

save(save_dir,'HGPScores', '-v7.3', '-append');

   
%moving averages to smooth all the values out, use gaussian kernel
%channelAverage = mean(config.trial{1}(chanelIndices,:),2);
%stanDev = std(config.trial{1}(channelIndices,:),0,2);
%channelAverage = mean(current.trial{1}(1:61,:),2);
%stanDev = std(current.trial{1}(1:61,:),0,2);





%save('/home/ishamie/Desktop/projects/24_7/NY442/442_day4_bpFeatures.mat', '*Scores', '*Zscore','*phase', '*envelope', '-v7.3');

end
