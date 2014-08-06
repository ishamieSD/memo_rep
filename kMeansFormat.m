
function [data_for_kmeans] =  kMeansFormat(name, varargin)
%created by ishamie in July, 2014

% name = file name, don't end it with .mat, it will be added automatically. file will be saved as kmeans_format_'name'.mat. example: HGP_alpha or alpha_theta_LFP
% varargin = [[binned data], w_type, path] where w_type is true or false whether a text
%   file is desired or not. default = false. path = path to save data, do not
%   include the file name, just the path, and end with /. default is current directory. example =/space/mdeh4/1/halgdev/projects/ishamie/   
%   Required in workspace: [binned data], a column vector that contains structs that contain the
%       binned waveforms. example: [HGPBinned] or [HGPBinned; alphaBinned] . these
%       structs must contain fields named after features, such as "mean", "variance", and "skewness"

% convert all the wave form bins into one long nXp matrix where n is num
% instances of time windows and p is 3 bandpassed wavesX 3 moments
% (mean var skew) X numChannels channels 


%% set up prerequisites
if nargin < 3 
    w_type = false;
    path = '' ;
elseif nargin == 3
    w_type = varargin{2};
    path = '' ;
else 
    w_type = varargin{2};
    path = varargin{3};

end

wavesBinned = varargin{1};
totalWaves = size(wavesBinned,1); %sees how many different bp'ed waves there are
n = length(wavesBinned(1,:));
numChannels = length(wavesBinned(1,1).mean);
numFeaturesPerWave = numChannels * 3; 

%%initialize the nXp matrix
data_for_kmeans = zeros(n,(numChannels*3)*(totalWaves)); %numChannels * 3 (mean, variance, skewness) * totalWaves (1,2,..)


%construct the matrix
for i = 1:n
    for wave = 1:(totalWaves)
        data_for_kmeans(i,numFeaturesPerWave*(wave-1)+1:numFeaturesPerWave*wave) = [wavesBinned(wave,i).mean', wavesBinned(wave,i).variance', wavesBinned(wave,i).skewness']  ;
    end 
end

%save to a .mat file
file = strcat(path, 'kmeans_format_', name, '.mat');
save(file,'data_for_kmeans');

%save to a .txt file if w_type == true
if w_type == true
    txtFile = strcat(path, 'kmeans_format_', name, '.txt');
    dlmwrite(txtFile, data_for_kmeans);
end

end

