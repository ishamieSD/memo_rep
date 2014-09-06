% reorganizes the output of align_rest_clin12.m into one fieldtrip struct
%
% Created by Xi Jiang, Aug. 6th, 2014
% Last edited by Xi Jiang, Aug. 8th, 2014
%
% Dependencies: N/A

function data_comb = clin_superglue(data)
   
% data: 1x2 cell with two aligned fieldtrip structures
    
    % concatenating and re-arranging channel labels while removing 
    % duplicates and preserving the indices of re-arranged channels in 
    % the variable "idx"
    data_comb.label = vertcat(data{1}.label,data{2}.label);
    [~,idx] = unique(strcat(data_comb.label(:,1)));
    data_comb.label = data_comb.label(idx,:);
    
    % sampling rate assumed to be the same
    data_comb.fsample = data{1}.fsample;
    
    % actual data (Nchan X Nsamples)
    data_comb.trial = vertcat(data{1}.trial{1},data{2}.trial{1});
    
    % again, assuming time vectors from clin1/2 to be the same
    data_comb.time = data{1}.time;
    
    % reorganizes data_comb.trial so that it matches the rearranged labels
    % It is necessary to use curly brackets for fieldtrip to recognize this
    % structure as raw data
    data_comb.trial = {data_comb.trial(idx,:)};
    
    data_comb.sampleinfo = [1 22116329];
    
end