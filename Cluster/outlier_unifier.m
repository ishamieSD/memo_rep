% Make time bins that contain HGP peak outliers across all channels: 
% only peak locations that are separated sufficiently are counted as 
% time bin centers.
%
% Created by Xi Jiang, Sept. 5th, 2014
% Last edited by Xi Jiang, Sept. 9th, 2014
%
% Dependencies: outlier_wavedecompC.m
%
% Outputs:
%  good_bins: logical indices of (allegedly) artifact-free time bins, with
%             respect to the elements in all_bins_seq
%  bad_bins: logical indices of (allegedly) artifact-containing time bins, 
%            with respect to the elements in all_bins_seq
%  new_bins: indices of time bin start/end points. To select only
%            the non-artifact bins, use the indices in "new_starts" and
%            "new_ends".
%  new_starts: starting point indices (in terms of raw data) for clean time
%              bins
%  new_ends: similar to the above, except containing end points instead
%        
% Inputs:
%  data_out: one of the outputs of outlier_cleaning.m, with required field
%            "newPOI", representing numerical indices of all HGP peak 
%            outliers
%  data_raw: raw data in matrix form
%  methods: a structure that contains at least the following fields:
%           -time_bin: integer representing size of time bins, in seconds
%           -sfreq: sampling frequency (in Hz)
%           -wavthresh: threshold (in z-score) for artifact detection
%

function [good_bins,bad_bins,new_bins,new_starts,new_ends] = ...
            outlier_unifier(data_out,data_raw,methods)
%% locate well-spaced HGP peak outliers
    
    step = methods.sfreq * methods.time_bin;
    test = data_out.newPOI;
    binarypeaks = zeros(size(data_raw));

% set distance for centering, i.e. anything within the separation will be 
% included in the new time bins, which will be centered around peaks that
% are sufficiently far apart

    for i = 1:size(data_raw,1)
        test_dist = [Inf diff(test{i})];
        test{i} = test{i}(test_dist > step);
        binarypeaks(i,test{i}) = 1;
    end
    test = logical(sum(binarypeaks,1));
    test_ind = find(test);
    test_dist = [Inf diff(test_ind)];
    test_ind = test_ind(test_dist > step);

    clear binarypeaks test test_dist
    
    
%% preassign space for outputs and obtain time bin indices

    good_bins = zeros(size(data_raw,1),length(test_ind));
    bad_bins = zeros(size(data_raw,1),length(test_ind));
    new_bins.start = [];
    new_bins.end = [];
    
    % obtain "universal" time bins
    for j = 1:length(test_ind)
        if test_ind(j) <= methods.sfreq-1
            new_bins.start = horzcat(new_bins.start,1);
            new_bins.end = horzcat(new_bins.end,step);
        elseif test_ind(j) > length(data_out.smoothed)-methods.sfreq
            new_bins.start = horzcat(new_bins.start,length(data_raw)-(step-1));
            new_bins.end = horzcat(new_bins.end,length(data_raw));
        else
            new_bins.start = horzcat(new_bins.start,test_ind(j)-(round(step/2)-1));
            new_bins.end = horzcat(new_bins.end,test_ind(j)+round(step/2));
        end
    end

    
%% remove artifact bins

    for i = 1:size(data_raw,1)

        num_bins = length(test_ind);       % number of bins
        bins_LFP = cell(1,num_bins);
        step = methods.time_bin * methods.sfreq;  % time bin size in sampling points

    % obtain LFP and convert data to matrix form
        for j = 1:num_bins
            if test_ind(j) <= methods.sfreq-1
                bins_LFP{j} = data_raw(i,1:step);
                new_bins.start = horzcat(new_bins.start,1);
                new_bins.end = horzcat(new_bins.end,step);
            elseif test_ind(j) > length(data_out.smoothed)-methods.sfreq
                bins_LFP{j} = data_raw(i,...
                    length(data_raw)-(step-1):length(data_raw));
                new_bins.start = horzcat(new_bins.start,length(data_raw)-(step-1));
                new_bins.end = horzcat(new_bins.end,length(data_raw));
            else
                bins_LFP{j} = data_raw(i,...
                    (test_ind(j)-(round(step/2)-1)):test_ind(j)+round(step/2));
                new_bins.start = horzcat(new_bins.start,test_ind(j)-(round(step/2)-1));
                new_bins.end = horzcat(new_bins.end,test_ind(j)+round(step/2));
            end
        end
 
        all_bins_LFP = cell2mat(bins_LFP);
    % new matrix: time bins by sampling points (note the transpose)
        all_bins_LFP = reshape(all_bins_LFP,step,num_bins)';
    
    % reject artifacts based on wavelet decomposition
        [good_bin_ind,bad_bin_ind] = ...
            outlier_wavedecompC(all_bins_LFP,methods);
        
        good_bins(i,:) = good_bin_ind;
        bad_bins(i,:) = bad_bin_ind;


        if any(new_bins.start{i}-new_bins.end{i} >= 0)
            disp('ERROR!')
        end
        
        disp([num2str(i) ' channel(s) processed!'])
    
    end
    
    % obtain start/end indices of "clean" time bins
    good_cumulative = sum(good_bins,1);
    good_ind = good_cumulative == size(data_raw,1);
    new_starts = new_bins.start(good_ind);
    new_ends = new_bins.end(good_ind);
    
end