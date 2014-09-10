% Make time bins that contain HGP peak outliers across all channels: 
% only peak locations that are separated sufficiently are counted as 
% time bin centers.
%
% Created by Xi Jiang, Sept. 5th, 2014
% Last edited by Xi Jiang, Sept. 10th, 2014
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
%  good_peaks: indices of HGP peak outliers in each channel that happen to
%              fall into one of the new time bins 
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
%           Optionally:
%           -savedir: string path of directory for saving outputs AND 
%                     filename, e.g. '/home/test/output.mat'
%

function [good_bins,bad_bins,new_bins,new_starts,new_ends,good_peaks] = ...
            outlier_unifier(data_out,data_raw,methods)
%% locate well-spaced HGP peak outliers
    
    step = methods.sfreq * methods.time_bin;
    test = data_out.newPOI;
    binarypeaks = zeros(size(data_raw));
    chan_num = size(data_raw,1);            % number of channels
    data_size = length(data_raw);           % number of sampling points

% set distance for centering, i.e. anything within the separation will be 
% included in the new time bins, which will be centered around peaks that
% are sufficiently far apart

    for i = 1:chan_num
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

    good_bins = zeros(chan_num,length(test_ind));
    bad_bins = zeros(chan_num,length(test_ind));
    new_bins.start = [];
    new_bins.end = [];
    
    % obtain "universal" time bins
    for j = 1:length(test_ind)
        if test_ind(j) <= methods.sfreq-1
            new_bins.start = horzcat(new_bins.start,1);
            new_bins.end = horzcat(new_bins.end,step);
        elseif test_ind(j) > data_size-methods.sfreq
            new_bins.start = horzcat(new_bins.start,data_size-(step-1));
            new_bins.end = horzcat(new_bins.end,data_size);
        else
            new_bins.start = horzcat(new_bins.start,test_ind(j)-(round(step/2)-1));
            new_bins.end = horzcat(new_bins.end,test_ind(j)+round(step/2));
        end
    end

    
%% remove artifact bins

    for i = 1:chan_num

        num_bins = length(test_ind);       % number of bins
        bins_LFP = cell(1,num_bins);
        step = methods.time_bin * methods.sfreq;  % time bin size in sampling points

    % obtain LFP and convert data to matrix form
        for j = 1:num_bins
            if test_ind(j) <= methods.sfreq-1
                bins_LFP{j} = data_raw(i,1:step);
                new_bins.start = horzcat(new_bins.start,1);
                new_bins.end = horzcat(new_bins.end,step);
            elseif test_ind(j) > data_size-methods.sfreq
                bins_LFP{j} = data_raw(i,...
                    data_size-(step-1):data_size);
                new_bins.start = horzcat(new_bins.start,data_size-(step-1));
                new_bins.end = horzcat(new_bins.end,data_size);
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
        
        disp([num2str(i) ' channel(s) processed!'])
    
    end
    
    % obtain start/end indices of "clean" time bins
    good_cumulative = sum(good_bins,1);
    good_ind = good_cumulative == chan_num;
    new_starts = new_bins.start(good_ind);
    new_ends = new_bins.end(good_ind);
    
%% obtain HGP peak indices that lie within the good bins

    % represet locations of all HGP peak outliers
    all_binary_peaks = zeros(chan_num,data_size);
    for i = 1:chan_num
        all_binary_peaks(i,data_out.newPOI{i}) = 1;
    end
    
    % represent locations of newly made time bins
    binary_bins = zeros(1,data_size);
    for i = 1:length(new_starts)
        binary_bins(new_starts(i):new_ends(i)) = 1;
    end
    
    % after addition of binary indices, '2's represent the indices where
    % HGP peak lands in one of the good bins
    all_binary_peaks = all_binary_peaks + repmat(binary_bins,chan_num,1);

    good_peaks = cell(chan_num,1);
    for i = 1:chan_num
        good_peaks{i} = find(all_binary_peaks(i,:) == 2);
    end
    
%% save outputs, if desired

    if isfield(methods,'savedir')
        save(methods.savedir,'good_bins','bad_bins','new_bins',...
                            'new_starts','new_ends','good_peaks')
    end
    
    
end