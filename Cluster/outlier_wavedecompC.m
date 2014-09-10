% Evaluate HGP peak outlier bins for jump/spike artifacts via wavelet
% decomposition (1-level haar and 2-level db4)
%
% Created by Xi Jiang, Sept. 5th, 2014
% !May be merged with outlier_cleaning after performance tests
%
% Dependency: N/A
%
% data_LFP:    binned raw data (LFP)
% chan:        numerical, denoting the channel being used for decomposition
% methods:     structure that contains at least the following:
%              -time_bin: integer representing size of time bins, in seconds
%              -sfreq: sampling frequency (in Hz)
%              -wavthresh: threshold (in z-score) for artifact detection
%              Optional:
%              -haar_only: can take values 0/1 (yes/no)
%              -db4_only: similar to the field above

function [good_bin_ind,bad_bin_ind] = outlier_wavedecompC(all_bins_LFP,methods)
%% obtain LFP bins corresponding to the HGP peak outliers

    num_bins = size(all_bins_LFP,1);
    
    % initialize bin indices
    good_bin_ind_haar = ones(1,num_bins);
    bad_bin_ind_haar = zeros(1,num_bins);
    good_bin_ind_db4 = ones(1,num_bins);
    bad_bin_ind_db4 = zeros(1,num_bins);
 
    
%% 1-level haar decomposition

    if ~isfield(methods,'db4_only')
        methods.db4_only = 0;
    end
    if methods.db4_only ~= 1
        for i = 1:num_bins
            test = all_bins_LFP(i,:);
            level = 1;
        
            [c,l] = wavedec(test,level,'haar');
            [d1] = detcoef(c,l,level);
        
            % convert component to absolute z-score
            d1 = abs(zscore(d1));
        
            % assign binary status to bin indices (for eventual conversion to
            % logical)
            if max(d1) > methods.wavthresh
                good_bin_ind_haar(i) = 0;
                bad_bin_ind_haar(i) = 1;
            end
        end
    end

%% 2-level db4 decomposition

    if ~isfield(methods,'haar')
        methods.haar_only = 0;
    end
    if methods.haar_only ~= 1
        for i = 1:num_bins
            test = all_bins_LFP(i,:);
            level = 2;
            
            [c,l] = wavedec(test,level,'db4');
            [d1,d2] = detcoef(c,l,1:level);
            % convert components to absolute z-scores
            d1 = abs(zscore(d1));
            d2 = abs(zscore(d2));
            
            % assign binary status to bin indices (for eventual conversion to
            % logical)
            if max(d1) > methods.wavthresh || max(d2) > methods.wavthresh
                good_bin_ind_db4(i) = 0;
                bad_bin_ind_db4(i) = 1;
            end
        end
    end
    
%% convert bin indices to logical
    
    % find common indices for both decompositions
%     good_bin_ind = good_bin_ind_db4.*good_bin_ind_haar;
%     bad_bin_ind = bad_bin_ind_db4.*bad_bin_ind_haar;

    % combine indices
    good_bin_ind = good_bin_ind_db4+good_bin_ind_haar;
    bad_bin_ind = bad_bin_ind_db4+bad_bin_ind_haar;    

    good_bin_ind = logical(good_bin_ind);
    bad_bin_ind = logical(bad_bin_ind);
    

end