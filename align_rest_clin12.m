

function [data, match_dur] = align_rest_clin12(params, data)

% Created by XJ on 07/21/2014
% Right now can only handle clin1&2 with only one data file (no parts)

% Edited by Xi Jiang, Aug. 6th, 2014 (commented out faulty trigger
% comparison segment)
% In order to use this script for continuous data, triggers need to be
% obtained from the modified kb_trig_mult_XJ_edit function

% Outputs:
%   data: cell containing two fieldtrip data structures (clin1/2)
%   match_dur: 4-element matrix marking the segments aligned
%       e.g. [start_of_clin1,end_of_clin1;start_of_clin2,end_of_clin2]
% Inputs:
%   params: a structure containing the fields "my_trl" and "clinsys"
%   "clinsys": [1 2] or [2 1], representing order of clin1/2 files
%   "my_trl": output from kb_trig_mult_XJ_edit.m

triggers = {};
trig_len = [];
last_samp = [];

for clin = 1:2
    
    fileid = find(params.clinsys==clin);
    triggers{clin} = params.my_trl{fileid};
    triggers{clin}(:,3) = [0; diff(triggers{clin}(:,1))];
    trig_len(clin) = size(triggers{clin},1);
    last_samp(clin) = size(data{fileid}.trial{1}, 2);
    
end

[max_len, moretrig_clin] = max(trig_len);   % The clinical system 
lesstrig_clin = find([1,2]~=moretrig_clin);


%% For each pair of triggers (m, m+1) in the clinical system with more triggers, loop through
% all trigger pairs in the other clinical system and find all trigger
% pairs (l, l+1) with matching trigger codes and inter-trigger interval,
% Then, find the entire matching sequences starting from (m, m+1) and (l, l+1) (until no longer matching).
% Finally, use the longest sequence among all M x L sequences
% This is probably quite redundant, but will not miss any matching sequences...

all_seeds = cell(1, max_len-1);
all_seeds_len = zeros(1, max_len-1);

for m = 1:max_len-1
    
    beg_m = m;
    sub_seeds =      cell(1, trig_len(lesstrig_clin)-1);
    sub_seeds_len = zeros(1, trig_len(lesstrig_clin)-1);
    
    for l = 1:trig_len(lesstrig_clin)-1
        
        beg_l = l;
        m = beg_m;
        match = 1;
        
        while l<=trig_len(lesstrig_clin)-1 && m<=max_len-1 && match==1
            
            if (triggers{moretrig_clin}(m,2) == triggers{lesstrig_clin}(l,2)) && ...
                    (triggers{moretrig_clin}(m+1,2) == triggers{lesstrig_clin}(l+1,2)) && ...
                    (abs(triggers{moretrig_clin}(m+1,3) - triggers{lesstrig_clin}(l+1,3)) <= 3)
                
                sub_seeds{beg_l} = [sub_seeds{beg_l}, [m, m+1; l, l+1]];
                sub_seeds_len(beg_l) = sub_seeds_len(beg_l)+1;
                
                m = m+1;
                
            else
                if(sub_seeds_len(beg_l)>0)
                    match = 0;
                end
            end
            
            l = l+1;
            
        end
        
    end
    
    [maxlen, maxlen_idx] = max(sub_seeds_len);
    all_seeds{beg_m} = sub_seeds{maxlen_idx};
    all_seeds_len(beg_m) = maxlen;
    
end

figure, plot(all_seeds_len, '.-');
title('all-seeds-len');

[~, maxlen_idx] = max(all_seeds_len);
seed = all_seeds{maxlen_idx};


%% Check if the intervals of matching trigger indices follow the 1 0 1 0 1 ... 1 pattern
% Then figure out during which samples to truncate the data so that clin1&2
% become aligned and also have (almost) the same numbers of samples

seed_diff1 = diff(seed(1,:));
seed_diff2 = diff(seed(2,:));

match_dur = [];

if seed_diff1(1)==1 && sum(abs(seed_diff1(2:end) - repmat([0,1], 1, length(seed_diff1(2:end))/2)))==0 ...
&& seed_diff2(1)==1 && sum(abs(seed_diff2(2:end) - repmat([0,1], 1, length(seed_diff2(2:end))/2)))==0

    if seed(1,1)==1 && seed(2,1)==1
        if triggers{moretrig_clin}(1,1) >= triggers{lesstrig_clin}(1,1)
            match_dur(lesstrig_clin,1) = 1;
            match_dur(moretrig_clin,1) = triggers{moretrig_clin}(1,1) - triggers{lesstrig_clin}(1,1) + 1;
        else
            match_dur(moretrig_clin,1) = 1;
            match_dur(lesstrig_clin,1) = triggers{lesstrig_clin}(1,1) - triggers{moretrig_clin}(1,1) + 1;
        end
    elseif seed(1,1)==1 && seed(2,1)>1
        match_dur(moretrig_clin,1) = 1;
        match_dur(lesstrig_clin,1) = triggers{lesstrig_clin}(seed(2,1),1) - triggers{moretrig_clin}(1,1) + 1;
    elseif seed(2,1)==1 && seed(1,1)>1
        match_dur(lesstrig_clin,1) = 1;
        match_dur(moretrig_clin,1) = triggers{moretrig_clin}(seed(1,1),1) - triggers{lesstrig_clin}(1,1) + 1;
    end
    
    
    if seed(1,end)==max_len && seed(2,end)==trig_len(lesstrig_clin)
        if (last_samp(moretrig_clin) - triggers{moretrig_clin}(end,1)) >= ...
           (last_samp(lesstrig_clin) - triggers{lesstrig_clin}(end,1))
            match_dur(lesstrig_clin,2) = last_samp(lesstrig_clin);
            match_dur(moretrig_clin,2) = triggers{moretrig_clin}(end,1) + ...
                                         (last_samp(lesstrig_clin) - triggers{lesstrig_clin}(end,1));
        else
            match_dur(moretrig_clin,2) = last_samp(moretrig_clin);
            match_dur(lesstrig_clin,2) = triggers{lesstrig_clin}(end,1) + ...
                                         (last_samp(moretrig_clin) - triggers{moretrig_clin}(end,1));
        end
    elseif seed(1,end)==max_len && seed(2,end)<trig_len(lesstrig_clin)
        match_dur(moretrig_clin,2) = last_samp(moretrig_clin);
        match_dur(lesstrig_clin,2) = triggers{lesstrig_clin}(seed(2,end),1) + ...
                                     (last_samp(moretrig_clin) - triggers{moretrig_clin}(end,1));
    elseif seed(2,end)==trig_len(lesstrig_clin) && seed(1,end)<max_len
        match_dur(lesstrig_clin,2) = last_samp(lesstrig_clin);
        match_dur(moretrig_clin,2) = triggers{moretrig_clin}(seed(1,end),1) + ...
                                     (last_samp(lesstrig_clin) - triggers{lesstrig_clin}(end,1));
    end

end


%% Extract the overlapping portions of the data from clin1&2

diff_totlen = length(match_dur(1,1):match_dur(1,2)) - length(match_dur(2,1):match_dur(2,2)) 
% There could be a difference because not all the inter-trigger intervals 
% are exactly the same between clin1&2. This number should be quite small 
% as there shouldn't be a cumulative effect across multiple triggers

if all(all(match_dur>0))
    
    % Do this simplest correction for now as the number of samples to be
    % corrected can't go over the index range of data
    if(diff_totlen>0)
        match_dur(1,2) = match_dur(1,2) - diff_totlen;
    elseif(diff_totlen<0)
        match_dur(2,2) = match_dur(2,2) + diff_totlen;
    end
    
end
    
for clin = 1:2
    fileid = find(params.clinsys==clin);
    data{fileid}.trial{1} = data{fileid}.trial{1}(:, match_dur(clin,1):match_dur(clin,2));
    data{fileid}.time{1}  = data{fileid}.time{1}(1:length(match_dur(clin,1):match_dur(clin,2)));
end
    

%% Detect trigger times and codes
% Check if the rest-period data have been correctly aligned

% triggers_method1 = kb_trig_mult_XJ_edit(params, match_dur); % Adapted from Kristen's script
% triggers_method2 = ft_NYU_iEEG_trialfun_XJedit(params, match_dur); % Adapted from UCSD script 

% compare trigger times & codes derived by the two methods
% for r = 1:params.filenum
%     
%     T1 = triggers_method1{r};
%     T2 = triggers_method2{r};
%     
%     if size(T1,1)~=size(T2,1) || sum(abs(T1(:,2)-T2(:,2)))~=0
%         error('Major error with trigger detection!');
%     else
%         mean(abs(T1(:,1)-T2(:,1)))
%     end
% end

% new_triggers = triggers_method2; % Use UCSD triggers since they are earlier by 1 sample...
% 
% disp('Compare trigger times and code for CLIN1 vs CLIN2:');
% [new_triggers{1}(:,1), new_triggers{2}(:,1)]
% [new_triggers{1}(:,2), new_triggers{2}(:,2)]
% max_trigtime_diff = max(abs(new_triggers{1}(:,1) - new_triggers{2}(:,1)))
% trigcode_diff =     sum(abs(new_triggers{1}(:,2) - new_triggers{2}(:,2)))











