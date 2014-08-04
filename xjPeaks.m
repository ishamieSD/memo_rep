peaks_x = {};
peaks_y = {};
for ch = 1:length(HGP_rest_data.label)
   
    HGP_cont = HGP_rest_data.trial{1}(ch,:);
   
    HGP_cont_mean = mean(HGP_cont);
    HGP_cont_std = std(HGP_cont);
    HGP_threshold1 = HGP_cont_mean + 2.5*HGP_cont_std;
   
    [peaks_y{ch}, peaks_x{ch}] = findpeaks(HGP_cont,'MINPEAKHEIGHT',HGP_threshold);

end

%For generating the raster plots:
%Here select_chan_idx has four cells, each containing selected electrodes that were responsive to sound, face, building, or face&building stimuli in our MS task.

    tick_color = {'g','b','r',[0.48 0.06 0.89]}; % match to {'Sound','Faces','Buildings','Faces+Buildings'}
   
    k=0;
    f_HGP_peaktimes = figure;
    chan_name_all = {};
    for i = 1:length(select_chan_idx)  % don't want to plot non-responsive channels
        for j = 1:length(select_chan_idx{i})
           
            k=k+1;
            chan_idx = select_chan_idx{i,1}(j);
            chan_name = select_chan{i}{j};
            chan_name = strrep(chan_name,'_','-');
            chan_name_all{k} = chan_name;
           
            peak_times = HGP_rest_data.time{1}(peaks_x{chan_idx});
           
            figure(f_HGP_peaktimes),
            hold on, plot([peak_times;peak_times],[k*ones(1,length(peak_times))-0.40;k*ones(1,length(peak_times))+0.40],...
                'Color',tick_color{i},'LineWidth',2);
        end
    end
    figure(f_HGP_peaktimes),
    set(gca,'YTick',1:length(chan_name_all),'YTickLabel',chan_name_all,'YLim',[0 length(chan_name_all)+1]);
    xlabel('Time (sec)');
