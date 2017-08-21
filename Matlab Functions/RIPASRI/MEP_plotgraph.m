function MEP_plotgraph

% PLOT ALL TRIALS - CONDITION x STEP

cd(fig_folder)
if exist('ALL_TRIALS', 'dir')==0; mkdir('ALL_TRIALS'); end
cd('ALL_TRIALS')
DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
cd(DateString)

s_rate=data_param(1,:);
points=data_param(2,:);

for m = 1:size(data_cond,3) % muscles
    h=figure; count=1;
    for j = 1:size(data_cond,5) % step - baseline, t1, t2
        duration_ms=points(j)/s_rate(j)*1000;
        step=duration_ms/points(j);
        time_axis=0:step:duration_ms-step; time_axis=time_axis';
        for c = 1:size(data_cond,4) % cond - single, licix3
            if m < 3 % FCR has different time windows than the other muscle
                curr_time_win = time_windows(:,c);
            else
                curr_time_win = time_windows_FCR(:,c);
            end
            
            ylim_all= abs(min(min(min(min(data_goodrecs(:,:,m,:,:))))));
            
            for i = 1:size(data_cond,2) % trials per condition
                
                %___________GRAPH PARAM__________
                curr_time_win_pt = round(curr_time_win*points(j)/duration_ms);
                max_start=curr_time_win_pt(1); max_end=curr_time_win_pt(2);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                hold on
                
                %___________ALL TRIALS____________
                %                 curr_data = squeeze(data_cond(:,:,m,c,j));
                curr_data = squeeze(data_goodrecs(:,:,m,c,j));
                subplot(3,8,count)
                plot(time_axis,curr_data); % data, trials, muscle, cond, step
                ylim([-ylim_all ylim_all]); xlim([curr_time_win(1)-15 curr_time_win(2)+15])
                set (gca,'Ydir','reverse')
                %                 fig_str=[cond_name{c}, ' - ',step_name{j}];
                %                 title(fig_str,'FontSize', 8)
                
                hold on; line([curr_time_win(1) curr_time_win(1)],[-3 3],'Color',[0 1 0]);
                hold on; line([curr_time_win(2) curr_time_win(2)],[-3 3],'Color',[0 1 0]);
                
            end
            %___________MEAN____________
            %                 mean_data = mean(squeeze(data_cond(:,:,m,c,j)),2);
            mean_data = mean(squeeze(data_goodrecs(:,:,m,c,j)),2);
            subplot(3,8,count+1)
            plot(time_axis,mean_data,'Color','b'); % data, trials, muscle, cond, step
            ylim([-ylim_all ylim_all]); xlim([curr_time_win(1)-15 curr_time_win(2)+15])
            set (gca,'Ydir','reverse')
            
            hold on; line([curr_time_win(1) curr_time_win(1)],[-3 3],'Color',[0 1 0]);
            hold on; line([curr_time_win(2) curr_time_win(2)],[-3 3],'Color',[0 1 0]);
            
            % Detected MAX
            [mean_max, mean_delay] = min(mean_data(max_start:max_end,:));
            h_peak_disp = round((mean_delay+curr_time_win_pt(1))*duration_ms/points(j));
            hold on; line([h_peak_disp h_peak_disp],[mean_max-0.5 mean_max+0.5],'Color',[1 0 1], 'Marker','.');
            % Detected MIN
            [mean_ref, mean_ref_delay] = max(mean_data(max_start+mean_delay-1:max_end,:));
            h_ref_disp = round((mean_ref_delay+curr_time_win_pt(1)+mean_delay)*duration_ms/points(j));
            hold on; line([h_ref_disp h_ref_disp],[mean_ref-0.5 mean_ref+0.5],'Color',[0 0 1], 'Marker','.');
            
            peak_str = sprintf('%.3f',mean_max);
            hold on; text(((curr_time_win(1)+curr_time_win(2))/2)-2.7,mean_max-1,peak_str)
            
            count=count+2;
        end
    end
    pause;
    print_str=m_name{m};
    disp(['Saving ',print_str]);
    print(h,print_str,'-dpng','-r900');
    save_fig=[print_str,'.fig'];
    savefig(save_fig);
    close
end


end