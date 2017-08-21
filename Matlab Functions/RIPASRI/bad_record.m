function resp_analysis(data_cond,rep)

mkdir('Checking Recs')
cd('Checking Recs')

DateString = datestr(clock); mkdir(DateString); 
cd(DateString);

fig_rep=input('Do you want to save the figures (y/n) ?:   ');
badrecs=zeros(max(n_trial),max(n_muscle), 4, 3);
% resp=zeros(max(n_trial), max(n_muscle), 4, 3,5);
for j = 1:size(data_cond,5) % step - baseline, t1, t2
    duration_ms=points(j)/s_rate(j)*1000;
    step=duration_ms/points(j);
    time_axis=0:step:duration_ms-step; time_axis=time_axis';
    for c = 1:size(data_cond,4) % cond - single, licix3
        for m = 1:size(data_cond,3) % muscles
            h=figure;
            for i = 1:size(data_cond,2) % trials per condition
                curr_data = squeeze(data_cond(1:points(j),i,m,c,j));
                subplot(4,3,i)
                plot(time_axis,curr_data); % data, trials, muscle, cond, step
                ylim_up=abs(min(curr_data))+0.5;
                %                 ylim_up=squeeze(min(min(min(min(min(min(data_cond,1)))))))-0.5;
                %                 ylim_down=squeeze(max(max(max(max(max(max(data_cond,1)))))))+0.5;
                %___________TIME WINDOWS VERTICAL LINES__________
                if m < 3 % FCR has different time windows than the other muscle
                    curr_time_win = time_windows(:,c);
                else
                    curr_time_win = time_windows_FCR(:,c);
                end
                hold on; line([curr_time_win(1) curr_time_win(1)],[-4 4],'Color',[1 0 0]);
                hold on; line([curr_time_win(2) curr_time_win(2)],[-4 4],'Color',[1 0 0]);
                
                %___________DETECT PEAKS__________
                start_pt=round(curr_time_win(1)*points(j)/duration_ms);
                end_pt=round(curr_time_win(2)*points(j)/duration_ms);
%                 start_pt=curr_time_win(1)*points(j)/duration_ms;
%                 end_pt=curr_time_win(2)*points(j)/duration_ms;
                [peak, delay] = min(curr_data(start_pt:end_pt));
                [ref, ref_delay] = max(curr_data(start_pt+delay-1:end_pt));
                amp = abs(peak - ref);
                resp(i,c,j,1) = peak;
                % TRYING TO REMOVE THE DEPENDENCY ON THE TIME WINDOW
                resp(i,m,c,j,2) = delay + start_pt;
                resp(i,m,c,j,3) = ref;
                resp(i,m,c,j,4) = ref_delay + start_pt + delay;
                resp(i,m,c,j,5) = amp;
                
                %___________PLOT DETECTED PEAKS_____________
                % H detected MAX
                peak_disp = (delay+start_pt)*duration_ms/points(j);
                hold on; line([peak_disp peak_disp],[peak-0.3 peak+0.3],'Color',[1 0 1], 'Marker','.');
                % H detected MIN
                ref_disp = (ref_delay + start_pt + delay)*duration_ms/points(j);
                hold on; line([ref_disp ref_disp],[ref-0.3 ref+0.3],'Color',[0 0 1], 'Marker','.');
                
                %___________GRAPH PARAM__________
%                 xlim([curr_time_win(1)-20 curr_time_win(2)+20])
                ylim([-ylim_up ylim_up])
                fig_str=[cond_name{c}, ' - ',m_name{m},' - ',step_name{j}, ' - trial ', int2str(i)];
                title(fig_str,'FontSize', 8)
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                set (gca,'Ydir','reverse')
                hold on
            end
            
            if strcmp(fig_rep,'y')==1
                print_str=[cond_name{c}, ' in ',m_name{m},' at ',step_name{j}];
                print(h,print_str,'-dpng','-r600');
                save_fig=[print_str,'.fig'];
                savefig(save_fig);
            end
            pause;
            close
            
            if strcmp(badrec_rep,'y') == 1
                rep=input('Please enter into [] separated with comas the bad channels:  ');
                if rep~= 0
                    for count = 1:length(rep)
                        badrecs(rep(count),m,c,j)=rep(count); % bad recordings, muscles, condition, step
                    end
                end
            end
            close all
        end
    end
end
cd(MLdata_folder)
DateString = datestr(clock);
save(['Resp_Amp',DateString],'resp')
save(['Bad_Recordings',DateString],'badrecs')

disp('End of BAD RECORDINGS IDENTIFICATION')

end