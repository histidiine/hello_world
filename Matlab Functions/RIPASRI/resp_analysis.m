function [resp, resp_singleExtra] = resp_analysis(data_goodrecs,data_param,time_windows,time_windows_FCR)

% ________GET PARAMETERS_______
s_rate=data_param(1,:);
points=data_param(2,:);

n_muscle=size(data_goodrecs,3);
n_step=size(data_goodrecs,5);
n_cond=size(data_goodrecs,4);
n_trial=size(data_goodrecs,2);

%_____LABELS INITIALIZATION_____
dataset = ['b', 't1', 't2'];
ds_folder = {'Baseline', 'T1', 'T2'};
m_name={'ADM','FDI','FCR','FPB'};
step_name={'baseline','T1','T2'};
cond_name={'single pulse', 'LICI 50', 'LICI 100', 'LICI 200'};

%______FOLDER TO SAVE VERFICATION FIGURES______
if exist('MEP_Detection', 'dir')==0; mkdir('MEP_Detection'); end
cd('MEP_Detection')
DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
cd(DateString)

for m = 1:n_muscle % muscles
    for j = 1:n_step % step - baseline, t1, t2
        testingEx=sum(sum(sum(sum(all(data_goodrecs(:,:,:,:,j),5)))));
        if testingEx ~= 0
            h=figure; count=1; % for figure
            duration_ms=points(j)/s_rate(j)*1000;
            step=duration_ms/points(j);
            time_axis=0:step:duration_ms-step; time_axis=time_axis'; % for figure
            start_ylim = round(20*points(j)/duration_ms);
            end_ylim = round(58*points(j)/duration_ms);
            for c = 1:n_cond % cond - single, licix3
                if m < 3 % FCR has different time windows than the other muscle
                    curr_time_win = time_windows(:,c);
                else
                    curr_time_win = time_windows_FCR(:,c);
                end
                
                if c~=1 % to extract single data in the first pulse of paired cond
                    single_timewin = [30;55];
                else
                    single_timewin = curr_time_win;
                end
                
                ylim_all= abs(min(min(min(data_goodrecs(start_ylim:end_ylim,:,m,:,j)))))+1;
                for i = 1:n_trial % trials per condition
                    curr_data = squeeze(data_goodrecs(1:points(j),i,m,c,j));
                    
                    %___________DETECT PEAKS - MATRIX__________
                    start_pt=round(curr_time_win(1)*points(j)/duration_ms);
                    end_pt=round(curr_time_win(2)*points(j)/duration_ms);
                    
                    [peak, delay] = min(curr_data(start_pt:end_pt));
                    [ref, ref_delay] = max(curr_data(start_pt+delay-1:end_pt));
                    amp = abs(peak - ref);
                    resp(i,m,c,j,1) = peak;
                    resp(i,m,c,j,2) = delay + start_pt;
                    resp(i,m,c,j,3) = ref;
                    resp(i,m,c,j,4) = ref_delay + start_pt + delay;
                    resp(i,m,c,j,5) = amp;
                    
                    if c~= 1
                        start_pt_single=round(single_timewin(1)*points(j)/duration_ms);
                        end_pt_single=round(single_timewin(2)*points(j)/duration_ms);
                        
                        [peak_single, delay_single] = min(curr_data(start_pt_single:end_pt_single));
                        [ref_single, ref_delay_single] = max(curr_data(start_pt_single+delay_single-1:end_pt_single));
                        amp_single = abs(peak_single - ref_single);
                        resp_singleExtra(i,m,c-1,j,1) = peak_single;
                        resp_singleExtra(i,m,c-1,j,2) = delay_single + start_pt_single;
                        resp_singleExtra(i,m,c-1,j,3) = ref_single;
                        resp_singleExtra(i,m,c-1,j,4) = ref_delay_single + start_pt_single + delay_single;
                        resp_singleExtra(i,m,c-1,j,5) = amp_single;
                    else
                        start_pt_single=start_pt;
                        end_pt_single=end_pt;
                    end
                    
                    %_________FIGURES SAVED FOR FURTHER VERIFICATIONS_______
                    curr_time_win_pt = round(curr_time_win*points(j)/duration_ms);
                    single_timewin_pt = round(single_timewin*points(j)/duration_ms);
                    max_start=curr_time_win_pt(1); max_end=curr_time_win_pt(2);
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                    
                    subplot(n_cond,n_trial,count)
                    plot(time_axis,curr_data,'LineWidth',0.5); % data, trials, muscle, cond, step
                    ylim([-ylim_all ylim_all]); xlim([single_timewin(1)-15 curr_time_win(2)+15])
                    set(gca,'Ydir','reverse')
                    set(gca,'FontSize',4)
                    
                    if c~=1
                        hold on; line([single_timewin(1) single_timewin(1)],[-3 3],'Color',[0 1 0],'LineWidth',0.3);
                        hold on; line([single_timewin(2) single_timewin(2)],[-3 3],'Color',[0 1 0],'LineWidth',0.3);
                        
                        % Detected MAX
                        h_peak_disp_single = round((delay_single+single_timewin_pt(1))*duration_ms/points(j));
                        hold on; line([h_peak_disp_single h_peak_disp_single],[peak_single-0.5 peak_single+0.5],'Color',[1 0 1], 'Marker','.','LineWidth',0.5,'MarkerSize',0.5);
                        % Detected MIN
                        h_ref_disp_single = round((ref_delay_single+single_timewin_pt(1)+delay_single)*duration_ms/points(j));
                        hold on; line([h_ref_disp_single h_ref_disp_single],[ref_single-0.5 ref_single+0.5],'Color',[0 0 1], 'Marker','.','LineWidth',0.5,'MarkerSize',0.5);
                        
                        peak_str_single = sprintf('%.3f',amp_single);
                        hold on; text(single_timewin(1)+2,peak_single-0.6,peak_str_single,'FontSize',6)
                    end
                    
                    hold on; line([curr_time_win(1) curr_time_win(1)],[-3 3],'Color',[0 1 0],'LineWidth',0.3);
                    hold on; line([curr_time_win(2) curr_time_win(2)],[-3 3],'Color',[0 1 0],'LineWidth',0.3);
                    
                    % Detected MAX
                    h_peak_disp = round((delay+curr_time_win_pt(1))*duration_ms/points(j));
                    hold on; line([h_peak_disp h_peak_disp],[peak-0.5 peak+0.5],'Color',[1 0 1], 'Marker','.','LineWidth',0.5,'MarkerSize',0.5);
                    % Detected MIN
                    h_ref_disp = round((ref_delay+curr_time_win_pt(1)+delay)*duration_ms/points(j));
                    hold on; line([h_ref_disp h_ref_disp],[ref-0.5 ref+0.5],'Color',[0 0 1], 'Marker','.','LineWidth',0.5,'MarkerSize',0.5);
                    
                    peak_str = sprintf('%.3f',amp);
                    hold on; text(curr_time_win(1)+2,peak-0.6,peak_str,'FontSize',6)
                    %                 hold on; text(((curr_time_win(1)+curr_time_win(2))/2)-2.7,peak-1,peak_str,'FontSize',6)
                    
                    count=count+1;
                end
            end
            print_str=[m_name{m},'-',step_name{j}];
            disp(['Saving ',print_str]);
            export_fig('print_str' ,'-bmp', '-nocrop', '-r160');
            %         print(h,print_str,'-dpng','-r900');
            close
        end
        count=1;
    end
end

disp('End of RESPONSE ANALYSIS')

end



%------------------------------------------------------------------------
% CIMETIERE
%------------------------------------------------------------------------

% for m = 1:n_muscle % muscles
%     for j = 1:n_step % step - baseline, t1, t2
%         h=figure; count=1; % for figure
%         duration_ms=points(j)/s_rate(j)*1000;
%         step=duration_ms/points(j);
%         time_axis=0:step:duration_ms-step; time_axis=time_axis'; % for figure
%         for c = 1:n_cond % cond - single, licix3
%             if m < 3 % FCR has different time windows than the other muscle
%                 curr_time_win = time_windows(:,c);
%             else
%                 curr_time_win = time_windows_FCR(:,c);
%             end
%
%             if c~=1
%                single_timewin = [30;55];
%             end
%
%             ylim_all= abs(min(min(min(data_goodrecs(:,:,m,:,j)))))+0.5;
%             for i = 1:n_trial % trials per condition
%                 curr_data = squeeze(data_goodrecs(1:points(j),i,m,c,j));
%
%                 %___________DETECT PEAKS - MATRIX__________
%                 start_pt=round(curr_time_win(1)*points(j)/duration_ms);
%                 end_pt=round(curr_time_win(2)*points(j)/duration_ms);
%
%                 [peak, delay] = min(curr_data(start_pt:end_pt));
%                 [ref, ref_delay] = max(curr_data(start_pt+delay-1:end_pt));
%                 amp = abs(peak - ref);
%                 resp(i,m,c,j,1) = peak;
%                 resp(i,m,c,j,2) = delay + start_pt;
%                 resp(i,m,c,j,3) = ref;
%                 resp(i,m,c,j,4) = ref_delay + start_pt + delay;
%                 resp(i,m,c,j,5) = amp;
%
%                 %_________FIGURES SAVED FOR FURTHER VERIFICATIONS_______
%                 curr_time_win_pt = round(curr_time_win*points(j)/duration_ms);
%                 max_start=curr_time_win_pt(1); max_end=curr_time_win_pt(2);
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%
%                 subplot(n_cond,n_trial,count)
%                 plot(time_axis,curr_data,'LineWidth',0.5); % data, trials, muscle, cond, step
%                 ylim([-ylim_all ylim_all]); xlim([curr_time_win(1)-15 curr_time_win(2)+15])
%                 set(gca,'Ydir','reverse')
%                 set(gca,'FontSize',4)
%
%                 hold on; line([curr_time_win(1) curr_time_win(1)],[-3 3],'Color',[0 1 0],'LineWidth',0.3);
%                 hold on; line([curr_time_win(2) curr_time_win(2)],[-3 3],'Color',[0 1 0],'LineWidth',0.3);
%
%                 % Detected MAX
%                 h_peak_disp = round((delay+curr_time_win_pt(1))*duration_ms/points(j));
%                 hold on; line([h_peak_disp h_peak_disp],[peak-0.5 peak+0.5],'Color',[1 0 1], 'Marker','.','LineWidth',0.5,'MarkerSize',0.5);
%                 % Detected MIN
%                 h_ref_disp = round((ref_delay+curr_time_win_pt(1)+delay)*duration_ms/points(j));
%                 hold on; line([h_ref_disp h_ref_disp],[ref-0.5 ref+0.5],'Color',[0 0 1], 'Marker','.','LineWidth',0.5,'MarkerSize',0.5);
%
%                 peak_str = sprintf('%.3f',amp);
%                 hold on; text(curr_time_win(1)+2,peak-1,peak_str,'FontSize',6)
% %                 hold on; text(((curr_time_win(1)+curr_time_win(2))/2)-2.7,peak-1,peak_str,'FontSize',6)
%
%                 count=count+1;
%             end
%         end
%         print_str=[m_name{m},'-',step_name{j}];
%         disp(['Saving ',print_str]);
%         export_fig('print_str' ,'-bmp', '-nocrop', '-r160');
% %         print(h,print_str,'-dpng','-r900');
%         close
%         count=1;
%     end
% end
