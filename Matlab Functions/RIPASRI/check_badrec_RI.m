function [badrecs, badrec_muscles] = check_badrec_RI(data_cond, data_param, hdelays)

% ________GET PARAMETERS_______
s_rate=data_param(1,:);
points=data_param(2,:);
n_muscle=data_param(3,:);
n_trial=data_param(4,:);

%_____LABELS INITIALIZATION_____
dataset = ['b', 't1', 't2'];
ds_folder = {'Baseline', 'T1', 'T2'};
m_name={'ADM','FDI','FCR','FPB'};
step_name={'baseline','T1','T2'};
cond_name={'single pulse', 'LICI 50', 'LICI 100', 'LICI 200'};

% _______FIGURES_______
fig_rep=input('Do you want to save the figures (y/n) ?:   ');
if strcmp(fig_rep,'y')==1
    % CREATE FOLDER TO SAVE FIGURES
    mkdir('Checking Recs')
    cd('Checking Recs')
    DateString = datestr(clock); mkdir(DateString);
    cd(DateString);
end

for m = 1:size(data_cond,3) % muscles
    h=figure; count=1;
    for j = 1:size(data_cond,5) % step - baseline, t1, t2
        duration_ms=points(j)/s_rate(j)*1000;
        step=duration_ms/points(j);
        time_axis=0:step:duration_ms-step; time_axis=time_axis';
        for c = 1:size(data_cond,4) % cond - single, licix3
            for i = 1:size(data_cond,2) % trials per condition
                curr_data = squeeze(data_cond(1:points(j),i,m,c,j));
                subplot(3,5,count)
                plot(time_axis,curr_data); % data, trials, muscle, cond, step
                ylim_up=abs(min(curr_data))+0.5;
                
                hold on; line([hdelays(1) hdelays(1)],[-4 4],'Color',[1 0 0]);
                hold on; line([hdelays(2) hdelays(2)],[-4 4],'Color',[1 0 0]);
                
                %___________DETECT PEAKS__________
                start_pt=round(hdelays(1)*points(j)/duration_ms);
                end_pt=round(hdelays(2)*points(j)/duration_ms);
                [peak, delay] = min(curr_data(start_pt:end_pt));
                [ref, ref_delay] = max(curr_data(start_pt+delay-1:end_pt));
                
                %___________PLOT DETECTED PEAKS_____________
%                 % Detected MAX
%                 peak_disp = (delay+start_pt)*duration_ms/points(j);
%                 hold on; line([peak_disp peak_disp],[peak-0.3 peak+0.3],'Color',[1 0 1], 'Marker','.');
%                 % Detected MIN
%                 ref_disp = (ref_delay + start_pt + delay)*duration_ms/points(j);
%                 hold on; line([ref_disp ref_disp],[ref-0.3 ref+0.3],'Color',[0 0 1], 'Marker','.');
                
                %___________GRAPH PARAM__________
                xlim([hdelays(1)-20 hdelays(2)+20])
                ylim([-ylim_up ylim_up])
                fig_str=[cond_name{c}, ' - ',m_name{m},' - ',step_name{j}];
                title(fig_str,'FontSize', 8)
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                set (gca,'Ydir','reverse')
                hold on
            end
            count=count+1;
        end
    end
    pause;
    if strcmp(fig_rep,'y')==1
        print_str=[cond_name{c}, ' in ',m_name{m},' at ',step_name{j}];
        print(h,print_str,'-dpng','-r600');
%         save_fig=[print_str,'.fig'];
%         savefig(save_fig);
    end
    close
end

badrec_rep = input('Did you see any bad recordings that need to be removed for proper MEP analysis? (y/n):   ');
badrec_muscles=zeros(max(n_muscle),1);
if strcmp(badrec_rep,'y')==1
    badrecs=zeros(max(n_trial),max(n_muscle), 4, 3);
    for m = 1:size(data_cond,3) % muscles
        check_muscle=input(['Where there any bad recordings in ', m_name{m}, ' ? (y/n):   ']);
        for j = 1:size(data_cond,5) % step - baseline, t1, t2
            duration_ms=points(j)/s_rate(j)*1000;
            step=duration_ms/points(j);
            time_axis=0:step:duration_ms-step; time_axis=time_axis';
            for c = 1:size(data_cond,4) % cond - single, licix3
                if strcmp(check_muscle,'y')==1
                    badrec_muscles(m)=1;
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
                            hdelays = time_windows(:,c);
                        else
                            hdelays = time_windows_FCR(:,c);
                        end
                        hold on; line([hdelays(1) hdelays(1)],[-4 4],'Color',[1 0 0]);
                        hold on; line([hdelays(2) hdelays(2)],[-4 4],'Color',[1 0 0]);
                        
                        %___________DETECT PEAKS__________
                        start_pt=round(hdelays(1)*points(j)/duration_ms);
                        end_pt=round(hdelays(2)*points(j)/duration_ms);
                        [peak, delay] = min(curr_data(start_pt:end_pt));
                        [ref, ref_delay] = max(curr_data(start_pt+delay-1:end_pt));
                        
                        %___________PLOT DETECTED PEAKS_____________
                        % Detected MAX
                        peak_disp = (delay+start_pt)*duration_ms/points(j);
                        hold on; line([peak_disp peak_disp],[peak-0.3 peak+0.3],'Color',[1 0 1], 'Marker','.');
                        % Detected MIN
                        ref_disp = (ref_delay + start_pt + delay)*duration_ms/points(j);
                        hold on; line([ref_disp ref_disp],[ref-0.3 ref+0.3],'Color',[0 0 1], 'Marker','.');
                        
                        %___________GRAPH PARAM__________
                        %                 xlim([hdelays(1)-20 hdelays(2)+20])
                        ylim([-ylim_up ylim_up])
                        fig_str=[cond_name{c}, ' - ',m_name{m},' - ',step_name{j}, ' - trial ', int2str(i)];
                        title(fig_str,'FontSize', 8)
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        set (gca,'Ydir','reverse')
                        hold on
                    end
                    rep=input('Please enter into [] separated with comas the bad channels (if non enter 0):  ');
                    if rep~= 0
                        for count = 1:length(rep)
                            badrecs(rep(count),m,c,j)=rep(count); % bad recordings, muscles, condition, step
                        end
                    end
                    close all
                end
            end
        end
    end
else
    badrecs = 0;
    disp('----------------------------------------');
    disp('No bad recordings to remove! Wouhouuuuu!');
    disp('----------------------------------------');
end

end