function [data_goodrecs] = remove_badrec(data_cond, badrecs, badrec_muscles, data_param, time_windows, time_windows_FCR)

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

data_goodrecs=zeros(size(data_cond,1),size(data_cond,2),size(data_cond,3),size(data_cond,4),size(data_cond,5));
for j = 1: size(badrecs,4)
    for c = 1:size(badrecs,3)
        for m = 1:size(badrecs,2)
            goodrecs=setdiff(1:size(data_cond,2),badrecs(:,m,c,j));
            for cc = 1 : length(goodrecs)
                data_goodrecs(:,cc,m,c,j)= data_cond(:,goodrecs(cc),m,c,j);
            end
        end
    end
end

check_rep=input('Do you want to check the removal of bad recordings (y/n) ?:  ');
if strcmp(check_rep,'y') == 1
    for m = 1:size(data_cond,3) % muscles
        if badrec_muscles(m)==1
            for j = 2:size(data_cond,5) % step - baseline, t1, t2
                duration_ms=points(j)/s_rate(j)*1000;
                step=duration_ms/points(j);
                time_axis=0:step:duration_ms-step; time_axis=time_axis';
                for c = 1:size(data_cond,4) % cond - single, licix3
                    
                    h=figure;
                    ylimup = abs(min(min(data_goodrecs(:,:,m,c,j))))+0.5;
                    for i = 1:size(data_cond,2) % trials per condition
                        subplot(4,3,i)
                        curr_data = squeeze(data_goodrecs(:,i,m,c,j));
                        if curr_data == 0
                            plot([0 length(time_axis)],[-5 5],'r')
                        else
                            plot(time_axis,curr_data); % data, trials, muscle, cond, step
                            if m < 4 % FCR has different time windows than the other muscle
                                curr_time_win = time_windows(:,c);
                            else
                                curr_time_win = time_windows_FCR(:,c);
                            end
                            hold on; line([curr_time_win(1) curr_time_win(1)],[-4 4],'Color',[1 0 0]);
                            hold on; line([curr_time_win(2) curr_time_win(2)],[-4 4],'Color',[1 0 0]);
                        end
                        
                        fig_str=[cond_name{c}, ' in ',m_name{m},' at ',step_name{j}, ' - trial ', int2str(i)];
                        title(['\fontsize{14}',fig_str])
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        set (gca,'Ydir','reverse')
                        ylim([-ylimup ylimup])
                        hold on
                    end
                    pause
                    print_str=['REMOVED',cond_name{c}, ' in ',m_name{m},' at ',step_name{j}];
                    print(h,print_str,'-dpng');
                    close all
                end
            end
        end
    end
end

end