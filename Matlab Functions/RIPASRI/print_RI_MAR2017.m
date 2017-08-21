function print_RI_MAR2017(data_cond,folder,fig_name)

data_cond = data_cond.*-1; % To get the data with minus up and positive down
cond_name={'-1ms ISI', '-0.5ms ISI', 'synchronous','0.5ms ISI', '1ms ISI',...
    '2ms ISI','5ms ISI','7.5ms ISI','10ms ISI','25ms ISI','50ms ISI','75ms ISI','100ms ISI','200ms ISI', 'single pulse'};
cd(folder)
for j = 1:size(data_cond,4) % STEP LOOP - baseline, post
    duration_ms=points(j)/s_rate(j)*1000;
    step=duration_ms/points(j);
    time_axis=0:step:duration_ms-step; time_axis=time_axis';
    for c = 1:size(data_cond,3) % COND LOOP
        h=figure;
        for i=1:size(data_cond,2) % TRIAL LOOP
              plot(time_axis,squeeze(data_cond(:,i,c,j)));
%             if c < 4
%                 plot_start = 1;
%                 plot_end = 650;
%             else
%                 [val_max ind_max]=max(data_cond(:,i,c));
%                 plot_start = ind_max-300;
%                 plot_end = size(data_cond,1);
%             end
%             data_zoom = data_cond(plot_start:plot_end,i,c);
            
%             subplot (2,5,i); plot(data_zoom);
            fig_str = [fig_name, cond_name{c}, ' trial ', int2str(i)];
            title(fig_str);
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            set (gca,'Ydir','reverse')
            ylim([-15 15]); xlim([0 size(data_zoom,1)]);
            hold on;
        end
        %print(h,fig_str,'-dpng');
        pause;
        close
    end
end

 for c = 1:size(data_cond,4) % cond - single, licix3
        for m = 1:size(data_cond,3) % muscles
            h=figure;
            for i = 1:size(data_cond,2) % trials per condition
                plot(time_axis,squeeze(data_cond(:,i,m,c,j))); % data, trials, muscle, cond, step
                %                 plot(squeeze(data_cond(:,i,m,c,j))); % data, trials, muscle, cond, step
                if m < 4 % FCR has different time windows than the other muscle
                    hold on; line([time_windows(1,c) time_windows(1,c)],[-4 4],'Color',[1 0 0]);
                    hold on; line([time_windows(2,c) time_windows(2,c)],[-4 4],'Color',[1 0 0]);
                else
                    hold on; line([time_windows_FCR(1,c) time_windows_FCR(1,c)],[-4 4],'Color',[1 0 0]);
                    hold on; line([time_windows_FCR(2,c) time_windows_FCR(2,c)],[-4 4],'Color',[1 0 0]);
                end
                fig_str=[cond_name{c},' in ',m_name{m},' at ',step_name{j}];
                title(['\fontsize{16}',fig_str])
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                set (gca,'Ydir','reverse')
                ylim([-12 4])
                hold on;
            end
%             cd(mlfig_folder); savefig(fig_str);
            cd(png_folder); print(h,fig_str,'-dpng','-r300');
            % mean of this muscle and condition
            mean_cond_m=mean(squeeze(data_cond(:,:,m,c,j)),2);
            h2=figure;
            plot(time_axis, mean_cond_m);
            if m < 4 % FCR has different time windows than the other muscle
                hold on; line([time_windows(1,c) time_windows(1,c)],[-4 4],'Color',[1 0 0]);
                hold on; line([time_windows(2,c) time_windows(2,c)],[-4 4],'Color',[1 0 0]);
            else
                hold on; line([time_windows_FCR(1,c) time_windows_FCR(1,c)],[-4 4],'Color',[1 0 0]);
                hold on; line([time_windows_FCR(2,c) time_windows_FCR(2,c)],[-4 4],'Color',[1 0 0]);
            end
            %             plot(mean_cond_m);
            fig_str2=['Mean of ',cond_name{c},' in ',m_name{m},' at ',step_name{j}];
            title(['\fontsize{16}',fig_str2])
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            set (gca,'Ydir','reverse')
            ylim([-12 4])
%             cd(mlfig_folder); savefig(fig_str2);
            cd(png_folder); print(h2,fig_str2,'-dpng','-r300');
%             pause
            close all
        end
    end