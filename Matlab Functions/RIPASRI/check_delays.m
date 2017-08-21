function [sub_delay2]= check_delays(RI_cond, sub_delays, data_param)

%% INI
step_name= {'Baseline', 'Post'};
cond_name={'single pulse','-1ms ISI', '-0,5ms ISI', 'synchronous','0,5ms ISI', '1ms ISI',...
    '2ms ISI','5ms ISI','7,5ms ISI','10ms ISI','25ms ISI','50ms ISI','75ms ISI','100ms ISI','200ms ISI'};

s_rate=data_param(1,:);
points=data_param(2,:);
sub_delay=sub_delays(1);
sub_delay2=sub_delays(2);

m=1; % for now I will concentrate on FCR
test=1; search=1; count = 1;

%% PROCESS
while search ~= 0
    duration_ms=points(1)/s_rate(1)*1000;
    step=duration_ms/points(1);
    time_axis=0:step:duration_ms-step; time_axis=time_axis';
    for j = 1:size(RI_cond,5) % STEP LOOP
        h=figure;
        for c = 1:size(RI_cond,4) % COND LOOP
            subplot(3,5,c)
            if c == 1
                plot_start = 199; plot_end = 240;
                sub_start = sub_delay+200; sub_end = sub_delay2+200;
            else if c < 5
                    plot_start = 9; plot_end = 50;
                    sub_start = sub_delay+10; sub_end = sub_delay2+10;
                else
                    plot_start = 204; plot_end = 245;
                    sub_start = sub_delay+205; sub_end = sub_delay2+205;
                end
            end
            
            sub_start_pt = round(sub_start*points(1)/duration_ms);
            sub_end_pt = round(sub_end*points(1)/duration_ms);
            plot_start_pt = round(plot_start*points(1)/duration_ms);
            plot_end_pt = round(plot_end*points(1)/duration_ms);
            ylim_down = squeeze(min(min(min(RI_cond(plot_start_pt+45:plot_end_pt,:,c,1)))))-0.2;
            ylim_up = squeeze(max(max(max(RI_cond(plot_start_pt+45:plot_end_pt,:,c,1)))))+0.2;
            
            for i=1:size(RI_cond,3) % TRIAL LOOP
                curr_data = squeeze(RI_cond(:,m,i,c,j)); % 5000 * 1
                plot(time_axis,curr_data);
                fig_str = [cond_name{c}, ' ', step_name{j}];
                title(fig_str);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                set (gca,'Ydir','reverse')
                ylim([ylim_down ylim_up]);
                xlim([plot_start plot_end]);
                hold on;
            end
            
            % Personal H-r delays from Viking
            if test==1
                sub_end_test = sub_end;
            else if c == 1
                    sub_end_test = rep_sub_end+190;
                else if c > 4
                        sub_end_test = rep_sub_end+195;
                    else
                        sub_end_test = rep_sub_end;
                    end
                end
            end
            hold on; line([sub_start sub_start],[ylim_down ylim_up],'Color',[0 1 0]);
            hold on; line([sub_end_test sub_end_test],[ylim_down ylim_up],'Color',[0 1 0])
            count = count +1;
        end
        
    end
    pause;
    rep_sub_end=input('Enter the value for the end of the H-reflex that would be most proper (if the value is good type 0):  ');
    
    if rep_sub_end ~= 0
        sub_end_test=rep_sub_end;
        save_str=['Check delays',step_name{j},DateString,'.png'];
    else
        sub_delay2 = sub_end_test-195-10;
        search=0;
        save_str=['Check delays_BEST_',step_name{j},DateString,'.png'];
    end
    
    print(h,save_str,'-dpng','-r600');
    close all
    test=test+1;
end
