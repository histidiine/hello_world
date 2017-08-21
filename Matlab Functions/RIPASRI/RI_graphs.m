function RI_graphs

step_name= {'Baseline', 'Post'};
cond_name={'single pulse','-1ms ISI', '-0,5ms ISI', 'synchronous','0,5ms ISI', '1ms ISI',...
    '2ms ISI','5ms ISI','7,5ms ISI','10ms ISI','25ms ISI','50ms ISI','75ms ISI','100ms ISI','200ms ISI'};


%% BOXPLOTS

rep_boxplot=input('Do you want to generate a box plot of h-reflexes per cond? (y/n):  ');


% H_resp - 1) peak, 2) peak delay, 3) ref, 4) ref delay, 5) amplitude

if strcmp(rep_boxplot,'y')==1
    cd(fig_folder);
    if exist('Boxplot','dir') == 0; mkdir('Boxplot'); else end
    cd('Boxplot')
    count=1;
    up_lim = (squeeze(max(max(max(max(H_resp(:,:,:,5))))))+0.1);
    h=figure;
    for j = 1:size(FCR_cond,4)
        curr_data=squeeze(H_resp(:,:,j,5));
        single_data=squeeze(H_resp(:,1,j,5));
        subplot(2,1,j);
        line([0 15],[median(single_data) median(single_data)], 'Color','r')
        hold on;
        boxplot(curr_data, 'labels', cond_name);
        count=count+1;
        ylim([0 up_lim])
        %         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    end
    DateString = datestr(clock);
    save_fig=[DateString,'.fig'];
    save_png=[DateString,'.png'];
    print(h,save_png,'-dpng','-r300');
    savefig(save_fig);
end

% In pourcentage of single
curr_data_pourcent=zeros(size(curr_data,1), size(curr_data,2));
if strcmp(rep_boxplot,'y')==1
    cd(fig_folder);
    if exist('Boxplot','dir') == 0; mkdir('Boxplot'); else end
    cd('Boxplot')
    h=figure;
    for j = 1:size(FCR_cond,4)
        curr_data=squeeze(H_resp(:,:,j,5));
        curr_data(curr_data==0)=0.001;
        single_data=squeeze(H_resp(:,1,j,5));
        single_data(single_data==0)=0.001; % to avoid division by 0
        subplot(2,1,j);
        line([0 15.5],[100 100], 'Color','r')
        hold on;
        for count = 1:size(curr_data,2) % cond
            curr_data_pourcent(:,count) = curr_data(:,count)./median(single_data)*100;
        end
        boxplot(curr_data_pourcent, 'labels', cond_name);
        for x_pt = 1:size(curr_data_pourcent,2)
            hold on; text(x_pt-0.1,105,int2str(median(curr_data_pourcent(:,x_pt))))
        end
        count=count+1;
        ylim([0 350])
        %         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    end
    DateString = datestr(clock);
    save_fig=[DateString,'_pourcent.fig'];
    save_png=[DateString,'_pourcent.png'];
    print(h,save_png,'-dpng','-r300');
    savefig(save_fig);
end

%% ALL COMPARE

cd(fig_folder)
if exist('PLOT','dir') == 0 mkdir('PLOT'); else end
plot_folder=[fig_folder,'/PLOT'];
cd(plot_folder)
if exist('ALL','dir') == 0 mkdir('ALL'); else end
% mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
% date2_folder = [plot_folder, '/ALL/', DateString];
% cd(date2_folder)

rep_big=input('Do you want to see h-reflexes big picture for Compare ALL?');
if strcmp(rep_big,'y')==1
    rep_save=input('Do you want to save the figures for "Compare ALL - h-reflex"? (y/n):   ');
    cd(plot_folder);cd('ALL')
    count=1; DateString = datestr(clock);
    for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
        h=figure;
        for c = 1:size(FCR_cond,3) % COND LOOP
            if c==1
                sub_start = sub_delay+200-2; sub_end = sub_delay2+200;
                plot_start = sub_start-5; plot_end = sub_end+5;
            else if  c < 5
                    sub_start = sub_delay+10-2; sub_end = sub_delay2+10;
                    plot_start = sub_start-5; plot_end = sub_end+5;
                else
                    sub_start = sub_delay+205-2; sub_end = sub_delay2+205;
                    plot_start = sub_start-5; plot_end = sub_end+5;
                end
            end
            
            duration_ms=points(j)/s_rate(j)*1000;
            step=duration_ms/points(j);
            time_axis=0:step:duration_ms-step; time_axis=time_axis';
            
            sub_start_pt = sub_start*points(j)/duration_ms;
            sub_end_pt = sub_end*points(j)/duration_ms;
            sub_single_start = sub_delay+200-2; sub_single_end = sub_delay2+200;
            sub_single_start_pt = sub_single_start*points(j)/duration_ms;
            sub_single_end_pt = sub_single_end*points(j)/duration_ms;
            
            plot_single_start = 190; plot_single_end = 245;
            
            plot_start_pt = round(plot_start*points(j)/duration_ms);
            plot_end_pt = round(plot_end*points(j)/duration_ms);
            
            %____________MEAAAAAAAN____________
            mean_FCR_cond = mean(squeeze(FCR_cond(:,:,c,j)),2);
            subplot(2,15,count)
            plot(time_axis,mean_FCR_cond);
            
            % Personal H-r delays from Viking
            hold on; line([sub_start-1 sub_start-1],[mean_FCR_cond(sub_start_pt)-2 mean_FCR_cond(sub_start_pt)+2],'Color',[0 1 0]);
            hold on; line([sub_end+1 sub_end+1],[mean_FCR_cond(sub_end_pt)-2 mean_FCR_cond(sub_end_pt)+2],'Color',[0 1 0]);
            % Detected MAX
            [mean_max, mean_delay] = min(mean_FCR_cond(sub_start_pt:sub_end_pt,:));
            h_peak_disp = round((mean_delay+sub_start_pt)*duration_ms/points(j));
            hold on; line([h_peak_disp h_peak_disp],[mean_max-0.5 mean_max+0.5],'Color',[1 0 1], 'Marker','.');
            % Detected MIN
            [mean_ref, mean_ref_delay] = max(mean_FCR_cond(sub_start_pt+mean_delay-1:sub_end_pt,:));
            h_ref_disp = round((mean_ref_delay+sub_start_pt+mean_delay)*duration_ms/points(j));
            hold on; line([h_ref_disp h_ref_disp],[mean_ref-0.5 mean_ref+0.5],'Color',[0 0 1], 'Marker','.');
            
            h_amp=abs(mean_max-mean_ref);
            peak_str = sprintf('%.3f',h_amp);
            hold on; text(sub_start,mean_ref-2,peak_str)
            
            mid_subD = mean_FCR_cond((sub_end_pt+sub_start_pt)/2);
            %         if mid_subD < 0
            %             ylim_d=mid_subD-5;
            %             ylim_u=mid_subD+5;
            %         else
            %             ylim_d=mid_subD-6;
            %             ylim_u=mid_subD+4;
            %         end
            %         ylim([ylim_d ylim_u])
            
            ylim_down = squeeze(min(min(min(FCR_cond(plot_start_pt:plot_end_pt,:,c,j)))));
            ylim_up = squeeze(max(max(max(FCR_cond(plot_start_pt:plot_end_pt,:,c,j)))));
            
            ylim([mean_max-1 mean_max+1]);
            %         ylim([ylim_down ylim_up])
            %         ylim([-1.5 1.5])
            xlim([plot_start plot_end]);
            set (gca,'Ydir','reverse')
            
            %____________ALL TRIALS - CURRENT COND____________
            for i=1:size(FCR_cond,2) % TRIAL LOOP
                % PLOT
                subplot(2,15,count+15)
                curr_data=FCR_cond(:,i,c,j);
                plot(time_axis,curr_data);
                ylim([mean_max-1 mean_max+1]);
                %             ylim([ylim_down ylim_up]);
                %             ylim([-1.5 1.5])
                xlim([plot_start plot_end]);
                set (gca,'Ydir','reverse')
                
                % Personal H-r delays from Viking
                hold on; line([sub_start-1 sub_start-1],[curr_data(sub_start_pt)-2 curr_data(sub_start_pt)+2],'Color',[0 1 0]);
                hold on; line([sub_end+1 sub_end+1],[curr_data(sub_end_pt)-2 curr_data(sub_end_pt)+2],'Color',[0 1 0]);
                
                hold on;
            end
            
            fig_str = [step_name{j}];
            %         save_str = [cond_name{c}, '_VS_SinglePulse_', step_name{j}];
            set(gcf,'NextPlot','add');
            axes;
            ht=title(fig_str);
            set(gca,'Visible','off');
            set(ht,'Visible','on');
            
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            hold on
            %         pause;
            count=count+1;
        end
        if strcmp(rep_save, 'y') == 1
            save_fig=[DateString, step_name{j},'_ALL_compare','.fig'];
            save_png=[DateString, step_name{j},'_ALL_compare','.png'];
            %         print(h,save_png,'-dpng','-r600');
            savefig(save_fig);
        end
        count=1;
    end
end
clear rep_save

% ZOOMED
rep_zoom=input('Do you want to see the zoomed on h-reflex ALL-Trials');
if strcmp(rep_zoom,'y')==1
    rep_save=input('Do you want to save the figures for "Compare ALL - h-reflex"? (y/n):   ');
    cd(plot_folder);cd('ALL')
    count=1; DateString = datestr(clock);
    for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
        h=figure;
        for c = 1:size(FCR_cond,3) % COND LOOP
            if c==1
                sub_start = sub_delay+200-2; sub_end = sub_delay2+200;
                plot_start = sub_start-5; plot_end = sub_end+5;
            else if  c < 5
                    sub_start = sub_delay+10-2; sub_end = sub_delay2+10;
                    plot_start = sub_start-5; plot_end = sub_end+5;
                else
                    sub_start = sub_delay+205-2; sub_end = sub_delay2+205;
                    plot_start = sub_start-5; plot_end = sub_end+5;
                end
            end
            
            duration_ms=points(j)/s_rate(j)*1000;
            step=duration_ms/points(j);
            time_axis=0:step:duration_ms-step; time_axis=time_axis';
            
            sub_start_pt = sub_start*points(j)/duration_ms;
            sub_end_pt = sub_end*points(j)/duration_ms;
            sub_single_start = sub_delay+200-2; sub_single_end = sub_delay2+200;
            sub_single_start_pt = sub_single_start*points(j)/duration_ms;
            sub_single_end_pt = sub_single_end*points(j)/duration_ms;
            
            plot_single_start = 190; plot_single_end = 245;
            
            plot_start_pt = round(plot_start*points(j)/duration_ms);
            plot_end_pt = round(plot_end*points(j)/duration_ms);
            
            %____________MEAAAAAAAN____________
            mean_FCR_cond = mean(squeeze(FCR_cond(:,:,c,j)),2);
            subplot(2,15,count)
            plot(time_axis,mean_FCR_cond);
            
            % Personal H-r delays from Viking
            hold on; line([sub_start-1 sub_start-1],[mean_FCR_cond(sub_start_pt)-2 mean_FCR_cond(sub_start_pt)+2],'Color',[0 1 0]);
            hold on; line([sub_end+1 sub_end+1],[mean_FCR_cond(sub_end_pt)-2 mean_FCR_cond(sub_end_pt)+2],'Color',[0 1 0]);
            % Detected MAX
            [mean_max, mean_delay] = min(mean_FCR_cond(sub_start_pt:sub_end_pt,:));
            h_peak_disp = round((mean_delay+sub_start_pt)*duration_ms/points(j));
            hold on; line([h_peak_disp h_peak_disp],[mean_max-0.5 mean_max+0.5],'Color',[1 0 1], 'Marker','.');
            % Detected MIN
            [mean_ref, mean_ref_delay] = max(mean_FCR_cond(sub_start_pt+mean_delay-1:sub_end_pt,:));
            h_ref_disp = round((mean_ref_delay+sub_start_pt+mean_delay)*duration_ms/points(j));
            hold on; line([h_ref_disp h_ref_disp],[mean_ref-0.5 mean_ref+0.5],'Color',[0 0 1], 'Marker','.');
            
            h_amp=abs(mean_max-mean_ref);
            peak_str = sprintf('%.3f',h_amp);
            hold on; text(sub_start,mean_ref-2,peak_str)
            
            mid_subD = mean_FCR_cond((sub_end_pt+sub_start_pt)/2);
            %         if mid_subD < 0
            %             ylim_d=mid_subD-5;
            %             ylim_u=mid_subD+5;
            %         else
            %             ylim_d=mid_subD-6;
            %             ylim_u=mid_subD+4;
            %         end
            %         ylim([ylim_d ylim_u])
            
            ylim_down = squeeze(min(min(min(FCR_cond(plot_start_pt:plot_end_pt,:,c,j)))));
            ylim_up = squeeze(max(max(max(FCR_cond(plot_start_pt:plot_end_pt,:,c,j)))));
            
            ylim([mean_max-1 mean_max+1]);
            %         ylim([ylim_down ylim_up])
            %         ylim([-1.5 1.5])
            xlim([plot_start plot_end]);
            set (gca,'Ydir','reverse')
            
            %____________ALL TRIALS - CURRENT COND____________
            for i=1:size(FCR_cond,2) % TRIAL LOOP
                % PLOT
                subplot(2,15,count+15)
                curr_data=FCR_cond(:,i,c,j);
                plot(time_axis,curr_data);
                ylim([mean_max-1 mean_max+1]);
                %             ylim([ylim_down ylim_up]);
                %             ylim([-1.5 1.5])
                xlim([plot_start plot_end]);
                set (gca,'Ydir','reverse')
                
                % Personal H-r delays from Viking
                hold on; line([sub_start-1 sub_start-1],[curr_data(sub_start_pt)-2 curr_data(sub_start_pt)+2],'Color',[0 1 0]);
                hold on; line([sub_end+1 sub_end+1],[curr_data(sub_end_pt)-2 curr_data(sub_end_pt)+2],'Color',[0 1 0]);
                
                hold on;
            end
            
            fig_str = [step_name{j}];
            %         save_str = [cond_name{c}, '_VS_SinglePulse_', step_name{j}];
            set(gcf,'NextPlot','add');
            axes;
            ht=title(fig_str);
            set(gca,'Visible','off');
            set(ht,'Visible','on');
            
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            hold on
            %         pause;
            count=count+1;
        end
        if strcmp(rep_save, 'y') == 1
            save_fig=[DateString, step_name{j},'_ALL_compare','.fig'];
            save_png=[DateString, step_name{j},'_ALL_compare','.png'];
            %         print(h,save_png,'-dpng','-r600');
            savefig(save_fig);
        end
        count=1;
    end
end

%% HISTOGRAMS - M AND H

figure;
for i = 1:size(H_ratio,2)
    plot(H_ratio(:,i))
    xlabel= cond_name;
    hold on;  line([0 15],[1 1],'Color',[1 0 0]);
    hold on
    %     pause
end

%% GRAPHIC CONDITIONS COMPARISON

cd(fig_folder)
if exist('PLOT','dir') == 0 mkdir('PLOT'); else end
plot_folder=[fig_folder,'/PLOT'];
cd(plot_folder)
if exist('Cond_VS_Single','dir') == 0 mkdir('Cond_VS_Single'); else end
cd('Cond_VS_Single')
DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
date2_folder = [plot_folder, '/Cond_VS_Single/', DateString];
cd(date2_folder)

rep_save=input('Do you want to save the figures? (y/n):   ');
for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
    for c = 1:size(FCR_cond,3)-1 % COND LOOP
        h=figure;
        if c < 5 && c~=1
            plot_start = 0; plot_end = 55;
            sub_start = sub_delay+10-2; sub_end = sub_delay2+10;
        else
            plot_start = 195; plot_end = 250;
            sub_start = sub_delay+205-2; sub_end = sub_delay2+205;
        end
        duration_ms=points(j)/s_rate(j)*1000;
        step=duration_ms/points(j);
        time_axis=0:step:duration_ms-step; time_axis=time_axis';
        
        sub_start_pt = sub_start*points(j)/duration_ms;
        sub_end_pt = sub_end*points(j)/duration_ms;
        sub_single_start = sub_delay+200-2; sub_single_end = sub_delay2+200;
        sub_single_start_pt = sub_single_start*points(j)/duration_ms;
        sub_single_end_pt = sub_single_end*points(j)/duration_ms;
        
        plot_single_start = 190; plot_single_end = 245;
        ylim_up=squeeze(abs(min(min(min(min(FCR_cond(:,:,c,j)),2)))));
        
        %____________ALL TRIALS - CURRENT COND____________
        for i=1:size(FCR_cond,2) % TRIAL LOOP
            % PLOT
            subplot(2,2,1)
            plot(time_axis,FCR_cond(:,i,c,j));
            ylim([-ylim_up ylim_up]); xlim([plot_start plot_end]);
            set (gca,'Ydir','reverse')
            
            % Personal H-r delays from Viking
            hold on; line([sub_start sub_start],[-4 4],'Color',[0 1 0]);
            hold on; line([sub_end sub_end],[-4 4],'Color',[0 1 0]);
            
            % H detected MAX
            h_peak_disp = H_resp(i,c,j,2)*duration_ms/points(j);
            hold on; line([h_peak_disp h_peak_disp],[(H_resp(i,c,j,1))-0.3 (H_resp(i,c,j,1))+0.3],'Color',[1 0 1], 'Marker','*');
            % H detected MIN
            h_ref_disp = H_resp(i,c,j,4)*duration_ms/points(j);
            hold on; line([h_ref_disp h_ref_disp],[(H_resp(i,c,j,3))-0.3 (H_resp(i,c,j,3))+0.3],'Color',[0 0 1], 'Marker','*');
            
            hold on;
        end
        %         pause;
        %____________ALL TRIALS - SINGLE PULSE___________
        for i=1:size(FCR_cond,2) % TRIAL LOOP
            % PLOT
            subplot(2,2,2)
            plot(time_axis,FCR_cond(:,i,15,j));
            ylim([-ylim_up ylim_up]); xlim([plot_single_start plot_single_end]);
            set (gca,'Ydir','reverse')
            
            % Personal H-r delays from Viking
            hold on; line([sub_single_start sub_single_start],[-4 4],'Color',[0 1 0]);
            hold on; line([sub_single_end sub_single_end],[-4 4],'Color',[0 1 0]);
            
            % H detected MAX
            h_single_disp = H_resp(i,15,j,2)*duration_ms/points(j);
            hold on; line([h_single_disp h_single_disp],[(H_resp(i,15,j,1))-0.3 (H_resp(i,15,j,1))+0.3],'Color',[1 0 1], 'Marker','*');
            % H detected MIN
            h_ref_single_disp = H_resp(i,15,j,4)*duration_ms/points(j);
            hold on; line([h_ref_single_disp h_ref_single_disp],[(H_resp(i,15,j,3))-0.3 (H_resp(i,15,j,3))+0.3],'Color',[0 0 1], 'Marker','*');
            
            hold on;
        end
        %         pause;
        %____________MEAAAAAAAN____________
        mean_FCR_cond = mean(squeeze(FCR_cond(:,:,c,j)),2);
        mean_FCR_single = mean(squeeze(FCR_cond(:,:,15,j)),2);
        
        % CURRENT CONDITION
        subplot(2,2,3)
        plot(time_axis,mean_FCR_cond);
        ylim([-ylim_up ylim_up]); xlim([plot_start plot_end]);
        set (gca,'Ydir','reverse')
        
        % Personal H-r delays from Viking
        hold on; line([sub_start sub_start],[-4 4],'Color',[0 1 0]);
        hold on; line([sub_end sub_end],[-4 4],'Color',[0 1 0]);
        
        % H detected MAX
        %         [mean_max, mean_delay] = min(mean_FCR_cond(sub_start_pt:sub_end_pt,:));
        %         h_peak_disp = mean_delay+sub_start_pt*duration_ms/points(j);
        %         hold on; line([h_peak_disp h_peak_disp],[mean_max-0.3 mean_max+0.3],'Color',[1 0 1], 'Marker','*');
        %         % H detected MIN
        %         [mean_ref, mean_ref_delay] = max(mean_FCR_cond(sub_start+mean_delay-1:sub_end,:));
        %         h_ref_disp = mean_ref_delay+sub_start*duration_ms/points(j);
        %         hold on; line([h_ref_disp h_ref_disp],[mean_ref-0.3 mean_ref+0.3],'Color',[0 0 1], 'Marker','*');
        % %         pause;
        
        % SINGLE PULSE
        subplot(2,2,4)
        plot(time_axis,mean_FCR_single);
        ylim([-ylim_up ylim_up]); xlim([plot_single_start plot_single_end]);
        set (gca,'Ydir','reverse')
        % Personal H-r delays from Viking
        hold on; line([sub_single_start sub_single_start],[-4 4],'Color',[0 1 0]);
        hold on; line([sub_single_end sub_single_end],[-4 4],'Color',[0 1 0]);
        
        % H detected MAX
        %         [mean_max_single, mean_delay_single] = min(mean_FCR_single(sub_single_start:sub_single_start,:));
        %         h_peak_single_disp = mean_delay+sub_single_start*duration_ms/points(j);
        %         hold on; line([h_peak_single_disp h_peak_single_disp],[mean_delay_single-0.3 mean_delay_single+0.3],'Color',[1 0 1], 'Marker','*');
        %         % H detected MIN
        %         [mean_ref_single, mean_ref_single_delay] = max(mean_FCR_single(sub_start+mean_delay-1:sub_end,:));
        %         h_ref_disp = mean_ref_single_delay+sub_single_start*duration_ms/points(j);
        %         hold on; line([h_ref_disp h_ref_disp],[mean_ref_single-0.3 mean_ref_single_delay+0.3],'Color',[0 0 1], 'Marker','*');
        %         pause;
        %title
        fig_str = [cond_name{c}, ' VS SinglePulse ', step_name{j}];
        %         save_str = [cond_name{c}, '_VS_SinglePulse_', step_name{j}];
        set(gcf,'NextPlot','add');
        axes;
        ht=title(fig_str);
        set(gca,'Visible','off');
        set(ht,'Visible','on');
        
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        hold on
        %         pause;
        
        if strcmp(rep_save, 'y') == 1
            save_fig=['Single_VS_',cond_name{c}, step_name{j},'.fig'];
            save_png=['Single_VS_',cond_name{c}, step_name{j},'.png'];
            print(h,save_png,'-dpng','-r600');
            savefig(save_fig);
        else
        end
        close;
    end
end



end