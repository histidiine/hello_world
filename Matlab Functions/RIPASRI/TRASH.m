% TRASH

%% graph_gen


% curr_ratio=zeros(size(curr_mean));
% for step = 1:size(pre_curr_mean,2)
%     for cond = 1:size(pre_curr_mean,1)
%         if cond == 1 % Single Pulse
%             curr_ratio(cond,step) = curr_mean(cond,step)./curr_mean(1,1)*100;
%         else % LICI 100
%             temp = curr_mean(cond,step)./curr_mean(1,step)*100;% absolute response
%             curr_ratio(cond,step) = -(100 - temp); % inhibition
%         end
%     end
% end




% curr_ratio_mean_INHIB = zeros(size(curr_ratio_mean));




% for step = 1:size(pre_curr_mean,2)
%     for cond = 1 :size(pre_curr_mean,1)
%         if cond == 1
%             curr_ratio_mean_INHIB(cond,step) = curr_ratio_mean(cond,step) ;
%         else
%             curr_ratio_mean_INHIB(cond,step) = -(100-curr_ratio_mean(cond,step));
%         end
%     end
% end


%% RI_process - previously treated 


% if strcmp(rep_proc, 'y')==1
%     rep_load = input('Do you want to load the data already processed?(y/n):  ');
%     if strcmp(rep_load, 'y')== 1
%         cd(ML_data_folder)
%         while ~exist('MandH.mat','file')==1
%             disp('MandH.cannot be found, please check in the folder')
%             pause
%         end
%         load MandH.mat
%         while ~exist('FCR_cond.mat','file')==1
%             disp('FCR_cond.cannot be found, please check in the folder')
%             pause
%         end
%         load FCR_cond.mat
%         rep_bd=input('Did you process bad detections? (y/n):  ');
%         if strcmp(rep_bd,'y')==1
%             while ~exist('Bad_Detect.mat','file')==1
%                 disp('Bad_Detect.cannot be found, please check in the folder')
%                 pause
%             end
%             load Bad_Detect.mat
%         end
%     end
% end

%% ILLUSTRATION OF RESPONSES

% graph_rep=input('Do you want to take a look at the amplitudes extracted (y/n) ?:  ');
% 
% if strcmpt(graph_rep,'y')== 1
%     
% end


%________________________END____________________________%

%%
%-------------------------------------------
% -------------- CEMETERY ------------------
%-------------------------------------------


% %% HISTOGRAMS OF AMPLITUDE
% 
% cd(fig_folder)
% if exist('Histograms','dir') == 0 mkdir('Histograms'); else end
% hist_folder=[fig_folder,'/Histograms'];
% cd(hist_folder)
% DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
% hist_date_folder = [hist_folder, '/', DateString];
% cd(hist_date_folder)
% 
% % SELECTED GOOD RECS
% up_lim = round(squeeze(max(max(max(max(resp_goodrecs(:,:,:,:,5))))))+0.5);
% figure; count=1;
% for m = 1:size(data_cond,3)    
%     for j = 1:3
%         curr_data=squeeze(resp_goodrecs(:,m,:,j,5));
%         subplot(size(data_cond,3),3,count);
%         boxplot(curr_data, 'labels', cond_name);
%         count=count+1;
%         ylim([0 up_lim])
% %         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%     end
% end
% 


% mean_data_cond = squeeze(mean(data_cond,4));
% mean_resp = squeeze(mean(resp,1)); % resp - trials, muscle, cond, step, param (peak, delay, ref, delayref, amp)

% for m = 1:size(data_cond,3)
%     for j = 1:3
%         for c = 1:size(data_cond,4)-1
%             mean_single = mean_data_cond(:,m,1,j);
%             curr_mean = mean_data_cond(:,m,c,j);
%             
%         end
%     end
% end
% 
% up_lim = round(squeeze(max(max(max(max(resp(:,:,:,:,5))))))+0.5);
% figure; count=1;
% for m = 1:size(data_cond,3)    
%     for j = 1:3
%         curr_data=squeeze(resp(:,m,:,j,5));
%         subplot(size(data_cond,3),3,count);
%         boxplot(curr_data, 'labels', cond_name);
%         count=count+1;
%         ylim([0 up_lim])
% %         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%     end
% end


% PLOT MEAN CONDITION x STEP
% 

% cd(fig_folder)
% if exist('MEAN', 'dir')==0; mkdir('MEAN'); end
% cd('MEAN')
% DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
% cd(DateString)
% 
% for m = 1:size(data_cond,3) % muscles
%     h=figure; count=1;
%     for j = 1:size(data_cond,5) % step - baseline, t1, t2
%         duration_ms=points(j)/s_rate(j)*1000;
%         step=duration_ms/points(j);
%         time_axis=0:step:duration_ms-step; time_axis=time_axis';
%         for c = 1:size(data_cond,4) % cond - single, licix3
%             curr_data = mean(squeeze(data_cond(:,:,m,c,j)),2);
%             subplot(3,4,count)
%             plot(time_axis,curr_data); % data, trials, muscle, cond, step
%             
%             %___________TIME WINDOWS VERTICAL LINES__________
%             if m < 4 % FCR has different time windows than the other muscle
%                 curr_time_win = time_windows(:,c);
%             else
%                 curr_time_win = time_windows_FCR(:,c);
%             end
%             hold on; line([curr_time_win(1) curr_time_win(1)],[-4 4],'Color',[1 0 0]);
%             hold on; line([curr_time_win(2) curr_time_win(2)],[-4 4],'Color',[1 0 0]);
%             
%             %___________GRAPH PARAM__________
%             xlim([curr_time_win(1)-20 curr_time_win(2)+50])
%             ylim([-8 8])
%             fig_str=[cond_name{c}, ' - ',m_name{m},' - ',step_name{j}];
%             title(fig_str,'FontSize', 8)
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%             set (gca,'Ydir','reverse')
%             count=count+1;
%             hold on
%         end
%     end
%     print_str=[cond_name{c}, ' - ',m_name{m},' - ',step_name{j}];
%     print(h,print_str,'-dpng','-r600');
%     save_fig=[print_str,'.fig'];
%     savefig(save_fig);
%     close
% end

% CREATE AND SAVE GRAPHS
% check that the data were properly transfered (since I cannot check 4D)
% cd(log_folder)
% 
% cd(fig_folder)
% if exist('PLOT','dir') == 0 mkdir('PLOT'); else end
% plot_folder=[fig_folder,'/PLOT'];
% cd(plot_folder)
% DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
% date_folder = [plot_folder, '/', DateString];
% cd(date_folder)
% mkdir('Matlab Fig'); mkdir('PNG');
% mlfig_folder=[date_folder,'/Matlab Fig']; png_folder=[date_folder,'/PNG'];
% for j = 1:size(data_cond,5) % step - baseline, t1, t2
%     duration_ms=points(j)/s_rate(j)*1000;
%     step=duration_ms/points(j);
%     time_axis=0:step:duration_ms-step; time_axis=time_axis';
%     for c = 1:size(data_cond,4) % cond - single, licix3
%         for m = 1:size(data_cond,3) % muscles
%             h=figure;
%             for i = 1:size(data_cond,2) % trials per condition
%                 plot(time_axis,squeeze(data_cond(:,i,m,c,j))); % data, trials, muscle, cond, step
%                 %                 plot(squeeze(data_cond(:,i,m,c,j))); % data, trials, muscle, cond, step
%                 if m < 4 % FCR has different time windows than the other muscle
%                     hold on; line([time_windows(1,c) time_windows(1,c)],[-4 4],'Color',[1 0 0]);
%                     hold on; line([time_windows(2,c) time_windows(2,c)],[-4 4],'Color',[1 0 0]);
%                 else
%                     hold on; line([time_windows_FCR(1,c) time_windows_FCR(1,c)],[-4 4],'Color',[1 0 0]);
%                     hold on; line([time_windows_FCR(2,c) time_windows_FCR(2,c)],[-4 4],'Color',[1 0 0]);
%                 end
%                 fig_str=[cond_name{c},' in ',m_name{m},' at ',step_name{j}];
%                 title(['\fontsize{16}',fig_str])
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                 set (gca,'Ydir','reverse')
%                 ylim([-12 4])
%                 hold on;
%             end
% %             cd(mlfig_folder); savefig(fig_str);
%             cd(png_folder); print(h,fig_str,'-dpng','-r300');
%             % mean of this muscle and condition
%             mean_cond_m=mean(squeeze(data_cond(:,:,m,c,j)),2);
%             h2=figure;
%             plot(time_axis, mean_cond_m);
%             if m < 4 % FCR has different time windows than the other muscle
%                 hold on; line([time_windows(1,c) time_windows(1,c)],[-4 4],'Color',[1 0 0]);
%                 hold on; line([time_windows(2,c) time_windows(2,c)],[-4 4],'Color',[1 0 0]);
%             else
%                 hold on; line([time_windows_FCR(1,c) time_windows_FCR(1,c)],[-4 4],'Color',[1 0 0]);
%                 hold on; line([time_windows_FCR(2,c) time_windows_FCR(2,c)],[-4 4],'Color',[1 0 0]);
%             end
%             %             plot(mean_cond_m);
%             fig_str2=['Mean of ',cond_name{c},' in ',m_name{m},' at ',step_name{j}];
%             title(['\fontsize{16}',fig_str2])
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%             set (gca,'Ydir','reverse')
%             ylim([-12 4])
% %             cd(mlfig_folder); savefig(fig_str2);
%             cd(png_folder); print(h2,fig_str2,'-dpng','-r300');
%             pause
%             close all
%         end
%     end
% end
% clear fig_str fig_str2 

% GRAPHS TO COMPARE STEPS AND CONDITIONS

% % generate random data - 32 channels
% t = linspace(0,2,points);
% sig = rand(10,points);
% 
% % calculate shift
% mi = min(sig,[],2);
% ma = max(sig,[],2);
% shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]);
% shift = repmat(shift,10,points);
% 
% %plot 'EMG' data of blocks
% for j=1:14 % 1 plot avec 10 trials pour chaque condition
%     figure;
%     h1=plot(t,FCR_cond(:,:,j)'+shift,'Color', [0 0 .5]);
%     ylim([-11 11])
%     %set(h1, 'Ydir', 'reverse')
% %     str = ['Block ', int2str(j)]
% %     print(h1,str, 'd-png')
% end

% if cond_num < 4 % if not all datasets have to be processed
%     cond_code = dataset(cond_num);
% else 
%     for i = 1:3 % if all datasets are to be processed
%        cond_code = dataset(i);
%     end
% end

% cd(fig_folder)
% if exist('PLOT3','dir')==0 mkdir('PLOT3'); else end
% plot3_folder=[fig_folder,'/PLOT3']; cd(plot3_folder)
% 
% duration_ms=points(j)/s_rate(j)*1000;
% step=duration_ms/points(j);
% time_axis=0:step:duration_ms-step; time_axis=time_axis';
% h=figure;
% for j = 1:size(data_cond,5) % step - baseline, t1, t2
%     mean_cond_m=squeeze(mean(squeeze(data_cond(:,:,2,:,j)),2)); % just FDI
%     subplot(1,3,j)
%     plot3(time_axis,mean_cond_m)
%     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%     set (gca,'Ydir','reverse')
%     fig_str2=['Mean of ',cond_name{c},' in ',m_name{m},' at ',step_name{j}];
%     title(['\fontsize{16}',fig_str2])
%     cd(mlfig_folder); savefig(fig_str);
%     cd(png_folder); print(h,fig_str,'-dpng');  
% end

%% RI_PROCESS

%% ---------------------------------------------------------------------
%                              CEMETERY
% ----------------------------------------------------------------------


% %% RATIO OVER SINGLE - H/M RATIO
% H_ratio = zeros(14,2); % Conditions, Steps
% for c = 1:size(FCR_cond,3)-1 % COND LOOP
%     for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
%         HoverM_resp = mean(H_resp(:,c,j,1),1)./mean(M_resp(:,c,j,1),1);
%         HoverM_base = mean(H_resp(:,15,j,1),1)./mean(M_resp(:,15,j,1),1);
%         H_ratio(c,j) = HoverM_resp./HoverM_base;
%     end
% end
%
% figure;
% for i = 1:size(H_ratio,2)
%     plot(H_ratio(:,i))
% %     xlabel= cond_name;
%     hold on;  line([0 15],[1 1],'Color',[1 0 0]);
%     hold on
% %     pause
% end
%
%
%
% % MEAN VALUES ACROSS CONDITIONS --> EXPORT TO EXCEL
% curr_data = zeros(15,6); % M before, H before, ratio before, M post, H post, ratio post
% for c = 1:size(FCR_cond,3) % COND LOOP
%     for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
%         curr_data(c,3*j-2) = mean(M_resp(:,c,j,1));
%         curr_data(c,3*j-1) = mean(H_resp(:,c,j,1));
%         curr_data(c,3*j) = curr_data(c,3*j-1) / curr_data(c,3*j-2);
%     end
%
% end
% filename = [sub_ID,'.csv'];
% csvwrite(filename,curr_data);
%
% % MEAN - BAR
% curr_data = zeros(15,2);
% figure;
% % for c = 1:size(FCR_cond,3) % COND LOOP
%     for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
%         curr_data(:,1) = mean(M_resp(:,:,j,1)); %pre
%         curr_data(:,2) = mean(H_resp(:,:,j,1)); %pre
%         subplot(2,1,j)
%         bar(curr_data)
%         ylim([0 8])
%         colormap(jet)
%         hold on
%     end
% % end
%
% figure;
% curr_data = zeros(15,1);
% % for c = 1:size(FCR_cond,3) % COND LOOP
%     for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
%         curr_data(:,1) = mean(H_resp(:,:,j,1))./mean(M_resp(:,:,j,1));
%         subplot(2,1,j)
%         bar(curr_data)
%         hold on
%     end
% % end
%
% % All trials
% curr_data = zeros(10,6);
% for c = 1:size(FCR_cond,3) % COND LOOP
%     %     filename = [sub_ID,'.csv']; count=2;
%     %     while exist(filename, 'file')==2
%     %         filename = [sub_ID,int2str(count),'.csv'];
%     %         count=count+1;
%     %     end
%     filename = [cond_name{c},'.csv'];
% %     for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
%         %         figure;
%         %         for i=1:size(FCR_cond,2) % TRIAL LOOP
%         curr_data(:,1) = M_resp(:,c,1,1); %pre
%         curr_data(:,2) = M_resp(:,c,2,1); %post
%         curr_data(:,3) = H_resp(:,c,1,1); %pre
%         curr_data(:,4) = H_resp(:,c,2,1); %post
%         curr_data(:,5) = curr_data(:,3)./curr_data(:,1); %pre
%         curr_data(:,6) = curr_data(:,4)./curr_data(:,2); %post
%         %         hist(curr_data, length(curr_data(:,2)))
%         %         pause;
%         %         end
%         %         sheet = int2str(c);
%         %         xlRange = ['B',int2str(10*j)];
%         %         xlswrite(filename,curr_data,sheet,xlRange)
%         csvwrite(filename,curr_data);
% %     end
% end
%
%
% %% Get a look at the data - plot all by condition
%
% % generate random data - 32 channels
% t = linspace(0,2,points);
% sig = rand(10,points);
%
% % calculate shift
% mi = min(sig,[],2);
% ma = max(sig,[],2);
% shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]);
% shift = repmat(shift,10,points);
%
% %plot 'EMG' data of blocks
% for j=1:14 % 1 plot avec 10 trials pour chaque condition
%     figure;
%     h1=plot(t,FCR_cond(:,:,j)'+shift,'Color', [0 0 .5]);
%     ylim([-11 11])
%     %set(h1, 'Ydir', 'reverse')
% %     str = ['Block ', int2str(j)]
% %     print(h1,str, 'd-png')
% end
%
% %% !!!EXTRA FIGURES!!! - ONE FIGURE PER CONDITION  - ALL TRIALS HOLD ON + MEAN - H RESPONSE
%
% cd(fig_folder)
% if exist('PLOT','dir') == 0; mkdir('PLOT'); else end
% plot_folder=[fig_folder,'/PLOT'];
% cd(plot_folder)
% if exist('DETECT','dir') == 0; mkdir('DETECT'); else end
% plot_detect_folder=[fig_folder,'/PLOT/DETECT'];
% cd(plot_detect_folder)
% DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
% date2_folder = [plot_detect_folder, '/', DateString];
% cd(date2_folder)
%
% for c = 1:size(FCR_cond,3) % COND LOOP
%     h=figure;
%     for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
%
%         %____________ALL TRIALS____________
%         for i=1:size(FCR_cond,2) % TRIAL LOOP
%             duration_ms=points(j)/s_rate(j)*1000;
%             step=duration_ms/points(j);
%             time_axis=0:step:duration_ms-step; time_axis=time_axis';
%             mean_FCR_cond = mean(squeeze(FCR_cond(:,:,c,j)),2);
%             if c < 4
%                 plot_start = 0; plot_end = 50;
%                 M_start = 2+10; M_end = 12+10;
%                 bar_start = 26; bar_end = 34;% delay + 10 ms preceeding the stimulus
%                 sub_start = sub_delay+10-2; sub_end = sub_delay2+10;
% %                 stim_time= 10;
%             else if c == 15
%                     plot_start = 190; plot_end = 240;
%                     M_start = 2+200; M_end = 12+200;
%                     bar_start = 216; bar_end = 224;
%                     sub_start = sub_delay+200-2; sub_end = sub_delay2+200;
% %                     stim_time= 200;
%                 else
%                     plot_start = 195; plot_end = 245;
%                     M_start = 2+205; M_end = 12+205;
%                     bar_start = 221; bar_end = 229;
%                     sub_start = sub_delay+205-2; sub_end = sub_delay2+205;
% %                     stim_time= 205;
%                 end
%             end
%             % CONVERT DELAYS IN TIMEPOINTS
%             M_start_pt = M_start*points(j)/duration_ms;
%             M_end_pt = M_end*points(j)/duration_ms;
%             sub_start_pt = sub_start*points(j)/duration_ms;
%             sub_end_pt = sub_end*points(j)/duration_ms;
%
%             % FIND M and H
%             M_onset = FCR_cond(M_start_pt,i,c,j);
%             [M_peak, m_peak_delay] = min(FCR_cond(M_start_pt:M_end_pt,i,c,j));
%             M_amp = abs(M_peak - M_onset);
%             M_resp(i,c,j,1) = M_amp;
%
%             H_onset = FCR_cond(sub_start_pt,i,c,j);
%             [H_peak, h_peak_delay] = min(FCR_cond(sub_start_pt:sub_end_pt,i,c,j));
%             H_amp = abs(H_peak - H_onset);
%             H_resp(i,c,j,1) = H_amp;
%
%             % PLOT
%             subplot(2,2,2*j-1)
%             plot(time_axis,FCR_cond(:,i,c,j));
%             % M-response delays
%             hold on; line([M_start M_start],[-4 4],'Color',[0 0 1]);
%             hold on; line([M_end M_end],[-4 4],'Color',[0 0 1]);
%             % M detected MAX
%             m_peak_disp = (m_peak_delay+M_start_pt)*duration_ms/points(j);
%             hold on; line([m_peak_disp m_peak_disp],[M_peak-1 M_peak+1],'Color',[1 0 1], 'Marker','*');
%
%             % Peter Fuhr H-r delays
%             hold on; line([bar_start bar_start],[-4 4],'Color',[1 0 0]);
%             hold on; line([bar_end bar_end],[-4 4],'Color',[1 0 0]);
%             % Personal H-r delays from Viking
%             hold on; line([sub_start sub_start],[-4 4],'Color',[0 1 0]);
%             hold on; line([sub_end sub_end],[-4 4],'Color',[0 1 0]);
%             % H detected MAX
%             h_peak_disp = (h_peak_delay+sub_start_pt)*duration_ms/points(j);
%             hold on; line([h_peak_disp h_peak_disp],[H_peak-1 H_peak+1],'Color',[1 0 1], 'Marker','*');
%
%             set (gca,'Ydir','reverse')
%             fig_str = [cond_name{c}, ' ALL ', step_name{j}];
%             title(fig_str);
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%             ylim([-8 8]); xlim([plot_start plot_end]);
%             hold on;
%         end
%
%         %____________MEAAAAAAAN____________
%         M_onset = mean_FCR_cond(M_start_pt,:);
%         [M_peak, m_peak_delay] = min(mean_FCR_cond(M_start_pt:M_end_pt,:));
%
%         H_onset = mean_FCR_cond(sub_start_pt,:);
%         [H_peak, h_peak_delay] = min(mean_FCR_cond(sub_start_pt:sub_end_pt,:));
%
%         subplot(2,2,2*j)
%         plot(time_axis,mean_FCR_cond);
%         hold on;
%         % M detected MAX
%         m_peak_disp = (m_peak_delay+M_start_pt)*duration_ms/points(j);
%         hold on; line([m_peak_disp m_peak_disp],[M_peak-1 M_peak+1],'Color',[1 0 1], 'Marker','*');
%         % H detected MAX
%         h_peak_disp = (h_peak_delay+sub_start_pt)*duration_ms/points(j);
%         hold on; line([h_peak_disp h_peak_disp],[H_peak-1 H_peak+1],'Color',[1 0 1], 'Marker','*');
%
%         set (gca,'Ydir','reverse')
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%         ylim([-8 8]); xlim([plot_start plot_end]);
%         hold on
%                 pause;
%     end
% %     save_fig=['DETECT','_',cond_name{c},'_', ' MEAN','.fig'];
% %     save_png=['DETECT','_',cond_name{c},'_', ' MEAN,''.png'];
% %     print(h,save_png,'-dpng','-r300');
% %     savefig(save_fig);
%     close all
% end

% PLOT BY CONDITION AND STEP SEPARATELY - EDC
% for j = 1:size(EDC_cond,4) % STEP LOOP - baseline, post
%     duration_ms=points(j)/s_rate(j)*1000;
%     step=duration_ms/points(j);
%     time_axis=0:step:duration_ms-step; time_axis=time_axis';
%     ylim_up_EDC = squeeze(abs(min(min(min(min(EDC_cond))))))+0.5;
%     for c = 1:size(EDC_cond,3) % COND LOOP
%         h=figure;
%         for i=1:size(EDC_cond,2) % TRIAL LOOP
%             plot(time_axis,squeeze(EDC_cond(:,i,c,j)));
% %             hold on; line([time_windows(1,c) time_windows(1,c)],[-4 4],'Color',[1 0 0]);
% %             hold on; line([time_windows(2,c) time_windows(2,c)],[-4 4],'Color',[1 0 0]);
%             fig_str = [fig_name, ' ', cond_name{c}, ' ', step_name{j}];
%             title(fig_str);
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%             set (gca,'Ydir','reverse')
%             ylim([-ylim_up_EDC ylim_up_EDC]); xlim([0 duration_ms]);
%             hold on;
%         end
%         save_str=['EDC_',cond_name{c},'_', step_name{j},'.fig'];
%         print(h,fig_str,'-dpng','-r300');
%         savefig(save_str);
%         close
%     end
% end

% [~, ~, test1matlab] = xlsread([nom_dossier,'\test_1_matlab.xlsx'],'Feuil1');
% test1matlab(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),test1matlab)) = {''};
% log=test1matlab; clear test1matlab

% %% In miliseconds
% % Organise data by condition (ISI)
%
%
% count1=1;count2=1;count3=1;count4=1;count5=1;count6=1;count7=1;count8=1;count9=1;count10=1;count11=1;count12=1;count13=1;count14=1;count15=1;
% for i = 1:num
%     if log(i,7)==10
%         if log(i,9)==0
%             FCR_cond_ms(:,count3,3) =  FCR_data_ms(:,i);
%             count3=count3+1;
%         else if log(i,9)==-0.5
%                 FCR_cond_ms(:,count2,2) =  FCR_data_ms(:,i);
%                 count2=count2+1;
%             else if log(i,9)==-1
%                     FCR_cond_ms(:,count1,1) =  FCR_data_ms(:,i);
%                     count1=count1+1;
%                 end
%             end
%         end
%     else if log(i,7)==200
%             FCR_cond_ms(:,count15,15) =  FCR_data_ms(:,i);
%             count15=count15+1;
%         else if log(i,9)==-0.5
%                 FCR_cond_ms(:,count4,4) =  FCR_data_ms(:,i);
%                 count4=count4+1;
%             else if log(i,9)==-1
%                     FCR_cond_ms(:,count5,5) =  FCR_data_ms(:,i);
%                     count5=count5+1;
%                 else if log(i,9)==-2
%                         FCR_cond_ms(:,count6,6) =  FCR_data_ms(:,i);
%                         count6=count6+1;
%                     else if log(i,9)==-5
%                             FCR_cond_ms(:,count7,7) =  FCR_data_ms(:,i);
%                             count7=count7+1;
%                         else if log(i,9)==-7.5
%                                 FCR_cond_ms(:,count8,8) =  FCR_data_ms(:,i);
%                                 count8=count8+1;
%                             else if log(i,9)==-10
%                                     FCR_cond_ms(:,count9,9) =  FCR_data_ms(:,i);
%                                     count9=count9+1;
%                                 else if log(i,9)==-25
%                                         FCR_cond_ms(:,count10,10) =  FCR_data_ms(:,i);
%                                         count10=count10+1;
%                                     else if log(i,9)==-50
%                                             FCR_cond_ms(:,count11,11) =  FCR_data_ms(:,i);
%                                             count11=count11+1;
%                                         else if log(i,9)==-75
%                                                 FCR_cond_ms(:,count12,12) =  FCR_data_ms(:,i);
%                                                 count12=count12+1;
%                                             else if log(i,9)==-100
%                                                     FCR_cond_ms(:,count13,13) =  FCR_data_ms(:,i);
%                                                     count13=count13+1;
%                                                 else if log(i,9)==-200
%                                                         FCR_cond_ms(:,count14,14) =  FCR_data_ms(:,i);
%                                                         count14=count14+1;
%
%                                                     end
%                                                 end
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%
%         end
%     end
% end
%     clear count1 count2 count3 count4 count5 count6 count7 count8 count9 count10 count11 count12 count13 count14 count15
%
% FCR_cond_ms = FCR_cond_ms.*-1; % To get the data with minus up and positive down
% for j = 1:15
%    h=figure;
%    fig_str = ['Condition number ', int2str(j), ' in ms'];
%    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%    for i=1:10
%        subplot (5,2,i); plot(FCR_cond_ms(:,i,j));
%        ylim([-13 13]); xlim([0 time_ms]);
%    end
%    title(fig_str);
%    print(h,fig_str,'-dpng');
%    close
% end
%
% %% In original timepoints
%
% FCR_cond = zeros(points,10,15);
% count1=1;count2=1;count3=1;count4=1;count5=1;count6=1;count7=1;count8=1;count9=1;count10=1;count11=1;count12=1;count13=1;count14=1;count15=1;
% for i = 1:num
%     if log(i,7)==10
%         if log(i,9)==0
%             FCR_cond(:,count3,3) =  FCR_data(:,i);
%             count3=count3+1;
%         else if log(i,9)==-0.5
%                 FCR_cond(:,count2,2) =  FCR_data(:,i);
%                 count2=count2+1;
%             else if log(i,9)==-1
%                     FCR_cond(:,count1,1) =  FCR_data(:,i);
%                     count1=count1+1;
%                 end
%             end
%         end
%     else if log(i,7)==200
%             FCR_cond_ms(:,count15,15) =  FCR_data_ms(:,i);
%             count15=count15+1;
%         else
%             if log(i,9)==-0.5
%                 FCR_cond(:,count4,4) =  FCR_data(:,i);
%                 count4=count4+1;
%             else if log(i,9)==-1
%                     FCR_cond(:,count5,5) =  FCR_data(:,i);
%                     count5=count5+1;
%                 else if log(i,9)==-2
%                         FCR_cond(:,count6,6) =  FCR_data(:,i);
%                         count6=count6+1;
%                     else if log(i,9)==-5
%                             FCR_cond(:,count7,7) =  FCR_data(:,i);
%                             count7=count7+1;
%                         else if log(i,9)==-7.5
%                                 FCR_cond(:,count8,8) =  FCR_data(:,i);
%                                 count8=count8+1;
%                             else if log(i,9)==-10
%                                     FCR_cond(:,count9,9) =  FCR_data(:,i);
%                                     count9=count9+1;
%                                 else if log(i,9)==-25
%                                         FCR_cond(:,count10,10) =  FCR_data(:,i);
%                                         count10=count10+1;
%                                     else if log(i,9)==-50
%                                             FCR_cond(:,count11,11) =  FCR_data(:,i);
%                                             count11=count11+1;
%                                         else if log(i,9)==-75
%                                                 FCR_cond(:,count12,12) =  FCR_data(:,i);
%                                                 count12=count12+1;
%                                             else if log(i,9)==-100
%                                                     FCR_cond(:,count13,13) =  FCR_data(:,i);
%                                                     count13=count13+1;
%                                                 else if log(i,9)==-200
%                                                         FCR_cond(:,count14,14) =  FCR_data(:,i);
%                                                         count14=count14+1;
%                                                     end
%                                                 end
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%
%         end
%     end
% end
% clear count1 count2 count3 count4 count5 count6 count7 count8 count9 count10 count11 count12 count13 count14 count15
%
%
% FCR_cond = FCR_cond.*-1; % To get the data with minus up and positive down
% for j = 1:14
%    h=figure;
%    fig_str = ['Condition number ', int2str(j),' in tp'];
%    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%    for i=1:10
%        subplot (5,2,i); plot(FCR_cond(:,i,j));
%        ylim([-13 13]); xlim([0 points]);
%    end
%    title(fig_str);
%    print(h,fig_str,'-dpng');
%    close
% end
%
% ntrial=size(FCR_cond,2);
% ncond=size(FCR_cond,3);
% max_min = zeros(ntrial,ncond,2) % max = max_min(:,:,1); min = max_min(:,:,2)
% all_ind_max = zeros(ntrial,ncond); all_ind_min = zeros(ntrial,ncond);
% for i = 1:ncond
%     for j = 1:ntrial
%         if i < 4
%             [max_min(j,i,1), ind_max]=max(FCR_cond(370:550,j,i));
%             [max_min(j,i,2), ind_min]=min(FCR_cond(370:550,j,i));
%             all_ind_max(j,i) = ind_max+370; clear ind_max
%             all_ind_min(j,i) = ind_min+370; clear ind_min
%         else
%             [max_min(j,i,1), ind_max]=max(FCR_cond(4300:4450,j,i));
%             [max_min(j,i,2), ind_min]=min(FCR_cond(4300:4450,j,i));
%             all_ind_max(j,i) = ind_max+4300; clear ind_max
%             all_ind_min(j,i) = ind_min+4300; clear ind_min
%         end
%     end
% end
%
% hr_amp = max_min(:,:,1)-max_min(:,:,2);
% mean_amp = mean(hr_amp,1);


%             [delay, amplitude]=ginput(4); % base and peak M-resp, base and peak H-resp
%             M_resp(i,c,j,1) = amplitude(2)-amplitude(1);
%             M_resp(i,c,j,2) = delay(1)-stim_time;
%             H_resp(i,c,j,1) = amplitude(4)-amplitude(3);
%             H_resp(i,c,j,2) = delay(3)-stim_time;


%curr_data = zeros(15,3);
% figure;
% for c = 1:size(FCR_cond,3) % COND LOOP
%     for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
%         curr_data(c,1) = mean(M_resp(:,c,j,1));
%         curr_data(c,2) = mean(H_resp(:,c,j,1));
%         curr_data(c,3) = curr_data(c,2)./curr_data(c,1); %post
%         subplot(2,2,j*2-1)
%         bar(curr_data(:,1:2))
%         subplot(2,2,2*j)
%         bar(curr_data(:,3))
%         hold on
%     end
% end
%
%
% % PLOT BY CONDITION MEAN AND STEP IN SUBPLOTS
% cd(fig_folder)
% if exist('PLOT','dir') == 0; mkdir('PLOT'); else end
% plot_folder=[fig_folder,'/PLOT'];
% cd(plot_folder)
% if exist('MEAN_COMPARE','dir') == 0; mkdir('MEAN_COMPARE'); else end
% plot2_folder=[fig_folder,'/PLOT/MEAN_COMPARE'];
% cd(plot2_folder)
% DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
% date2_folder = [plot2_folder, '/', DateString];
% cd(date2_folder)
%
% % Without zoom
% for c = 1:size(FCR_cond,3) % COND LOOP
%     h=figure;
%     for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
%         duration_ms=points(j)/s_rate(j)*1000;
%         step=duration_ms/points(j);
%         time_axis=0:step:duration_ms-step; time_axis=time_axis';
%         ylim_up = squeeze(abs(min(min(min(min(FCR_cond))))))+0.5;
%         mean_FCR_cond = mean(squeeze(FCR_cond(:,:,c,j)),2);
%         if c < 4
%             bar_start = 26; bar_end = 34;
%         else if c == 15
%                 bar_start = 216; bar_end = 224;
%             else
%                 bar_start = 221; bar_end = 229;
%             end
%         end
%
%         subplot(2,1,j)
%         plot(time_axis,mean_FCR_cond);
%         hold on; line([bar_start bar_start],[-4 4],'Color',[1 0 0]);
%         hold on; line([bar_end bar_end],[-4 4],'Color',[1 0 0]);
%         fig_str = [fig_name, ' ', cond_name{c}, ' ', step_name{j}];
%         title(fig_str);
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%         set (gca,'Ydir','reverse')
%         ylim([-ylim_up ylim_up]); xlim([0 duration_ms]);
%         hold on;
%     end
% %     pause;
%             save_str=[fig_name,'_',cond_name{c},'_', step_name{j},'.fig'];
%             print(h,fig_str,'-dpng','-r300');
%             savefig(save_str);
%     close
% end
%
% % With zoom
% for c = 1:size(FCR_cond,3) % COND LOOP
%     h=figure;
%     for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
%         duration_ms=points(j)/s_rate(j)*1000;
%         step=duration_ms/points(j);
%         time_axis=0:step:duration_ms-step; time_axis=time_axis';
%         ylim_up_mean = squeeze(abs(min(min(min(min(mean_FCR_cond))))))+0.5;
%         mean_FCR_cond = mean(squeeze(FCR_cond(:,:,c,j)),2);
%         if c < 4
%             plot_start = 0; plot_end = 50;
%             bar_start = 26; bar_end = 34;
%         else if c == 15
%                 plot_start = 190; plot_end = 240;
%                 bar_start = 216; bar_end = 224;
%             else
%                 plot_start = 195; plot_end = 245;
%                 bar_start = 221; bar_end = 229;
%             end
%         end
%         subplot(2,1,j)
%         plot(time_axis,mean_FCR_cond);
%         hold on; line([bar_start bar_start],[-4 4],'Color',[1 0 0]);
%         hold on; line([bar_end bar_end],[-4 4],'Color',[1 0 0]);
%         fig_str = ['ZOOM ',fig_name, ' ', cond_name{c}, ' ', step_name{j}];
%         title(fig_str);
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%         set (gca,'Ydir','reverse')
%         ylim([-ylim_up_mean ylim_up_mean]); xlim([plot_start plot_end]);
%         hold on;
%     end
% %     pause;
%     save_fig=['ZOOM_',fig_name,'_',cond_name{c},'_', step_name{j},'.fig'];
%     save_png=['ZOOM_',fig_name,'_',cond_name{c},'_', step_name{j},'.png'];
%     print(h,save_png,'-dpng','-r300');
%     savefig(save_str);
%     close
% end

%% CREATE AND SAVE GRAPHS OF EMG SIGNAL - With zoom - MEAN and ALL trials
% for c = 1:size(FCR_cond,3) % COND LOOP
%     h=figure;
%     for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
%         duration_ms=points(j)/s_rate(j)*1000;
%         step=duration_ms/points(j);
%         time_axis=0:step:duration_ms-step; time_axis=time_axis';
%         mean_FCR_cond = mean(squeeze(FCR_cond(:,:,c,j)),2);
%         if c < 4
%             plot_start = 0; plot_end = 50;
%             bar_start = 26; bar_end = 34;% dely + 10 ms preceeding the stimulus
%         else if c == 15
%                 plot_start = 190; plot_end = 240;
%                 bar_start = 216; bar_end = 224;
%             else
%                 plot_start = 195; plot_end = 245;
%                 bar_start = 221; bar_end = 229;
%             end
%         end
%         for i=1:size(FCR_cond,2) % TRIAL LOOP
%             subplot(2,2,2*j-1)
%             plot(time_axis,FCR_cond(:,i,c,j));
%             hold on; line([bar_start bar_start],[-4 4],'Color',[1 0 0]);
%             hold on; line([bar_end bar_end],[-4 4],'Color',[1 0 0]);
%             set (gca,'Ydir','reverse')
%             ylim([-10 10]); xlim([plot_start plot_end]);
%             hold on;
%         end
%         subplot(2,2,2*j)
%         plot(time_axis,mean_FCR_cond);
%         hold on; line([bar_start bar_start],[-4 4],'Color',[1 0 0]);
%         hold on; line([bar_end bar_end],[-4 4],'Color',[1 0 0]);
%         fig_str = ['ZOOM ',fig_name, ' ', cond_name{c}, ' ', step_name{j}];
%         title(fig_str);
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%         set (gca,'Ydir','reverse')
%         ylim([-10 10]); xlim([plot_start plot_end]);
%         hold on;
%     end
%     pause;
% %     save_fig=['2ZOOM_',fig_name,'_',cond_name{c},'_', step_name{j},'.fig'];
% %     save_png=['2ZOOM_',fig_name,'_',cond_name{c},'_', step_name{j},'.png'];
% %     print(h,save_png,'-dpng','-r300');
% %     savefig(save_fig);
%     close
% end


% ONE FIGURE PER TRIALS - PRE/POST

% cd(fig_folder)
% if exist('PLOT','dir') == 0 mkdir('PLOT'); else end
% plot_folder=[fig_folder,'/PLOT'];
% cd(plot_folder)
% if exist('DETECT','dir') == 0 mkdir('DETECT'); else end
% plot_detect_folder=[fig_folder,'/PLOT/DETECT'];
% cd(plot_detect_folder)
% DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
% date2_folder = [plot_detect_folder, '/', DateString];
% cd(date2_folder)
% %
% M_resp = zeros(size(FCR_cond,2), size(FCR_cond,3), size(FCR_cond,4), 2); % trials, condition,  pre/post, amplitude / delay
% H_resp = zeros(size(FCR_cond,2), size(FCR_cond,3), size(FCR_cond,4), 2); % trials, condition,  pre/post, amplitude / delay
% rep=input('Do you want the graphs too? (y/n) - into simple quotes :    ');
% for c = 1:size(FCR_cond,3) % COND LOOP
%     for i=1:size(FCR_cond,2) % TRIAL LOOP
%         if strcmp(rep,'y')==1
%             h=figure;
%         end
%         for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
%             duration_ms=points(j)/s_rate(j)*1000;
%             step=duration_ms/points(j);
%             time_axis=0:step:duration_ms-step; time_axis=time_axis';
%             mean_FCR_cond = mean(squeeze(FCR_cond(:,:,c,j)),2);
%             if c < 4
%                 plot_start = 0; plot_end = 50;
%                 M_start = 2+10; M_end = 12+10;
%                 bar_start = 26; bar_end = 34;% dely + 10 ms preceeding the stimulus
%                 sub_start = sub_delay+10-2; sub_end = sub_delay2+10;
%                 stim_time= 10;
%             else if c == 15
%                     plot_start = 190; plot_end = 240;
%                     M_start = 2+200; M_end = 12+200;
%                     bar_start = 216; bar_end = 224;
%                     sub_start = sub_delay+200-2; sub_end = sub_delay2+200;
%                     stim_time= 200;
%                 else
%                     plot_start = 195; plot_end = 245;
%                     M_start = 2+205; M_end = 12+205;
%                     bar_start = 221; bar_end = 229;
%                     sub_start = sub_delay+205-2; sub_end = sub_delay2+205;
%                     stim_time= 205;
%                 end
%             end
%             % CONVERT DELAYS IN TIMEPOINTS
%             M_start_pt = M_start*points(j)/duration_ms;
%             M_end_pt = M_end*points(j)/duration_ms;
%             sub_start_pt = sub_start*points(j)/duration_ms;
%             sub_end_pt = sub_end*points(j)/duration_ms;
%
%             % FIND M and H
%             M_onset = FCR_cond(M_start_pt,i,c,j);
%             [M_peak, m_peak_delay] = min(FCR_cond(M_start_pt:M_end_pt,i,c,j));
%             M_amp = abs(M_peak - M_onset);
%             M_resp(i,c,j,1) = M_amp;
%
%             H_onset = FCR_cond(sub_start_pt,i,c,j);
%             [H_peak, h_peak_delay] = min(FCR_cond(sub_start_pt:sub_end_pt,i,c,j));
%             H_amp = abs(H_peak - H_onset);
%             H_resp(i,c,j,1) = H_peak;
%             % TRYING TO REMOVED THE DEPENDENCY ON THE TIME WINDOW
%             H_resp(i,c,j,2) = h_peak_delay + sub_start_pt;
%             H_resp(i,c,j,3) = H_onset;
%
%             if strcmp(rep,'y')==1 % ONLY GENERATE GRAPHS WHEN WANTED - SEE LINE 276
%                 % PLOT
%                 subplot(1,2,j)
%                 plot(time_axis,FCR_cond(:,i,c,j));
%                 % M-response delays
%                 hold on; line([M_start M_start],[-4 4],'Color',[0 0 1]);
%                 hold on; line([M_end M_end],[-4 4],'Color',[0 0 1]);
%                 % M detected MAX
%                 m_peak_disp = (m_peak_delay+M_start_pt)*duration_ms/points(j);
%                 hold on; line([m_peak_disp m_peak_disp],[M_peak-1 M_peak+1],'Color',[1 0 1], 'Marker','*');
%
%                 % Peter Fuhr H-r delays
%                 hold on; line([bar_start bar_start],[-4 4],'Color',[1 0 0]);
%                 hold on; line([bar_end bar_end],[-4 4],'Color',[1 0 0]);
%                 % Personal H-r delays from Viking
%                 hold on; line([sub_start sub_start],[-4 4],'Color',[0 1 0]);
%                 hold on; line([sub_end sub_end],[-4 4],'Color',[0 1 0]);
%                 % H detected MAX
%                 h_peak_disp = (h_peak_delay+sub_start_pt)*duration_ms/points(j);
%                 hold on; line([h_peak_disp h_peak_disp],[H_peak-1 H_peak+1],'Color',[1 0 1], 'Marker','*');
%
%                 set (gca,'Ydir','reverse')
%                 fig_str = [cond_name{c}, ' - ', int2str(i),' ', step_name{j}];
%                 title(fig_str);
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                 ylim([-10 10]); xlim([plot_start plot_end]);
%                 hold on;
%             end
%             %             pause;
%         end
%         if strcmp(rep, 'y')==1 % ONLY SAVE GRAPHS WHEN WANTED - SEE LINE 276
%             save_fig=['DETECT','_',cond_name{c},'_', ' trial ', int2str(i),'.fig'];
%             save_png=['DETECT','_',cond_name{c},'_', ' trial ', int2str(i),'.png'];
%             print(h,save_png,'-dpng','-r300');
%             savefig(save_fig);
%         end
%     end
%     close all
% end
%
% cd(ML_data_folder)
% save('MandH2.mat','M_resp','H_resp')