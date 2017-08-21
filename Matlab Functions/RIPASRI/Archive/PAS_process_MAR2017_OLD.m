function [data_goodrecs, resp] = PAS_process_MAR2017(PAS_ID,gp_ID, sub_ID)

% NOTES:
% 1- Ajouter les barres de la fen?tre de temps que le code va
% consid?rer dans sa recherche de l'amplitude des r?ponses.
% 2 - Ajouter pour chaque ?tape un input qui demande si l'utilisateur veut
% faire cette ?tape (mettre un if autour de tout le processus en question qui d?pend de la r?ponse positive)
% !! RESP SHOULD BE CALCULATED ONLY ON THE GOOD RECORDINGS


%% _______INITIALIZATION_____
gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
log_folder = [gp_folder,'/', sub_ID];
data_folder = [log_folder,'/Extract/'];
cd(log_folder)
if exist('Figures', 'dir')==0; mkdir('Figures'); end
fig_folder = [log_folder,'/Figures'];
cd(log_folder) 
if exist('Data', 'dir')==0; mkdir('Data'); end
MLdata_folder = [log_folder, '/Data']; 

cd(gp_folder);load time_windows; load time_windows_FCR

%% ____DATA PRE-PROCESSING_____
[data_cond, data_param] = PAS_preproc(PAS_ID,gp_ID,sub_ID);

cd(MLdata_folder)
if exist('data_cond.mat','file')==2
%     rep_replace=input('A file with the name data_cond.mat already exists, do you want to replace or keep both? (r/b):   ');
%     if strcmp(rep_replace,'b')
        copyfile('data_cond.mat','data_cond_OLD.mat');
%     end
end
save('data_cond', 'data_cond')

%% CHECKING PRESENCE OF BAD RECORDINGS

prev_badrec_rep=input('Did you already check for bad recordings for this subject (y/n) ?:  ');
if strcmp(prev_badrec_rep,'y') == 1 % checking if bad rec definition was already done 
    cd(MLdata_folder)
    load data_goodrecs.mat % if yes load the data of previously defined good recrodings
else
    badrec_rep=input('Do you want to check for bad recordings (y/n) ?:  ');
    badrecs = 0;
    if strcmp(badrec_rep,'y')==1
        cd(fig_folder)
        [badrecs, badrec_muscles] = check_badrec(data_cond, data_param, time_windows, time_windows_FCR); % else proceed to bad rec checking
    end
end
%% REMOVING BAD RECORDINGS
DateStamp=datestr(clock);
if strcmp(prev_badrec_rep,'y') == 0
    if badrecs==0
        data_goodrecs = data_cond;
    else
        if exist('Bad_Recordings.mat','file')==2
%             rep_replace=input('A file with the name Bad_Recordings.mat already exists, do you want to replace or keep both? (r/b):   ');
%             if strcmp(rep_replace,'b')
                copyfile('Bad_Recordings.mat',['Bad_Recordings_',DateStamp,'.mat']);
%             end
        end
        save('Bad_Recordings.mat','badrecs')
        
        [data_goodrecs] = remove_badrec(data_cond, badrecs, badrec_muscles, data_param, time_windows, time_windows_FCR);
        
        cd(MLdata_folder)
        if exist('data_goodrecs.mat','file')==2
%             rep_replace=input('A file with the same name already exists, do you want to replace or keep both? (r/b):   ');
%             if strcmp(rep_replace,'b')
                copyfile('data_goodrecs.mat',['data_goodrecs_',DateStamp,'.mat']);
%             end
        end
        save('data_goodrecs','data_goodrecs');
        
        disp('End of BAD REC REMOVAL CHECK');
    end
end
%% MEASURING RESPONSES

[resp] = resp_analysis(data_goodrecs,data_param,time_windows, time_windows_FCR);

cd(MLdata_folder); 
if exist('Resp_Amp.mat','file')==2
%     rep_replace=input('A file with the name Resp_Amp.mat already exists, do you want to replace or keep both? (r/b):   ');
%     if strcmp(rep_replace,'b')
        copyfile('Resp_Amp.mat',['Resp_Amp_',DateStamp,'.mat']);    
%     end
end
save('Resp_Amp','resp');

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
% %% PLOT ALL TRIALS - CONDITION x STEP
% 
% cd(fig_folder)
% if exist('ALL_TRIALS', 'dir')==0; mkdir('ALL_TRIALS'); end
% cd('ALL_TRIALS')
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
%             if m < 3 % FCR has different time windows than the other muscle
%                 curr_time_win = time_windows(:,c);
%             else
%                 curr_time_win = time_windows_FCR(:,c);
%             end
%             
%             ylim_all= abs(min(min(min(min(data_goodrecs(:,:,m,:,:))))));
%             
%             for i = 1:size(data_cond,2) % trials per condition
%                 
%                 %___________GRAPH PARAM__________
%                 max_start=curr_time_win_pt(1); max_end=curr_time_win_pt(2);
%                 curr_time_win_pt = round(curr_time_win*points(j)/duration_ms);
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                 hold on
%                 
%                 %___________ALL TRIALS____________
%                 %                 curr_data = squeeze(data_cond(:,:,m,c,j));
%                 curr_data = squeeze(data_goodrecs(:,:,m,c,j));
%                 subplot(3,8,count)
%                 plot(time_axis,curr_data); % data, trials, muscle, cond, step
%                 ylim([-ylim_all ylim_all]); xlim([curr_time_win(1)-15 curr_time_win(2)+15])
%                 set (gca,'Ydir','reverse')
%                 %                 fig_str=[cond_name{c}, ' - ',step_name{j}];
%                 %                 title(fig_str,'FontSize', 8)
%                 
%                 hold on; line([curr_time_win(1) curr_time_win(1)],[-3 3],'Color',[0 1 0]);
%                 hold on; line([curr_time_win(2) curr_time_win(2)],[-3 3],'Color',[0 1 0]);
%                 
%             end
%             %___________MEAN____________
%             %                 mean_data = mean(squeeze(data_cond(:,:,m,c,j)),2);
%             mean_data = mean(squeeze(data_goodrecs(:,:,m,c,j)),2);
%             subplot(3,8,count+1)
%             plot(time_axis,mean_data,'Color','b'); % data, trials, muscle, cond, step
%             ylim([-ylim_all ylim_all]); xlim([curr_time_win(1)-15 curr_time_win(2)+15])
%             set (gca,'Ydir','reverse')
%             
%             hold on; line([curr_time_win(1) curr_time_win(1)],[-3 3],'Color',[0 1 0]);
%             hold on; line([curr_time_win(2) curr_time_win(2)],[-3 3],'Color',[0 1 0]);
%             
%             % Detected MAX
%             [mean_max, mean_delay] = min(mean_data(max_start:max_end,:));
%             h_peak_disp = round((mean_delay+curr_time_win_pt(1))*duration_ms/points(j));
%             hold on; line([h_peak_disp h_peak_disp],[mean_max-0.5 mean_max+0.5],'Color',[1 0 1], 'Marker','.');
%             % Detected MIN
%             [mean_ref, mean_ref_delay] = max(mean_data(max_start+mean_delay-1:max_end,:));
%             h_ref_disp = round((mean_ref_delay+curr_time_win_pt(1)+mean_delay)*duration_ms/points(j));
%             hold on; line([h_ref_disp h_ref_disp],[mean_ref-0.5 mean_ref+0.5],'Color',[0 0 1], 'Marker','.');
%             
%             peak_str = sprintf('%.3f',mean_max);
%             hold on; text(((curr_time_win(1)+curr_time_win(2))/2)-2.7,mean_max-1,peak_str)
%             
%             count=count+2;
%         end
%     end
%     pause;
%     print_str=m_name{m};
%     disp(['Saving ',print_str]);
%     print(h,print_str,'-dpng','-r900');
%     save_fig=[print_str,'.fig'];
%     savefig(save_fig);
%     close
% end

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
