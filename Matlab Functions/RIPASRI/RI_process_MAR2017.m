function [H_resp, FCR_cond]=RI_process_MAR2017(PAS_ID, gp_ID,sub_ID, subChoice)

% _____________NOTES_______________

% Once the RI_process_ALL finished and this code finalized,
% RI_process_MAR2017 will take as a parameter the string of the group to be
% analyzed so that it can be called from RI_process_ALL

% Once everything works properly, maybe transform main steps into nested
% functions so that one can easily decide to do EDC without having to
% manually change every lines and just call the functin with either FDC or
% EDC in parameter

% Change line colors and use the rgb function

%% INITIALIZATION
gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',gp_ID];
log_folder = ['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/HV/', sub_ID];
data_folder = [log_folder,'/Extract/'];
fig_folder = [log_folder,'/Figures/RI'];
ML_data_folder = [log_folder,'/Data'];

disp('End of INITIALIZATION')

%% INITIALIZATION FOR GRAPHS
step_name= {'Baseline', 'Post'};
cond_name={'single pulse','-1ms ISI', '-0,5ms ISI', 'synchronous','0,5ms ISI', '1ms ISI',...
    '2ms ISI','5ms ISI','7,5ms ISI','10ms ISI','25ms ISI','50ms ISI','75ms ISI','100ms ISI','200ms ISI'};
fig_name='FCR';
sub_ID_all={'BAM','ALS','ETB','CAM','ELB','CAR','GIB','LUM','MAS','ISH'};

cd(log_folder)
if ~exist('Figures','dir'); mkdir('Figures'); else end
cd('Figures')
if ~exist('RI','dir'); mkdir('RI'); else end
fig_folder=[log_folder,'/Figures/RI'];
cd(fig_folder)

cd(gp_folder)
list_folders=dir;
for k = length(list_folders):-1:1
    if ~list_folders(k).isdir
        list_folders(k) = [ ];
        continue
    end
    fname = list_folders(k).name;
    if fname(1) == '.' || fname(1) == '*'
        list_folders(k) = [ ];
    end
end
for k = 1:length(list_folders)
    disp(list_folders(k).name)
    temp(k)=strcmp(list_folders(k).name,sub_ID);
end

data_xls=xlsread('data_collection.xlsx');
[osef,col]=max(temp); clear osef
sub_row=col+2;

sub_delay=data_xls(sub_row,4);
sub_delay2=data_xls(sub_row,5);


disp('end of GRAPH INITIALIZATION')

%% CHECK IF SUBJECT WAS ALREADY PROCESSED

rep_proc = input(['Did you already process this subject - ', sub_ID,' ? (y/n):  ']);
if strcmp(rep_proc, 'y')==1
    rep_load = input('Do you want to load the data already processed?(y/n):  ');
    if strcmp(rep_load, 'y')== 1
        cd(ML_data_folder)
        while ~exist('MandH.mat','file')==1
            disp('MandH.cannot be found, please check in the folder')
            pause
        end
        load MandH.mat
        while ~exist('FCR_cond.mat','file')==1
            disp('FCR_cond.cannot be found, please check in the folder')
            pause
        end
        load FCR_cond.mat
        rep_bd=input('Did you process bad detections? (y/n):  ');
        if strcmp(rep_bd,'y')==1
            while ~exist('Bad_Detect.mat','file')==1
                disp('Bad_Detect.cannot be found, please check in the folder')
                pause
            end
            load Bad_Detect.mat
        end
    end
end

%% EXTRACTING INFO FROM LOG
% if strcmp(rep_proc, 'y')==0
cd(log_folder)
log_base_temp=importdata('RI_base.log');
log_base = log_base_temp.data;
log_post_temp=importdata('RI_post.log');
log_post = log_post_temp.data;

s_rate=[log_base(1),log_post(1)]; % I could only take one as they should all be the same, but jic to avoid an error
points=[log_base(2),log_post(2)];
n_muscle=[log_base(3),log_post(3)];
n_trial=[log_base(4),log_post(4)];

disp('Please verify that RI_base.log and RI_post.log were copied to excel files and the first three rows deleted and press any key')
pause;

% open RI_base.xlsx; pause;
% open RI_post.xlsx; pause;

RI_base = xlsread('RI_base.xlsx');
RI_post = xlsread('RI_post.xlsx');

conditions = zeros(max(n_trial), 2, 2); %trials, delas, step
conditions(:,1,1)= squeeze(RI_base(:,7));conditions(:,2,1)= squeeze(RI_base(:,9)); % conditions (:,:,1) = RIbase
conditions(:,1,2)= squeeze(RI_post(:,7));conditions(:,2,2)= squeeze(RI_post(:,9)); % conditions (:,:,2) = RIpost

clear log_base_temp log_base log_post_temp log_post

disp('End of LOG EXTRACTION')
% pause;
%% GETTING DATA FROM NGUYET EXTRACT
cd(data_folder)

% Initialisation
ds_folder = {'RI_base', 'RI_post'};
ri_data=zeros(max(points),2,max(n_trial),2); % timeopints, muscle, trials, step(baseline; post)
% Extraction loop - loads file 1 buy 1 and add the data in the ri_data matrix
for j = 1:2 % baseline / post
    cd((ds_folder{j}))
    for i = 1:max(n_trial)
        while exist(int2str(i),'file') == 0
            disp(['Programmed paused because file ', int2str(i),' is missing in ', ds_folder{j}]);
            pause
        end
        temp=importdata(int2str(i));
        ri_data(:,1,i,j) = temp(1:points(j),1); % I put row 1:points in case there was a mistake during the nguyet and a file was recorded twice
        ri_data(:,2,i,j) = temp(1:points(j),2); clear temp;
    end
    cd(data_folder)
end

% Separating FCR from EDC (we are mostly interested in FCR)
FCR_data = squeeze(ri_data(:,1,:,:)); % FCR - [timepoints x trials x step]
EDC_data = squeeze(ri_data(:,2,:,:));

cd(log_folder)
if exist('Data','dir') == 0; mkdir('Data'); else end
ML_data_folder=[log_folder,'/Data'];
cd(ML_data_folder)
DateString = datestr(clock);

name_file_RI = [DateString, 'RI_data'];
save(name_file_RI, 'ri_data', 'FCR_data','EDC_data')

disp('End of DATA EXTRACTION')
% pause;
%% EXTRACT DATA IN CONDITIONS

FCR_cond = order_by_cond_MAY2017(FCR_data,conditions,max(n_trial));
EDC_cond = order_by_cond_MAY2017(EDC_data,conditions,max(n_trial));

% name_file = 'FCR_cond'; count=1;
% while exist([name_file, '.mat'], 'file') == 2
%     count=count+1;
%     name_file = [name_file, '_', int2str(count)];
% end
name_file = [DateString, 'FCR_cond'];
save(name_file, 'FCR_cond');

% name_file_EDC = 'EDC_cond'; count=1;
% while exist([name_file_EDC, '.mat'], 'file') == 2
%     count=count+1;
%     name_file_EDC = [name_file_EDC, '_', int2str(count)];
% end
% save(name_file_EDC, 'EDC_cond');

disp('End of COND REARRANGING')
% pause;

%% CHECK SUBJECT H-REFLEX DELAYS

test=1; search=1; count = 1;
rep_fig_save = input('Do you want to save the H-reflex DELAYS figures?(y/n):  ');

cd(fig_folder)
if ~exist('Delays','dir'); mkdir('Delays');end
cd('Delays')

DateString = datestr(clock);
while search ~= 0
    duration_ms=points(1)/s_rate(1)*1000;
    step=duration_ms/points(1);
    time_axis=0:step:duration_ms-step; time_axis=time_axis';
    for j = 1:size(FCR_cond,4) % STEP LOOP
        h=figure;
        for c = 1:size(FCR_cond,3) % COND LOOP
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
            ylim_down = squeeze(min(min(min(FCR_cond(plot_start_pt+45:plot_end_pt,:,c,1)))))-0.2;
            ylim_up = squeeze(max(max(max(FCR_cond(plot_start_pt+45:plot_end_pt,:,c,1)))))+0.2;
            
            for i=1:size(FCR_cond,2) % TRIAL LOOP
                plot(time_axis,squeeze(FCR_cond(:,i,c,j)));
                %             hold on; line([time_windows(1,c) time_windows(1,c)],[-4 4],'Color',[1 0 0]);
                %             hold on; line([time_windows(2,c) time_windows(2,c)],[-4 4],'Color',[1 0 0]);
                fig_str = [fig_name, ' ', cond_name{c}, ' ', step_name{j}];
                title(fig_str);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                set (gca,'Ydir','reverse')
                ylim([ylim_down ylim_up]);
                %                 ylim([-1.5 1.5])
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
            %             hold on; line([sub_start sub_start],[-3 3],'Color',[0 1 0]);
            %             hold on; line([sub_end_test sub_end_test],[-3 3],'Color',[0 1 0])
            count = count +1;
        end
        if strcmp(rep_fig_save,'y')==1
            save_str=[DateString,step_name{j},'All_Sup_BEST_D.fig'];
            %         print(h,fig_str,'-dpng','-r300');
            savefig(save_str);
        end
        pause;
    end
    rep_sub_end=input('Enter the value for the end of the H-reflex that would be most proper (if the value is good type 0):  ');
    close all
    
    if rep_sub_end ~= 0
        sub_end_test=rep_sub_end;
    else
        sub_delay2 = sub_end_test-195-10;
        search=0;
    end
    test=test+1;
end



%% CREATE AND SAVE GRAPHS OF RAW EMG SIGNAL

cd(fig_folder)
if exist('Raw_Signal','dir') == 0; mkdir('Raw_Signal'); else end
cd('Raw_signal')

cd(fig_folder);cd('Raw_signal')

rep_sup=input('Do you want to se all trials super-imposed? (y/n):  ');

if exist('SuperImposed','dir') == 0; mkdir('SuperImposed'); else end
cd('SuperImposed')
% mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
% cd(DateString)

DateString = datestr(clock);
if strcmp(rep_sup,'y')==1
    rep_save_figS=input('Do you ant to save the raw separate figures? (y/n):  ');
    % PLOT BY CONDITION AND STEP SEPARATELY
    for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
        duration_ms=points(j)/s_rate(j)*1000;
        step=duration_ms/points(j);
        time_axis=0:step:duration_ms-step; time_axis=time_axis';
        for c = 1:size(FCR_cond,3) % COND LOOP
            h=figure;
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
            sub_start_pt = round(sub_start*points(j)/duration_ms);
            sub_end_pt = round(sub_end*points(j)/duration_ms);
            plot_start_pt = round(plot_start*points(j)/duration_ms);
            plot_end_pt = round(plot_end*points(j)/duration_ms);
            ylim_down = squeeze(min(min(min(FCR_cond(plot_start_pt+25:plot_end_pt,:,c,j)))))-0.2;
            ylim_up = squeeze(max(max(max(FCR_cond(plot_start_pt+25:plot_end_pt,:,c,j)))))+0.2;
            
            for i=1:size(FCR_cond,2) % TRIAL LOOP
                plot(time_axis,squeeze(FCR_cond(:,i,c,j)));
                %             hold on; line([time_windows(1,c) time_windows(1,c)],[-4 4],'Color',[1 0 0]);
                %             hold on; line([time_windows(2,c) time_windows(2,c)],[-4 4],'Color',[1 0 0]);
                fig_str = [fig_name, ' ', cond_name{c}, ' ', step_name{j}];
                title(fig_str);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                set (gca,'Ydir','reverse')
                ylim([ylim_down ylim_up]); xlim([plot_start plot_end]);
                hold on;
            end
            % Personal H-r delays from Viking
            hold on; line([sub_start sub_start],[-ylim_up+0.7 ylim_up-0.7],'Color',[0 1 0]);
            hold on; line([sub_end sub_end],[-ylim_up+0.7 ylim_up-0.7],'Color',[0 1 0]);
            
            if strcmp(rep_save_figS,'y')==1
                save_str=[DateString,'_',cond_name{c},'_', step_name{j},'.fig'];
                print(h,fig_str,'-dpng','-r300');
                savefig(save_str);
            else
                pause;
            end
            close
        end
    end
end

rep_sep=input('Do you want to se all trials in separate subplots? (y/n):  ');

cd(fig_folder);cd('Raw_signal')
if exist('Separate','dir') == 0; mkdir('Separate'); else end
cd('Separate')

DateString = datestr(clock);
if strcmp(rep_sep,'y')==1
    rep_save_fig=input('Do you ant to save the raw separate figures? (y/n):  ');
    for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
        duration_ms=points(j)/s_rate(j)*1000;
        step=duration_ms/points(j);
        time_axis=0:step:duration_ms-step; time_axis=time_axis';
        for c = 1:size(FCR_cond,3) % COND LOOP
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
            
            sub_start_pt = round(sub_start*points(j)/duration_ms);
            sub_end_pt = round(sub_end*points(j)/duration_ms);
            plot_start_pt = round(plot_start*points(j)/duration_ms);
            plot_end_pt = round(plot_end*points(j)/duration_ms);
            
            ylim_down = squeeze(min(min(min(FCR_cond(plot_start_pt+45:plot_end_pt,:,c,j)))))-0.2;
            ylim_up = squeeze(max(max(max(FCR_cond(plot_start_pt+45:plot_end_pt,:,c,j)))))+0.2;
            
            h=figure;
            for i=1:size(FCR_cond,2) % TRIAL LOOP
                subplot(2,5,i)
                curr_data =squeeze(FCR_cond(:,i,c,j));
                plot(time_axis,curr_data);
                %             hold on; line([time_windows(1,c) time_windows(1,c)],[-4 4],'Color',[1 0 0]);
                %             hold on; line([time_windows(2,c) time_windows(2,c)],[-4 4],'Color',[1 0 0]);
                fig_str = [sub_ID, ' ', cond_name{c}, ' ', step_name{j}];
                %title(fig_str);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                set (gca,'Ydir','reverse')
                ylim([ylim_down ylim_up]);
                xlim([plot_start plot_end]);
                % Personal H-r delays from Viking
                hold on; line([sub_start-1 sub_start-1],[curr_data(sub_start_pt)-2 curr_data(sub_start_pt)+2],'Color',[0 1 0]);
                hold on; line([sub_end+1 sub_end+1],[curr_data(sub_end_pt)-2 curr_data(sub_end_pt)+2],'Color',[0 1 0]);
                hold on;
            end
            set(gcf,'NextPlot','add');
            axes;
            ht=title(fig_str);
            set(gca,'Visible','off');
            set(ht,'Visible','on');
            if strcmp(rep_save_fig,'y')==1
                save_str=[DateString, sub_ID,'_',cond_name{c},'_', step_name{j},'.fig'];
                print(h,fig_str,'-dpng','-r600');
                %savefig(save_str);
            else
                pause;
            end
            close
        end
    end
end
clear rep_save_fig

rep1_bad_rec = input('Did you see any bad recordings that need to be removed? (y/n):  ');

if strcmp(rep1_bad_rec,'y')==1
    %% PLOT ALL TRIALS (SUBPLOT) TO SEE BAD RECORDINGS
    
    rep2_bad_rec=input('Do you want to enter bad recordings? (y/n):   ');
    
    if strcmp(rep2_bad_rec,'y')==1
        cd(fig_folder);
        if exist('Bad_recs','dir') == 0; mkdir('Bad_recs'); else end
        cd('Bad_recs')
        DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
        cd(DateString)
        
        Bad_Recs=zeros(size(FCR_cond,2), size(FCR_cond,3), size(FCR_cond,4));
        for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
            for c = 1:size(FCR_cond,3) % COND LOOP
                h=figure;
                duration_ms=points(j)/s_rate(j)*1000;
                step=duration_ms/points(j);
                time_axis=0:step:duration_ms-step; time_axis=time_axis';
                ylim_up = squeeze(abs(min(min(min(min(FCR_cond))))))+0.5;
                if c == 1
                    plot_start = 190; plot_end = 240;
                else if c < 5
                        plot_start = 0; plot_end = 50;
                    else
                        plot_start = 195; plot_end = 245;
                    end
                end
                for i=1:size(FCR_cond,2) % TRIAL LOOP
                    subplot(2,5,i)
                    plot(time_axis,FCR_cond(:,i,c,j));
                    %             title(fig_str);
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                    set (gca,'Ydir','reverse')
                    ylim([-ylim_up ylim_up]); xlim([plot_start plot_end]);
                    hold on;
                end
                fig_str = ['WRITE DOWN THE BAD TRIALS FOR THIS CONDITION AND STEP --->', cond_name{c}, ' ', step_name{j}];
                set(gcf,'NextPlot','add');
                axes;
                ht=title(fig_str);
                set(gca,'Visible','off');
                set(ht,'Visible','on');
                
                save_fig=['Check_',cond_name{c},'_', step_name{j},'.fig'];
                save_png=['Check_',cond_name{c},'_', step_name{j},'.png'];
                print(h,save_png,'-dpng','-r300');
                savefig(save_fig);
                pause; close;
                bad_recs=input('Please enter in square brackets and separated by comas the bad trials for this condition (0 if none):  ');
                Bad_Recs(1:length(bad_recs),c,j)=bad_recs';
            end
        end
        cd(ML_data_folder)
        save('Bad_Recs.mat','Bad_Recs')
        
        %% REMOVE BAD RECORDINGS
        
        rep_remove=input('Is there any bad recordings to completely remove ? (y/n):    ');
        
        if strcmp(rep_remove, 'y')
            
            data_goodrecs=zeros(size(FCR_cond,1),size(FCR_cond,2),size(FCR_cond,3),size(FCR_cond,4));
            for j = 1: size(bad_recs,4)
                for c = 1:size(bad_recs,3)
                    for i = 1:size(bad_recs,2)
                        goodrecs=setdiff(1:size(data_cond,2),bad_recs(:,i,c,j));
                        for cc = 1 : length(goodrecs)
                            data_goodrecs(:,cc,i,c,j)= FCR_cond(:,goodrecs(cc),i,c,j);
                        end
                    end
                end
            end
            
        end
    end
end

if rep1_bad_rec == 1
    if rep2_bad_rec == 1
        FCR_cond_original = FCR_cond;
        clear FCR_cond;
        FCR_cond = data_goodrecs;
    end
end
%% CHECK AUTOMATIC DETECTION OF M AND H (M_resp and H_resp) - Define those to be corrected

% ALL TRIALS (1 FIGURE 10 SUBPLOT PER COND)


rep_bad_detect = input('Do you want to enter bad H-r detections? (y/n):   ');

M_resp = zeros(size(FCR_cond,2), size(FCR_cond,3), size(FCR_cond,4), 2); % trials, condition,  pre/post, amplitude / delay
H_resp = zeros(size(FCR_cond,2), size(FCR_cond,3), size(FCR_cond,4), 5); % t, cond, step, Hpeak / delay / Href / Href delay Hamp
Bad_Detect=zeros(size(FCR_cond,2), size(FCR_cond,3), size(FCR_cond,4));

rep_fig=input('Do you want to save the figures of H-r detections ? (y/n):     ');
rep_fig2=input('Do you want to see the figures of H-r detections ? (y/n):     ');

if strcmp(rep_fig,'y')==1
    cd(fig_folder);
    if exist('H_detect','dir') == 0; mkdir('H_detect'); else end
    cd('H_detect')
    DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
    cd(DateString)
end
for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
    for c = 1:size(FCR_cond,3) % COND LOOP
        h=figure;
        for i=1:size(FCR_cond,2) % TRIAL LOOP
            duration_ms=points(j)/s_rate(j)*1000;
            step=duration_ms/points(j);
            time_axis=0:step:duration_ms-step; time_axis=time_axis';
            mean_FCR_cond = mean(squeeze(FCR_cond(:,:,c,j)),2);
            ylim_up=squeeze(abs(min(min(min(min(mean_FCR_cond))))));
            if c == 1
                plot_start = 199; plot_end = 235;
                M_start = 2+200; M_end = 12+200;
                bar_start = 216; bar_end = 224;
                % TEMPORARILY REMOVED THE -2 TO SEE IF THE DETECTION IS BETTER
                sub_start = sub_delay+200; sub_end = sub_delay2+200;
                stim_time= 200;
            else if c < 5
                    plot_start = 9; plot_end = 45;
                    M_start = 2+10; M_end = 12+10;
                    bar_start = 26; bar_end = 34;% dely + 10 ms preceeding the stimulus
                    % TEMPORARILY REMOVED THE -2 TO SEE IF THE DETECTION IS BETTER
                    sub_start = sub_delay+10; sub_end = sub_delay2+10;
                    stim_time= 10;
                else
                    plot_start = 204; plot_end = 240;
                    M_start = 2+205; M_end = 12+205;
                    bar_start = 221; bar_end = 229;
                    % TEMPORARILY REMOVED THE -2 TO SEE IF THE DETECTION IS BETTER
                    sub_start = sub_delay+205; sub_end = sub_delay2+205;
                    stim_time= 205;
                end
            end
            % CONVERT DELAYS IN TIMEPOINTS
            M_start_pt = M_start*points(j)/duration_ms;
            M_end_pt = M_end*points(j)/duration_ms;
            sub_start_pt = sub_start*points(j)/duration_ms;
            sub_end_pt = sub_end*points(j)/duration_ms;
            
            % FIND M and H
            M_onset = FCR_cond(M_start_pt,i,c,j);
            [M_peak, m_peak_delay] = min(FCR_cond(M_start_pt:M_end_pt,i,c,j));
            M_amp = abs(M_peak - M_onset);
            M_resp(i,c,j,1) = M_amp;
            
            [H_peak, h_peak_delay] = min(FCR_cond(sub_start_pt:sub_end_pt,i,c,j));
            [H_ref, h_ref_delay] = max(FCR_cond(sub_start_pt+h_peak_delay-1:sub_end_pt,i,c,j));
            H_amp = abs(H_peak - H_ref);
            H_resp(i,c,j,1) = H_peak;
            % TRYING TO REMOVE THE DEPENDENCY ON THE TIME WINDOW
            H_resp(i,c,j,2) = h_peak_delay + sub_start_pt;
            H_resp(i,c,j,3) = H_ref;
            H_resp(i,c,j,4) = h_ref_delay + sub_start_pt + h_peak_delay;
            H_resp(i,c,j,5) = H_amp;
            
            % PLOT
            subplot(2,5,i)
            plot(time_axis,FCR_cond(:,i,c,j));
            % M-response delays
            hold on; line([M_start M_start],[-ylim_up+0.7 ylim_up-0.7],'Color',[0 0 1]);
            hold on; line([M_end M_end],[-ylim_up+0.7 ylim_up-0.7],'Color',[0 0 1]);
            % M detected MAX
            m_peak_disp = (m_peak_delay+M_start_pt)*duration_ms/points(j);
            hold on; line([m_peak_disp m_peak_disp],[M_peak-0.3 M_peak+0.3],'Color',[1 0 1], 'Marker','.');
            
            %             % Peter Fuhr H-r delays
            %             hold on; line([bar_start bar_start],[-ylim_up+0.5 ylim_up-0.5],'Color',[1 0 0]);
            %             hold on; line([bar_end bar_end],[-ylim_up+0.5 ylim_up-0.5],'Color',[1 0 0]);
            % Personal H-r delays from Viking
            hold on; line([sub_start sub_start],[-ylim_up+0.7 ylim_up-0.7],'Color',[0 1 0]);
            hold on; line([sub_end sub_end],[-ylim_up+0.7 ylim_up-0.7],'Color',[0 1 0]);
            % H detected MAX
            h_peak_disp = (h_peak_delay+sub_start_pt)*duration_ms/points(j);
            hold on; line([h_peak_disp h_peak_disp],[H_peak-0.3 H_peak+0.3],'Color',[1 0 1], 'Marker','.');
            % H detected MIN
            h_ref_disp = (h_ref_delay + sub_start_pt + h_peak_delay)*duration_ms/points(j);
            hold on; line([h_ref_disp h_ref_disp],[H_ref-0.3 H_ref+0.3],'Color',[0 0 1], 'Marker','.');
            
            % Display peak-to-peak amplitude
            peak_str = sprintf('%.3f',H_amp);
            hold on; text(((sub_start+sub_end)/2)-2.7,H_peak-1,peak_str)
            
            set (gca,'Ydir','reverse')
            %             title(fig_str);
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            ylim([-ylim_up ylim_up]); xlim([plot_start plot_end]);
            % Title
            set(gcf,'NextPlot','add');
            axes;
            fig_str = [cond_name{c}, ' - ', step_name{j}];
            ht = title(fig_str);
            set(gca,'Visible','off');
            set(ht,'Visible','on');
            hold on;
        end
        %         pause;
        if strcmp(rep_fig,'y')==1
            save_fig=['DETECT','_',cond_name{c},'_', step_name{j},'.fig'];
            save_png=['DETECT','_',cond_name{c},'_', step_name{j},'.png'];
            print(h,save_png,'-dpng','-r600');
            %             savefig(save_fig);
        end
        if strcmp(rep_fig2,'y')==1
            pause;
        end
        close;
        if strcmp(rep_bad_detect,'y')==1
            bad_detect=input('Please enter in square brackets and separated by comas the bad trials for this condition (0 if none):  ');
            Bad_Detect(1:length(bad_detect),c,j)=bad_detect';
            clear bad_detect
        end
    end
end
cd(ML_data_folder)
DateString = datestr(clock);
save(['MandH_', DateString,'.mat'],'M_resp','H_resp')
if strcmp(rep_bad_detect,'y')==1
    save('Bad_Detect.mat','Bad_Detect')
end

disp('End of BAD H-R DETECTION')
%% MANUALLY DEFINE H-R MAX FOR WRONG DETECTIONS

if strcmp(rep_bad_detect,'y')==1
    rep_detect=input('Is there any wrong detections for the h-r? (y/n):   ');
    Corrected_H=zeros(size(FCR_cond,2), size(FCR_cond,3),size(FCR_cond,4),2);
    if strcmp(rep_detect,'y')
        for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
            Err_Col=any(Bad_Detect(:,:,j)); % finds non-zero columns
            for c = 1:size(FCR_cond,3) % COND LOOP
                duration_ms=points(j)/s_rate(j)*1000;
                step=duration_ms/points(j);
                time_axis=0:step:duration_ms-step; time_axis=time_axis';
                if c < 5 && c~=1
                    plot_start = 0; plot_end = 50;
                    sub_start = sub_delay+10; sub_end = sub_delay2+10;
                else if c == 1
                        plot_start = 190; plot_end = 240;
                        sub_start = sub_delay+200; sub_end = sub_delay2+200;
                    else
                        plot_start = 195; plot_end = 245;
                        sub_start = sub_delay+205; sub_end = sub_delay2+205;
                    end
                end
                if Err_Col(c)==1
                    for i=1:size(FCR_cond,2) % TRIAL LOOP
                        if Bad_Detect(i,c,j) ~= 0
                            
                            sub_start_pt = sub_start*points(j)/duration_ms;
                            sub_end_pt = sub_end*points(j)/duration_ms;
                            
                            % Detected H-max
                            H_onset = FCR_cond(sub_start_pt,i,c,j);
                            [H_peak, h_peak_delay] = min(FCR_cond(sub_start_pt:sub_end_pt,trial,c,j));
                            
                            trial = Bad_Detect(i,c,j);
                            h=figure;
                            plot(time_axis,FCR_cond(:,trial,c,j));
                            % Personal H-r delays from Viking
                            hold on; line([sub_start sub_start],[-ylim_up+0.5 ylim_up-0.5],'Color',[0 1 0]);
                            hold on; line([sub_end sub_end],[-ylim_up+0.5 ylim_up-0.5],'Color',[0 1 0]);
                            % H detected MAX
                            h_peak_disp = (h_peak_delay+sub_start_pt)*duration_ms/points(j);
                            hold on; line([h_peak_disp h_peak_disp],[H_peak-1 H_peak+1],'Color',[1 0 1], 'Marker','*');
                            
                            set (gca,'Ydir','reverse')
                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                            ylim_up = squeeze(abs(min(min(min(min(FCR_cond(:,trial,c,j)))))))+0.5;
                            ylim([-ylim_up ylim_up]); xlim([plot_start plot_end]);
                            
                            % New manually defined H-max
                            [delay, amplitude]=ginput(1); % base and peak M-resp, base and peak H-resp
                            Corrected_H(i,c,j,1) = amplitude;
                            Corrected_H(i,c,j,2) = delay;
                            clear amplitude
                        else
                            break
                        end
                    end
                end
            end
            close all
        end
    end
    cd(ML_data_folder)
    save('Corrected_H.mat','Corrected_H')
    
    
    
    %% REPLACED WRONG MAX BY CORRECTED
    
    corr_idx = Corrected_H ~=0 ;
    H_resp_corr = H_resp;
    H_resp_corr(corr_idx)=Corrected_H(corr_idx);
    
    % _____________NOTES______________
    % problem with the new detected H-r
    % probably a problem of indexes, check how they are saved in the previous
    % section
    
    % CHECK CORRECTED H
    
    cd(fig_folder);
    if exist('H_correct','dir') == 0; mkdir('H_correct'); else end
    cd('H_correct')
    DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
    cd(DateString)
    
    for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
        for c = 1:size(FCR_cond,3) % COND LOOP
            h=figure;
            for i=1:size(FCR_cond,2) % TRIAL LOOP
                duration_ms=points(j)/s_rate(j)*1000;
                step=duration_ms/points(j);
                time_axis=0:step:duration_ms-step; time_axis=time_axis';
                ylim_up=squeeze(abs(min(min(min(min(mean_FCR_cond))))))+0.5;
                if c < 5 && c~=1
                    plot_start = 0; plot_end = 50;
                    M_start = 2+10; M_end = 12+10;
                    bar_start = 26; bar_end = 34;% dely + 10 ms preceeding the stimulus
                    % TEMPORARILY REMOVED THE -2 TO SEE IF THE DETECTION IS BETTER
                    sub_start = sub_delay+10; sub_end = sub_delay2+10;
                    stim_time= 10;
                else if c == 1
                        plot_start = 190; plot_end = 240;
                        M_start = 2+200; M_end = 12+200;
                        bar_start = 216; bar_end = 224;
                        % TEMPORARILY REMOVED THE -2 TO SEE IF THE DETECTION IS BETTER
                        sub_start = sub_delay+200; sub_end = sub_delay2+200;
                        stim_time= 200;
                    else
                        plot_start = 195; plot_end = 245;
                        M_start = 2+205; M_end = 12+205;
                        bar_start = 221; bar_end = 229;
                        % TEMPORARILY REMOVED THE -2 TO SEE IF THE DETECTION IS BETTER
                        sub_start = sub_delay+205; sub_end = sub_delay2+205;
                        stim_time= 205;
                    end
                end
                
                % PLOT
                subplot(2,5,i)
                plot(time_axis,FCR_cond(:,i,c,j));
                % M detected MAX
                m_peak_disp = (m_peak_delay+M_start_pt)*duration_ms/points(j);
                hold on; line([m_peak_disp m_peak_disp],[M_peak-1 M_peak+1],'Color',[1 0 1], 'Marker','*');
                
                % Personal H-r delays from Viking
                hold on; line([sub_start sub_start],[-ylim_up-1 ylim_up-1],'Color',[0 1 0]);
                hold on; line([sub_end sub_end],[-ylim_up-1 ylim_up-1],'Color',[0 1 0]);
                % H detected MAX
                h_peak_delay=H_resp_corr(i,c,j,2);
                H_peak=H_resp_corr(i,c,j,1);
                h_peak_disp = h_peak_delay*duration_ms/points(j);
                hold on; line([h_peak_disp h_peak_disp],[H_peak-0.5 H_peak+0.5],'Color',[1 0 1], 'Marker','*');
                peak_str=sprintf('%.3f',H_peak);
                hold on; text(h_peak_disp,H_peak-0.6,peak_str)
                
                set (gca,'Ydir','reverse')
                fig_str = [cond_name{c}, ' - ', int2str(i),' ', step_name{j}];
                title(fig_str);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                ylim([-ylim_up ylim_up]); xlim([plot_start plot_end]);
                hold on;
                %             pause;
            end
            save_fig=['DETECT','_',cond_name{c},'_', step_name{j},'.fig'];
            save_png=['DETECT','_',cond_name{c},'_', step_name{j},'.png'];
            print(h,save_png,'-dpng','-r300');
            savefig(save_fig);
            
            pause;
        end
        close all
    end
end

%% CALCULATE INHIBITION IN PERCENTAGE

% RATIO OVER SINGLE - H-RESPONSES
H_ratio = zeros(14,3,2); % Conditions, Steps
for c = 1:size(FCR_cond,3)-1 % COND LOOP
    for j = 1:size(FCR_cond,4) % STEP LOOP - baseline, post
        mean_H = mean(H_resp(:,c,j,5));
        mean_H_base = mean(H_resp(:,15,j,5));
        H_ratio(c,1,j) = mean_H;
        H_ratio(c,2,j) = mean_H_base;
        H_ratio(c,3,j) = abs(mean_H./mean_H_base);
    end
end
cd(ML_data_folder)
DateString = datestr(clock);
save(['H_ratio_',DateString,'.mat'],'H_ratio')

H_percent = squeeze(H_ratio(:,3,:))*100-100; % conditions, steps

Mean_Phase1 = mean(H_percent(1:5,:));
Mean_Phase2 = mean(H_percent(7:10,:));
Mean_Phase3 = mean(H_percent(11:14,:));

save('H_inhib', 'H_percent');
save('Mean_inhib_phases','Mean_Phase1','Mean_Phase2','Mean_Phase3'); 


%%------------------------------------------------------------------------
%%------------------------------------------------------------------------
%%------------------- ILLUSTRATIONS / GRAPHS -----------------------------
%%------------------------------------------------------------------------
%%------------------------------------------------------------------------

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


% end


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
