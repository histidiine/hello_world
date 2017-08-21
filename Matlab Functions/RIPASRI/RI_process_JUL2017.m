function RI_process_JUL2017(PAS_ID, gp_ID,sub_ID, subChoice)

%% INITIALIZATION
disp('1- 5 Hz');
disp('2- 0.2 Hz');
disp('3- 5 Hz only');
PAS_choice=0;
while PAS_choice~=1 && PAS_choice~=2 && PAS_choice~=3
    PAS_choice = input('What PAS protocole do you want to process (see choice above)?:   ');
end
if PAS_choice == 1
    PAS_ID = '5hz';
else if PAS_choice == 2
        PAS_ID = '02hz';
    else if PAS_choice ==3
            PAS_ID = '5hzonly';
        end
    end
end
disp(['You chose: ', PAS_ID]);

disp('1- HV');
disp('2- CRPS');
disp('3- CRPSD');
disp('4- FHD');
disp('5- FUNCT');
gp_choice=0;
while gp_choice~=1 && gp_choice~=2 && gp_choice~=3 && gp_choice~=4 && gp_choice~=5
    gp_choice = input('What group do you want to process (see choice above)?:   ');
end
if gp_choice == 1;
    gp_ID = 'HV';
else if gp_choice == 2;
        gp_ID = 'CRPS';
    else if gp_choice == 3;
            gp_ID = 'CRPSD';
        else if gp_choice == 4;
                gp_ID = 'FHD';
            else if gp_choice == 5;
                    gp_ID = 'FUNCT';
                end
            end
        end
    end
end
disp(['You chose: ', gp_ID]);

gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
log_folder = [gp_folder,'/', sub_ID];
data_folder = [log_folder,'/Extract/'];
fig_folder = [log_folder,'/Figures/RI'];
ML_data_folder = [log_folder,'/Data'];

cd(gp_folder)
list_folders=dir; % create a variable (struc) containing list of subjects to process
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
for k = 1:length(list_folders) % to print the list of subjects to remind the IDs to the user
    disp(list_folders(k).name)
end
sub_ID='none';
while ~exist(sub_ID,'dir') % Choosing subject ID = when in the main menu (RI_routine) single subject processing was chosen
    sub_ID = input('Please enter the ID of the subject you want to process into single quotes (see choices above):   ');
end

load hdelays; 
for k = 1:length(hdelays) % To get the H-reflex delays of the current subject
    temp(k)=strcmp(hdelays(k).names,sub_ID);
end
[osef,col]=max(temp); clear osef
sub_delays(1)=hdelays(col).delays(1);
sub_delays(2)=hdelays(col).delays(2);

disp('End of INITIALIZATION')

%% RI Pre-Processing: Log extraction, Data re-ordering by conditions
[RI_data,RI_cond,data_param] = RI_preproc(PAS_ID,gp_ID,sub_ID); % timepoints, muscle, trials (per cond), cond, step [5000, 2, 10, 15, 2]

%% CHECK SUBJECT H-REFLEX DELAYS

cd(fig_folder)

[subdelay2]= check_delays(RI_cond, sub_delays, data_param);
sub_delays(2)=sub_delay2;

%% CREATE AND SAVE GRAPHS OF RAW EMG SIGNAL

[badrecs, badrec_muscles] = check_badrec_RI(data_cond, data_param, hdelays);

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


end




