function EBCC_process(sub_id)

disp('Check the line 18 to see if folder is in the right group of subject')
pause;
sub_id=input('Please enter the subject ID used to name the folder in between simple quotes: ');


%% EXTRACTING INFO FROM LOG
log_folder=['/Users/histidiine/Google Drive/Capsule de sauvegarde/PhD/Matlab/Matlab Data/FUNCT/', sub_id];
nguyet_folger= [log_folder, '/Extract'];
cond_folder=[nguyet_folger,'/Cond'];
ext_folder=[nguyet_folger, '/Decond'];

cd(log_folder)
if exist('Fig','dir') == 0; mkdir('Fig'); else end
fig_folder = [log_folder, '/Fig'];
if exist('Data','dir') == 0; mkdir('Data'); else end
ML_data_folder = [log_folder, '/Data'];

cd(log_folder) 
log_cond_temp=importdata('cond.log');
log_cond = log_cond_temp.data;
log_decond_temp=importdata('decond.log');
log_decond = log_decond_temp.data;

s_rate=[log_cond(1),log_decond(1)]; % I could only take one as they should all be the same, but jic to avoid an error
points=[log_cond(2),log_decond(2)]; 
n_muscle=[log_cond(3),log_decond(3)];
n_trial=[log_cond(4),log_decond(4)];

duration_ms(1)=points(1)/s_rate(1)*1000;
duration_ms(2)=points(2)/s_rate(2)*1000;
step(1)=duration_ms(1)/points(1);
step(2)=duration_ms(2)/points(1);
time_axis(1,:)=0:step(1):duration_ms(1)-step(1); time_axis=time_axis';
temp=0:step(2):duration_ms(2)-step(2); time_axis(:,2)=temp;

% disp('Please verify that cond.log and decond.log were copied to excel files and the first three rows deleted and press any key')
% pause;
% 
% open cond.xlsx; pause;
% open decond.xlsx; pause;

clear log_cond_temp log_cond log_decond_temp log_decond

%% CONDITIONNING - Getting data from both eyes ---> cond_data [data, eye, trials]


% getting data from conditioning blocks
cd(cond_folder)
cond_data=zeros(points(1),n_muscle(1),n_trial(1));

for i = 1:n_trial(1)
   temp=importdata(int2str(i));
   cond_data(:,1,i) = temp(1:points(1),1); % I put row 1:8750 (points) in case there was a mistake during the nguyet and a file was recorded twice
   cond_data(:,2,i) = temp(1:points(1),2); 
   clear temp
end

% oeil_droit=squeeze(cond_data(:,1,:)); oeil_gauche=squeeze(cond_data(:,2,:)); 
% ms_OD = ricampiona(oeil_droit', time_ms); ms_OD = ms_OD';


%% ------------------------------------------------------------------------
%            Number of CR per block new technique (8.12.16) - tp
% -------------------------------------------------------------------------

cd(log_folder)
oeil_droit=squeeze(cond_data(:,1,:));
oeil_gauche=squeeze(cond_data(:,2,:));

[CR_block, extra_data] = CRwDuration(oeil_droit);
[CR_block, extra_data] = CRwDuration(oeil_gauche);

cd(fig_folder) 
DateString = datestr(clock); mkdir(DateString); % Creates a folder to not overwrite previous graphs saved
date_folder = [fig_folder, '/', DateString];
cd(date_folder)
% check_CR_detection(oeil_droit,extra_data,time_axis(:,1));

%% CHECK CR DETECTION

noise1=1237; noise1_ms=noise1*duration_ms(1)/points(1);
noise2=2487; noise2_ms=noise2*duration_ms(1)/points(1);
cr1=3750; cr1_ms=cr1*duration_ms(1)/points(1);
cr2=5000; cr2_ms=cr2*duration_ms(1)/points(1);
dur=306; dur_ms=dur*duration_ms(1)/points(1);
dur_pre=94; dur_pre_ms=dur_pre*duration_ms(1)/points(1);
       
resp=input('Do you want to see all trials ? (char: y/n)  ');
if strcmp(resp, 'y')
    
    
    %% Visualizer all trials per block
    for j = 1:7 % LOOP FOR BLOCKS
        h=figure; % One figure per block
        set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
        tlstr=['Block number ', int2str(j)];
        for i=1:10 % LOOP FOR TRIALS
            
            indice_take = j*10+i; % in order not to take first 10 trials
            indice_save = (j-1)*10+i; % in order to save starting at 1
            ylim_up = max(max(oeil_droit(1:cr2)))+0.05;
            
            if size(oeil_droit,1)==8750
                dur_value=extra_data(indice_save,6)./6.25;
            else
                dur_value=extra_data(indice_save,6);
            end
            str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
            subplot (5,2,i); plot(time_axis,oeil_droit(:,indice_take), 'Color','b'); title(str);
            set (gca,'Ydir','reverse')
            hold on; line([cr1_ms cr1_ms],[-0.5 0.5],'Color',[1 0 0]);
            hold on; line([cr2_ms cr2_ms],[-0.5 0.5],'Color',[1 0 0]);
            hold on;
            if extra_data(indice_save,3)~=0 % only for trials where a CR was detected
                dur_txt=num2str(dur_value);
                line([extra_data(indice_save,4)*duration_ms(1)/points(1) extra_data(indice_save,5)*duration_ms(1)/points(1)],[-0.09 -0.09],'Color',[0 1 0]); % put a red line under the CR
                descr = [dur_txt, ' ms'];
                text(extra_data(indice_save,4)*duration_ms(1)/points(1),-0.12,descr)
            end
            hold on;
            ylim([-ylim_up ylim_up]);
        end
        kk=1;fig_str = tlstr;
        while 7 == exist(tlstr, 'file'); fig_str = [tlstr,int2str(kk)]; end
        print(h,fig_str,'-dpng','-r600');
        save_fig=[fig_str,'.fig'];
        savefig(save_fig);
%         pause
        close
    end
else
    disp('Moving to next step without visualizing all trials')
end

resp=input('Do you want to visualize the trials with detected CR only ? (char: y/n)  ');
if strcmp(resp, 'y')
    %% Visualize only the trials with detected CR
    for j = 1:7 % LOOP FOR BLOCKS
        for i=1:10 % LOOP FOR TRIALS
            indice_take = j*10+i; % in order not to take first 10 trials
            indice_save = (j-1)*10+i; % in order to save starting at 1
            ylim_up = max(oeil_droit(:,i))+0.2;
            
            detect_start = extra_data(indice_save,4)*duration_ms(1)/points(1);
            detect_end = extra_data(indice_save,5)*duration_ms(1)/points(1); 
            cr_duration = extra_data(indice_save,6)*duration_ms(1)/points(1); ;
            
            if extra_data(indice_save,3)==0 % only for trials where a CR was detected
            else
                if size(raw_data,1)==8750
                    dur_value=extra_data(indice_save,6)./6.25;
                else
                    dur_value=extra_data(indice_save,6);
                end
                h=figure; % One figure per trial where a CR was detected
                set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
                tlstr=['Block number ', int2str(j),' - Trial number ', int2str(i)];
                % ------ NOISE -------
                subplot (1,2,1); plot(time_axis,oeil_droit(:,indice_take));
                ylim([-ylim_up ylim_up]); xlim([noise1_ms-50 noise2_ms+50])
                descr=[num2str(extra_data(indice_save,1)*1000), ' uV']; title('noise');
                text(900, -0.11, descr);
                hold on;
                line([noise1_ms noise1_ms],[-0.05 0.05],'Color',[1 0 0]); % Vertical Line at start
                hold on;
                line([noise2_ms noise2_ms],[-0.05 0.05],'Color',[1 0 0]); % Vertical Line at end
                
                % ------ ZOOM CR --------
                subplot (1,2,2); plot(time_axis,oeil_droit(:,indice_take)); title(tlstr);
                ylim([-ylim_up ylim_up]); xlim([detect_start-50 detect_end+50]);
                line_length=cr_duration;
                % ----- CR lines ----
                hold on;
                line([detect_start detect_start],[-0.05 0.05],'Color',[1 0 0]); % Vertical Line at start
                hold on;
                line([line_length line_length],[-0.05 0.05],'Color',[1 0 0]); % Vertical Line at end
                
                hold on;
                line([detect_start line_length],[-0.07 -0.07],'Color',[0 1 0]); % put a red line under the CR
                descr = [num2str(dur_value), ' ms'];
                text(500,-0.09,descr);
                descr2=[num2str(extra_data(indice_save,2)*1000), ' uV'];
                text(500,-0.11,descr2);
                kk=1;fig_str = tlstr;
                while 7 == exist(tlstr, 'file'); fig_str = [tlstr,int2str(kk)]; end
                print(h,fig_str,'-dpng');
                pause
                close
            end
        end
    end
else
    disp('End of check Cr Detection protocol')
end


cd(main_dossier)
if 7 == exist('ML_data','dir'); else mkdir('ML_data');end
cd(data_folder)
save CR_block; save extra_data
xlswrite([sub_id, '_ExtraData'],extra_data); xlswrite([sub_id, '_CR_Block'],CR_block);
xlswrite([sub_id, '_OD_data'],oeil_droit);
% Add somewhere to save with xlswrite


h2=figure; plot(CR_block); title('Learning Curve'); ylim([0 10]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
cd(ML_data_folder)
print(h2,'learning curve_10then40','-dpng');


%% CONDITIONNING - Counting alpha blinks and CRs

%-------------------------------------------------
% Number of CR per block with 50ms - 10 then 40
%-------------------------------------------------

% CR_block = zeros(7,1); CR=0; noise_resp=zeros(70,3);
% for j = 1:7 % without the first 10 trials / LOOP OF BLOCKS
%     for i=  1:10 % LOOP OF TRIALS
%         indice_take = j*10+i; % in order not to take first 10 trials
%         indice_save = (j-1)*10+i; % in order to save starting at 1 
%         for k = 598:748 % LOOP OF TIME CHUNKS
%             ref= mean(abs(ms_OD(198:248,indice_take)),1); % ref for noise 200 ms in the first 400
%             test=mean(abs(ms_OD(k:k+9,indice_take)),1); % Chunks of 50 ms
%             test2=mean(abs(ms_OD(k+10:k+49,indice_take)),1);
%             noise_resp(indice_save,1) = ref; noise_resp(indice_save,2) = test;
%             if test > 2*ref % if 200 ms before US activity is bigger than 1.5*noise = CR
%                 if test2 > 2*ref
%                 CR=CR+1;
%                 noise_resp(indice_save,3) = k;
%                 double_check(i,j)=k;
%                 break
%                 end
%             else
%                 noise_resp(indice_save,3) = 0;
%             end
%         end
%         for k = 598:748
%             ref= mean(abs(ms_OD(198:298,indice_take)),1); % ref for noise 200 ms in the first 400
%             alpha_test=mean(abs(ms_OD(k-200:k+9-200,indice_take)),1); 
%             alpha_test2=mean(abs(ms_OD(k+10-200:k+19-200,indice_take)),1); 
%             if alpha_test > 2*ref %% Alpha blinks
%                 if alpha_test2 > 2*ref
%                 noise_resp(indice_save,4) = k-200;
%                 break
%                 end
%             else
%                 noise_resp(indice_save,4) = 0;
%             end
%         end
%     end
%     CR_block(j,1)=CR;
%     CR=0;
% end

% cd('/Users/aureliestephan/Documents/MATLAB/Data/Archive')
% load oeil_droit

%% ------------------------------------------------------------------------
%            Number of CR per block new technique (8.12.16) - ms
% -------------------------------------------------------------------------


CR_block = zeros(7,1); % = nb of CR
CR=0;
noise_resp=zeros(70,5);%1-noise, 2-resp, 3-start, 4-end, 5-duration 
Peak=zeros(10,7); Onset=zeros(10,7)


for j = 1:7 % without the first 10 trials / LOOP OF BLOCKS
    for i=  1:10 % LOOP OF TRIALS
        indice_take = j*10+i; % in order not to take first 10 trials
        indice_save = (j-1)*10+i; % in order to save starting at 1
        search='on';k=598;
        while strcmp(search, 'on') && k < 749 % LOOP OF TIME CHUNKS
            
            %___________INITIALIZATION_________________
            ref= mean(abs(ms_OD(198:398,indice_take)),1); % ref for noise 200 ms in the first 400
            test=mean(abs(ms_OD(k:k+49,indice_take)),1); % Chunks of 50 ms
            noise_resp(indice_save,1) = ref; noise_resp(indice_save,2) = test;
            
            %-------------DETECTING CR---------------------
            if test > 2*ref % if 200 ms before US activity is bigger than 2*noise = CR
                start_search='on'; kk = k;
                while strcmp(start_search,'on') && kk < (k+50) % finding beginning
                    test_start=mean(abs(ms_OD(kk:kk+4,indice_take)),1);
                    if test_start > 2*ref
                        timepoint_start=kk;
                        start_search='off';
                    end
                    kk=kk+1;
                end
                end_search='on';
                if start_search == 'off'
                    while strcmp(end_search, 'off') && ll < (k+50) % finding end
                        test_end=mean(abs(ms_OD(k+49-(ll+4):k+49-(ll),indice_take)),1);
                        if test_end > 2*ref
                            timepoint_end=ll;
                            end_search='off';
                            search='off';
                            CR_duration=timepoint_end-timepoint_start;
                        end
                        ll=ll+1;
                    end
                end
            end
        end
        noise_resp(indice_save,3)=timepoint_start;
        noise_resp(indice_save,4)=timepoint_end;
        noise_resp(indice_save,5)=CR_duration;
        CR_block(j,1)=CR;
        CR=0;
        k=k+1;
    end
end





%%

%-------------------------------------------------
% Number of CR per block with 50ms - at once
%-------------------------------------------------

CR_block = zeros(7,1); CR=0; noise_resp=zeros(70,3); onset=0; 
Peak=zeros(10,7); Onset=zeros(10,7)
for j = 1:7 % without the first 10 trials / LOOP OF BLOCKS
    for i=  1:10 % LOOP OF TRIALS
        indice_take = j*10+i; % in order not to take first 10 trials
        indice_save = (j-1)*10+i; % in order to save starting at 1 
        for k = 598:748 % LOOP OF TIME BINS
            
            %___________INITIALIZATION_________________
            ref= mean(abs(ms_OD(198:398,indice_take)),1); % ref for noise 200 ms in the first 400
            test=mean(abs(ms_OD(k:k+49,indice_take)),1); % Chunks of 50 ms
            test2=mean(abs(ms_OD(k:k+9,indice_take)),1);
            test3=mean(abs(ms_OD(k+10:k+19,indice_take)),1);
            test4=mean(abs(ms_OD(k+20:k+29,indice_take)),1);
            test5=mean(abs(ms_OD(k+30:k+39,indice_take)),1);
            noise_resp(indice_save,1) = ref; noise_resp(indice_save,2) = test;
            
            %-------------DETECTING CR---------------------
            if test > 2*ref % if 200 ms before US activity is bigger than 2*noise = CR
                CR=CR+1;
                noise_resp(indice_save,3) = k;
                
                %________TESTS for time to onset_______
                if test2 > 2*ref
                    onset=k;
                    Onset(i,j)=onset; clear onset
                else if test3 > 2*ref
                        onset=k+10;
                        Onset(i,j)=onset; clear onset
                    else if test4 > 2*ref
                            onset=k+20;
                            Onset(i,j)=onset; clear onset
                        else if test5 > 2*ref
                                   onset=k+30;
                                   Onset(i,j)=onset; clear onset
                            else
                                onset=k;
                                Onset(i,j)=onset; clear onset
                            end
                        end
                    end
                end
                
                %____________TIME TO PEAK___________
                [max_CR(i,j), ind_max]=max(abs(ms_OD(k:798,indice_take)));
                Peak(i,j) = ind_max+k; clear ind_max
                break
               % end
            else
                noise_resp(indice_save,3) = 0;
            end
        end
        for k = 598:748
            ref= mean(abs(ms_OD(198:298,indice_take)),1); % ref for noise 200 ms in the first 400
            alpha_test=mean(abs(ms_OD(k-200:k+9-200,indice_take)),1); 
            alpha_test2=mean(abs(ms_OD(k+10-200:k+19-200,indice_take)),1); 
            if alpha_test > 2*ref %% Alpha blinks
                if alpha_test2 > 2*ref
                noise_resp(indice_save,4) = k-200;
                break
                end
            else
                noise_resp(indice_save,4) = 0;
            end
        end
    end
    CR_block(j,1)=CR;
    CR=0;
end


%---------------------------------------------
%  Raw data conditioning plots of 50 ms loop
%---------------------------------------------

% cd(fig_folder)
for j = 1:7 % LOOP FOR BLOCKS
   h=figure; % One figure per block
   set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
   tlstr=['Block number ', int2str(j), '_withpeakandonset'];
   for i=1:10 % LOOP FOR TRIALS
       indice_take = j*10+i; % in order not to take first 10 trials
       indice_save = (j-1)*10+i; % in order to save starting at 1 
       str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
       subplot (5,2,i); plot(ms_OD(:,indice_take)); title(str);
       hold on; 
       if noise_resp(indice_save,3)==0 % only for trials where a CR was detected
       else
           line([noise_resp(indice_save,3) noise_resp(indice_save,3)+50],[-0.1 -0.1],'Color',[1 0 0]); % put a red line under the CR
           %hold on; line([Onset(i,j) Onset(i,j)],[-0.2 0.2],'Color',[0 0 0]); % TIME TO ONSET
%            if Peak(i,j) == 0
%            else
%            hold on; line([Peak(i,j) Peak(i,j)],[-0.2 0.2],'Color',[1 0 0]); % TIME TO PEAK
%            end
       end
       hold on;
       if noise_resp(indice_save,4)==0 % only for trials where a alpha blink detected
       else
           line([noise_resp(indice_save,4) noise_resp(indice_save,4)+50],[-0.1 -0.1],'Color',[0 1 0]); % put a green line under the alpha blink
       end
       ylim([-0.2 0.2]);
   end
   if 7 == exist(tlstr, 'file'); fig_str = [tlstr,'_2']; else fig_str = tlstr; end
   print(h,fig_str,'-dpng');
   pause
end
close all


%____________ TIME TO ONSET AND PEAK____________________
for j = 1:7 % LOOP FOR BLOCKS
   for i=1:10 % LOOP FOR TRIALS
       indice_take = j*10+i; % in order not to take first 10 trials
       indice_save = (j-1)*10+i; % in order to save starting at 1
       str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
       if noise_resp(indice_save,3)==0 % only for trials where a CR was detected
       else
           h=figure; % One figure per block
           set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
           tlstr=['Block', int2str(j), ' trial ', int2str(i), '_peakandonset'];
           plot(ms_OD(:,indice_take)); title(str);
           hold on; line([Onset(i,j) Onset(i,j)],[-0.12 0.12],'Color',[0 0 0]); % TIME TO ONSET
           if Peak(i,j) == 0
           else
               hold on; line([Peak(i,j) Peak(i,j)],[-0.12 0.12],'Color',[1 0 0]); % TIME TO PEAK
           end
           ylim([-0.2 0.2]);
           pause
           if 7 == exist(tlstr, 'file'); fig_str = [tlstr,'_2']; else fig_str = tlstr; end
           print(h,fig_str,'-dpng');
       end
       close
   end
end
close all

%-------------------------------------------------------------------------
% Number of CR per block with 50ms - at once + check onset and end then
%                               duration
%--------------------------------------------------------------------------

CR_block = zeros(7,1); CR=0; noise_resp=zeros(70,3); onset=0; 
Peak=zeros(10,7); Onset=zeros(10,7)
for j = 1:7 % without the first 10 trials / LOOP OF BLOCKS
    for i=  1:10 % LOOP OF TRIALS
        indice_take = j*10+i; % in order not to take first 10 trials
        indice_save = (j-1)*10+i; % in order to save starting at 1 
        for k = 598:748 % LOOP OF TIME CHUNKS
            
            %___________INITIALIZATION_________________
            ref= mean(abs(ms_OD(198:398,indice_take)),1); % ref for noise 200 ms in the first 400
            test=mean(abs(ms_OD(k:k+49,indice_take)),1); % Chunks of 50 ms
            test2=mean(abs(ms_OD(k:k+9,indice_take)),1);
            test3=mean(abs(ms_OD(k+10:k+19,indice_take)),1);
            test4=mean(abs(ms_OD(k+20:k+29,indice_take)),1);
            test5=mean(abs(ms_OD(k+30:k+39,indice_take)),1);
            noise_resp(indice_save,1) = ref; noise_resp(indice_save,2) = test;
            
            %-------------DETECTING CR---------------------
            if test > 2*ref % if 200 ms before US activity is bigger than 2*noise = CR
                CR=CR+1;
                noise_resp(indice_save,3) = k;
                
                %________TESTS for time to onset_______
                if test2 > 2*ref
                    onset=k;
                    Onset(i,j)=onset; clear onset
                else if test3 > 2*ref
                        onset=k+10;
                        Onset(i,j)=onset; clear onset
                    else if test4 > 2*ref
                            onset=k+20;
                            Onset(i,j)=onset; clear onset
                        else if test5 > 2*ref
                                   onset=k+30;
                                   Onset(i,j)=onset; clear onset
                            else
                                onset=k;
                                Onset(i,j)=onset; clear onset
                            end
                        end
                    end
                end
                
                %____________TIME TO PEAK___________
                [max_CR(i,j), ind_max]=max(abs(ms_OD(k:798,indice_take)));
                Peak(i,j) = ind_max+k; clear ind_max
                break
               % end
            else
                noise_resp(indice_save,3) = 0;
            end
        end
        for k = 598:748
            ref= mean(abs(ms_OD(198:298,indice_take)),1); % ref for noise 200 ms in the first 400
            alpha_test=mean(abs(ms_OD(k-200:k+9-200,indice_take)),1); 
            alpha_test2=mean(abs(ms_OD(k+10-200:k+19-200,indice_take)),1); 
            if alpha_test > 2*ref %% Alpha blinks
                if alpha_test2 > 2*ref
                noise_resp(indice_save,4) = k-200;
                break
                end
            else
                noise_resp(indice_save,4) = 0;
            end
        end
    end
    CR_block(j,1)=CR;
    CR=0;
end


%---------------------------------------------
%  Raw data conditioning plots of 50 ms loop
%---------------------------------------------

% cd(fig_folder)
for j = 1:7 % LOOP FOR BLOCKS
   h=figure; % One figure per block
   set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
   tlstr=['Block number ', int2str(j), '_withpeakandonset'];
   for i=1:10 % LOOP FOR TRIALS
       indice_take = j*10+i; % in order not to take first 10 trials
       indice_save = (j-1)*10+i; % in order to save starting at 1 
       str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
       subplot (5,2,i); plot(ms_OD(:,indice_take)); title(str);
       hold on; 
       if noise_resp(indice_save,3)==0 % only for trials where a CR was detected
       else
           line([noise_resp(indice_save,3) noise_resp(indice_save,3)+50],[-0.1 -0.1],'Color',[1 0 0]); % put a red line under the CR
           %hold on; line([Onset(i,j) Onset(i,j)],[-0.2 0.2],'Color',[0 0 0]); % TIME TO ONSET
%            if Peak(i,j) == 0
%            else
%            hold on; line([Peak(i,j) Peak(i,j)],[-0.2 0.2],'Color',[1 0 0]); % TIME TO PEAK
%            end
       end
       hold on;
       if noise_resp(indice_save,4)==0 % only for trials where a alpha blink detected
       else
           line([noise_resp(indice_save,4) noise_resp(indice_save,4)+50],[-0.1 -0.1],'Color',[0 1 0]); % put a green line under the alpha blink
       end
       ylim([-0.2 0.2]);
   end
   if 7 == exist(tlstr, 'file'); fig_str = [tlstr,'_2']; else fig_str = tlstr; end
   print(h,fig_str,'-dpng');
   pause
end
close all


%____________ TIME TO ONSET AND PEAK PLOTS____________________
for j = 1:7 % LOOP FOR BLOCKS
   for i=1:10 % LOOP FOR TRIALS
       indice_take = j*10+i; % in order not to take first 10 trials
       indice_save = (j-1)*10+i; % in order to save starting at 1
       str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
       if noise_resp(indice_save,3)==0 % only for trials where a CR was detected
       else
           h=figure; % One figure per block
           set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
           tlstr=['Block', int2str(j), ' trial ', int2str(i), '_peakandonset'];
           plot(ms_OD(:,indice_take)); title(str);
           hold on; line([Onset(i,j) Onset(i,j)],[-0.12 0.12],'Color',[0 0 0]); % TIME TO ONSET
           if Peak(i,j) == 0
           else
               hold on; line([Peak(i,j) Peak(i,j)],[-0.12 0.12],'Color',[1 0 0]); % TIME TO PEAK
           end
           ylim([-0.2 0.2]);
           pause
           if 7 == exist(tlstr, 'file'); fig_str = [tlstr,'_2']; else fig_str = tlstr; end
           print(h,fig_str,'-dpng');
       end
       close
   end
end
close all

h2=figure; plot(CR_block); title('Learning Curve'); ylim([0 10]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
print(h2,'learning curve_10then40','-dpng');


%% EXTINCTION

cd(ext_folder)
ext_trial = 20; ext_data=zeros(points,2,20);

for i = 1:ext_trial
   temp=importdata(int2str(i));
   ext_data(:,1,i) = temp(1:points,1); 
   ext_data(:,2,i) = temp(1:points,2); clear temp
end
ext_OD=squeeze(ext_data(:,1,:)); ext_OG=squeeze(ext_data(:,2,:)); 
ms_ext_OD = ricampiona(ext_OD', time_ms); ms_ext_OD = ms_ext_OD';

for j = 1:2
   h=figure; 
   for i=1:10
       str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
       subplot (5,2,i); plot(ms_ext_OD(:,j*10+i-10)); title(str);
       ylim([-0.2 0.2]);
   end
   fig_str = ['Extinction Block number ', int2str(j)];
   cd(ML_data_folder); print(h,fig_str,'-dpng');
end

%-----------------------------------------------------------------------
% Counting CR - 200 to 600 ms after tone (600-1000 ms) - Teo et al, 2008
%-----------------------------------------------------------------------
CR=0; CR_ext_block = zeros(2,1);
for j = 1:2
    for i=  1:10
        ref= mean(abs(ms_ext_OD(1:399,(j-1)*10+i)),1); % ref for noise = first 400 ms
        test=mean(abs(ms_ext_OD(601:999,(j-1)*10+i)),1); % 200 to 600 ms after sound
        noise_resp_ext((j-1)*10+i,1) = ref; noise_resp_ext((j-1)*10+i,2) = test;
        if test > 1.5*ref 
            CR=CR+1;
            noise_resp_ext((j-1)*10+i,3) = 1;
        else
            noise_resp_ext((j-1)*10+i,3) = 0;
        end
    end
    CR_ext_block(j,1)=CR; CR=0;
end
h2=figure; plot(CR_ext_block); title('Extinction Curve'); ylim([0 10]);
cd(ML_data_folder); print(h2,'extinction curve','-dpng');

%% SAVING 
% Saving matlab matrices in ML_data's subject folder
cd(main_dossier)
if 7 == exist('ML_data','dir'); cd(data_folder);else mkdir('ML_data'); cd(data_folder);end

xlswrite('CR_block',CR_block); xlswrite('noise_resp',noise_resp); 
xlswrite('CR_ext_block',CR_ext_block); xlswrite('noise_resp_ext',noise_resp_ext); 

% Save all matlab data
cd(data_folder); 
save('oeil_droit', 'oeil_droit'); save('oeil_gauche', 'oeil_gauche'); save('CR_block', 'CR_block'); save('noise_resp', 'noise_resp'); save('ms_OD', 'ms_OD');
save('CR_ext_block','CR_ext_block');save('noise_resp_ext','noise_resp_ext');save('ms_ext_OD','ms_ext_OD');
disp(['Finished ', sub_id]);




%% ----------------------------------------------------------------------
% -----------------------------------------------------------------------
%                                  CEMETERY
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------


%---------------------------------------
%  Number of CRs per block - classic 2*
%---------------------------------------

% CR_block = zeros(7,1); CR=0; noise_resp=zeros(70,3);
% for j = 1:7 % without the first 10 trials / LOOP OF BLOCKS
%     for i=  1:10 % LOOP OF TRIALS
%         indice_take = j*10+i; % in order not to take first 10 trials
%         indice_save = (j-1)*10+i; % in order to save starting at 1 
%         ref= mean(abs(ms_OD(198:398,indice_take)),1); % ref for noise = first 400 ms
%         test=mean(abs(ms_OD(598:798,indice_take)),1); % 200 ms before US
%         noise_resp(indice_save,1) = ref; noise_resp(indice_save,2) = test;
%         if test > 2*ref % if 200 ms before US activity is bigger than 2*noise = CR
%             CR=CR+1;
%             noise_resp(indice_save,3) = 1;
%         else
%             noise_resp(indice_save,3) = 0;
%         end
%     end
%     CR_block(j,1)=CR;
%     CR=0;
% end
% 
% %-------------------------------------------------
% %  Raw data conditioning plots without 50 ms loop
% %-------------------------------------------------
% 
% % cd(fig_folder)
% for j = 1:7 % LOOP FOR BLOCKS
%    h=figure; % One figure per block
%    set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
%    tl_str=['Block number ', int2str(j), '_2fois'];
%    for i=1:10 % LOOP FOR TRIALS
%        indice_take = j*10+i; % in order not to take first 10 trials
%        indice_save = (j-1)*10+i; % in order to save starting at 1 
%        str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
%        subplot (5,2,i); plot(ms_OD(:,indice_take)); title(str);
%        hold on; 
%        if noise_resp(indice_save,3)==0 % only for trials where a CR was detected
%        else
%            line([599 799],[-0.1 -0.1],'Color',[1 0 0]); % put a red line under the CR
%        end
%        ylim([-0.2 0.2]);
%    end
%    if 7 == exist(tl_str, 'file'); fig_str = [tl_str,'_2']; else fig_str = tl_str; end
%    print(h,fig_str,'-dpng');
% end
% close all

%---------------------------------------
%       Number of CRs per block - 1.5*
%---------------------------------------

% CR_block = zeros(7,1); CR=0; noise_resp=zeros(70,3);
% for j = 1:7 % without the first 10 trials / LOOP OF BLOCKS
%     for i=  1:10 % LOOP OF TRIALS
%         indice_take = j*10+i; % in order not to take first 10 trials
%         indice_save = (j-1)*10+i; % in order to save starting at 1 
%         ref= mean(abs(ms_OD(198:398,indice_take)),1); % ref for noise = first 400 ms
%         test=mean(abs(ms_OD(598:798,indice_take)),1); % 200 ms before US
%         noise_resp(indice_save,1) = ref; noise_resp(indice_save,2) = test;
%         if test > 1.5*ref % if 200 ms before US activity is bigger than 2*noise = CR
%             CR=CR+1;
%             noise_resp(indice_save,3) = 1;
%         else
%             noise_resp(indice_save,3) = 0;
%         end
%     end
%     CR_block(j,1)=CR;
%     CR=0;
% end
% 
% %-------------------------------------------------
% %  Raw data conditioning plots without 50 ms loop
% %-------------------------------------------------
% 
% % cd(fig_folder)
% for j = 1:7 % LOOP FOR BLOCKS
%    h=figure; % One figure per block
%    set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
%    tl_str=['Block number ', int2str(j),'1point5rms'];
%    for i=1:10 % LOOP FOR TRIALS
%        indice_take = j*10+i; % in order not to take first 10 trials
%        indice_save = (j-1)*10+i; % in order to save starting at 1 
%        str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
%        subplot (5,2,i); plot(ms_OD(:,indice_take)); title(str);
%        hold on; 
%        if noise_resp(indice_save,3)==0 % only for trials where a CR was detected
%        else
%            line([599 799],[-0.1 -0.1],'Color',[1 0 0]); % put a red line under the CR
%        end
%        ylim([-0.2 0.2]);
%    end
%    if 7 == exist(tl_str, 'file'); fig_str = [tl_str,'_2']; else fig_str = tl_str; end
%    print(h,fig_str,'-dpng');
% end
% close all

%---------------------------------------
%       Number of CRs per block - sqrt
%---------------------------------------

% CR_block = zeros(7,1); CR=0; noise_resp=zeros(70,3);
% for j = 1:7 % without the first 10 trials / LOOP OF BLOCKS
%     for i=  1:10 % LOOP OF TRIALS
%         indice_take = j*10+i; % in order not to take first 10 trials
%         indice_save = (j-1)*10+i; % in order to save starting at 1 
%         ref= sqrt(mean(ms_OD(198:398,indice_take).^2,1)); % ref for noise = first 400 ms
%         test=sqrt(mean(ms_OD(598:798,indice_take).^2,1)); % 200 ms before US
%         noise_resp(indice_save,1) = ref; noise_resp(indice_save,2) = test;
%         if test > 2*ref % if 200 ms before US activity is bigger than 2*noise = CR
%             CR=CR+1;
%             noise_resp(indice_save,3) = 1;
%         else
%             noise_resp(indice_save,3) = 0;
%         end
%     end
%     CR_block(j,1)=CR;
%     CR=0;
% end

%-------------------------------------------------------
%  Raw data conditioning plots without 50 ms loop - sqrt
% %-------------------------------------------------------
% 
% % cd(fig_folder)
% for j = 1:7 % LOOP FOR BLOCKS
%    h=figure; % One figure per block
%    set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
%    tl_str=['Block number ', int2str(j), '_sqrt'];
%    for i=1:10 % LOOP FOR TRIALS
%        indice_take = j*10+i; % in order not to take first 10 trials
%        indice_save = (j-1)*10+i; % in order to save starting at 1 
%        str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
%        subplot (5,2,i); plot(ms_OD(:,indice_take)); title(str);
%        hold on; 
%        if noise_resp(indice_save,3)==0 % only for trials where a CR was detected
%        else
%            line([599 799],[-0.1 -0.1],'Color',[1 0 0]); % put a red line under the CR
%        end
%        ylim([-0.2 0.2]);
%    end
%    if 7 == exist(tl_str, 'file'); fig_str = [tl_str,'_2']; else fig_str = tl_str; end
%    print(h,fig_str,'-dpng');
% end
% close all



%-------------------------------------------------
% Number of CR per block with 50ms by 10 ms chunck
%-------------------------------------------------
% 
% CR_block = zeros(7,1); CR=0; noise_resp=zeros(70,3);
% for j = 1:7 % without the first 10 trials / LOOP OF BLOCKS
%     for i=  1:10 % LOOP OF TRIALS
%         indice_take = j*10+i; % in order not to take first 10 trials
%         indice_save = (j-1)*10+i; % in order to save starting at 1 
%         for k = 598:748 % LOOP OF TIME CHUNKS
%             ref= mean(abs(ms_OD(198:298,indice_take)),1); % ref for noise 200 ms in the first 400
%             test=mean(abs(ms_OD(k:k+9,indice_take)),1); % Chunks of 50 ms
%             noise_resp(indice_save,1) = ref; noise_resp(indice_save,2) = test;
%             if test > 2*ref % if 200 ms before US activity is bigger than 1.5*noise = CR
%                 if mean(abs(ms_OD(k+10:k+19,indice_take)),1) > 2*ref & mean(abs(ms_OD(k+20:k+29,indice_take)),1) > 2*ref...
%                         & mean(abs(ms_OD(k+30:k+39,indice_take)),1) > 2*ref & mean(abs(ms_OD(k+40:k+49,indice_take)),1) > 2*ref
%                 CR=CR+1;
%                 noise_resp(indice_save,3) = 1;
%                 double_check(i,j)=k;
%                 break
%                 end
%             else
%                 noise_resp(indice_save,3) = 0;
%             end
%         end
%     end
%     CR_block(j,1)=CR;
%     CR=0;
% end
% 
% %---------------------------------------------
% %  Raw data conditioning plots of 50 ms loop
% %---------------------------------------------
% 
% % cd(fig_folder)
% for j = 1:7 % LOOP FOR BLOCKS
%    h=figure; % One figure per block
%    set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
%    tlstr=['Block number ', int2str(j), '_with50msin10'];
%    for i=1:10 % LOOP FOR TRIALS
%        indice_take = j*10+i; % in order not to take first 10 trials
%        indice_save = (j-1)*10+i; % in order to save starting at 1 
%        str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
%        subplot (5,2,i); plot(ms_OD(:,indice_take)); title(str);
%        hold on; 
%        if double_check(i,j)==0 % only for trials where a CR was detected
%        else
%            line([double_check(i,j) double_check(i,j)+50],[-0.1 -0.1],'Color',[1 0 0]); % put a red line under the CR
%        end
%        ylim([-0.2 0.2]);
%    end
%    if 7 == exist(tlstr, 'file'); fig_str = [tlstr,'_2']; else fig_str = tlstr; end
%    print(h,fig_str,'-dpng');
% end
% close all

%-------------------------------------------------
% Number of CR per block with 50ms by 10 ms chunck
%-------------------------------------------------
% load ms_OD
% CR_block = zeros(7,1); CR=0; noise_resp=zeros(70,7);
% for j = 1:7 % without the first 10 trials / LOOP OF BLOCKS
%     for i=  1:10 % LOOP OF TRIALS
%         indice_take = j*10+i; % in order not to take first 10 trials
%         indice_save = (j-1)*10+i; % in order to save starting at 1 
%         for k = 598:748 % LOOP OF TIME CHUNKS
%             ref= mean(abs(ms_OD(198:298,indice_take)),1); % ref for noise 200 ms in the first 400
%             test=mean(abs(ms_OD(k:k+9,indice_take)),1); % Chunks of 50 ms
%             noise_resp(indice_save,1) = ref; noise_resp(indice_save,2) = test;
%             if test > 2*ref % if 200 ms before US activity is bigger than 2*noise = CR
%                 if mean(abs(ms_OD(k+10:k+19,indice_take)),1) > 2*ref & mean(abs(ms_OD(k+20:k+29,indice_take)),1) > 2*ref...
%                         & mean(abs(ms_OD(k+30:k+39,indice_take)),1) > 2*ref & mean(abs(ms_OD(k+40:k+49,indice_take)),1) > 2*ref
%                     CR=CR+1;
%                     noise_resp(indice_save,3) = k;
%                     noise_resp(indice_save,5) = mean(abs(ms_OD(k+10:k+19,indice_take)),1);
%                     noise_resp(indice_save,6) = mean(abs(ms_OD(k+20:k+29,indice_take)),1);
%                     noise_resp(indice_save,7) = mean(abs(ms_OD(k+30:k+39,indice_take)),1);
%                     double_check(i,j)=k;
%                     break
%                 end
%             else
%                 noise_resp(indice_save,3) = 0;
%                 noise_resp(indice_save,5) = mean(abs(ms_OD(k+10:k+19,indice_take)),1);
%                 noise_resp(indice_save,6) = mean(abs(ms_OD(k+20:k+29,indice_take)),1);
%                 noise_resp(indice_save,7) = mean(abs(ms_OD(k+30:k+39,indice_take)),1);
%             end
%         end
%         for k = 598:748
%             ref= mean(abs(ms_OD(198:298,indice_take)),1); % ref for noise 200 ms in the first 400
%             alpha_test=mean(abs(ms_OD(k-200:k+9-200,indice_take)),1); % Chunks of 50 ms
%             if alpha_test > 2*ref %% Alpha blinks
%                 if mean(abs(ms_OD(k+10-200:k+19-200,indice_take)),1) > 2*ref & mean(abs(ms_OD(k+20-200:k+29-200,indice_take)),1) > 2*ref
%                 noise_resp(indice_save,4) = k-200;
%                 break
%                 end
%             else
%                 noise_resp(indice_save,4) = 0;
%             end
%         end
%     end
%     CR_block(j,1)=CR;
%     CR=0;
% end
% 
% %-------------------------------------------------
% %  Raw data conditioning plots of 50 ms loop - 10
% %-------------------------------------------------
% 
% % cd(fig_folder)
% for j = 1:7 % LOOP FOR BLOCKS
%    h=figure; % One figure per block
%    set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
%    tlstr=['Block number ', int2str(j), '_with50msin10_3'];
%    for i=1:10 % LOOP FOR TRIALS
%        indice_take = j*10+i; % in order not to take first 10 trials
%        indice_save = (j-1)*10+i; % in order to save starting at 1 
%        str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
%        subplot (5,2,i); plot(ms_OD(:,indice_take)); title(str);
%        hold on; 
%        if noise_resp(indice_save,3)==0 % only for trials where a CR was detected
%        else
%            line([noise_resp(indice_save,3) noise_resp(indice_save,3)+50],[-0.1 -0.1],'Color',[1 0 0]); % put a red line under the CR
%        end
%        hold on;
%        if noise_resp(indice_save,4)==0 % only for trials where a alpha blink detected
%        else
%            line([noise_resp(indice_save,4) noise_resp(indice_save,4)+50],[-0.1 -0.1],'Color',[0 1 0]); % put a green line under the alpha blink
%        end
%        ylim([-0.2 0.2]);
%    end
%    if 7 == exist(tlstr, 'file'); fig_str = [tlstr,'_2']; else fig_str = tlstr; end
%    print(h,fig_str,'-dpng');
%    pause
% end
% close all

%-----------------------------------------
%        Counting alpha blinks 
%-----------------------------------------
% cd(data_folder);
% save('AB_block','AB_block');save('noise_resp_a','noise_resp_a'); xlswrite('AB_block',AB_block); xlswrite('noise_resp_a',noise_resp_a); 
% Learning curve figure




%---------------------------------------
%       With the first 10 trials
%---------------------------------------

% for j = 1:8
%     for i=  1:10
%         ref= mean(abs(ms_OD(198:398,(j-1)*10+i)),1); % ref for noise = first 400 ms
%         test=mean(abs(ms_OD(598:798,(j-1)*10+i)),1); % 200 ms before US
%         noise_resp((j-1)*10+i,1) = ref; noise_resp((j-1)*10+i,2) = test;
%         if test > 1.5*ref % if 200 ms before US activity is bigger than 1.5*noise = CR
%             CR=CR+1;
%             noise_resp((j-1)*10+i,3) = 1;
%         else
%             noise_resp((j-1)*10+i,3) = 0;
%         end
%     end
%     CR_block(j,1)=CR;
%     CR=0;
% end

%-------------------------------------------
%  Raw data conditioning plots w/ first 10
%-------------------------------------------

% for j = 1:8
%    h=figure; 
%    for i=1:10
%        str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
%        subplot (5,2,i); plot(ms_OD(:,(j-1)*10+i)); title(str);
%        ylim([-0.2 0.2]);
%    end
%    if 7 == exist(['Block number ', int2str(j)], 'file'); fig_str = ['Block number ', int2str(j),'_2']; else fig_str = ['Block number ', int2str(j)]; end
%    print(h,fig_str,'-dpng');
% end


%---------------------------------------
%         Total number of CRs
%---------------------------------------

% for i= 11:num
%     ref= mean(ms_OD(198:398,i),1);
%     test=mean(ms_OD(598:798,i),1);
%     if abs(test) > 2*abs(ref)
%         CR=CR+1;
%     end
% end

end