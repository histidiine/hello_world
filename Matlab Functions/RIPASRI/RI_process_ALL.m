function RI_process_ALL()

%_________________NOTES_________________

% This script should take the data from one group at a time, generate
% matrices for each subjects by calling RI_process_MAR2017 in a loop,
% generate then matrices per group and finally generate a matrix with all
% groups in it.

% Left to decide - should it reprocess all subjects everytime, or give an
% option to not reprocess those who were already done before and just
% "append" new data to the already existing

%% DATA PROCESSING
test=1;

while test == 1;
    
    gp_ID=0;ini=0;
    while strcmp(gp_ID,'HV')==0 && strcmp(gp_ID,'CRPS')==0 && strcmp(gp_ID,'CRPSD')==0 && strcmp(gp_ID,'FHD')==0 && strcmp(gp_ID,'FUNCT')==0
        if ini == 1
            disp('SORRY ! That is not a valid group name, check out what I suggested in the parenthesis down here v')
        end
        gp_ID = input('Please enter the code for the group you want to process into single quotes (HV,CRPS,CRPSD,FHD,FUNCT):   ');
        ini=1;
    end
    
    % CHANGE ALL IDs IN FOLDERS TO NUMBERS!!!!! SO THEN IT CAN JUST BE AN
    % ITERATION WHATEVER THE GROUP (AND GP_FOLDER)
    %sub_ID_all={'BAM','ALS','ETB','CAM','ELB','CAR','GIB','LUM','MAS','ISH'};
    step_name= {'Baseline', 'Post'};
    cond_name={'single pulse','-1ms ISI', '-0,5ms ISI', 'synchronous','0,5ms ISI', '1ms ISI',...
    '2ms ISI','5ms ISI','7,5ms ISI','10ms ISI','25ms ISI','50ms ISI','75ms ISI','100ms ISI','200ms ISI'}

    
    disp('End of INITIALIZATION');
    
    rep_all_proc = input('Did you already process the data of some subjects? (y/n):  ');
    if strcmp(rep_all_proc,'y')==1
        rep_numS=input('What was the ID of the last subject you processed? (BAM, ALS, ETB, CAM, ELB; CAR; GIB, LUM, MAS, ISH - if none 0) : ');
        if rep_numS == 0
            start_ID =sub_ID_all{3};
        else
            [osef,col]=max(strcmp(sub_ID_all,rep_numS)); clear osef
            if col == length(sub_ID_all)
                start_ID = 3;
                launch_RI_process = 0;
            else
                start_ID = col+1;
                launch_RI_process = 1;
            end
        end
    end
    
    if launch_RI_process ==1
        ALL_H_resp=zeros(10,15,2,5,length(sub_ID_all)); %trial, cond, step, param, subject
        ALL_FCR_cond=zeros(5000,10,15,2,length(sub_ID_all)); %points, trial, cond, step, subject
        for i = start_ID:length(sub_ID_all)
            disp('------------------------')
            disp(['Start of ', sub_ID_all{i}])
            disp('------------------------')
            
            [curr_H, curr_FCR]=RI_process_MAR2017(gp_ID,sub_ID_all{i});
            ALL_H_resp(:,:,:,:,i)=curr_H;
            ALL_FCR_cond(:,:,:,:,i)=curr_FCR;
            disp('------------------------')
            disp(['End of ', sub_ID_all{i}])
            disp('------------------------')
        end
    else
        rep_calculate = input('Do yo want to proceed to group calculation (group mean H) with the subjects already processed? (y/n):  ');
        if strcmp(rep_calculate,'y')
            gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',gp_ID];
            cd(gp_folder)
            disp('Are there any subjects to discard from the analysis?');
            id_discard=input('If yes please enter their ID into curly brackets separated by comas: ');
            sub_ID_proc = setdiff(sub_ID_all,id_discard);
            ALL_H_resp=zeros(10,15,2,5,length(sub_ID_proc)); %trial, cond, step, param, subject
            for i = 1:length(sub_ID_proc)
                disp('------------------------')
                disp(['Start of ', sub_ID_proc{i}])
                disp('------------------------')
                
                sub_ID=sub_ID_proc{i};
                log_folder = ['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/HV/', sub_ID];
                ML_data_folder=[log_folder,'/Data'];
                cd(ML_data_folder)
                disp('Check that you put in archive other MandH you want to keep and left only the last corect one in the folder without DateStamp')
                while ~exist('MandH.mat','file')==1
                    disp('MandH.cannot be found, please check in the folder')
                    pause
                end
                load MandH.mat H_resp
                ALL_H_resp(:,:,:,:,i)=H_resp;
                disp('------------------------')
                disp(['End of ', sub_ID_proc{i}])
                disp('------------------------')
            end
        end
        
        for j = 1:size(ALL_H_resp,3)
            figure;
            for s = 1:size(ALL_H_resp,5)
                subplot(size(ALL_H_resp,5),1,s)
                curr_mean=squeeze(mean(ALL_H_resp(:,:,j,5,s),1));
                plot(curr_mean)
                title(sub_ID_proc{s})
                line([0 15],[curr_mean(:,1) curr_mean(:,1)], 'Color','r');
                hold on
            end
        end
        
        for j = 1:size(ALL_H_resp,3)
            figure;
            for s = 1:size(ALL_H_resp,5)-1
                curr_data=squeeze(ALL_H_resp(:,:,j,5,s));
                curr_data(curr_data==0)=0.001;
                single_data=squeeze(ALL_H_resp(:,1,j,5,s));
                single_data(single_data==0)=0.001; % to avoid division by 0
                subplot(size(ALL_H_resp,5)-1,1,s)
                line([0 15],[100 100], 'Color','r')
                hold on;
                for count = 1:size(curr_data,2) % cond
                    curr_data_pourcent(:,count) = curr_data(:,count)./median(single_data)*100;
                end
                boxplot(curr_data_pourcent, 'labels', cond_name);
                for x_pt = 1:size(curr_data_pourcent,2)
                    hold on; text(x_pt-0.1,120,int2str(median(curr_data_pourcent(:,x_pt))))
                end
                title(sub_ID_proc{s})
                ylim([0 max(max(curr_data_pourcent))])
            end
        end
        
        mean_H=squeeze(mean(mean(ALL_H_resp(:,:,:,5,:),1),5)); % cond, step
        figure;
        for j = 1:size(mean_H,2)
            subplot(10,1,j)
           plot(mean_H(:,j))
           line([0 15],[mean_H(1,j) mean_H(1,j)], 'Color','r');
           hold on
        end
        
      
        
    end
    
    test = input('Do you want to process another group ? (y/n - in single quotes):   ');
end
