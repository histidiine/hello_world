% Routine to process MEP amplitudes of all subject of one group

sub_ID_all={'BAM','ALS','ETB','CAM','ELB','CAR','GIB','LUM','MAS','ISH'};

PAS_ID=0;ini=0;
while strcmp(PAS_ID,'5hz')==0 && strcmp(PAS_ID,'02hz')==0
    if ini == 1
        disp('SORRY ! That is not a valid group name, check out what I suggested in the parenthesis down here v')
    end
    PAS_ID = input('Please enter the code for the PAS protocole you want to process into single quotes (5hz, 02hz):   ');
    ini=1;
end

gp_ID=0;ini=0;
while strcmp(gp_ID,'HV')==0 && strcmp(gp_ID,'CRPS')==0 && strcmp(gp_ID,'CRPSD')==0 && strcmp(gp_ID,'FHD')==0 && strcmp(gp_ID,'FUNCT')==0
    if ini == 1
        disp('SORRY ! That is not a valid group name, check out what I suggested in the parenthesis down here v')
    end
    gp_ID = input('Please enter the code for the group you want to process into single quotes (HV,CRPS,CRPSD,FHD,FUNCT):   ');
    ini=1;
end

gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
cd(gp_folder)

all_files = dir;
all_dir = all_files([all_files(:).isdir]);
num_dir = numel(all_dir)-3;

processChoice = input('Do you want to process subjects data or create figures from already treated data (d/f)?:   ');
% dataProc_rep = input('Do you want to process subjects data (y/n)?  ');
% figProc_rep = input('Do you want to process figures of already treated data (y/n)?  ');

%% DATA PROCESSING

switch processChoice
    
    case 'd'
        
%         if ~exist(['MEP_amp_',gp_ID,'.mat'],'file')
%         else
%             DateString = datestr(clock);
%             copyfile(['MEP_amp_',gp_ID,'.mat'],['MEP_amp_PREV',DateString,'.mat'])
%         end
        
        already_rep=input('Did you already process some subjects (y/n)?:  ');
        
        if strcmp(already_rep,'y') == 1
            load(['MEP_amp_',gp_ID,'.mat']);
            load(['data_goodrecs_',gp_ID,'.mat'],'data_goodrecs_ALL');
            load('alreadyProc');
            lastSub=max(processedSub);
            start_ID=lastSub+1;
            count=start_ID-2;
            disp(['It appears the last subject you processed was ', int2str(lastSub)]);
        else
            MEP_amp=zeros(12,4,4,3,5,num_dir); % trial, cond, muscle, step, subjects
            data_goodrecs_ALL=zeros(14999,12,4,4,3,num_dir);
            processedSub=zeros(num_dir,1);
            start_ID=1;
            count = 1;
        end
        
        for sub = start_ID:num_dir
            check_sub=input(['Do you want to process subject ', int2str(sub), ' ?  ']);
            if strcmp(check_sub,'y')==1
                
                disp(['BEGINNING OF PROCESSING OF SUBJECT ', int2str(sub)]);
                sub_ID = sub_ID_all{sub};
                [data_goodrecs, resp] = PAS_process_MAR2017(PAS_ID,gp_ID,sub_ID);
                disp(['END OF PROCESSING OF SUBJECT ', int2str(sub)]);
                
                % _____SAVE______
                %             MEP_amp(:,:,:,:,:,count)=resp;% trial,cond, muscle, step
                resp_amp=squeeze(resp(:,:,:,:,5));
                MEP_amp(:,:,:,:,count)=resp_amp;% amplitude -->trial,cond, muscle, step,subject
                curr_ntrial=size(data_goodrecs,1);
                curr_nmuscle=size(data_goodrecs,3);
                data_goodrecs_ALL(1:curr_ntrial,:,1:curr_nmuscle,:,:,count)=data_goodrecs;
                cd(gp_folder)
                save(['MEP_amp_',gp_ID,'.mat'],'MEP_amp')
                save(['data_goodrecs_',gp_ID,'.mat'],'data_goodrecs_ALL')
                processedSub(count,1)=sub;
                save('alreadyProc.mat','processedSub')
                
                count = count+1;
            end
        end
        
        %% FIGURE PROCESSING
        
    case 'f'
        
        cd(gp_folder)
        gp_data=['MEP_amp_', gp_ID, '.mat'];
        load(gp_data);
        cd('Fig_Save')
        MEP_graph(MEP_amp);
        
end

