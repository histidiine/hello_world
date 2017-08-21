function group_analysis

% NOTES
% This function takes the analyzed MEP amplitudes of processed subject with
% PAS_process_MAR2017 and creates a matrix of the amplitudes of all
% subjects of a same group and saves it in the group folder

% ______________INITIALIZATION_______________
dataset = ['b', 't1', 't2'];
ds_folder = {'Baseline', 'T1', 'T2'};
m_name={'ADM','FDI','FCR','FPB'};
step_name={'baseline','T1','T2'};
cond_name={'single', 'LICI 50', 'LICI 100', 'LICI 200'};
sub_ID_all={'BAM','ETB','CAM','ELB','CAR','GIB','LUM','ISH', 'MAS'};

% __________________________ID & DIRECTORY___________________________
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

n_subjects = length(list_folders);

%% __________________GROUP MATRIX CREATION__________________
totalSub=1;
MEP_amp_GP = zeros(12,4,4,3,5,n_subjects);
data_cond_GP = zeros(14999, 12, 4, 4, 3, n_subjects);
all_single_GP = zeros(48, 4, 3, 5,n_subjects);
for sub = 1:n_subjects
    cd(gp_folder)
    sub_ID = list_folders(sub).name;
    cd(sub_ID); cd('Data');
    if ~exist('OK.mat','file')
        disp([sub_ID, ' was not yet processed with PAS_process_MAR2017 and was thus ignored in this group calculation']);
    else
        load Resp_Amp % Amplitude and other parameters treated with resp_analysis
%         load Resp_SingleExt % Amplitude and other parameters treated with resp_analysis for responses after first burst of paired pulse 
        load Data_Cond % All datapoints of muscle response
%         resp_amp=squeeze(resp(:,:,:,:,5));
        MEP_amp_GP(:,1:size(resp,2),:,:,:,totalSub)=resp;
        data_cond_GP(1:size(data_cond,1),:,1:size(data_cond,3),:,:,totalSub)=data_cond;
        
        %______________SINGLE PULSE PROCESSING______________
%         all_single_GP(1:12,1:size(resp,2),:,:,totalSub)=squeeze(resp(:,:,1,:,:));
%         for curr_cond = 1:3
%             end_trial = 12*(curr_cond+1);
%             start_trial = end_trial-11;
%             all_single_GP(start_trial:end_trial,1:size(resp,2),:,:,totalSub)=squeeze(resp_singleExtra(:,:,curr_cond,:,:));
%         end
        
        totalSub = totalSub +1;
        disp(['Just finished adding the data of ', sub_ID]);
    end
end
disp([int2str(totalSub-1), ' subjects were processed for this group calculation!']);

cd(gp_folder)
if exist('Resp_Amp.mat','file')==2
    DateStamp=datestr(clock);
    copyfile('Resp_Amp.mat',['Resp_Amp_',DateStamp,'.mat']);
end
save('Resp_Amp','MEP_amp_GP');

if exist('data_cond_GP.mat','file')==2
    DateStamp=datestr(clock);
    copyfile('data_cond_GP.mat',['data_cond_GP',DateStamp,'.mat']);
end
save('data_cond_GP','data_cond_GP');

% if exist('all_single_GP.mat','file')==2
%     DateStamp=datestr(clock);
%     copyfile('all_single_GP.mat',['all_single_GP_',DateStamp,'.mat']);
% end
% save('all_single_GP','all_single_GP');


