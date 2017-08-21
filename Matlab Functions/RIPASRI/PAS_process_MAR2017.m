function PAS_process_MAR2017(PAS_ID, gp_ID,sub_ID, subChoice)

% NOTES:
% 1- This code can be called from PAS_routine
% 2- It uses the following nested functions: PAS_preproc; check_badrec;
% remove_badrec; resp_analysis


%% _______INITIALIZATION_____

if subChoice ~= 0 % When the option 1 (single subject) was chosen in PAS_routine
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
    
    % Get list of subjects from existing folders
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
    end
    
    sub_ID='none';
    while ~exist(sub_ID,'dir')
        sub_ID = input('Please enter the ID of the subject you want to process into single quotes (see choices above):   ');
    end
end

gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
cd(gp_folder)

log_folder = [gp_folder,'/', sub_ID];
data_folder = [log_folder,'/Extract/'];
cd(log_folder)
mkdir('Figures');
fig_folder = [log_folder,'/Figures'];
cd(log_folder)
mkdir('Data');
MLdata_folder = [log_folder, '/Data'];

cd(gp_folder);load time_windows; load time_windows_FCR

%% ____DATA PRE-PROCESSING_____
[data_cond, data_param] = PAS_preproc(PAS_ID,gp_ID,sub_ID);

cd(MLdata_folder)
if exist('data_cond.mat','file')==2
    copyfile('data_cond.mat','data_cond_OLD.mat');
end
save('data_cond', 'data_cond')

%% CHECKING PRESENCE OF BAD RECORDINGS

disp('1- Yes');
disp('2- No');
prev_badrec_rep=0;
while prev_badrec_rep ~=1 && prev_badrec_rep ~=2
    prev_badrec_rep=input('Did you already check for bad recordings for this subject (y/n) ?:  ');
end
if prev_badrec_rep== 1 % checking if bad rec definition was already done
    cd(MLdata_folder)
    load data_goodrecs.mat % if yes load the data of previously defined good recrodings
else
    disp('1- Yes');
    disp('2- No');
    badrec_rep=0;
    while badrec_rep ~=1 && badrec_rep ~=2
        badrec_rep=input('Do you want to check for bad recordings (y/n) ?:  ');
    end
    badrecs = 0;
    if badrec_rep==1
        cd(fig_folder)
        [badrecs, badrec_muscles] = check_badrec(data_cond, data_param, time_windows, time_windows_FCR); % else proceed to bad rec checking
    end
end
%% REMOVING BAD RECORDINGS
DateStamp=datestr(clock);
if prev_badrec_rep==2
    if badrecs==0
        data_goodrecs = data_cond;
    else
        if exist('Bad_Recordings.mat','file')==2
            copyfile('Bad_Recordings.mat',['Bad_Recordings_',DateStamp,'.mat']);
        end
        save('Bad_Recordings.mat','badrecs')
        
        [data_goodrecs] = remove_badrec(data_cond, badrecs, badrec_muscles, data_param, time_windows, time_windows_FCR);
        
        cd(MLdata_folder)
        if exist('data_goodrecs.mat','file')==2
            copyfile('data_goodrecs.mat',['data_goodrecs_',DateStamp,'.mat']);
        end
        save('data_goodrecs','data_goodrecs');
        
        disp('End of BAD REC REMOVAL CHECK');
    end
end

%% MEASURING RESPONSES
cd(fig_folder)
[resp, resp_singleExtra] = resp_analysis(data_goodrecs,data_param,time_windows,time_windows_FCR);

cd(MLdata_folder);
if exist('Resp_Amp.mat','file')==2
    %     rep_replace=input('A file with the name Resp_Amp.mat already exists, do you want to replace or keep both? (r/b):   ');
    %     if strcmp(rep_replace,'b')
    copyfile('Resp_Amp.mat',['Resp_Amp_',DateStamp,'.mat']);
    %     end
end
save('Resp_Amp','resp');

if exist('Resp_SingleExt.mat','file')==2
    %     rep_replace=input('A file with the name Resp_Amp.mat already exists, do you want to replace or keep both? (r/b):   ');
    %     if strcmp(rep_replace,'b')
    copyfile('Resp_SingleExt.mat',['Resp_SingleExt_',DateStamp,'.mat']);
    %     end
end
save('Resp_SingleExt','resp_singleExtra');

OK = '';
save OK

disp('End of PAS_process ! Hourray!')


