function [data_cond, data_param] = PAS_preproc(PAS_ID,gp_ID,sub_ID)
% This function will extract the data from all the trial files created with
% Nguyet, orders the data by conditions in the data_cond matrix

% ____FOLDERS INITIALIZATION_____
gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
log_folder = [gp_folder,'/', sub_ID];
data_folder = [log_folder,'/Extract/'];
cd(log_folder)
if exist('Figures', 'dir')==0; mkdir('Figures'); end
fig_folder = [log_folder,'/Figures'];

%_____LABELS INITIALIZATION_____
dataset = ['b', 't1', 't2'];
ds_folder = {'Baseline', 'T1', 'T2'};
m_name={'ADM','FDI','FCR','FPB'};
step_name={'baseline','T1','T2'};
cond_name={'single pulse', 'LICI 50', 'LICI 100', 'LICI 200'};

cd(log_folder) 
if exist('Data', 'dir')==0; mkdir('Data'); else end
MLdata_folder = [log_folder, '/Data']; 

disp('End of INITIALIZATION')
%% EXTRACTING INFO FROM LOG

cd(log_folder) 
if ~exist('Baseline.log')
    log_b=zeros(4,1);
else
    log_b_temp=importdata('Baseline.log');
    log_b = log_b_temp.data;
end
if ~exist('T1.log')
    log_t1=zeros(4,1);
else
    log_t1_temp=importdata('T1.log');
    log_t1 = log_t1_temp.data;
end
if ~exist('T2.log')
    log_t2=zeros(4,1);
else
    log_t2_temp=importdata('T2.log');
    log_t2 = log_t2_temp.data;
end

s_rate=[log_b(1),log_t1(1),log_t2(1)]; % I could only take one as they should all be the same, but jic to avoid an error
points=[log_b(2),log_t1(2),log_t2(2)]; 
n_muscle=[log_b(3),log_t1(3),log_t2(3)];
n_trial=[log_b(4),log_t1(4),log_t2(4)];

data_param=zeros(4,3);
data_param(1,:)=s_rate;data_param(2,:)=points;data_param(3,:)=n_muscle;data_param(4,:)=n_trial;

% disp('Please verify that Baseline.log, T1.log and T2.log were copied to excel files and the first three rows deleted and press any key'); pause;
conditions = zeros(max(n_trial), 2, 3); %trials, delas, step
if exist('Baseline.xlsx') == 2
Baseline = xlsread('Baseline.xlsx');
conditions(:,1,1)= squeeze(Baseline(:,3));conditions(:,2,1)= squeeze(Baseline(:,5));
end
if exist('T1.xlsx') == 2
T1 = xlsread('T1.xlsx');
conditions(:,1,2)= squeeze(T1(:,3));conditions(:,2,2)= squeeze(T1(:,5));
end
if exist('T2.xlsx') == 2
T2 = xlsread('T2.xlsx');
conditions(:,1,3)= squeeze(T2(:,3));conditions(:,2,3)= squeeze(T2(:,5));
end

% clear log_b_temp log_b log_t1_temp log_t1 log_t2_temp log_t2

disp('End of LOG EXTRACTION')
%% EXTRACTING DATA
% Extraction loop - loads file 1 by 1 and adds the data in the pre_post_data
cd(data_folder)
pre_post_data=zeros(max(points),max(n_trial),max(n_muscle),3); % data, trials, muscles and step (baseline, t1 or t2)
for j = 1:size(pre_post_data,4) % pre, t1 and t2
    if ~exist(ds_folder{j})
    else
        cd(ds_folder{j})
        for i = 1:n_trial(j) % each trial
            while exist(int2str(i),'file') == 0
                disp(['Programmed paused because file ', int2str(i), ' of ', step_name{j},' is missing'])
                pause
            end
            temp=importdata(int2str(i));
            pre_post_data(1:points(j),i,1:n_muscle(j),j) = temp(1:points(j),:); % I put row 1:points in case  a file was recorded twice in Nguyet
        end
        cd(data_folder)
    end
end

disp('End of DATA EXTRACTION')
%% ORDER BY CONDITIONS
cd(log_folder) 
count1=1;count2=1;count3=1;count4=1;
for j = 1:3 %step (baseline, t1, t2)
    temp_data = squeeze(pre_post_data(:,:,:,j)); % data, trials, muscle
    for i = 1:n_trial(j) %trials
        if conditions(i,2,j)==0  % single pulse
            data_cond(:,count1,:,1,j) =  temp_data(:,i,:); % data_cond - data, trials, muscle, cond, step
            count1=count1+1;
        else if conditions(i,2,j)==-50  % LICI50
                data_cond(:,count2,:,2,j) =  temp_data(:,i,:);
                count2=count2+1;
            else if conditions(i,2,j)==-100  % LICI100
                    data_cond(:,count3,:,3,j) =  temp_data(:,i,:);
                    count3=count3+1;
                else if conditions(i,2,j)==-200  % LICI200
                        data_cond(:,count4,:,4,j) =  temp_data(:,i,:);
                        count4=count4+1;
                    else while isempty(conditions(i,2,j))==1
                            disp('it appears the condition is empty, check that they were imported properly');
                        end
                    end
                end
            end
            
        end
    end
    count1=1;count2=1;count3=1;count4=1;
end

clear count1 count2 count3 count4

disp('End of COND RE-ORDERING')


end