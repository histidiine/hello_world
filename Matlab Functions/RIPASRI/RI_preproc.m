function [RI_data,RI_cond,data_param] = RI_preproc(PAS_ID,gp_ID,sub_ID)

%% FOLDER INITIATION
gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
log_folder = [gp_folder,'/', sub_ID];
data_folder = [log_folder,'/Extract/'];

%% EXTRACTING INFO FROM LOG
cd(log_folder)
log_base_temp=importdata('RI_base.log');
log_base = log_base_temp.data;
log_post_temp=importdata('RI_post.log');
log_post = log_post_temp.data;

s_rate=[log_base(1),log_post(1)]; % I could only take one as they should all be the same, but jic to avoid an error
points=[log_base(2),log_post(2)];
n_muscle=[log_base(3),log_post(3)];
n_trial=[log_base(4),log_post(4)];

data_param=zeros(4,3);
data_param(1,:)=s_rate;data_param(2,:)=points;data_param(3,:)=n_muscle;data_param(4,:)=n_trial;

RI_base = xlsread('RI_base.xlsx');
RI_post = xlsread('RI_post.xlsx');

conditions = zeros(max(n_trial), 2, 2); %trials, delas, step
conditions(:,1,1)= squeeze(RI_base(:,7));conditions(:,2,1)= squeeze(RI_base(:,9)); % conditions (:,:,1) = RIbase
conditions(:,1,2)= squeeze(RI_post(:,7));conditions(:,2,2)= squeeze(RI_post(:,9)); % conditions (:,:,2) = RIpost

clear log_base_temp log_base log_post_temp log_post
disp('End of LOG EXTRACTION')

%% GETTING DATA FROM NGUYET EXTRACT
cd(data_folder)

% Initialisation
ds_folder = {'RI_base', 'RI_post'};
RI_data=zeros(max(points),2,max(n_trial),2); % timepoints, muscle, trials, step(baseline; post)
% Extraction loop - loads file 1 by 1 and add the data in the ri_data matrix
for j = 1:2 % baseline / post
    cd((ds_folder{j}))
    for i = 1:max(n_trial)
        while exist(int2str(i),'file') == 0
            disp(['Programmed paused because file ', int2str(i),' is missing in ', ds_folder{j}]);
            pause
        end
        temp=importdata(int2str(i));
        RI_data(:,1,i,j) = temp(1:points(j),1); % I put row 1:points in case there was a mistake during the nguyet and a file was recorded twice
        RI_data(:,2,i,j) = temp(1:points(j),2); clear temp;
    end
    cd(data_folder)
end

cd(log_folder)
if ~exist('Data','dir'); mkdir('Data'); end
ML_data_folder=[log_folder,'/Data'];
cd(ML_data_folder)

save('RI_data','RI_data')

disp('End of DATA EXTRACTION')
%% EXTRACT DATA IN CONDITIONS

RI_cond = order_by_cond_JUL2017(RI_data,conditions,max(n_trial));
save('RI_cond', 'RI_cond');

disp('End of COND REARRANGING')

end