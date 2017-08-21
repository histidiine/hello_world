function process_mean_Hinhib

% INITIALIZATION
gp_ID = input('Please enter the ID of the group as named in the folders containing his/her data - in between single quotes:  ');
gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',gp_ID];

step_name= {'Baseline', 'Post'};
cond_name={'single pulse','-1ms ISI', '-0,5ms ISI', 'synchronous','0,5ms ISI', '1ms ISI',...
    '2ms ISI','5ms ISI','7,5ms ISI','10ms ISI','25ms ISI','50ms ISI','75ms ISI','100ms ISI','200ms ISI'};
fig_name='FCR';
sub_ID_all={'BAM','ALS','ETB','CAM','CAR','GIB','LUM','MAS','ISH'}; % ,'ELB'

for i = 3:length(sub_ID_all)
    
    sub_ID=sub_ID_all{i};
    
    log_folder = ['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/HV/', sub_ID];
    ML_data_folder = [log_folder,'/Data'];
    
    cd(ML_data_folder)
    
        while ~exist('MandH.mat','file')==1
            disp(['MandH.cannot be found, please check in ', sub_ID ,' folder'])
            pause
        end
        load MandH.mat
    
        while ~exist('FCR_cond.mat','file')==1
            disp('FCR_cond.cannot be found, please check in the folder')
            pause
        end
        load FCR_cond.mat
    
%     if ~exist('MandH.mat','file')==1 || ~exist('FCR_cond.mat','file')==1
%         break
%     else
%         load MandH.mat
%         load FCR_cond.mat
%     end
    
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
    
    H_percent = squeeze(H_ratio(:,3,:))*100-100; % conditions, steps
    
    Mean_Phase1 = mean(H_percent(1:5,:)); 
    Mean_Phase2 = mean(H_percent(7:10,:)); 
    Mean_Phase3 = mean(H_percent(11:14,:)); 
    
    ALL_H_percent(:,:,i)=H_percent;
    ALL_MEAN(i,:,1)=Mean_Phase1;
    ALL_MEAN(i,:,2)=Mean_Phase2;
    ALL_MEAN(i,:,3)=Mean_Phase3;
    
end

% ALL_MEAN(:,:,1)=ALL_Mean_Phase1;
% ALL_MEAN(:,:,2)=ALL_Mean_Phase2;
% ALL_MEAN(:,:,3)=ALL_Mean_Phase3;

cd(gp_folder);save('ALL_MEAN','ALL_MEAN');save('ALL_H_percent','ALL_H_percent');

end