clear all
clc
close all

% disp('*****************************************************');
% disp('EMG signals EEG-EMG-glove');
% disp('*****************************************************');
% disp(' ');

directory_iniziale=cd;

dir_sub=input('Folder path of the subject to analyze: ','s');
name_file={'block_01.mat','block_02.mat','block_03.mat','block_04.mat','block_05.mat','block_06.mat','block_07.mat','block_08.mat','block_09.mat','block_10.mat','block_11.mat','block_12.mat','block_13.mat'};
n_file1=input('How many files? ');

cd(dir_sub)
load('MVC_proc');

for n=1:n_file1 
    
%     for c=1:size(name_conf{n},2)
%      MVCval(c)=MVCb(find(strcmp(name_conf{n}{c},conf_2)==1));
%     end

     
cd(dir_sub);
load(name_file{n});
cd(directory_iniziale);
   
for i=1:noChans-1  
EMG(i,:)=Data{i};
end
clear Activities Data length_sec samplingRate

[EMG_proc]=preprocessing(EMG,1500,MVCb,noChans-1); % MVC and the other EMG are preprocessed in the same way

name=[name_file{n}(1:length(name_file{n})-4),'_p']; 
cd(dir_sub);  
eval(['  save ' name ' EMG_proc MVCb']);
clear EMG MVCval EMG_proc i 
end
  