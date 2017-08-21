function graph_gen 

%_____________________INITIALIZATION______________________

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

gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
cd(gp_folder)

fileName='Resp_Amp.mat';
load(fileName)
fileName='data_cond_GP.mat';
load(fileName)
load time_windows;

points=14999;
s_rate=40000;

duration_ms=points/s_rate*1000;
step=duration_ms/points;
time_axis=0:step:duration_ms-step; time_axis=time_axis';


%% (1) INDIVIDUAL MEAN MEP AMPLITUDE + GROUP MEAN (Leo style) 

curr_data = squeeze(MEP_amp_GP(:,1,1,:,5,:)); % trials, step, subject - first 1 = ADM, second 1 = single pulse, 5 = amplitude

% This loop allows to not count in the calculation people who didn't get a
% T1 (otherwise it brings down the mean) + I added nanmean everywhere
% instead of mean 
for sub = 1:size(curr_data,3)
    for step = 1:size(curr_data,2)
        testingEx=sum(all(curr_data(:,step,sub),2));
        if testingEx == 0
            curr_data(:,step,sub)=NaN;
        end
    end
end

curr_mean = squeeze(nanmean(curr_data,1)); % step, subject
all_mean = squeeze(nanmean(curr_mean,2));

curr_mean_ratio = zeros(size(curr_mean));
for step = 1:size(curr_mean,1)
    curr_mean_ratio(step,:) = curr_mean(step,:)./curr_mean(1,:);
end
curr_mean_ratio = curr_mean_ratio*100;
maxlim=max(max(curr_mean_ratio));

all_mean_ratio = squeeze(nanmean(curr_mean_ratio,2));

% RAW
figure;
for step = 1:size(curr_mean,2)
    line([1,2,3],curr_mean(:,step),'LineWidth', 2);
    hold on
end
line([1,2,3],all_mean,'Color', 'r', 'LineWidth', 3);
hold on;
mean_str=num2str(all_mean);
text([1,2,3]-0.01, all_mean+0.1, mean_str, 'Color', 'r', 'FontSize', 14);
ylim([0 3.5]);

% RATIO - baseline, t1, t2
figure;
for sub = 1:size(curr_mean_ratio,2)
    line([1,2,3],curr_mean_ratio(:,sub), 'Color','k', 'LineWidth', 2);
    hold on;
end
line([1,2,3],all_mean_ratio, 'Color', 'r', 'LineWidth', 3);
hold on;
mean_str=int2str(round(all_mean_ratio));
text([1,2,3]-0.01, all_mean_ratio+5, mean_str, 'Color', 'r', 'FontSize', 14);
ylim([0 maxlim]);

% RATIO - baseline, t2
curr_mean_ratio_2=curr_mean_ratio;
curr_mean_ratio_2(2,:)=[];
all_mean_ratio_2 = squeeze(nanmean(curr_mean_ratio_2,2));
figure;
for sub = 1:size(curr_mean_ratio,2)
    line([1,2],curr_mean_ratio_2(:,sub), 'Color','k', 'LineWidth', 2);
    hold on;
    mean_str=int2str(round(curr_mean_ratio_2(:,sub)));
    text([1,2]-0.01, curr_mean_ratio_2(:,sub)+5, mean_str, 'Color', 'k', 'FontSize', 14);
    hold on;
end
line([1,2],all_mean_ratio_2, 'Color', 'r', 'LineWidth', 3);
hold on;
mean_str=int2str(round(all_mean_ratio_2));
text([1,2]-0.01, all_mean_ratio_2+3, mean_str, 'Color', 'r', 'FontSize', 14);
ylim([0 maxlim]);

% RATIO - baseline, t1
curr_mean_ratio_2=curr_mean_ratio;
curr_mean_ratio_2(3,:)=[];
all_mean_ratio_2 = squeeze(nanmean(curr_mean_ratio_2,2));
figure;
for sub = 1:size(curr_mean_ratio,2)
    line([1,2],curr_mean_ratio_2(:,sub), 'Color','k', 'LineWidth', 2);
    hold on;
    mean_str=int2str(round(curr_mean_ratio_2(:,sub)));
    text([1,2]-0.01, curr_mean_ratio_2(:,sub)+5, mean_str, 'Color', 'k', 'FontSize', 14);
    hold on;
end
line([1,2],all_mean_ratio_2, 'Color', 'r', 'LineWidth', 3);
hold on;
mean_str=int2str(round(all_mean_ratio_2));
text([1,2]-0.01, all_mean_ratio_2+3, mean_str, 'Color', 'r', 'FontSize', 14);
ylim([0 maxlim]);


%% (2) All subject of a group - PAS EFFECT ON SINGLE AND LICI 100 - SIMPLE SINGLE ONLY


gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
cd(gp_folder)

fileName='Resp_Amp.mat';
load(fileName)
% singleFile = 'all_single_GP.mat';
% load(singleFile);

n_subjects = size(MEP_amp_GP,6);
ind = sqrt(n_subjects);

if ~mod(ind,1) % if number of subject is an integer
    row = ind;
    col = row;
else
    if ind < floor(ind)+0.5
        row = floor(ind);
        col = row +1;
    else
        row = floor(ind)+1;
        col = row;
    end
end

m=1;
for MyUnit =1:2
    hist=figure;count=1;
    group_data=squeeze(MEP_amp_GP(:,m,:,:,5,:));% trials, cond, step, sub
    gp_mean=squeeze(mean(group_data,1));
    gp_ratio=zeros(size(gp_mean));
    for step = 1:size(gp_mean,2)
        for cond = 1:size(gp_mean,1)
            gp_ratio(cond,step,:) = gp_mean(cond,step,:)./gp_mean(1,1,:);
        end
    end
    max_lim = squeeze(max(max(max(max(gp_mean(1,:,:))))))+1;
    max_lim_r=squeeze(max(max(max(max(gp_ratio(1,:,:))))))+1;
    for sub = 1:size(MEP_amp_GP,6)
        subplot(row,col,count)
        %_________________CALCULATING__________________
        curr_data = squeeze(MEP_amp_GP(:,m,:,:,5,sub)); % trials, cond, step
%         curr_single = squeeze(all_single_GP(:,m,:,5,sub)); % trials, step
        
        curr_mean=squeeze(mean(curr_data,1));
        curr_mean(2,:) = []; curr_mean(3,:) = [];
        
        curr_dev = squeeze(std(curr_data,1)); % cond, step
        curr_dev(2,:)=[];curr_dev(3,:)=[];
        
        curr_ratio=zeros(size(curr_mean));
        for step = 1:size(curr_mean,2)
            for cond = 1:size(curr_mean,1)
                curr_ratio(cond,step) = curr_mean(cond,step)./curr_mean(1,1);
            end
        end
        
        %________________CREATING FIGURE_________________
        
        x=[0.78,1.78,2.78,3.78 ;...
            1,2,3,4; 1.22,2.22,3.22,4.22];
        
        if MyUnit == 1 % PLOT IN MV WITH ERROR
            bar(curr_mean, 1,'edgecolor','k', 'linewidth', 2);
            hold on;
            for col = 1:size(curr_mean,2)
                %                 std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
                errorbar(x(col,1:2),curr_mean(:,col),curr_dev(1:2,col),'.','LineWidth',2,'Color','k');
            end
             ylim([0 max_lim]);
            unitStr='mV';
        else if MyUnit == 2 % PLOT IN RATIO
                bar(curr_ratio, 1,'edgecolor','k', 'linewidth', 2);
                hold on;
                for col = 1:size(curr_mean,2)
                    %                     std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
                    errorbar(x(col,1:2),curr_ratio(:,col),curr_dev(1:2,col),'.','LineWidth',2,'Color','k');
                end
                ylim([0 max_lim_r]);
                unitStr='ratio';
            end
        end
        
        %_____________GRAPH PARAMETERS______________
        colormap(bone);
        hold on
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        count = count +1;
    end
end


%% (2 BIS) All subject of a group - PAS EFFECT ON SINGLE AND LICI 100 - SIMPLE SINGLE ONLY - LICI AS % REDUCTION


gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
cd(gp_folder)

fileName='Resp_Amp.mat';
load(fileName)
% singleFile = 'all_single_GP.mat';
% load(singleFile);

n_subjects = size(MEP_amp_GP,6);
ind = sqrt(n_subjects);

if ~mod(ind,1) % if number of subject is an integer
    row = ind;
    col = row;
else
    if ind < floor(ind)+0.5
        row = floor(ind);
        col = row +1;
    else
        row = floor(ind)+1;
        col = row;
    end
end

m=1;
for MyUnit =1:2
    hist=figure;count=1;
    group_data=squeeze(MEP_amp_GP(:,m,:,:,5,:));% trials, cond, step, sub
    gp_mean=squeeze(mean(group_data,1));
    gp_ratio=zeros(size(gp_mean));
    for step = 1:size(gp_mean,2)
        for cond = 1:size(gp_mean,1)
            gp_ratio(cond,step,:) = gp_mean(cond,step,:)./gp_mean(1,1,:);
        end
    end
    max_lim = squeeze(max(max(max(max(gp_mean(1,:,:))))))+1;
    max_lim_r=squeeze(max(max(max(max(gp_ratio(1,:,:))))))+1;
    for sub = 1:size(MEP_amp_GP,6)
        subplot(row,col,count)
        %_________________CALCULATING__________________
        curr_data = squeeze(MEP_amp_GP(:,m,:,:,5,sub)); % trials, cond, step
%         curr_single = squeeze(all_single_GP(:,m,:,5,sub)); % trials, step
        
        curr_mean=squeeze(mean(curr_data,1));
        curr_mean(2,:) = []; curr_mean(3,:) = [];
        
        curr_dev = squeeze(std(curr_data,1)); % cond, step
        curr_dev(2,:)=[];curr_dev(3,:)=[];
        
        curr_ratio=zeros(size(curr_mean));
        for step = 1:size(curr_mean,2)
            for cond = 1:size(curr_mean,1)
                curr_ratio(cond,step) = curr_mean(cond,step)./curr_mean(1,1);
            end
        end
        
        %________________CREATING FIGURE_________________
        
        x=[0.78,1.78,2.78,3.78 ;...
            1,2,3,4; 1.22,2.22,3.22,4.22];
        
        if MyUnit == 1 % PLOT IN MV WITH ERROR
            bar(curr_mean, 1,'edgecolor','k', 'linewidth', 2);
            hold on;
            for col = 1:size(curr_mean,2)
                %                 std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
                errorbar(x(col,1:2),curr_mean(:,col),curr_dev(1:2,col),'.','LineWidth',2,'Color','k');
            end
             ylim([0 max_lim]);
            unitStr='mV';
        else if MyUnit == 2 % PLOT IN RATIO
                bar(curr_ratio, 1,'edgecolor','k', 'linewidth', 2);
                hold on;
                for col = 1:size(curr_mean,2)
                    %                     std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
                    errorbar(x(col,1:2),curr_ratio(:,col),curr_dev(1:2,col),'.','LineWidth',2,'Color','k');
                end
                ylim([0 max_lim_r]);
                unitStr='ratio';
            end
        end
        
        %_____________GRAPH PARAMETERS______________
        colormap(bone);
        hold on
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        count = count +1;
    end
end



%% (3) MEAN - PAS EFFECT ON SINGLE AND LICI 100


%_________________CALCULATING__________________
m=1;
curr_data = squeeze(MEP_amp_GP(:,m,:,:,5,:)); % trials, cond, step, subject
% curr_single = squeeze(all_single_GP(:,m,:,5,:)); % trials, step, subject

%     curr_mean=zeros(4,3);
%     curr_mean(1,:)=squeeze(mean(mean(curr_single,1),3));
%     curr_mean(2:4,:) = squeeze(mean(mean(curr_data(:,2:4,:,:),1),4)); % cond, step
%     curr_mean(4,:) = []; curr_mean(2,:) = [];

% This loop allows to not count in the calculation people who didn't get a
% T1 (otherwise it brings down the mean) + I added nanmean everywhere
% instead of mean 
for sub = 1:size(curr_data,3)
    for step = 1:size(curr_data,2)
        testingEx=sum(all(curr_data(:,step,sub),2));
        if testingEx == 0
            curr_data(:,step,sub)=NaN;
        end
    end
end

pre_curr_mean=squeeze(nanmean(curr_data,1)); % 4, 3, 9 = cond, step, sub
curr_mean = squeeze(nanmean(pre_curr_mean,3)); % cond x step
curr_data_conc=curr_data(:,:,:,1);
for k = 2:size(curr_data,4)
    addData=curr_data(:,:,:,k);
    curr_data_conc=cat(1,curr_data_conc,addData);
end


curr_dev = squeeze(std(curr_data_conc,1)); % cond, step
curr_ratio=zeros(size(pre_curr_mean));
curr_ratio_conc=zeros(size(curr_data_conc)); % cond, ste, sub
for step = 1:size(pre_curr_mean,2)
    for cond = 1:size(pre_curr_mean,1)
        curr_ratio(cond,step,:) = curr_mean(cond,step,:)./curr_mean(1,1,:);
        curr_ratio_conc(:,cond,step) = curr_data_conc(:,cond,step)./mean(curr_data_conc(:,1,1),1);
    end
end
curr_ratio_dev = squeeze(std(curr_ratio_conc,1));

curr_mean(4,:)=[];curr_mean(2,:)=[];
curr_dev(4,:)=[];curr_dev(2,:)=[];
curr_ratio_conc(:,4,:)=[];curr_ratio_conc(:,2,:)=[];
curr_ratio(4,:,:)=[];curr_ratio(2,:,:)=[];
curr_ratio_dev(4,:)=[];curr_ratio_dev(2,:)=[];
curr_ratio_mean = squeeze(nanmean(curr_ratio,3));

%________________CREATING FIGURE_________________
%     max_V=max(max(max(curr_dev)))+ max(max(max(curr_mean)));
%     max_R=max(max(max(curr_ratio))) + max(max(max(curr_dev)));

hist=figure;
MyUnit = 2;

x=[0.78,1.78,2.78,3.78 ;...
    1,2,3,4; 1.22,2.22,3.22,4.22];

if MyUnit == 1 % PLOT IN MV WITH ERROR
    bar(curr_mean, 1,'edgecolor','k', 'linewidth', 2);
    hold on;
    for col = 1:size(curr_mean,2)
%         std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
        std_err = curr_dev(:,col)/sqrt(size(curr_data,4));
        errorbar(x(col,1:2),curr_mean(:,col),std_err,'.','LineWidth',2,'Color','k');
    end
    %  ylim([0 (max_V+0.1)]);
    unitStr='mV';
else if MyUnit == 2 % PLOT IN RATIO
        bar(curr_ratio_mean, 1,'edgecolor','k', 'linewidth', 2);
        hold on;
        for col = 1:size(curr_mean,2)
%             std_err = curr_dev(:,col)/sqrt(size(curr_ratio_conc,1));
            std_err = curr_ratio_dev(:,col)/sqrt(size(curr_data,4));
            errorbar(x(col,1:2),curr_ratio_mean(:,col),std_err,'.','LineWidth',2,'Color','k');
        end
        unitStr='ratio';
    end
end
colormap(bone);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% count = count +1;
% end

% ________________________RATIO + LICI AS % INHIBITION - BASELINE SINGLE ________________________

curr_data = squeeze(MEP_amp_GP(:,m,:,:,5,:)); % trials, cond, step, subject
pre_curr_mean=squeeze(nanmean(curr_data,1)); % 4, 3, 9 = cond, step, sub
curr_mean = squeeze(nanmean(pre_curr_mean,3)); % cond x step

% This loop allows to not count in the calculation people who didn't get a
% T1 (otherwise it brings down the mean) + I added nanmean everywhere
% instead of mean 
for sub = 1:size(curr_data,3)
    for step = 1:size(curr_data,2)
        testingEx=sum(all(curr_data(:,step,sub),2));
        if testingEx == 0
            curr_data(:,step,sub)=NaN;
        end
    end
end

curr_ratio=zeros(size(pre_curr_mean));
for step = 1:size(pre_curr_mean,2)
    for cond = 1:size(pre_curr_mean,1)
        if cond == 1 % Single Pulse
            curr_ratio(cond,step,:) = pre_curr_mean(cond,step,:)./pre_curr_mean(1,1,:)*100;
        else % LICI 100
            temp = pre_curr_mean(cond,step,:)./pre_curr_mean(1,step,:)*100;% absolute response
            curr_ratio(cond,step,:) = -(100 - temp); % inhibition
        end
    end
end

curr_mean(4,:)=[];curr_mean(2,:)=[];
% curr_dev(4,:)=[];curr_dev(2,:)=[];
% curr_ratio_conc(:,4,:)=[];curr_ratio_conc(:,2,:)=[];
curr_ratio(4,:,:)=[];curr_ratio(2,:,:)=[];
% curr_ratio_dev(4,:)=[];curr_ratio_dev(2,:)=[];
curr_ratio_mean = squeeze(nanmean(curr_ratio,3));

figure;
bar(curr_ratio_mean, 1,'edgecolor','k', 'linewidth', 2);
hold on;
for col = 1:size(curr_mean,2)
    %             std_err = curr_dev(:,col)/sqrt(size(curr_ratio_conc,1));
%     std_err = curr_ratio_dev(:,col)/sqrt(size(curr_data,4));
%     errorbar(x(col,1:2),curr_ratio_mean(:,col),std_err,'.','LineWidth',2,'Color','k');
end

colormap(bone);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);





%__________________________________________________________________
for sub = 1:size(curr_data_all,4)
    h=figure;
    for PAS = 1:2
        %_________________CALCULATING__________________
        curr_data=squeeze(curr_data_all(:,:,:,sub,PAS));
        alter_mean = squeeze(mean(curr_data_all(:,:,:,sub,:),1));
        curr_mean=squeeze(mean(curr_data,1)); % cond, step
        curr_mean(2,:) = []; curr_mean(3,:) = []; % keeping only single and LICI100
        
        curr_dev = squeeze(std(curr_data,1));
        curr_dev(2,:)=[];curr_dev(3,:)=[];
        
        curr_ratio=zeros(size(curr_mean));
        for step = 1:size(curr_mean,2)
            for cond = 1:size(curr_mean,1)
                if cond == 1 % Single Pulse
                    curr_ratio(cond,step) = curr_mean(cond,step)./curr_mean(1,1)*100;
                else % LICI 100 
                    temp = curr_mean(cond,step)./curr_mean(1,1)*100; % absolute response
                    curr_ratio(cond,step) = -(100 - temp); % inhibition
                end
            end
        end
        subplot(2,1,PAS)
        bar(curr_ratio, 1,'edgecolor','k', 'linewidth', 2);
        hold on;
%         for col = 1:size(curr_ratio,2)
%             %                 std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
%             errorbar(x(col,1:2),curr_ratio(:,col),curr_dev(1:2,col),'.','LineWidth',2,'Color','k');
%         end
        max_sub=max(max(max(alter_mean)));
        lim_curr= max_sub+(max_sub*0.5);
%         ylim([0 lim_curr]);
        colormap(bone)
    end
end

%% MEAN PAS EFFECT ON SINGLE AND LICI - LICI AS % INHIBITION

% ____RATIO + LICI AS % INHIBITION - SINGLE - % calculated from mean ____

curr_data = squeeze(MEP_amp_GP(:,m,:,:,5,:)); % trials, cond, step, subject
pre_curr_mean=squeeze(nanmean(curr_data,1)); % 4, 3, 9 = cond, step, sub
curr_mean = squeeze(nanmean(pre_curr_mean,3)); % cond x step

% This loop allows to not count in the calculation people who didn't get a
% T1 (otherwise it brings down the mean) + I added nanmean everywhere
% instead of mean 
for sub = 1:size(curr_data,3)
    for step = 1:size(curr_data,2)
        testingEx=sum(all(curr_data(:,step,sub),2));
        if testingEx == 0
            curr_data(:,step,sub)=NaN;
        end
    end
end

curr_ratio=zeros(size(pre_curr_mean));
for step = 1:size(pre_curr_mean,2)
    for cond = 1:size(pre_curr_mean,1)
%         if cond == 1 % Single Pulse
            curr_ratio(cond,step,:) = pre_curr_mean(cond,step,:)./pre_curr_mean(1,1,:)*100; % absolute response
%         end
    end
end
curr_ratio(:,:,1)=[];
curr_ratio_mean = squeeze(nanmean(curr_ratio,3));
curr_ratio_mean(4,:)=[];curr_ratio_mean(2,:)=[]; % 2x3
curr_ratio_mean_INHIB= zeros(3,3) ; % Single LICI absolute LICI inhibition
curr_ratio_mean_INHIB(1:2,:) = curr_ratio_mean;
curr_ratio_mean_INHIB(3,:) =  -(100-curr_ratio_mean(2,:));

figure;
bar(curr_ratio_mean_INHIB, 1,'edgecolor','k', 'linewidth', 2);
hold on;

colormap(bone);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);


% ____RATIO + LICI AS % INHIBITION - STEP-WISE - % calculated from mean ____

curr_data = squeeze(MEP_amp_GP(:,m,:,:,5,:)); % trials, cond, step, subject
pre_curr_mean=squeeze(nanmean(curr_data,1)); % 4, 3, 9 = cond, step, sub
curr_mean = squeeze(nanmean(pre_curr_mean,3)); % cond x step

% This loop allows to not count in the calculation people who didn't get a
% T1 (otherwise it brings down the mean) + I added nanmean everywhere
% instead of mean 
for sub = 1:size(curr_data,3)
    for step = 1:size(curr_data,2)
        testingEx=sum(all(curr_data(:,step,sub),2));
        if testingEx == 0
            curr_data(:,step,sub)=NaN;
        end
    end
end

curr_ratio=zeros(size(pre_curr_mean));
for step = 1:size(pre_curr_mean,2)
    for cond = 1:size(pre_curr_mean,1)
        if cond == 1 % Single Pulse
            curr_ratio(cond,step,:) = pre_curr_mean(cond,step,:)./pre_curr_mean(1,1,:)*100; % absolute response
        else
            curr_ratio(cond,step,:) = pre_curr_mean(cond,step,:)./pre_curr_mean(1,step,:)*100; % absolute response
        end
    end
end
curr_ratio(:,:,1)=[];
curr_ratio_mean = squeeze(nanmean(curr_ratio,3));
curr_ratio_mean(4,:)=[];curr_ratio_mean(2,:)=[]; % 2x3
curr_ratio_mean_INHIB= zeros(3,3) ; % Single LICI absolute LICI inhibition
curr_ratio_mean_INHIB(1:2,:) = curr_ratio_mean;
curr_ratio_mean_INHIB(3,:) =  -(100-curr_ratio_mean(2,:));

figure;
bar(curr_ratio_mean_INHIB, 1,'edgecolor','k', 'linewidth', 2);
hold on;

colormap(bone);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);



%% (4) SAME SUBJECT - 0.2 to 5 HZ comparison

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

gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/5hz/',gp_ID];
cd(gp_folder)
list_folders_5=dir;
for k = length(list_folders_5):-1:1
    if ~list_folders_5(k).isdir
        list_folders_5(k) = [ ];
        continue
    end
    fname = list_folders_5(k).name;
    if fname(1) == '.' || fname(1) == '*'
        list_folders_5(k) = [ ];
    end
end

gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/02hz/',gp_ID];
cd(gp_folder)
list_folders_02=dir;
for k = length(list_folders_02):-1:1
    if ~list_folders_02(k).isdir
        list_folders_02(k) = [ ];
        continue
    end
    fname = list_folders_02(k).name;
    if fname(1) == '.' || fname(1) == '*'
        list_folders_02(k) = [ ];
    end
end

for k = 1:length(list_folders_5)
   names_5{k}= list_folders_5(k).name;
end

for k = 1:length(list_folders_02)
   names_02{k}= list_folders_02(k).name;
end


col_com = ismember(names_5,names_02);
common=names_5(col_com==1);

common_data=zeros(12,4,4,3,5,length(common),2); % to store for each common sub in the first column 02 and in the second 5 hz
for sub = 1:length(common)
    % __________________LOADING 02 HZ DATA__________________
    gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/02hz/',gp_ID];
    cd(gp_folder)
    sub_ID = common{sub};
    cd(sub_ID); cd('Data');
    load Resp_Amp
    common_data(:,1:size(resp,2),1:size(resp,3),:,:,sub,1)=resp;
    
    % __________________LOADING 5 HZ DATA__________________
    gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/5hz/',gp_ID];
    cd(gp_folder)
    sub_ID = common{sub};
    cd(sub_ID); cd('Data');
    load Resp_Amp
    common_data(:,1:size(resp,2),1:size(resp,3),:,:,sub,2)=resp; % trials, muscle, condition, step, param, sub, PASprot
end

curr_data_all = squeeze(common_data(:,m,:,:,5,:,:));%trials, condition, step, sub, PASprot

for sub = 1:size(curr_data_all,4)
    h=figure;
    for PAS = 1:2
        %_________________CALCULATING__________________
        curr_data=squeeze(curr_data_all(:,:,:,sub,PAS));
        alter_mean = squeeze(mean(curr_data_all(:,:,:,sub,:),1));
        curr_mean=squeeze(mean(curr_data,1)); % cond, step
        curr_mean(2,:) = []; curr_mean(3,:) = []; % keeping only single and LICI100
        
        curr_dev = squeeze(std(curr_data,1));
        curr_dev(2,:)=[];curr_dev(3,:)=[];
        
        %         curr_ratio=zeros(size(curr_mean));
        %         for step = 1:size(curr_mean,2)
        %             for cond = 1:size(curr_mean,1)
        %                 curr_ratio(cond,step) = curr_mean(cond,step)./curr_mean(1,1);
        %             end
        %         end
        subplot(2,1,PAS)
        bar(curr_mean, 1,'edgecolor','k', 'linewidth', 2);
        hold on;
        for col = 1:size(curr_mean,2)
            %                 std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
            errorbar(x(col,1:2),curr_mean(:,col),curr_dev(1:2,col),'.','LineWidth',2,'Color','k');
        end
        max_sub=max(max(max(alter_mean)));
        lim_curr= max_sub+(max_sub*0.5);
%         ylim([0 lim_curr]);
        colormap(bone)
    end
end



% ________________________RATIO________________________


for sub = 1:size(curr_data_all,4)
    h=figure;
    for PAS = 1:2
        %_________________CALCULATING__________________
        curr_data=squeeze(curr_data_all(:,:,:,sub,PAS));
        alter_mean = squeeze(mean(curr_data_all(:,:,:,sub,:),1));
        curr_mean=squeeze(mean(curr_data,1)); % cond, step
        curr_mean(2,:) = []; curr_mean(3,:) = []; % keeping only single and LICI100
        
        curr_dev = squeeze(std(curr_data,1));
        curr_dev(2,:)=[];curr_dev(3,:)=[];
        
        curr_ratio=zeros(size(curr_mean));
        for step = 1:size(curr_mean,2)
            for cond = 1:size(curr_mean,1)
                curr_ratio(cond,step) = curr_mean(cond,step)./curr_mean(1,1);
            end
        end
        subplot(2,1,PAS)
        bar(curr_ratio, 1,'edgecolor','k', 'linewidth', 2);
        hold on;
%         for col = 1:size(curr_ratio,2)
%             %                 std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
%             errorbar(x(col,1:2),curr_ratio(:,col),curr_dev(1:2,col),'.','LineWidth',2,'Color','k');
%         end
        max_sub=max(max(max(alter_mean)));
        lim_curr= max_sub+(max_sub*0.5);
%         ylim([0 lim_curr]);
        colormap(bone)
    end
end

% ________________________RATIO + LICI AS % INHIBITION ________________________


for sub = 1:size(curr_data_all,4)
    h=figure;
    for PAS = 1:2
        %_________________CALCULATING__________________
        curr_data=squeeze(curr_data_all(:,:,:,sub,PAS));
        alter_mean = squeeze(mean(curr_data_all(:,:,:,sub,:),1));
        curr_mean=squeeze(mean(curr_data,1)); % cond, step
        curr_mean(2,:) = []; curr_mean(3,:) = []; % keeping only single and LICI100
        
        curr_dev = squeeze(std(curr_data,1));
        curr_dev(2,:)=[];curr_dev(3,:)=[];
        
        curr_ratio=zeros(size(curr_mean));
        for step = 1:size(curr_mean,2)
            for cond = 1:size(curr_mean,1)
                if cond == 1 % Single Pulse
                    curr_ratio(cond,step) = curr_mean(cond,step)./curr_mean(1,1)*100;
                else % LICI 100 
                    temp = curr_mean(cond,step)./curr_mean(1,step)*100; % absolute response
                    curr_ratio(cond,step) = -(100 - temp); % inhibition
                end
            end
        end
        
        curr_plot = zeros(3,3);
        curr_plot(1,:)
        %___________________FIGURE____________________
        subplot(2,1,PAS)
        bar(curr_ratio, 1,'edgecolor','k', 'linewidth', 2);
        hold on;
%         for col = 1:size(curr_ratio,2)
%             %                 std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
%             errorbar(x(col,1:2),curr_ratio(:,col),curr_dev(1:2,col),'.','LineWidth',2,'Color','k');
%         end
        max_sub=max(max(max(alter_mean)));
        lim_curr= max_sub+(max_sub*0.5);
%         ylim([0 lim_curr]);
        colormap(bone)
    end
end

% ________________________RATIO + LICI AS % INHIBITION of baseline single ________________________


for sub = 1:size(curr_data_all,4)
    h=figure;
    for PAS = 1:2
        %_________________CALCULATING__________________
        curr_data=squeeze(curr_data_all(:,:,:,sub,PAS));
        alter_mean = squeeze(mean(curr_data_all(:,:,:,sub,:),1));
        curr_mean=squeeze(mean(curr_data,1)); % cond, step
        curr_mean(2,:) = []; curr_mean(3,:) = []; % keeping only single and LICI100
        
        curr_dev = squeeze(std(curr_data,1));
        curr_dev(2,:)=[];curr_dev(3,:)=[];
        
        curr_ratio=zeros(size(curr_mean));
        for step = 1:size(curr_mean,2)
            for cond = 1:size(curr_mean,1)
                if cond == 1 % Single Pulse
                    curr_ratio(cond,step) = curr_mean(cond,step)./curr_mean(1,1)*100;
                else % LICI 100 
                    temp = curr_mean(cond,step)./curr_mean(1,1)*100; % absolute response
                    curr_ratio(cond,step) = -(100 - temp); % inhibition
                end
            end
        end
        subplot(2,1,PAS)
        bar(curr_ratio, 1,'edgecolor','k', 'linewidth', 2);
        hold on;
%         for col = 1:size(curr_ratio,2)
%             %                 std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
%             errorbar(x(col,1:2),curr_ratio(:,col),curr_dev(1:2,col),'.','LineWidth',2,'Color','k');
%         end
        max_sub=max(max(max(alter_mean)));
        lim_curr= max_sub+(max_sub*0.5);
%         ylim([0 lim_curr]);
        colormap(bone)
    end
end

%% (5) STATISTICS

m=1;
curr_data = squeeze(MEP_amp_GP(:,m,:,:,5,:)); % trials, cond, step, subject
% curr_single = squeeze(all_single_GP(:,m,:,5,sub)); % trials, step, subject

for sub = 1:size(curr_data,3)
    for step = 1:size(curr_data,2)
        testingEx=sum(all(curr_data(:,step,sub),2));
        if testingEx == 0
            curr_data(:,step,sub)=NaN;
        end
    end
end

% Two-tail t test
single=squeeze(curr_ratio_conc(:,1,:));

testRes=size(2,3);
[testRes(1,1), testRes(2,1)] = ttest2(single(:,1), single(:,2), 'Tail','both');
[testRes(1,2), testRes(2,2)] = ttest2(single(:,2), single(:,3), 'Tail','both');
[testRes(1,3), testRes(2,3)] = ttest2(single(:,1), single(:,3), 'Tail','both');


% trials_step = nanmean(curr_single,3);

% KRUSKAL-WALLIS - One Way - 2 or more
p=kruskalwallis(trials_step);
p=kruskalwallis(curr_ratio(1,:));

% TWO WAY ANOVA - Anova2 - only suitabke for same size samples 
temp_mean = nanmean(curr_data,4); % trials, cond, step

cond_mean=zeros(36,4);
cond_mean(1:12,:)=squeeze(temp_mean(:,:,1));
cond_mean(13:24,:)=squeeze(temp_mean(:,:,2));
cond_mean(25:36,:)=squeeze(temp_mean(:,:,3));

p=anova2(cond_mean,12);

n_sub = size(curr_data,4);
n_trials = size(curr_data,1);
n_cond = size(curr_data,2);

% temp_mean = squeeze(mean(curr_data,1));
% 
% cond_mean=zeros(n_sub*n_trials*n_cond,3);
% count = 0;
% for cond = 1:n_cond
%     for sub = 1:n_sub
%         end_t=n_trials*sub+count;
%         start_t=end_t-11;
%         cond_mean(start_t:end_t,:)=squeeze(curr_data(:,cond,:,sub));
%     end
%     count=count+84;
% end
% 
% p=anova2(cond_mean,n_trials*n_sub);

curr_data = squeeze(MEP_amp_GP(:,m,:,:,5,:)); % trials, cond, step, subject
temp_mean = squeeze(nanmean(curr_data,1)); % cond, step, subject

n_sub = size(temp_mean,3);
n_cond = size(temp_mean,1);

cond_mean=zeros(n_sub*n_cond,3);
for cond = 1:n_cond
    end_t=cond*n_sub;
    start_t=end_t-(n_sub-1);
    cond_mean(start_t:end_t,:)=squeeze(temp_mean(cond,:,:))';
end

[pa1,tbla1,statsa1]=anova2(cond_mean,n_sub);


% without non-resp ELB - only single and LICI 100

temp_mean = squeeze(nanmean(curr_data,1)); % cond, step, subject
temp_mean(:,:,2)=[];
temp_mean(4,:,:)=[];
temp_mean(2,:,:)=[];

n_sub = size(temp_mean,3);
n_cond = size(temp_mean,1);

cond_mean=zeros(n_sub*n_cond,3);
for cond = 1:n_cond
    end_t=cond*n_sub;
    start_t=end_t-(n_sub-1);
    cond_mean(start_t:end_t,:)=squeeze(temp_mean(cond,:,:))';
end

[pa1,tbla1,statsa1]=anova2(cond_mean,n_sub);
[c,m] = multcompare(statsa1,'Dimension',[1 2]);

% Two-tail t test
single=squeeze(curr_ratio(1,:));
for step = 1:length(single)
[h,p,ci,stats] = ttest2(single(), y, 'Tail','both');

end



%%  (6) _________________ALL SUBJECTS SINGLE PULSE_____________________


single_timewin = time_windows(:,1);
n_subjects = size(data_cond_GP,6);

% all in one figure one plot

h=figure;
for n = 1:n_subjects
n=4;
    for trial = 1:size(data_cond_GP,2)-1
        if trial == 3 %|| trial == 11 || trial == 12
        else
        curr_data=data_cond_GP(:,trial,1,1,1,n);
        curr_time=time_axis(1:size(curr_data,1));
        plot(time_axis,curr_data, 'LineWidth', 2);
        set(gca,'FontSize',22)
        hold on;
                end
        set (gca,'Ydir','reverse');
        xlim([single_timewin(1)-5 single_timewin(2)]);
        xlim([-5 200]);
        ylim([-2 1]);
    end
end

h=figure;
for n = 1:n_subjects
n=5;
    for trial = 1:size(data_cond_GP,2)
        if trial == 3 %|| trial == 11 || trial == 12
        else
        curr_data=data_cond_GP(:,trial,1,3,1,n);
        curr_time=time_axis(1:size(curr_data,1));
        plot(time_axis,curr_data, 'LineWidth', 2);
        set(gca,'FontSize',22)
        hold on;
                end
        set (gca,'Ydir','reverse');
        xlim([single_timewin(1)-5 single_timewin(2)]);
        xlim([35 65]);
        ylim([-2 1]);
    end
end

% ______________All in one figure each sub a plot__________________
ind = sqrt(n_subjects);

if ~mod(ind,1) % if number of subject is an integer
    row = ind;
    col = row;
else
    if ind < floor(ind)+0.5
        row = floor(ind);
        col = row +1;
    else
        row = floor(ind)+1;
        col = row;
    end
end

h=figure;
for n = 1:n_subjects
    subplot(row,col,n)
    curr_sub=data_cond_GP(:,:,1,1,1,n);
    for trial = 1:size(curr_sub,2)
        curr_data=curr_sub(:,trial);
        curr_time=time_axis(1:size(curr_data,1));
        plot(time_axis,curr_data);
        hold on;
    end
    
    set (gca,'Ydir','reverse');
    xlim([single_timewin(1)-15 single_timewin(2)+15]);
    %     ylim([-3 2]);
end


%% _____________CHECK DIFFERENT SINGLE RESPONSES___________

cd(gp_folder);
load('all_single_GP.mat');

all_singleAmp_GP = all_single_GP(:,:,:,5,:);

singleSep=zeros(12,4,4,3,5,7);
for curr_cond = 1:4
    end_trial = 12*curr_cond;
    start_trial = end_trial-11;
    singleSep(1:12,:,curr_cond,:,:,:)=all_single_GP(start_trial:end_trial,:,:,:,:);
    
end

hist=figure;count=1;
for m = 1:2
    %_________________CALCULATING__________________
    curr_data = squeeze(singleSep(:,m,:,:,5,:)); % trials, cond, step, subject
    curr_mean = squeeze(mean(mean(curr_data,1),4)); % cond, step
    
    curr_data_conc=curr_data(:,:,:,1);
    for k = 2:size(curr_data,4) % concatenating trials of all subjects in one same column
        addData=curr_data(:,:,:,k);
        curr_data_conc=cat(1,curr_data_conc,addData);
    end
    curr_dev = squeeze(std(curr_data_conc,1)); % cond, step
    
    %________________CREATING FIGURE_________________
    max_V=max(max(max(curr_dev)))+ max(max(max(curr_mean)));
    subplot(1,2,count)
    %     x=[0.78,1.78,2.78,3.78 ;...
    %         1,2,3,4; 1.22,2.22,3.22,4.22];
    x=[0.78,1.78,2.78,3.78 ;...
        1,2,3,4; 1.22,2.22,3.22,4.22];
    
    bar(curr_mean', 1,'edgecolor','k', 'linewidth', 2);
    hold on;
%     for row = 1:size(curr_mean,1)
%         std_err = curr_dev(row,:)/sqrt(size(all_single_GP,1));
%         errorbar(x(row,:),curr_mean(row,:),std_err,'.','LineWidth',2,'Color','k');
%     end
    %  ylim([0 (max_V+0.1)]);
        
    %_____________GRAPH PARAMETERS______________
    colormap(bone);
    hold on
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'NextPlot','add');
    axes;
    count = count +1;
end

%% ONLY SINGLE PULSE - effect of PAS (all single counted)

cd(gp_folder);
load('all_single_GP.mat');
all_singleAmp_GP = squeeze(all_single_GP(:,:,:,5,:)); % 48     4     3     5     7
meanAmp = squeeze(mean(all_singleAmp_GP,1)); % muscle, step, subjects

n_subjects = size(meanAmp,3);
ind = sqrt(n_subjects);

if ~mod(ind,1) % if number of subject is an integer
    row = ind;
    col = row;
else
    if ind < floor(ind)+0.5
        row = floor(ind);
        col = row +1;
    else
        row = floor(ind)+1;
        col = row;
    end
end

figure; count=1; 
ylimup=max(max(max(meanAmp(1:2,:,:))));
for sub = 1:n_subjects
    subplot(row,col,count)
    bar(meanAmp(1:2,:,sub), 1,'edgecolor','k', 'linewidth', 2);
    hold on;
    colormap(bone);
    hold on
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'NextPlot','add');
    ylim([0 ylimup]);
    count=count+1;
end

% group mean
figure; count=1; 
ylimup=max(max(max(meanAmp(1:2,:,:))));
curr = squeeze(all_single_GP(:,1,:,5,:));
curr_mean= squeeze(mean(mean(all_single_GP(:,1,:,5,:),1),5));
curr_data_conc=curr(:,:,1);
for k = 2:size(curr,3) % concatenating trials of all subjects in one same column
    addData=curr(:,:,k);
    curr_data_conc=cat(1,curr_data_conc,addData);
end
curr_dev = squeeze(std(curr_data_conc,1)); % cond, step

subplot(1,2,1);
figure;
bar(curr_mean, 1,'edgecolor','k', 'linewidth', 2);
hold on;
x=[0.78,1.78,2.78,3.78 ;...
        1,2,3,4; 1.22,2.22,3.22,4.22];
% for col = 1:size(curr_mean,2)
    std_err = curr_dev/sqrt(size(all_single_GP,1));
    errorbar(x(:,1),curr_mean,std_err,'.','LineWidth',2,'Color','k');
% end
colormap(bone);
hold on
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'NextPlot','add');
ylim([0 ylimup]);




% RATIO

%all_singleAmp_GP [48,4,3,7]
ratio=zeros(size(meanAmp));
dev = squeeze(std(all_singleAmp_GP,1)); % cond, step,sub

for step = 1:size(meanAmp,2)
    ratio(:,step,:) = meanAmp(:,step,:)./meanAmp(:,1,:);
    ratio_dev(:,step,:) = dev(:,step,:)./dev(:,1,:);
end

x=[0.78,1.78,2.78,3.78 ;1,2,3,4 ; 1.22,2.22,3.22,4.22];

figure; count=1; 
ylimup=max(max(max(ratio(1:2,:,:))))+max(max(max(ratio_dev)));
ampStr={};
for sub = 1:n_subjects
    subplot(row,col,count)
    bar(ratio(1:2,:,sub), 1,'edgecolor','k', 'linewidth', 2);
    hold on;
    for col = 1:size(ratio,2)
        errorbar(x(col,1:2),ratio(1:2,col,sub),ratio_dev(1:2,col,sub),'.','LineWidth',2,'Color','k');
        ampStr{1} = sprintf('%.3f',ratio(1,col,sub)); ampStr{2} = sprintf('%.3f',ratio(2,col,sub));
        hold on;
        text(x(col,1:2)'-0.05, (ratio(1:2,col,sub)+ratio_dev(1:2,col,sub))+0.5 , ampStr )
    end
    colormap(bone);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'NextPlot','add');
    ylim([0 ylimup]);
    count=count+1;
end

% peak_str_single = sprintf('%.3f',amp_single);
% hold on; text(single_timewin(1)+2,peak_single-0.6,peak_str_single,'FontSize',6)

%% ILLUSTRATION OF PAS EFFECT

m=1; %ADM
sub = 6; %ETB

curr_data = squeeze(data_cond_GP(:,:,m,1,:,sub)); % [datapt,trials, step]
points=size(curr_data,1);
duration_ms=points/s_rate*1000;
step=duration_ms/points;
time_axis=0:step:duration_ms-step; time_axis=time_axis';


ylimup = abs(squeeze(min(min(min(curr_data)))));

for step = 1:size(curr_data,3)
    %     subplot(3,1,step)
    figure;
    for trial = 1:size(curr_data,2)
        plot(time_axis, curr_data(:,trial,step), 'LineWidth', 2);
        hold on;
    end
    set (gca,'Ydir','reverse');
    ylim([-3.4 2.]);
    xlim([-5 200]);
end


%% HISTOGRAM, ONLY SINGLE, PAS EFFECT - DOESNT WOOOOORK 

gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
cd(gp_folder)

fileName='Resp_Amp.mat';
load(fileName)
% singleFile = 'all_single_GP.mat';
% load(singleFile);

m=1;

hist=figure;count=1;
for MyUnit = 1:2 % 1 - in mV 2- ratio
    
    %_________________CALCULATING__________________
    curr_data = squeeze(MEP_amp_GP(:,m,1,:,5,:)); % trials, step, subject
%     curr_single = squeeze(all_single_GP(:,m,:,5,:)); % trials, step, subject
    %
    %         curr_mean=zeros(4,3);
    %         curr_mean(1,:)=squeeze(mean(mean(curr_single,1),3));
    %         curr_mean(2:4,:) = squeeze(mean(mean(curr_data(:,2:4,:,:),1),4)); % cond, step
    
    curr_mean = squeeze(mean(mean(curr_data,1),3));
    
    curr_data_conc=curr_data(:,:,1);
    for k = 2:size(curr_data,3)
        addData=curr_data(:,:,k);
        curr_data_conc=cat(1,curr_data_conc,addData);
    end
    
    curr_dev = squeeze(std(curr_data_conc,1)); % cond, step
    curr_ratio=zeros(size(curr_mean));
    for step = 1:size(curr_mean,2)
        for cond = 1:size(curr_mean,1)
            curr_ratio(cond,step,:) = curr_mean(cond,step,:)./curr_mean(1,1,:);
        end
    end
    
    %________________CREATING FIGURE_________________
    max_V=max(max(max(curr_dev)))+ max(max(max(curr_mean)));
    max_R=max(max(max(curr_ratio)));
    subplot(1,2,count)
    
    x=[1,2,3];
    
    if MyUnit == 1 % PLOT IN MV WITH ERROR
        bar(curr_mean', 1,'edgecolor','k', 'linewidth', 2);
        hold on;
        for col = 1:size(curr_mean,2)
            std_err = curr_dev/sqrt(size(curr_data_conc,1));
            errorbar(x',curr_mean',std_err','.','LineWidth',2,'Color','k');
        end
        %  ylim([0 (max_V+0.1)]);
        unitStr='mV';
    else if MyUnit == 2 % PLOT IN RATIO
            bar(curr_ratio, 1,'edgecolor','k', 'linewidth', 2);
            ylim([0 (max_R+0.1)]);
            unitStr='ratio';
        end
    end
    
    %_____________GRAPH PARAMETERS______________
    colormap(bone);
    hold on
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%     set(gcf,'NextPlot','add');
%     axes;
    count = count +1;
    
end


%% All subject of a group - PAS EFFECT ON SINGLE AND LICI 100 - ALL SINGLE INCLUDED


gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
cd(gp_folder)

fileName='Resp_Amp.mat';
load(fileName)
singleFile = 'all_single_GP.mat';
load(singleFile);

m=1;
hist=figure;count=1;
MyUnit = 2;
for sub = 1:size(MEP_amp_GP,6)
    if sub == 2
    else
        subplot(1,6,count)
        %_________________CALCULATING__________________
        curr_data = squeeze(MEP_amp_GP(:,m,:,:,5,sub)); % trials, cond, step, subject
        curr_single = squeeze(all_single_GP(:,m,:,5,sub)); % trials, step, subject
        
        curr_mean=zeros(4,3);
        curr_mean(1,:)=squeeze(mean(mean(curr_single,1),3));
        curr_mean(2:4,:) = squeeze(mean(mean(curr_data(:,2:4,:,:),1),4)); % cond, step
        
        curr_mean(2,:) = []; curr_mean(3,:) = [];
        
        curr_data_conc=curr_data(:,:,:,1);
        for k = 2:size(curr_data,4)
            addData=curr_data(:,:,:,k);
            curr_data_conc=cat(1,curr_data_conc,addData);
        end
        
        curr_dev = squeeze(std(curr_data_conc,1)); % cond, step
        curr_ratio=zeros(size(curr_mean));
        for step = 1:size(curr_mean,2)
            for cond = 1:size(curr_mean,1)
                curr_ratio(cond,step,:) = curr_mean(cond,step,:)./curr_mean(1,1,:);
            end
        end
        
        %________________CREATING FIGURE_________________
        %     max_V=max(max(max(curr_dev)))+ max(max(max(curr_mean)));
        %     max_R=max(max(max(curr_ratio))) + max(max(max(curr_dev)));
        %     subplot(1,2,count)
        
        x=[0.78,1.78,2.78,3.78 ;...
            1,2,3,4; 1.22,2.22,3.22,4.22];
        
        if MyUnit == 1 % PLOT IN MV WITH ERROR
            bar(curr_mean, 1,'edgecolor','k', 'linewidth', 2);
            hold on;
            for col = 1:size(curr_mean,1)
%                 std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
                errorbar(x(col,:),curr_mean(:,col),std_err,'.','LineWidth',2,'Color','k');
            end
            %  ylim([0 (max_V+0.1)]);
            unitStr='mV';
        else if MyUnit == 2 % PLOT IN RATIO
                bar(curr_ratio, 1,'edgecolor','k', 'linewidth', 2);
                hold on;
                for col = 1:size(curr_mean,2)
%                     std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
                    errorbar(x(col,1:2),curr_ratio(:,col),curr_dev(1:2,col),'.','LineWidth',2,'Color','k');
                end
                ylim([0 3]);
                unitStr='ratio';
            end
        end
        
        %_____________GRAPH PARAMETERS______________
        colormap(bone);
        hold on
        %     set(gca,'XTickLabel',cond_name);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        %     set(gcf,'NextPlot','add');
        %     axes;
        count = count +1;
        % end
    end
end





%% INDIVIDUAL MEAN MEP AMPLITUDE + GROUP MEAN (Leo style) - 0.2 Hz

curr_single = squeeze(all_single_GP(:,1,:,5,:)); % trials, step, subject
curr_mean = squeeze(mean(curr_single,1)); % step, subject
all_mean = squeeze(mean(curr_mean,2));

curr_mean_ratio = zeros(size(curr_mean));
for step = 1:size(curr_mean,1)
    curr_mean_ratio(step,:) = curr_mean(step,:)./curr_mean(1,:);
end
curr_mean_ratio = curr_mean_ratio*100;

all_mean_ratio = squeeze(mean(curr_mean_ratio,2));

% RAW
figure;
for step = 1:size(curr_mean,2)
    line([1,2,3],curr_mean(:,step));
    hold on
end
line([1,2,3],all_mean);

% RATIO - baseline, t1, t2
figure;
for sub = 1:size(curr_mean_ratio,2)
    line([1,2,3],curr_mean_ratio(:,sub), 'Color','k', 'LineWidth', 2);
    hold on;
end
line([1,2,3],all_mean_ratio, 'Color', 'r', 'LineWidth', 3);
hold on;
mean_str=int2str(round(all_mean_ratio));
text([1,2,3]-0.01, all_mean_ratio+5, mean_str, 'Color', 'r', 'FontSize', 14);

% RATIO - baseline, t2
curr_mean_ratio_2=curr_mean_ratio;
curr_mean_ratio_2(2,:)=[];
all_mean_ratio_2 = squeeze(mean(curr_mean_ratio_2,2));
figure;
for sub = 1:size(curr_mean_ratio,2)
    line([1,2],curr_mean_ratio_2(:,sub), 'Color','k', 'LineWidth', 2);
    hold on;
    mean_str=int2str(round(curr_mean_ratio_2(:,sub)));
    text([1,2]-0.01, curr_mean_ratio_2(:,sub)+5, mean_str, 'Color', 'k', 'FontSize', 14);
    hold on;
end
line([1,2],all_mean_ratio_2, 'Color', 'r', 'LineWidth', 3);
hold on;
mean_str=int2str(round(all_mean_ratio_2));
text([1,2]-0.01, all_mean_ratio_2+3, mean_str, 'Color', 'r', 'FontSize', 14);

% RATIO - baseline, t1
curr_mean_ratio_2=curr_mean_ratio;
curr_mean_ratio_2(3,:)=[];
all_mean_ratio_2 = squeeze(mean(curr_mean_ratio_2,2));
figure;
for sub = 1:size(curr_mean_ratio,2)
    line([1,2],curr_mean_ratio_2(:,sub), 'Color','k', 'LineWidth', 2);
    hold on;
    mean_str=int2str(round(curr_mean_ratio_2(:,sub)));
    text([1,2]-0.01, curr_mean_ratio_2(:,sub)+5, mean_str, 'Color', 'k', 'FontSize', 14);
    hold on;
end
line([1,2],all_mean_ratio_2, 'Color', 'r', 'LineWidth', 3);
hold on;
mean_str=int2str(round(all_mean_ratio_2));
text([1,2]-0.01, all_mean_ratio_2+3, mean_str, 'Color', 'r', 'FontSize', 14);
ylim([0 500]);


%% BOXPLOTS

% SINGLE PULSE RAW - Baseline, T1, T2
curr_data = squeeze(MEP_amp_GP(:,1,1,:,5,:)); % trials, steps, subjects - first 1 = ADM, second 1 = single pulse, 5 = amplitude

% All subjects
figure; count=1;
for sub = 1:size(curr_data,3)
    subplot(1,size(curr_data,3),count)
    curr_data_sub = squeeze(curr_data(:,:,sub));
    boxplot(curr_data_sub);
    count=count+1;
end

% Group mean
curr_data = squeeze(MEP_amp_GP(:,1,1,:,5,:));% trials, steps, subjects - first 1 = ADM, second 1 = single pulse, 5 = amplitude
mean_data = squeeze(mean(curr_data,3));
figure;
boxplot(mean_data);



%% NOTES
curr_data = squeeze(MEP_amp_GP(:,m,:,:,5,:)); % trials, cond, step, subject
curr_mean=squeeze(mean(mean(curr_data,1),4));

curr_ratio=zeros(size(curr_mean));
for step = 1:size(curr_mean,2)
    for cond = 1:size(curr_mean,1)
        curr_ratio(cond,step,:) = curr_mean(cond,step,:)./curr_mean(1,1,:);
    end
end

curr_data = squeeze(MEP_amp_GP(:,1,1,:,5,:)); % trials, step, subject - first 1 = ADM, second 1 = single pulse, 5 = amplitude
curr_mean = squeeze(mean(curr_data,1)); % step, subject
all_mean = squeeze(mean(curr_mean,2));
curr_mean_ratio = zeros(size(curr_mean));
for step = 1:size(curr_mean,1)
    curr_mean_ratio(step,:) = curr_mean(step,:)./curr_mean(1,:);
end
all_mean_ratio = squeeze(mean(curr_mean_ratio,2));

% curr_mean_ratio = curr_mean_ratio*100;




end











