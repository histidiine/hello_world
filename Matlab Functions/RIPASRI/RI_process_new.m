function RI_process_new

%% Translating timepoints into milliseconds

% rep=input('Is the frequency of acquisition 20000 (y/n - into simple quotes) ?) ')
% if rep == 'y' || 'yes'; freq = 20000; else freq=input('what is the frequency of acquisition? '); end
% rep2=input('Is the number of points 5000 (y/n - into simple quotes) ?) ')
% if rep2 == 'y' || 'yes'; points = 5000; else points=input('what is the frequency of acquisition? '); end

freq=20000; points=5000;
time_ms = points/freq*1000;
%% Defining Subject and Folders

sub_id='AuS';
main_dossier=['M:\NLG\STIM_LAB\16_CRPS-Dys\Data\RI\', sub_id];
nguyet_dossier= ['M:\NLG\STIM_LAB\16_CRPS-Dys\Data\RI\', sub_id, '\Nguyet_extract'];
cond_dossier=['M:\NLG\STIM_LAB\16_CRPS-Dys\Data\RI\', sub_id, '\Nguyet_extract\Cond'];
ext_dossier=['M:\NLG\STIM_LAB\16_CRPS-Dys\Data\RI\', sub_id, '\Nguyet_extract\Ext'];
fig_folder = [main_dossier, '\Fig'];
data_folder = [main_dossier, '\ML_data'];

% nguyet_dossier = '/Users/aureliestephan/Documents/MATLAB/Data/Nguyet extract';

%% Getting data from files extracted with Nguyet
cd(nguyet_dossier)

% Initialisation
num = 150; % number of trials
ri_data=zeros(points,2,num);

% Extraction loop - loads file 1 buy 1 and add the data in the ri_data matrix
for i = 1:num
    %if exist(int2str(i),'file') == 7
   temp=importdata(int2str(i));
    %else
      % disp(['Programmed paused because file ', int2str(i),' is missing'])
        %pause
    %end
   ri_data(:,1,i) = temp(1:points,1); % I put row 1:points in case there was a mistake during the nguyet and a file was recorded twice
   ri_data(:,2,i) = temp(1:points,2); clear temp;
end

% Separating FCR from EDC (we are mostly interested in FCR)
FCR_data = squeeze(ri_data(:,1,:)); % FCR - [data x cond]
EDC_data = squeeze(ri_data(:,2,:));

% Get information on condition - from log
% disp('First check that the log file was copied in a .xls file and the first lines removed')
% pause;
% [~, ~, test1matlab] = xlsread([nom_dossier,'\test_1_matlab.xlsx'],'Feuil1');
% test1matlab(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),test1matlab)) = {''};
% log=test1matlab; clear test1matlab
% cd('/Users/aureliestephan/Downloads')
cd(main_dossier)
log=xlsread('log.xlsx');
FCR_data_ms = ricampiona(EDC_data', time_ms);
FCR_data_ms = FCR_data_ms(:,1:250)';
FCR_cond_ms = zeros(time_ms,10,15);

EDC_data_ms = ricampiona(EDC_data', time_ms);
EDC_data_ms = EDC_data_ms(:,1:250)';
EDC_cond_ms = zeros(time_ms,10,15);

%% Extract data in conditions

FCR_cond = order_by_cond(FCR_data,log,num);
FCR_cond_ms = order_by_cond(FCR_data_ms,log,num);
EDC_cond_ms = order_by_cond(EDC_data,log,num);
EDC_cond = order_by_cond(EDC_data,log,num);

%% Plot data and remove bad recordings
% ML_data = '/Users/aureliestephan/Desktop';
% fig_name='EDC_';
% print_RI(EDC_cond,ML_data,fig_name)

fig_name='FCR_';
print_RI(FCR_cond,ML_data,fig_name)
cd('M:\NLG\STIM_LAB\16_CRPS-Dys\Data\RI\AuS\Test1');
for j = 1:15
    disp(['Processing condition', int2str(j), ' ...'])
    if j < 4
        plot_start = 1;
        plot_end = 45;
        data=FCR_cond_ms(plot_start:plot_end,:,j);
    else
        [val_max ind_max]=max(FCR_cond_ms(:,:,j));
        plot_start = 200;
        plot_end = 240;
        if plot_end < length(FCR_cond_ms)
            data=FCR_cond_ms(plot_start:plot_end,:,j);
        else
            data=FCR_cond_ms(plot_start:length(FCR_cond_ms),:,j);
        end
    end
    data=data.*-1;
    %xlswrite('Test_RI2_ms', data, j);
    for i = 1: 10
        h1=figure; plot(data(:,i));
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        fig_str=['RI cond ', int2str(j),' trial ', int2str(i)];
        title(fig_str); ylim([-10 15])
        print(h1,fig_str,'-dpng');
    end
    h2=figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i = 1:10
        plot(data(:,i));
        hold on;
    end
    ylim([-10 15])
    fig_str2=['RI_cond ', int2str(j)]; title(fig_str2);
    print(h2,fig_str2,'-dpng');
    %pause;
    clear plot_start plot_end data
    close all
end

%%
Bad_Trials={};
for j = 1:size(FCR_cond,3)
    bad_trials=input(['Bad trials for condition ', int2str(j), ' (vector):  '])
    Bad_Trials{j}=bad_trials;
end
FCR_cond_gtrials=zeros(points,10,15);
for j = 1:size(FCR_cond,3)
    temp_btrials=Bad_Trials{j};
    for i = 1:size(temp_btrials,2)
        good_trials=setdiff(1:10,temp_btrials);
        n_gtrials=size(good_trials,2)
        FCR_cond_gtrials(:,1:n_gtrials,j) = FCR_cond(:,good_trials,j);
    end
    clear temp_btrials good_trials n_gtrials
end

FCR_cond_gtrials = FCR_cond_gtrials.*-1; % To get the data with minus up and positive down
temp_data=zeros(5000,1);
for j = 1:size(FCR_cond_gtrials,3) % COND LOOP
    for i = 1:size(FCR_cond_gtrials,2)
        data_cond = FCR_cond_gtrials(:,i,j);
        test=zeros(5000,1);
        if data_cond ~= test
            temp_data(:,i)=data_cond;
        end
    end
    h=figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i=1:size(temp_data,2) % TRIAL LOOP
        fig_str = [fig_name, 'Condition number ', int2str(j), ' trial ', int2str(i)];
        if j < 4
            plot_start = 1;
            plot_end = 650;
        else
            [val_max ind_max]=max(temp_data(:,i));
            plot_start = ind_max-300;
            plot_end = size(temp_data,1);
        end
        data_zoom = temp_data(plot_start:plot_end,i,j);
        subplot (2,5,i); plot(data_zoom);
        title(fig_str);
        ylim([-15 15]); xlim([0 size(data_zoom,1)]);
        clear plot_start plot_end ind_max data_zoom
    end
    %print(h,fig_str,'-dpng');
    pause;
    close
end

%% Plot data grouped by blocks 



%% Find h-response
[EDC_hr, EDC_hr_mean] = hr_amp(EDC_cond);
[FCR_hr, FCR_hr_mean, ind_max_EDC, ind_min_FCR, maxmin_ECR] = hr_amp(FCR_cond);

%% Mean conditions


%% Double-checking the selected maximum
count=0;
for i = 1:ncond
    for j = 1:ntrial
        figure; 
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        plot(FCR_cond(:,j,i)); title(['Condition ', int2str(i), 'Trial number ', int2str(j)])
        hold on; line([all_ind_max(j,i) all_ind_max(j,i)],[-13 13],'Color',[1 0 0]);
        hold on; line([all_ind_min(j,i) all_ind_min(j,i)],[-13 13],'Color',[1 0 0]);
        descr = ' Verify that the maximum and minimum is properly detected and click right arrow' ;
        text(1500,9,descr)
        descr2 = ' If it is wrong, write the condition and trial number somewhere' ;
        text(1500,8,descr2)
        pause;
%         rep=input(' If it is right type y in between simple quotes else type n in simple quotes :   ');
%         if rep == 'y'
%         else
%             error(count,1)=i;error(count,2)=j;
%             count=count+1;
%         end
        close;
    end
    hold on 
end

for i = 1:ncond
        figure; 
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for j = 1:10
        plot(FCR_cond(:,j,i)); title(['Condition ', int2str(i), 'Trial number ', int2str(j)])
        hold on; line([all_ind_max(j,i) all_ind_max(j,i)],[-13 13],'Color',[1 0 0]);
        hold on; line([all_ind_min(j,i) all_ind_min(j,i)],[-13 13],'Color',[1 0 0]);
        descr = ' Verify that the maximum and minimum is properly detected and click right arrow' ;
        text(1500,9,descr)
        descr2 = ' If it is wrong, write the condition and trial number somewhere' ;
        text(1500,8,descr2)
        
%         rep=input(' If it is right type y in between simple quotes else type n in simple quotes :   ');
%         if rep == 'y'
%         else
%             error(count,1)=i;error(count,2)=j;
%             count=count+1;
%         end
        
    end
    hold on 
    pause;
    close;
end

%% Get a look at the data - plot all by condition

% generate random data - 32 channels
t = linspace(0,2,points);
sig = rand(10,points);

% calculate shift
mi = min(sig,[],2);
ma = max(sig,[],2);
shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]);
shift = repmat(shift,10,points);

%plot 'EMG' data of blocks
for j=1:14 % 1 plot avec 10 trials pour chaque condition
    figure;
    h1=plot(t,FCR_cond(:,:,j)'+shift,'Color', [0 0 .5]);
    ylim([-11 11])
    %set(h1, 'Ydir', 'reverse')
%     str = ['Block ', int2str(j)]
%     print(h1,str, 'd-png')
end



%% ---------------------------------------------------------------------
%                              CEMETERY
% ----------------------------------------------------------------------

% %% In miliseconds
% % Organise data by condition (ISI)
% 
% 
% count1=1;count2=1;count3=1;count4=1;count5=1;count6=1;count7=1;count8=1;count9=1;count10=1;count11=1;count12=1;count13=1;count14=1;count15=1;
% for i = 1:num
%     if log(i,7)==10
%         if log(i,9)==0
%             FCR_cond_ms(:,count3,3) =  FCR_data_ms(:,i);
%             count3=count3+1;
%         else if log(i,9)==-0.5
%                 FCR_cond_ms(:,count2,2) =  FCR_data_ms(:,i);
%                 count2=count2+1;
%             else if log(i,9)==-1
%                     FCR_cond_ms(:,count1,1) =  FCR_data_ms(:,i);
%                     count1=count1+1;
%                 end
%             end
%         end
%     else if log(i,7)==200
%             FCR_cond_ms(:,count15,15) =  FCR_data_ms(:,i);
%             count15=count15+1;
%         else if log(i,9)==-0.5
%                 FCR_cond_ms(:,count4,4) =  FCR_data_ms(:,i);
%                 count4=count4+1;
%             else if log(i,9)==-1
%                     FCR_cond_ms(:,count5,5) =  FCR_data_ms(:,i);
%                     count5=count5+1;
%                 else if log(i,9)==-2
%                         FCR_cond_ms(:,count6,6) =  FCR_data_ms(:,i);
%                         count6=count6+1;
%                     else if log(i,9)==-5
%                             FCR_cond_ms(:,count7,7) =  FCR_data_ms(:,i);
%                             count7=count7+1;
%                         else if log(i,9)==-7.5
%                                 FCR_cond_ms(:,count8,8) =  FCR_data_ms(:,i);
%                                 count8=count8+1;
%                             else if log(i,9)==-10
%                                     FCR_cond_ms(:,count9,9) =  FCR_data_ms(:,i);
%                                     count9=count9+1;
%                                 else if log(i,9)==-25
%                                         FCR_cond_ms(:,count10,10) =  FCR_data_ms(:,i);
%                                         count10=count10+1;
%                                     else if log(i,9)==-50
%                                             FCR_cond_ms(:,count11,11) =  FCR_data_ms(:,i);
%                                             count11=count11+1;
%                                         else if log(i,9)==-75
%                                                 FCR_cond_ms(:,count12,12) =  FCR_data_ms(:,i);
%                                                 count12=count12+1;
%                                             else if log(i,9)==-100
%                                                     FCR_cond_ms(:,count13,13) =  FCR_data_ms(:,i);
%                                                     count13=count13+1;
%                                                 else if log(i,9)==-200
%                                                         FCR_cond_ms(:,count14,14) =  FCR_data_ms(:,i);
%                                                         count14=count14+1;
%                                                         
%                                                     end
%                                                 end
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%             
%         end
%     end
% end
%     clear count1 count2 count3 count4 count5 count6 count7 count8 count9 count10 count11 count12 count13 count14 count15
% 
% FCR_cond_ms = FCR_cond_ms.*-1; % To get the data with minus up and positive down
% for j = 1:15
%    h=figure;
%    fig_str = ['Condition number ', int2str(j), ' in ms'];
%    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%    for i=1:10       
%        subplot (5,2,i); plot(FCR_cond_ms(:,i,j));
%        ylim([-13 13]); xlim([0 time_ms]);
%    end
%    title(fig_str);
%    print(h,fig_str,'-dpng');
%    close
% end
% 
% %% In original timepoints
% 
% FCR_cond = zeros(points,10,15);
% count1=1;count2=1;count3=1;count4=1;count5=1;count6=1;count7=1;count8=1;count9=1;count10=1;count11=1;count12=1;count13=1;count14=1;count15=1;
% for i = 1:num
%     if log(i,7)==10
%         if log(i,9)==0
%             FCR_cond(:,count3,3) =  FCR_data(:,i);
%             count3=count3+1;
%         else if log(i,9)==-0.5
%                 FCR_cond(:,count2,2) =  FCR_data(:,i);
%                 count2=count2+1;
%             else if log(i,9)==-1
%                     FCR_cond(:,count1,1) =  FCR_data(:,i);
%                     count1=count1+1;
%                 end
%             end
%         end
%     else if log(i,7)==200
%             FCR_cond_ms(:,count15,15) =  FCR_data_ms(:,i);
%             count15=count15+1;
%         else
%             if log(i,9)==-0.5
%                 FCR_cond(:,count4,4) =  FCR_data(:,i);
%                 count4=count4+1;
%             else if log(i,9)==-1
%                     FCR_cond(:,count5,5) =  FCR_data(:,i);
%                     count5=count5+1;
%                 else if log(i,9)==-2
%                         FCR_cond(:,count6,6) =  FCR_data(:,i);
%                         count6=count6+1;
%                     else if log(i,9)==-5
%                             FCR_cond(:,count7,7) =  FCR_data(:,i);
%                             count7=count7+1;
%                         else if log(i,9)==-7.5
%                                 FCR_cond(:,count8,8) =  FCR_data(:,i);
%                                 count8=count8+1;
%                             else if log(i,9)==-10
%                                     FCR_cond(:,count9,9) =  FCR_data(:,i);
%                                     count9=count9+1;
%                                 else if log(i,9)==-25
%                                         FCR_cond(:,count10,10) =  FCR_data(:,i);
%                                         count10=count10+1;
%                                     else if log(i,9)==-50
%                                             FCR_cond(:,count11,11) =  FCR_data(:,i);
%                                             count11=count11+1;
%                                         else if log(i,9)==-75
%                                                 FCR_cond(:,count12,12) =  FCR_data(:,i);
%                                                 count12=count12+1;
%                                             else if log(i,9)==-100
%                                                     FCR_cond(:,count13,13) =  FCR_data(:,i);
%                                                     count13=count13+1;
%                                                 else if log(i,9)==-200
%                                                         FCR_cond(:,count14,14) =  FCR_data(:,i);
%                                                         count14=count14+1;
%                                                     end
%                                                 end
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%             
%         end
%     end
% end
% clear count1 count2 count3 count4 count5 count6 count7 count8 count9 count10 count11 count12 count13 count14 count15
% 
% 
% FCR_cond = FCR_cond.*-1; % To get the data with minus up and positive down
% for j = 1:14
%    h=figure;
%    fig_str = ['Condition number ', int2str(j),' in tp'];
%    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%    for i=1:10       
%        subplot (5,2,i); plot(FCR_cond(:,i,j));
%        ylim([-13 13]); xlim([0 points]);
%    end
%    title(fig_str);
%    print(h,fig_str,'-dpng');
%    close
% end
% 
% ntrial=size(FCR_cond,2);
% ncond=size(FCR_cond,3);
% max_min = zeros(ntrial,ncond,2) % max = max_min(:,:,1); min = max_min(:,:,2)
% all_ind_max = zeros(ntrial,ncond); all_ind_min = zeros(ntrial,ncond);
% for i = 1:ncond
%     for j = 1:ntrial
%         if i < 4
%             [max_min(j,i,1), ind_max]=max(FCR_cond(370:550,j,i));
%             [max_min(j,i,2), ind_min]=min(FCR_cond(370:550,j,i));
%             all_ind_max(j,i) = ind_max+370; clear ind_max
%             all_ind_min(j,i) = ind_min+370; clear ind_min
%         else
%             [max_min(j,i,1), ind_max]=max(FCR_cond(4300:4450,j,i));
%             [max_min(j,i,2), ind_min]=min(FCR_cond(4300:4450,j,i));
%             all_ind_max(j,i) = ind_max+4300; clear ind_max
%             all_ind_min(j,i) = ind_min+4300; clear ind_min
%         end
%     end
% end
% 
% hr_amp = max_min(:,:,1)-max_min(:,:,2);
% mean_amp = mean(hr_amp,1);

end