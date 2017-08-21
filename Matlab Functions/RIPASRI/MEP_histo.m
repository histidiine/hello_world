function MEP_histo

% Leave only histograms --> create another for plot graphs
% Add a choice at the beginning for individual subject or group histograms
%

%% ______ INITIALIZATION__________
dataset = ['b', 't1', 't2'];
ds_folder = {'Baseline', 'T1', 'T2'};
m_name={'ADM','FDI','FCR','FPB'};
step_name={'baseline','T1','T2'};
cond_name={'single', 'LICI 50', 'LICI 100', 'LICI 200'};


%% _______ HISTOGRAM ________

histProc='on';

while strcmp(histProc,'on')
    
    %__________________CHOICE SUBJECT OR GROUP__________________
    disp('1- One single subject');
    disp('2- All subjects of a group');
    disp('3- Mean of a group');
    switch_pop=0;
    while switch_pop~=1 && switch_pop~=2 && switch_pop~=3
        switch_pop = input('What population do you want to process (see choice above)?:   ');
    end
    
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
    
    %______________________LAUNCHING CHOSEN PROCESS______________________
    switch switch_pop
        
        case 1 %______________________ONE SINGLE SUBJECT______________________
            
            gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
            cd(gp_folder)
            
            for k = 1:length(list_folders)
                disp(list_folders(k).name)
            end
            sub_ID='none';
            while ~exist(sub_ID,'dir')
                sub_ID = input('Please enter the ID of the subject you want to process into single quotes (see above):   ');
            end
            
            cd(sub_ID); cd('Data');
            load Resp_Amp
            
            hist=figure;count=1;
            for MyUnit = 1:2 % 1 - in mV 2- ratio
                for m = 1:3
                    
                    %_________________CALCULATING__________________
                    curr_data = squeeze(resp(:,m,:,:,5)); % trials, cond, step
                    curr_mean = squeeze(mean(curr_data,1)); % cond, step
                    curr_dev = squeeze(std(curr_data,1)); % cond, step
                    curr_ratio=zeros(size(curr_mean));
                    for step = 1:size(curr_mean,2)
                        for cond = 1:size(curr_mean,1)
                            curr_ratio(cond,step,:) = curr_mean(cond,step,:)./curr_mean(1,1,:);
                        end
                    end
                    
                    %________________CREATING FIGURE_________________
                    max_V=max(max(max(curr_dev)))+ max(max(max(curr_mean)));
                    max_R=max(max(max(curr_ratio)));
                    subplot(2,3,count)
                    
                    x=[0.78,1.78,2.78,3.78 ;...
                        1,2,3,4; 1.22,2.22,3.22,4.22];
                    
                    if MyUnit == 1 % PLOT IN MV WITH ERROR
                        bar(curr_mean, 1,'edgecolor','k', 'linewidth', 2);
                        hold on;
                        for col = 1:size(curr_mean,2)
                            %std_err = dev(:,col)/sqrt(size(curr_data,1));
                            errorbar(x(col,:),curr_mean(:,col),curr_dev(:,col),'.','LineWidth',2,'Color','k');
                        end
                        %  ylim([0 (max_V+0.1)]);
                    else if MyUnit == 2 % PLOT IN RATIO
                            bar(curr_ratio, 1,'edgecolor','k', 'linewidth', 2);
                            ylim([0 (max_R+0.1)]);
                        end
                    end
                    
                    %_____________GRAPH PARAMETERS______________
                    colormap(bone);
                    hold on
                    set(gca,'XTickLabel',cond_name);
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                    set(gcf,'NextPlot','add');
                    axes;
                    %xlabel('TMS Paradigm (Single/Paired Pulse)', 'fontsize',14);
                    ylabel('MEP Amplitude (in mV)', 'fontsize',8);
                    %legend(step_name, 'location', 'best', 'fontsize',11);
                    title(['Subject ', int2str(count)], 'fontsize',8);
                    count = count +1;
                end
            end
            ht=title(['MEP amplitude for ', m_name{m}]);
            set(gca,'Visible','off');
            set(ht,'Visible','on');
            
            cd(gp_folder);cd(sub_ID);cd('Figures')
            save_bmp='HIST_MuscleXTrials.bmp';
            export_fig(save_bmp, '-bmp', '-nocrop', '-r160');
            
            
            
            %________________________________________________________________________
            %________________________________________________________________________
            
            
            
        case 2 %______________________ALL SUBJECTS OF A GROUP______________________
            
            disp(['Starting - Process all subjects of ', gp_ID]);
            cd(gp_folder)
            fileName='Resp_Amp.mat';
            load(fileName)
            %__________________PARAMETERS FOR PLOT STRUCTURE_______________
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
            
            for m = 1:3 % two fig per muscle - one in mV one in ratio
                for MyUnit = 1:2
                    hist=figure;count = 1;
                    for sub = 1:n_subjects
                        
                        %_________________CALCULATING__________________
                        curr_data= MEP_amp_GP(:,m,:,:,5,sub);
                        curr_mean = squeeze(mean(curr_data,1));
                        curr_dev = squeeze(std(curr_data,1));
                        curr_ratio=zeros(size(curr_mean));
                        for step = 1:size(curr_mean,3)
                            for cond = 1:size(curr_mean,2)
                                curr_ratio(:,cond,step) = curr_mean(:,cond,step)./curr_mean(:,1,1);
                            end
                        end
                        
                        %________________CREATING FIGURE_________________
                        subplot(row,col,count)
                        x=[0.78,1.78,2.78,3.78 ;...
                            1,2,3,4; 1.22,2.22,3.22,4.22];
                        
                        if MyUnit == 1 % PLOT IN MV WITH ERROR
                            bar(curr_mean, 1,'edgecolor','k', 'linewidth', 2);
                            hold on;
                            for col = 1:size(curr_mean,2)
                                errorbar(x(col,:),curr_mean(:,col),curr_dev(:,col),'.','LineWidth',2,'Color','k');
                            end
                            unitStr='mV';
                        else if MyUnit == 2 % PLOT IN RATIO
                                bar(curr_ratio, 1,'edgecolor','k', 'linewidth', 2);
                                unitStr='ratio';
                            end
                        end
                        
                        %_____________GRAPH PARAMETERS______________
                        colormap(bone);
                        hold on
                        set(gca,'XTickLabel',cond_name);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        set(gcf,'NextPlot','add');
                        axes;
                        %xlabel('TMS Paradigm (Single/Paired Pulse)', 'fontsize',14);
                        ylabel('MEP Amplitude (in mV)', 'fontsize',14);
                        %legend(step_name, 'location', 'best', 'fontsize',12);
                        title(m_name{m}, 'fontsize',14);
                        
                        count = count+1;
                    end
                    cd(gp_folder);
                    cd('*Fig_Save');
                    titleStr = ['Hist ',gp_ID,' ' ,m_name{m},' ',unitStr];
                    ht=title(titleStr);
                    set(gca,'Visible','off');
                    set(ht,'Visible','on');
                    save_bmp=['Hist_',gp_ID,'_' ,m_name{m},'_',unitStr,'.bmp'];
                    export_fig('save_bmp', '-bmp', '-nocrop', '-r160');
%                     print(hist,save_png,'-dpng','-r900');
                end
            end
            
            
            %________________________________________________________________________
            %________________________________________________________________________
            
            
            
        case 3 %______________________GROUP MEAN_________________________
            
            gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
            cd(gp_folder)
            
            fileName='Resp_Amp.mat';
            load(fileName)
            singleFile = 'all_single_GP.mat';
            load(singleFile);
            
            hist=figure;count=1;
            for MyUnit = 1:2 % 1 - in mV 2- ratio
                for m = 1:3
                    
                    %_________________CALCULATING__________________
                    curr_data = squeeze(MEP_amp_GP(:,m,:,:,5,:)); % trials, cond, step, subject
                    curr_single = squeeze(all_single_GP(:,m,:,5,:)); % trials, step, subject
                    
                    curr_mean=zeros(4,3);
                    curr_mean(1,:)=squeeze(mean(mean(curr_single,1),3));
                    curr_mean(2:4,:) = squeeze(mean(mean(curr_data(:,2:4,:,:),1),4)); % cond, step
                    
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
                    max_V=max(max(max(curr_dev)))+ max(max(max(curr_mean)));
                    max_R=max(max(max(curr_ratio)));
                    subplot(2,3,count)
                    
                    x=[0.78,1.78,2.78,3.78 ;...
                        1,2,3,4; 1.22,2.22,3.22,4.22];
                    
                    if MyUnit == 1 % PLOT IN MV WITH ERROR
                        bar(curr_mean, 1,'edgecolor','k', 'linewidth', 2);
                        hold on;
                        for col = 1:size(curr_mean,2)
                            std_err = curr_dev(:,col)/sqrt(size(curr_data_conc,1));
                            errorbar(x(col,:),curr_mean(:,col),std_err,'.','LineWidth',2,'Color','k');
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
                    set(gca,'XTickLabel',cond_name);
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                    set(gcf,'NextPlot','add');
                    axes;
                    count = count +1;
                end
            end
            cd(gp_folder); cd('*Fig_Save');
            titleStr = ['MeanHist ',gp_ID];
            ht=title(titleStr);
            set(gca,'Visible','off');
            set(ht,'Visible','on');
            save_bmp=['MeanHist_',gp_ID,'.bmp'];
            export_fig(save_bmp, '-bmp', '-nocrop', '-r160')
%             print(hist,save_png,'-dbmp16m');
            
            
            %________________________________________________________________________
            %________________________________________________________________________
    end
    
    
    %_________________CHECK IF USER WANTS TO DO ANOTHER HIST_______________
    continueRep=0;
    while strcmp(continueRep,'y') == 0 && strcmp(continueRep,'n') ==0
        continueRep=input('Do you want to do another histogram? (y/n):   ');
        if strcmp(continueRep,'y') == 1
            histProc = 'on';
        else
            histProc = 'off';
            disp('End of FIGURE CREATION - Thank you for using PhDLife Airline !');
        end
    end
end


end


