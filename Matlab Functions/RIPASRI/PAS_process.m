function PAS_process(cond_code,sub_num)

%_Folders
ID = {'1_AuS','2_NnNa', '3_LeC', '4_RAA'};
dataset_folder = ['\\HOME\austepha\Mes Documents\Austepha_PhD\1_MainProject\Data\PAS Data\Healthy\Raw\', ID{sub_num}];
extract_folder = ['\\HOME\austepha\Mes Documents\Austepha_PhD\1_MainProject\Data\PAS Data\Healthy\Extract\Healthy', int2str(sub_num)];

cd(dataset_folder) 

%% Determining what recordings we are processing

if strcmp(cond_code, 'b') == 1
    Extract_lici = 'b_lici.txt';
    Extract_single = 'b_single.txt';
    Save_Graph = 'b_graph.jpeg';
    Save_Matrices = 'b_mat.mat';
    open 'baseline.xlsx'
    pause
end

if strcmp(cond_code, 't1') == 1
    Extract_lici = 't1_lici.txt';
    Extract_single = 't1_single.txt';
    Save_Graph = 't1_graph.jpeg';
    Save_Matrices = 't1_mat.mat';
    open 'baseline2.xlsx'
    pause
    baseline = baseline2;
end

if strcmp(cond_code, 't2') == 1
%     Extract_lici = 't2_lici.txt';
    Extract_lici = 't2_lici2.txt';
%     Extract_single = 't2_single.txt';
    Extract_single = 't2_single2.txt';
    Save_Graph = 't2_graph.jpeg';
    Save_Matrices = 't2_mat.mat';
    open 'baseline3.xlsx'
    pause
    baseline = baseline3;
end

conditions = baseline(:,3);

%% Loading data
cd(extract_folder)
sici_all=load(Extract_lici);
single=load(Extract_single);

% Initializing
j=1; k=1; l=1; m=1; % row_size=size(sici_all,1)/4;

% If we have the 2 muscles
LICI_data = zeros(12,4,5,2);% Data,Condition, Parameters(1-Area, 2-RMS, 3-Max, 4-Min, 5-Max-Min), Muscle

%% ___________ Filling LICI_data 2 MUSCLES____________________
for i = 1:size(conditions,1)
    curr_cond = conditions(i,1); 
    if curr_cond == 0
        % Abductor Digiti Mini
        LICI_data(m,1,1,1)= single(i,5);
        LICI_data(m,1,2,1)= single(i,6);
        LICI_data(m,1,3,1)= single(i,7);
        LICI_data(m,1,4,1)= single(i,8);
        LICI_data(m,1,5,1)= single(i,9);
        
        % First Dorsal Interosseous
        LICI_data(m,1,1,2)= single(i,29);
        LICI_data(m,1,2,2)= single(i,30);
        LICI_data(m,1,3,2)= single(i,31);
        LICI_data(m,1,4,2)= single(i,32);
        LICI_data(m,1,5,2)= single(i,33);
        m=m+1;
    end
    if curr_cond == 60 %|| 50 
        % Abductor Digiti Mini
        LICI_data(j,2,1,1)= sici_all(i,5);
        LICI_data(j,2,2,1)= sici_all(i,6);
        LICI_data(j,2,3,1)= sici_all(i,7);
        LICI_data(j,2,4,1)= sici_all(i,8);
        LICI_data(j,2,5,1)= sici_all(i,9);
        
        % First Dorsal Interosseous
        LICI_data(j,2,1,2)= sici_all(i,29);
        LICI_data(j,2,2,2)= sici_all(i,30);
        LICI_data(j,2,3,2)= sici_all(i,31);
        LICI_data(j,2,4,2)= sici_all(i,32);
        LICI_data(j,2,5,2)= sici_all(i,33);
        j=j+1;
    end
    if curr_cond == 110 %|| 100
        % Abductor Digiti Mini
        LICI_data(k,3,1,1)= sici_all(i,13);
        LICI_data(k,3,2,1)= sici_all(i,14);
        LICI_data(k,3,3,1)= sici_all(i,15);
        LICI_data(k,3,4,1)= sici_all(i,16);
        LICI_data(k,3,5,1)= sici_all(i,17);
       
        % First Dorsal Interosseous
        LICI_data(k,3,1,2)= sici_all(i,37);
        LICI_data(k,3,2,2)= sici_all(i,38);
        LICI_data(k,3,3,2)= sici_all(i,39);
        LICI_data(k,3,4,2)= sici_all(i,40);
        LICI_data(k,3,5,2)= sici_all(i,41);
        k=k+1;
    end
    if curr_cond == 210 % || 200
        % Abductor Digiti Mini
        LICI_data(l,4,1,1)= sici_all(i,21);
        LICI_data(l,4,2,1)= sici_all(i,22);
        LICI_data(l,4,3,1)= sici_all(i,23);
        LICI_data(l,4,4,1)= sici_all(i,24);
        LICI_data(l,4,5,1)= sici_all(i,25);
        
        % First Dorsal Interosseous
        LICI_data(l,4,1,2)= sici_all(i,45);
        LICI_data(l,4,2,2)= sici_all(i,46);
        LICI_data(l,4,3,2)= sici_all(i,47);
        LICI_data(l,4,4,2)= sici_all(i,48);
        LICI_data(l,4,5,2)= sici_all(i,49);
        l=l+1;
    end
end
All_Data_Amp = squeeze(LICI_data(:,:,5,:))

%% ___________ PARAMETERS MEAN ______________

Mean_LICI_data = squeeze(mean(LICI_data,1));
% Rows = Single, LICI50, LICI100, LICI200
% Column = Area, RMS, Max, Min, Max-Min

%% Saving
figure_fold = [dataset_folder,'\Figures'];
if ~exist(figure_fold, 'dir')
  mkdir(figure_fold);
end
cd(figure_fold)
save(Save_Matrices, 'LICI_data', 'Mean_LICI_data', 'All_Data_Amp')
% saveas(graph, Save_Graph)
end




%-------------------------------------------
% -------------- CEMETERY ------------------
%-------------------------------------------

%FID = fopen('baseline.txt')
%C = textscan(FID,'string')
%logfile=dlmread('baseline.txt'); 
%logfile =logfile(4:end,:);
% fid = fopen('baseline.txt');
% tline = fgetl(fid);
% size_file=size(tline,1);
% i=0;
% while ischar(tline)
%     disp(tline)
%     tline = fgetl(fid);
%     i=i+1
%     if i>2
%         disp(tline(:,1))
%         disp(tline(:,5))
%     end
% end
% fclose(fid);

% if cond_code == 'b'
%     Extract_lici = 'b_lici.txt';
%     Extract_single = 'b_single.txt';
%     Save_Graph = 'b_graph';
%     Save_Matrices = 'b_mat.mat';
%     open 'baseline.xlsx'
% end

% If we have the 1 muscle
% LICI_data = zeros(12,4,5);

% %% ___________ Filling LICI_data 1 MUSCLES____________________
% for i = 1:size(conditions,1)
%     curr_cond = conditions(i,1); 
%     if curr_cond == 0
%         % Abductor Digiti Mini
%         LICI_data(m,1,1)= single(i,5);
%         LICI_data(m,1,2)= single(i,6);
%         LICI_data(m,1,3)= single(i,7);
%         LICI_data(m,1,4)= single(i,8);
%         LICI_data(m,1,5)= single(i,9);
%         m=m+1;
%     end
%     if curr_cond == 50 
%         % Abductor Digiti Mini
%         LICI_data(j,2,1)= sici_all(i,5);
%         LICI_data(j,2,2)= sici_all(i,6);
%         LICI_data(j,2,3)= sici_all(i,7);
%         LICI_data(j,2,4)= sici_all(i,8);
%         LICI_data(j,2,5)= sici_all(i,9);
%         j=j+1;
%     end
%     if curr_cond == 100
%         % Abductor Digiti Mini
%         LICI_data(k,3,1)= sici_all(i,13);
%         LICI_data(k,3,2)= sici_all(i,14);
%         LICI_data(k,3,3)= sici_all(i,15);
%         LICI_data(k,3,4)= sici_all(i,16);
%         LICI_data(k,3,5)= sici_all(i,17);
%         k=k+1;
%     end
%     if curr_cond == 200
%         % Abductor Digiti Mini
%         LICI_data(l,4,1)= sici_all(i,21);
%         LICI_data(l,4,2)= sici_all(i,22);
%         LICI_data(l,4,3)= sici_all(i,23);
%         LICI_data(l,4,4)= sici_all(i,24);
%         LICI_data(l,4,5)= sici_all(i,25);
%         l=l+1;
%     end
% end


%% _______________ SPECIAL - SINGLE ALL TRIALS _________________

% Single_AllTrials(1,1)= single(i,5);
%         LICI_data(1,2)= single(i,6);
%         LICI_data(1,3)= single(i,7);
%         LICI_data(1,4)= single(i,8);
%         LICI_data(1,5)= single(i,9);

% %% ______________ GRAPHS ___________________
% 
% t_str = {'Area', 'RMS', 'Max', 'Min', 'Max-Min'};
% col = {'RoyalBlue', 'Salmon', 'FireBrick', 'DarkRed'};
% 
% graph = figure;
% count=1;
% for j = 1:size(Mean_LICI_data,3)
%     for i = 1:size(Mean_LICI_data,2)%Condition, Parameters(1-Area, 2-RMS,...
%         ...3-Max, 4-Min, 5-Max-Min), Muscle
%             subplot(2,5,count);
%         
%         numberOfBars = size(Mean_LICI_data,1);
%         for b = 1 : numberOfBars
%             h = bar(Mean_LICI_data(:,i,j));
%             %set(h, 'FaceColor', rgb(col(1,b)))
%         end
% %         ylim([])
%         ax = gca;
%         ax.XTickLabel = {'Single', 'LICI50', 'LICI100', 'LICI200'};
%         title(t_str{1,i})
%         count = count+1;
%         hold on
%     end
%     hold on 
% end

%\\HOME\austepha\Mes Documents\Austepha_PhD\1_MainProject\Data\PAS Data\Healthy\Raw\1_NnNa\Figures
%\\HOME\austepha\Mes Documents\Austepha_PhD\1_MainProject\Data\PAS Data\Healthy\Raw\2_LeC\Figures
%\\HOME\austepha\Mes Documents\Austepha_PhD\1_MainProject\Data\PAS Data\Healthy\Raw\3_AuS\Figures

