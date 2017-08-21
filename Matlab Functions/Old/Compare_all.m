function Compare_all(subject)

% sub= input('Subject ? : ');
ID = {'1_AuS','2_NnNa', '3_LeC', '4_RAA'};
cd(['\\HOME\austepha\Mes Documents\Austepha_PhD\1_MainProject\Data\PAS Data\Healthy\Raw\', ID{subject}]);
ALL = zeros(4,5,2,3);
%ALL = zeros(4,5,3);

%(COND, PARAM, BASELINE)
for i = 1:3
    if i == 1
        PAS_process('b',subject)
        load b_mat.mat
    end
    if i == 2
        PAS_process('t1',subject)
        load t1_mat.mat
    end
    if i == 3
        PAS_process('t2',subject)
        load t2_mat.mat
    end
    ALL(:,:,:,i) = Mean_LICI_data;
    %     ALL(:,:,i) = Mean_LICI_data;
end

subt_str = {'Baseline', 'Post 0-5 min', 'Post 20 min'};
t_str = {'Area', 'RMS', 'Max', 'Min', 'Max-Min'};
col = {'RoyalBlue', 'Salmon', 'FireBrick', 'DarkRed'};

for j = 1 : size(ALL,2) % Param
    figure; % One figure per parameter
    count=1;
    for k = 1:size(ALL,3) % Muscle
        for i = 1:size(ALL,4) % Time
            subplot(2,3,count);
            h = bar(ALL(:,j,k,i));
            ax = gca;
            ax.XTickLabel = {'Single', 'LICI50', 'LICI100', 'LICI200'};
            title(subt_str{1,i})
            ylim([0 1])
            count=count+1;
            hold on
        end
    end
    %title(t_str{1,j})
end

%% Only MEP Peak to Peak

ALL_Amplitude = squeeze(ALL(:,5,:,:));
% ALL_Amplitude = squeeze(ALL(:,5,:));
Amplitude_ADM = squeeze(ALL_Amplitude(:,1,:));
Amplitude_FDI = squeeze(ALL_Amplitude(:,2,:));

%% _____________________ RATIOS _______________________

% MEP Peak to Peak - Ratio over single
ratios = zeros(3,3,2);
for i = 1:3
    curr_cond = squeeze(ALL(i+1,5,:,:))'
    baseline_cond = squeeze(ALL(1,5,:,:))'
    ratios(i,:,:)= curr_cond./baseline_cond*100
end
% OR
RatioOverSingle_ADM = zeros(3,3);
for i = 1:3
    curr_cond = ALL_Amplitude(i+1,:)
    baseline_cond = ALL_Amplitude(1,:)
    RatioOverSingle_ADM(i,:)= curr_cond./baseline_cond*100
end
% ratios = squeeze(ratios);

RatioOverSingle_ADM = squeeze(ratios(:,:,1));
RatioOverSingle_FDI = squeeze(ratios(:,:,2));


figure;
for i=1:size(ratios,3)
    subplot(2,1,i)
    bar(ratios(:,:,i))
    hold on
    plot(xlim,[100 100])
    ax = gca;
    ax.XTickLabel = {'LICI50', 'LICI100', 'LICI200'};
    if i==1
        title('Ratios (over single)- ADM')
    else
        title('Ratios (over single) - FDI')
    end
    ylim([0 350])
end

% MEP Peak to Peak - Ratio over baseline
RatioOverBaseline = zeros(4,2,2);
for i = 1:2
    curr_cond = squeeze(ALL(:,5,:,i+1))
    baseline_cond = squeeze(ALL(:,5,:,1))
    RatioOverBaseline(:,:,i)= curr_cond./baseline_cond*100
end
% OR
RatioOverBaseline_ADM = zeros(4,2);
for i = 1:2
    %curr_cond = Amplitude_ADM(:,i+1)
    curr_cond = ALL_Amplitude(:,i+1)
    %baseline_cond = Amplitude_ADM(:,1)
    baseline_cond = ALL_Amplitude(:,1)
    RatioOverBaseline_ADM(:,i)= curr_cond./baseline_cond*100
end
% ratios = squeeze(ratios);

RatioOverBaseline_ADM = squeeze(RatioOverBaseline(:,:,1))';
RatioOverBaseline_FDI = squeeze(RatioOverBaseline(:,:,2))';

figure;
for i=1:size(ratios,3)
    subplot(2,1,i)
    bar(ratios(:,:,i))
    hold on
    plot(xlim,[100 100])
    ax = gca;
    ax.XTickLabel = {'LICI50', 'LICI100', 'LICI200'};
    if i==1
        title('Ratios (over single)- ADM')
    else
        title('Ratios (over single) - FDI')
    end
    ylim([0 350])
end
%%

figure;
count = 1;
for j = 1:size(ALL,3)
    for i = 1:size(ALL,4)
        subplot(2,3,count);
        h = bar(ALL(:,5,j,i));
        ax = gca;
        ax.XTickLabel = {'Single', 'LICI50', 'LICI100', 'LICI200'};
        title(subt_str{1,i})
        %ylim([-2 8])
        count = count +1;
        hold on
    end
end

ratios = zeros(2,4,2);
for i = 1:2
    curr_cond = squeeze(ALL(:,5,:,i+1))'
    baseline_cond = squeeze(ALL(:,5,:,1))'
    ratios(i,:,:)= curr_cond./baseline_cond*100
end
ratios = squeeze(ratios);

figure;
for i=1:size(ratios,3)
    subplot(2,1,i)
    bar(ratios(:,:,i))
    hold on
    plot(xlim,[100 100])
    ax = gca;
    ax.XTickLabel = {'LICI50', 'LICI100', 'LICI200'};
    if i==1
        title('Ratios (over single)- ADM')
    else
        title('Ratios (over single) - FDI')
    end
    ylim([0 350])
end
end


%-------------------------------------------
% -------------- CEMETERY ------------------
%-------------------------------------------


%ALL = zeros(4,5,2); For Nat because we don't have t2
% for i = 1:2
%     if i == 1
%         load b_mat.mat 
%     end
%     if i == 2
%         load t1_mat.mat 
%     end
%     ALL(:,:,i) = Mean_LICI_data;
% end

% FOR NAT
% for j = 1 : size(ALL,2) % Param
%     figure; % One figure per parameter
%     count=1;
%         for i = 1:size(ALL,4) % Time
%             subplot(2,3,count);
%             h = bar(ALL(:,j,i));
%             ax = gca;
%             ax.XTickLabel = {'Single', 'LICI50', 'LICI100', 'LICI200'};
%             title(subt_str{1,i})
%             ylim([0 1])
%             count=count+1;
%             hold on
%         end
%     %title(t_str{1,j})
% end

% Area under the curve
% ratios = zeros(3,3);
% for i = 1:3
%     curr_cond = squeeze(ALL(i+1,1,:))'
%     baseline_cond = squeeze(ALL(1,1,:))'
%     ratios(i,:)= curr_cond./baseline_cond*100
% end
% ratios = squeeze(ratios);
