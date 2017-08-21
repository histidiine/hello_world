function MEP_boxplot

%% _______ BOX PLOT _________

load ALL_H_percent;

figure; count=1;
for sub = 1:size(ALL_H_percent,3)
    subplot(size(ALL_H_percent,3),2,2*sub)
    curr_data=squeeze(ALL_H_percent(:,1,sub)');
    boxplot(curr_data);
    
    subplot(size(ALL_H_percent,3),2,2*sub-1)
    curr_data=squeeze(ALL_H_percent(:,2,sub)');
    boxplot(curr_data);
    
end

%         ylim([0 2])
    %         set(gcf,'units','normalized','outerposition',[0 0 1 1]);

end