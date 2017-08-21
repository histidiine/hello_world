function print_RI(data_cond,folder,fig_name)

data_cond = data_cond.*-1; % To get the data with minus up and positive down
cd(folder)
for j = 1:size(data_cond,3) % COND LOOP
    h=figure;
    if size(data_cond,1)==5000; data_scale = '_in_tp'; else data_scale = '_in_ms'; end
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for i=1:size(data_cond,2) % TRIAL LOOP
        fig_str = [fig_name, 'Condition number ', int2str(j), ' trial ', int2str(i)];
        if j < 4
            plot_start = 1;
            plot_end = 650;
        else
            [val_max ind_max]=max(data_cond(:,i,j));
            plot_start = ind_max-300;
            plot_end = size(data_cond,1);
        end
        data_zoom = data_cond(plot_start:plot_end,i,j);
        subplot (2,5,i); plot(data_zoom);
       title([fig_str, data_scale]);
       ylim([-15 15]); xlim([0 size(data_zoom,1)]);
       clear plot_start plot_end ind_max data_zoom
   end
   %print(h,fig_str,'-dpng');
   pause;
   close
end