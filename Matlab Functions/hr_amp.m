function [output1, output2, output3, output4, output5] = hr_amp(data_matrix)

ntrial=size(data_matrix,2);
ncond=size(data_matrix,3);
max_min = zeros(ntrial,ncond,2); % max = max_min(:,:,1); min = max_min(:,:,2)
all_ind_max = zeros(ntrial,ncond); all_ind_min = zeros(ntrial,ncond);
for i = 1:ncond
    for j = 1:ntrial
        if i == 15
            [max_min(j,i,1), ind_max]=max(data_matrix(4210:4300,j,i));
            [max_min(j,i,2), ind_min]=min(data_matrix(ind_max+4210:4310,j,i));
            all_ind_max(j,i) = ind_max+4210;  
            all_ind_min(j,i) = ind_min+ind_max+4210; 
            clear ind_min ind_max
        else if i < 4
                [max_min(j,i,1), ind_max]=max(data_matrix(380:560,j,i));
                [max_min(j,i,2), ind_min]=min(data_matrix(ind_max+380:570,j,i));
                all_ind_max(j,i) = ind_max+380;  
                all_ind_min(j,i) = ind_min+ind_max+380; 
                clear ind_min ind_max
            else
                [max_min(j,i,1), ind_max]=max(data_matrix(4310:4415,j,i));
                [max_min(j,i,2), ind_min]=min(data_matrix(ind_max+4310:4425,j,i));
                all_ind_max(j,i) = ind_max+4310; 
                all_ind_min(j,i) = ind_min+ind_max+4310; 
                clear ind_min ind_max
            end
        end
        disp(['Finished trial number ', int2str(j), ' of condition ', int2str(i)]);
    end
end

hr_amp = max_min(:,:,1)-max_min(:,:,2);
mean_amp = mean(hr_amp,1);
output1 = hr_amp; output2=(mean_amp); output3=all_ind_max; output4=all_ind_min; output5=max_min;
resp=input('Do you want to see the figures ? (y/n) ');

if resp == 'n'
else
    count=0;
    for i = 1:ncond
        for j = 1:ntrial
            figure;
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            plot(data_matrix(:,j,i)); title(['Condition ', int2str(i), 'Trial number ', int2str(j)])
            hold on; line([all_ind_max(j,i) all_ind_max(j,i)],[-15 15],'Color',[1 0 0]);
            hold on; line([all_ind_min(j,i) all_ind_min(j,i)],[-15 15],'Color',[1 0 0]);
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
end
disp('End of hr_amp')
end