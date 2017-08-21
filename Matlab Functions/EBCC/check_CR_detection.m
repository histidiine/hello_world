function check_CR_detection(raw_data,extra_data,time_axis)

resp=input('Do you want to see all trials ? (char: y/n)  ');
if size(raw_data,1)==8750
        noise1=1237; 
        noise2=2487;
        cr1=3750;
        cr2=5000;
        dur=306;
        dur_pre=94;
    else
        noise1=198;
        noise2=398;
        cr1=600;
        cr2=800;
        dur=49;
        dur_pre=15;
end
    
if strcmp(resp, 'y')
    
    
%% Visualizer all trials per block 
for j = 1:7 % LOOP FOR BLOCKS
   h=figure; % One figure per block
   set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
   tlstr=['Block number ', int2str(j)];
   for i=1:10 % LOOP FOR TRIALS
       
       indice_take = j*10+i; % in order not to take first 10 trials
       indice_save = (j-1)*10+i; % in order to save starting at 1
       ylim_up = max(raw_data(:,i))+0.2;
       
       if size(raw_data,1)==8750
           dur_value=extra_data(indice_save,6)./6.25;
       else
           dur_value=extra_data(indice_save,6);
       end
       str = ['Block number ', int2str(j), ' - Trial number ', int2str(i)];
       %t=0:(1/35):250;t=t(:,1:8750);t=t';
       %subplot (5,2,i); plot(t,raw_data(:,indice_take)); title(str);
       subplot (5,2,i); plot(time_axis,raw_data(:,indice_take)); title(str);
       hold on;
       line([cr1 cr1],[-0.15 0.15],'Color',[1 0 0]);
       hold on;
       line([cr2 cr2],[-0.15 0.15],'Color',[1 0 0]);
       hold on; 
       if extra_data(indice_save,3)==0 % only for trials where a CR was detected
       else
           dur_txt=num2str(dur_value);
           line([extra_data(indice_save,4) extra_data(indice_save,5)],[-0.09 -0.09],'Color',[0 1 0]); % put a red line under the CR
           descr = [dur_txt, ' ms'];
           text(extra_data(indice_save,4),-0.12,descr)
       end
       hold on;
%        if noise_resp(indice_save,4)==0 % only for trials where a alpha blink detected
%        else
%            line([noise_resp(indice_save,4) noise_resp(indice_save,4)+50],[-0.1 -0.1],'Color',[0 1 0]); % put a green line under the alpha blink
%        end
       ylim([-ylim_up ylim_up]);
   end
   kk=1;fig_str = tlstr; 
   while 7 == exist(tlstr, 'file'); fig_str = [tlstr,int2str(kk)]; end
   print(h,fig_str,'-dpng');
   pause
   close
end
else
    disp('Moving to next step without visualizing all trials')
end
resp=input('Do you want to visualize the trials with detected CR only ? (char: y/n)  ');

if strcmp(resp, 'y')
    %% Visualize only the trials with detected CR
    for j = 1:7 % LOOP FOR BLOCKS
        for i=1:10 % LOOP FOR TRIALS
            indice_take = j*10+i; % in order not to take first 10 trials
            indice_save = (j-1)*10+i; % in order to save starting at 1
            if extra_data(indice_save,3)==0 % only for trials where a CR was detected
            else
                if size(raw_data,1)==8750
                    dur_value=extra_data(indice_save,6)./6.25;
                else
                    dur_value=extra_data(indice_save,6);
                end
                h=figure; % One figure per trial where a CR was detected
                set(gcf,'units','normalized','outerposition',[0 0 1 1]); % Full screen
                tlstr=['Block number ', int2str(j),' - Trial number ', int2str(i)];
                % ------ NOISE -------
                subplot (1,2,1); plot(time_axis,raw_data(noise1-375:noise2+375,indice_take)); 
                ylim([-0.15 0.15]); xlim([0 2000]) 
                descr=[num2str(extra_data(indice_save,1)*1000), ' uV']; title('noise');
                text(900, -0.11, descr);
                hold on;
                line([375 375],[-0.05 0.05],'Color',[1 0 0]); % Vertical Line at start
                hold on;
                line([1625 1625],[-0.05 0.05],'Color',[1 0 0]); % Vertical Line at end

                % ------ ZOOM CR --------
                subplot (1,2,2); plot(time_axis,raw_data(extra_data(indice_save,4)-500:extra_data(indice_save,4)+1500,indice_take)); title(tlstr);
                ylim([-0.15 0.15]); xlim([0 2000]); 
                line_length=(500+extra_data(indice_save,6));
                    % ----- CR lines ----
                    hold on;
                    line([500 500],[-0.05 0.05],'Color',[1 0 0]); % Vertical Line at start
                    hold on;
                    line([line_length line_length],[-0.05 0.05],'Color',[1 0 0]); % Vertical Line at end

                    hold on;
                    line([500 line_length],[-0.07 -0.07],'Color',[0 1 0]); % put a red line under the CR
                    descr = [num2str(dur_value), ' ms'];
                    text(500,-0.09,descr);
                    descr2=[num2str(extra_data(indice_save,2)*1000), ' uV'];
                    text(500,-0.11,descr2);
                kk=1;fig_str = tlstr;
                while 7 == exist(tlstr, 'file'); fig_str = [tlstr,int2str(kk)]; end
                print(h,fig_str,'-dpng');
                pause
                close
            end
        end
    end
else
    disp('End of check Cr Detection protocol')
end
end