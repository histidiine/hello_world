function [output1, output2]= CRwDuration(data_matrix)


% l'indice k est la référence de détection générale du CR
% l'indice kk est la référence de détection du début du CR
% l'indice ll est la référence de détection de la fin du CR dans le timebin
% défini avec l'indice k à la première détection
% l'indice mm est la référence de détection de la fin du CR en dehors du
% timebin (dans le cas où le CR serait plus long que les premières 50 ms détectée avec l'indice k)
% ce qui sera sûrement souvent le cas!!

CR_block = zeros(7,1); % = nb of CR
CR=0;coef=1.4;
noise_resp=zeros(70,5);%1-noise, 2-resp, 3-result, 4-start, 5-end, 6-duration

if size(data_matrix,1)==8750
    test1=1237;
    test2=2487;
    cr1=3737;
    cr2=4681;
    dur=306;
    dur_pre=94;
else
    test1=198;
    test2=398;
    cr1=598;
    cr2=749;
    dur=49;
    dur_pre=15;
end

for j = 1:7 % without the first 10 trials / LOOP OF BLOCKS
    for i=  1:10 % LOOP OF TRIALS
        indice_take = j*10+i; % in order not to take first 10 trials
        indice_save = (j-1)*10+i; % in order to save starting at 1
        search='on';k=cr1;
        timepoint_start=0; timepoint_end=0; CR_duration=0;
        noise_resp(indice_save,3)=0;
        while strcmp(search, 'on') && k < (cr2) % LOOP OF TIME CHUNKS
       
            
            %___________INITIALIZATION_________________
            ref= mean(abs(data_matrix(test1:test2,indice_take)),1); % ref for noise 200 ms in the first 400
            test=mean(abs(data_matrix(k:k+dur,indice_take)),1); % Chunks of 50 ms
            noise_resp(indice_save,1) = ref; noise_resp(indice_save,2) = test;
            
            %-------------DETECTING CR---------------------
            if test > 2*ref % if 50 ms timebin activity is bigger than 2*noise = CR
                start_search='on'; kk = k;
                while strcmp(start_search,'on') && kk < (k+dur+1) % finding beginning inside the selected timebin
                    test_start=mean(abs(data_matrix(kk:kk+10,indice_take)),1); % =0.5 ms timebin
                    if test_start > 2*ref
                        timepoint_start=kk;
                        start_search='off';
                        %end
                        %kk=kk+1;
                        %end
                        end_search='on';
                        %--------- IF Start was found > Looking for end ---------
                        if strcmp(start_search, 'off')
                            ll=1;
                            while strcmp(end_search, 'on') && ll < (dur+1) % finding end inside the selected timebin
                                test_end=mean(abs(data_matrix(k+dur-(ll+10):k+dur-(ll),indice_take)),1);
                                if test_end > 2*ref
                                    mm = dur-(ll+10); % starting at the timepoint of the end inside selected timebin
                                    mm_past = dur-(ll+10);
                                    while strcmp(end_search, 'on') && k+mm < 4990
                                        % checking if the response is longer than the selected timebin
                                        % and looking for the total duration of the CR
                                        test_end=mean(abs(data_matrix(k+mm:k+mm+35,indice_take)),1);
                                        if test_end > 2*ref
                                            mm_past=mm;
                                            if k+mm == 4989
                                                timepoint_end=4999;
                                                CR_duration=timepoint_end-timepoint_start;
                                                if CR_duration > dur-dur_pre
                                                    end_search='off';
                                                    search='off';
                                                    CR=CR+1;
                                                    noise_resp(indice_save,3)=1;
                                                    noise_resp(indice_save,4)=timepoint_start;
                                                    noise_resp(indice_save,5)=timepoint_end;
                                                    noise_resp(indice_save,6)=CR_duration;
                                                    disp('Ouh! Definitely got one!')
                                                end
                                            end
                                        else
                                            timepoint_end=k+mm_past;
                                            CR_duration=timepoint_end-timepoint_start;
                                            if CR_duration > dur
                                                end_search='off';
                                                search='off';
                                                CR=CR+1;
                                                noise_resp(indice_save,3)=1;
                                                noise_resp(indice_save,4)=timepoint_start;
                                                noise_resp(indice_save,5)=timepoint_end;
                                                noise_resp(indice_save,6)=CR_duration;
                                                disp('Ouh! Definitely got one!')
                                            else
                                                end_search='off';
                                            end
                                        end
                                        mm=mm+1;
                                    end
                                end
                                ll=ll+1;
                            end
                        end
                    end
                    kk=kk+1;
                end
                
            end
%             prog=k/cr2;
%             fprintf('Trial crunching progress %d \r', prog);
            k=k+1;
        end
        disp(['finished trial ', int2str(i), ' condition ', int2str(j)])
    end
    CR_block(j,1)=CR;
    CR=0;
    total_CR=sum(CR_block);
end
disp(['Aaaaand we got a total of ', int2str(total_CR), ' CRs!!'])
output1=CR_block; output2=noise_resp;
end