% Routine to process data and create graphs for single subjects or whole
% groups for PAS experiments (SBSU - 2017)

% NOTES: 
% 1) I have to add a report file in which all disp outputs are written for future checking 

routineProc = 'on';
while strcmp(routineProc,'on')
    
    disp('1- Process MEP responses of a single subject');
    disp('2- Process MEP responses of all subjects of a group');
    disp('3- Process group response from already analyzed subjects data');
    disp('4- Create histograms of single subject or group MEPs');
    disp('5- Create plot graphs of single subject or group recordings');
    disp('6- Create box plot of results');
    
    choice = input('What do you want to do? (see choices above):   ');
    
    switch choice
        
        case 1 %_______Process MEP responses of a single subject
            PAS_process_MAR2017('','','',1);
            
        case 2 %_______Process MEP responses of all subjects of a group
            %______________Type of intervention______________
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
            %___________________Group choice___________________
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
            
            n_subjects = length(list_folders);
            
            %___________________Process loop___________________
            
            subChoice = 0;
            for sub = 1:n_subjects
                disp(['Starting ', list_folders(sub).name, '...']);
                PAS_process_MAR2017(PAS_ID, gp_ID,list_folders(sub).name, subChoice);
                disp(['Finished ', list_folders(sub).name, ' !']);
            end
            subChoice = 1;
            
        case 3 %___Process group response from already analyzed subjects data
            group_analysis;
            
        case 4 %_______Create histograms of single subject or group MEPs
            MEP_histo;
            
        case 5 %_____Create plot graphs of single subject or group recordings
            MEP_plotgraph;
            
        case 6 %_______Create box plot of results
            MEP_boxplot;
            
    end
    
    continueRep=0;
    while continueRep ~= 1 && continueRep~= 2
        disp('1- Yes, pleaaase');
        disp('2- No, no more I beg you')
        continueRep=input('Do you want to do another processing / result illustration (see above)?   ');
        if continueRep == 1
            routineProc = 'on';
        else if continueRep == 2
                routineProc = 'off';
                disp('End of FIGURE CREATION - Thank you for using PhDLife Airline !');
            end
        end
    end
end