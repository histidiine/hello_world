% Routine to process data and create graphs for single subjects or whole
% groups for Reciprocal Inhibition experiments (SBSU - 2017)

routineProc = 'on';
while strcmp(routineProc,'on')
    
    disp('1- Process MEP responses and H-reflex of a single subject');
    disp('2- Process MEP responses and H-reflex  of all subjects of a group');
    disp('3- Process group response from already analyzed subjects data');
    disp('4- Create histograms of single subject or group MEPs and H-reflexes ');
    disp('5- Create plot graphs of single subject or group recordings');
    disp('6- Create box plot of results');
    
    choice = input('What do you want to do? (see choices above):   ');
    
    switch choice
        
        case 1 %_______Process MEP responses of a single subject
            RI_process_MAR2017('','','',1);
            
        case 2 %_______Process MEP responses of all subjects of a group
           
            
        case 3 %___Process group response from already analyzed subjects data
            
            
        case 4 %_______Create histograms of single subject or group MEPs
            
            
        case 5 %_____Create plot graphs of single subject or group recordings
            
            
        case 6 %_______Create box plot of results
            
            
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