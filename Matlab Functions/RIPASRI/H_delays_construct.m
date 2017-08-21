% H_delays_construct
% Create a function to enter the H-reflex delays of subjects newly recorded
% in a matlab matrix in which we already have the subjects ID

%% INITIALIZATION
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


%% MATRIX CREATION
hdelays.names={};
hdelays.delays=zeros(1,2);

for i = 1:7
    hdelays(i).names=input('Subject ID:  ');
    hdelays(i).delays=input('H-reflex delays (into square braquets):  ');
end

%% SAVING

cd(gp_folder);
save('hdelays','hdelays');


