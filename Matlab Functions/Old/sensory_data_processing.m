function sensory_data_processing

%Extracting data from sensory excel file
cd('\\HOME\austepha\Mes Documents\Austepha_PhD\1_MainProject\Data\Sensory Data\2PD')
open 'sensory_matlabExtract.xlsx' % 1:2:16 = distance; 2:2:16 = response >> 4 sujets / 2 hands
pause;
sensory_raw = sensorymatlabExtract; clear sensorymatlabExtract;
%pilot1= sensory_raw(:,1:4);

sum = zeros(8,2);
distances= zeros(8,1);distances(2,1)=0.7; distances(3,1)=1; distances(4,1)=1.3; distances(5,1)=1.6; distances(6,1)=1.9;distances(7,1)=2.2; distances(8,1)=2.5;  
testing = zeros(8,2); % to check that there is indeed seven trials for each distance

for k=1:4 % sujets
    for j=1:size(sensory_raw,1) % trials
        for test=1:size(sum,1)
            %Main Droite
            if sensory_raw(j,4*k-3) == distances(test,1)
                sum(test,1)= sum(test,1) + sensory_raw(j,4*k-2);
                testing(test,1) = testing(test,1)+1;
                break
            end
            %Main Gauche
        end
        for test2=1:size(sum,1)
            if sensory_raw(j,4*k-1) == distances(test2,1)
                sum(test2,2)= sum(test2,2) + sensory_raw(j,4*k);
                testing(test2,2) = testing(test2,2)+1;
                break
            end
        end
    end
    pilot1= zeros(8,3); pilot1(:,1)=distances; pilot1(:,2:3)=sum; pilot1(1,2:3)=pilot1(1,2:3)/7;
    pilot1=flipud(pilot1);
    all_subjects{k}=pilot1;
    
    % Remise à 0
    testing = zeros(8,2); 
    sum = zeros(8,2);
end

save all_subjects

%% GRAPH

figure;
for k = 1:4
    curr=all{k};
    subplot(4,2,2*k-1)
    scatter(curr(:,1),curr(:,2),'filled'); 
    hold on
    subplot(4,2,2*k)
    scatter(curr(:,1),curr(:,3),'filled');
    clear curr
end
%clear
end
