% .     MEP_graph_IND

% ______________INITIALIZATION_______________
dataset = ['b', 't1', 't2'];
ds_folder = {'Baseline', 'T1', 'T2'};
m_name={'ADM','FDI','FCR','FPB'};
step_name={'baseline','T1','T2'};
cond_name={'single', 'LICI 50', 'LICI 100', 'LICI 200'};
sub_ID_all={'BAM','ETB','CAM','ELB','CAR','LUM','ISH'};

% __________________________ID & DIRECTORY___________________________
PAS_ID=0;ini=0;
while strcmp(PAS_ID,'5hz')==0 && strcmp(PAS_ID,'02hz')==0 
    if ini == 1
        disp('SORRY ! That is not a valid group name, check out what I suggested in the parenthesis down here v')
    end
    PAS_ID = input('Please enter the code for the PAS protocole you want to process into single quotes (5hz, 02hz):   ');
    ini=1;
end

gp_ID=0;ini=0;
while strcmp(gp_ID,'HV')==0 && strcmp(gp_ID,'CRPS')==0 && strcmp(gp_ID,'CRPSD')==0 && strcmp(gp_ID,'FHD')==0 && strcmp(gp_ID,'FUNCT')==0
    if ini == 1
        disp('SORRY ! That is not a valid group name, check out what I suggested in the parenthesis down here v')
    end
    gp_ID = input('Please enter the code for the group you want to process into single quotes (HV,CRPS,CRPSD,FHD,FUNCT):   ');
    ini=1;
end

gp_folder=['/Users/histidiine/AMS Drive/PhD/Matlab/Matlab Data/',PAS_ID,'/',gp_ID];
cd(gp_folder)


%_________NUMBER OF SUBJECTS________
n_subjects=input('How many subjects do you want to process?:     ');

figure; count=1;
for sub = 1:n_subjects
    cd(gp_folder)
    sub_ID = sub_ID_all{sub};
    cd(sub_ID); cd('Data');
    load Resp_Amp
    resp_amp_mean=squeeze(mean(resp_amp,1));
    %     figure;
    %________HIST_________
    for muscle = 1:3
        subplot(n_subjects,3,count)
        curr_data=squeeze(resp_amp_mean(muscle,:,:));
        bar(curr_data,1,'edgecolor','k', 'linewidth', 2);
        colormap(bone);
        hold on
        set(gca,'XTickLabel',cond_name);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'NextPlot','add');
        axes;
        count = count+1;
    end
    resp_amp=squeeze(resp(:,:,:,:,5));
    MEP_amp(:,:,:,:,sub)=resp_amp;
end
ht=title(['MEP amplitude for ', sub_ID_all{sub}]);
set(gca,'Visible','off');
set(ht,'Visible','on');

%_________________GROUP___________



