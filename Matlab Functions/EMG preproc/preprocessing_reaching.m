function preprocessing_reaching

c_dir=cd;


%% Future input arguments and global variables

w_dir1='E:\Jesus\New subjects\EEG';
w_dir2='E:\Jesus\New subjects\Events\EEG';
ch_file=[c_dir '\channels_location(64).ced'];
dip_path_hdm='E:\Jesus\eeglab11_0_5_4b\plugins\dipfit2.2\standard_BESA\standard_BESA.mat';
dip_path_mri='E:\Jesus\eeglab11_0_5_4b\plugins\dipfit2.2\standard_BESA\avg152t1.mat';
dip_path_chan='E:\Jesus\eeglab11_0_5_4b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp';
dip_coordformat='Spherical';
% dip_path_hdm='E:\Jesus\eeglab11_0_5_4b\plugins\dipfit2.2\standard_BEM\standard_vol.mat';
% dip_path_mri='E:\Jesus\eeglab11_0_5_4b\plugins\dipfit2.2\standard_BEM\standard_mri.mat';
% dip_path_chan='E:\Jesus\eeglab11_0_5_4b\plugins\dipfit2.2\standard_BEM\elec\standard_1005.elc';
% dip_coordformat='MNI';
dip_threshold=15;
sr=2048;
sr_new=128;
low_cf=1;
high_cf=40;
filt_order=8;
rule_ICA=20;
prestimulus=0.4;
subject='8';
blocks='reaching';
sufix='_event';
n_frames=av_n_frames_reach(w_dir2,{'1','2','3','4','6','8','9'},...
    [blocks sufix],sr,sr_new);
preframes=round(prestimulus*128);
namefolder=['Reaching preprocessing subject ' subject];
mkdir(namefolder);


%% Preprocessing on continuous data

% Data loading
cd(w_dir1);
cd(['subject_' subject]);
disp('Loading data... ');
reaching_original=pop_biosig(blocks);
reaching_original.setname=['Sub' subject '_reaching_original'];
cd(c_dir);
cd(namefolder);
save('reaching_original.mat','reaching_original','-v7.3');

% Extra channels removing and channel locations reading
disp('Removing channels from 65 to 88...');
reaching_ch_removed1=reaching_original;
reaching_ch_removed1.data=reaching_ch_removed1.data(1:64,:);
reaching_ch_removed1.nbchan=64;
reaching_ch_removed1.chanlocs=readlocs(ch_file);
reaching_ch_removed1.setname=['Sub' subject '_reaching_ch_removed1'];
save('reaching_ch_removed1.mat','reaching_ch_removed1');

% Data resampling
disp('Resampling data...');
reaching_resampled=pop_resample(reaching_ch_removed1,sr_new);        
reaching_resampled.setname=['Sub' subject '_reaching_resampled'];          
save('reaching_resampled.mat','reaching_resampled','-v7.3');

% Data filtering
disp('Filtering data...');
reaching_filtered=reaching_resampled;
reaching_filtered.data=reaching_filtered.data-repmat(mean(...
    reaching_filtered.data,2),[1 size(reaching_filtered.data,2)]);
h1=fdesign.highpass('N,F3dB',filt_order,low_cf,sr_new);
d1=design(h1,'Butter');
h2=fdesign.lowpass('N,F3dB',filt_order,high_cf,sr_new);
d2=design(h2,'Butter');
for nn=1:size(reaching_filtered.data,1)
    reaching_filtered.data(nn,:)=filtfilt(d1.sosMatrix,d1.ScaleValues,...
        reaching_filtered.data(nn,:));
    reaching_filtered.data(nn,:)=filtfilt(d2.sosMatrix,d2.ScaleValues,...
        reaching_filtered.data(nn,:));
end
reaching_filtered.setname=['Sub' subject '_reaching_filtered'];
save('reaching_filtered.mat','reaching_filtered','-v7.3');

% Channels visual removal
disp('Channels visual removal...');
reaching_ch_removed2=reaching_filtered;
eegplot(reaching_ch_removed2.data,'srate',sr_new);
bad_chans=input('Bad channels (vector): ');
good_chans=setdiff(1:64,bad_chans);
n_channels=length(good_chans);
reaching_ch_removed2.data=reaching_ch_removed2.data(good_chans,:);
reaching_ch_removed2.nbchan=n_channels;
reaching_ch_removed2.chanlocs=reaching_ch_removed2.chanlocs(good_chans);
reaching_ch_removed2.setname=['Sub' subject '_reaching_ch_removed2'];
save('reaching_ch_removed2.mat','reaching_ch_removed2');

% Common average re-referencing
disp('Re-referencing data...');
reaching_rereferenced=reaching_ch_removed2;
reaching_rereferenced.data=reaching_rereferenced.data-repmat(mean...
    (reaching_rereferenced.data,1),[size(reaching_rereferenced.data,1) 1]);
reaching_rereferenced.setname=['Sub' subject '_reaching_rereferenced'];
save('reaching_rereferenced.mat','reaching_rereferenced');

% Independent Component Analisys (ICA): Infomax
disp('Running Infomax ICA...');
tmpdata=reaching_rereferenced.data-repmat(mean(...
    reaching_rereferenced.data,2),[1 size(reaching_rereferenced.data,2)]);
n_comp=min([rank(tmpdata),floor(sqrt((size(reaching_rereferenced.data,...
    2)/rule_ICA)))]);
reaching_ICA1=pop_runica(reaching_rereferenced,'icatype','runica',...
    'extended',1,'pca',n_comp,'stop',1e-7);
reaching_ICA1.setname=['Sub' subject '_reaching_ICA1'];
save('reaching_ICA1.mat','reaching_ICA1');

% Condition number revision
cond_num=cond(reaching_ICA1.icaweights);
disp(['Cond. number -> ' num2str(cond_num)]);

% Visual inspection and artifacts removal
disp('Independent Components Inspection...');
reaching_IC_visual_removal1=reaching_ICA1;
pop_topoplot(reaching_IC_visual_removal1,0,1:n_comp,...
    'Independent Components');
eegplot(eeg_getdatact(reaching_IC_visual_removal1,'component',1:n_comp,...
    'reshape','2d'));
bad_IC=input('Bad componets (vector): ');
if(~isempty(bad_IC))
    projection=reaching_IC_visual_removal1.icawinv(:,bad_IC)*...
        eeg_getdatact(reaching_IC_visual_removal1,'component',bad_IC,...
        'reshape','2d');
    reaching_IC_visual_removal1.data=reaching_IC_visual_removal1.data-...
        projection;
    good_IC=setdiff(1:size(reaching_IC_visual_removal1.icaweights,1),...
        bad_IC);
    reaching_IC_visual_removal1.icawinv=...
        reaching_IC_visual_removal1.icawinv(:,good_IC);
    reaching_IC_visual_removal1.icaweights=...
        reaching_IC_visual_removal1.icaweights(good_IC,:);
end
reaching_IC_visual_removal1.reject=[];
reaching_IC_visual_removal1.setname=['Sub' subject ...
    '_reaching_IC_visual_removal1'];
save('reaching_IC_visual_removal1.mat','reaching_IC_visual_removal1');


%% Trials concatenating

disp('Extraction and trials concatenating... ');
data=[];
events=struct([]);
urevents=struct([]);
cd(w_dir2);
cd(['subject_' subject]);
load([blocks sufix]);
n_trials=length(hand_up_EEG);
trials=1:n_trials;
cd(c_dir);
for ii=1:n_trials
    block_start=round(hand_up_EEG(trials(ii))*(sr_new/sr));
    block_end=round(object_up_EEG(trials(ii))*(sr_new/sr));
    data_res=zeros(n_channels,n_frames);
    for nn=1:n_channels
        data_res(nn,:)=ricampiona(reaching_IC_visual_removal1.data(nn,...
            block_start:block_end),n_frames);
    end
    data=[data reaching_IC_visual_removal1.data(:,block_start-...
    preframes:block_start-1) data_res];
    event=struct('type',['trial_' num2str(trials(ii))],'latency',...
        (n_frames+preframes)*(trials(ii)-1)+1,'duration',n_frames+...
        preframes,'init_index',trials(ii),'init_time',NaN,'urevent',...
        trials(ii));
    urevent=struct('type',['trial_' num2str(trials(ii))],'latency',...
        (n_frames+preframes)*(trials(ii)-1)+1,'duration',n_frames+...
        preframes,'init_index',trials(ii),'init_time',NaN);
    events=struct([events event]);
    urevents=struct([urevents urevent]);
end
reaching_concatenated=pop_importdata('setname',['Sub' subject ...
    '_reaching_concatenated'],'data',data,'dataformat','array',...
    'chanlocs',reaching_IC_visual_removal1.chanlocs,'nbchan',...
    reaching_IC_visual_removal1.nbchan,'srate',...
    reaching_IC_visual_removal1.srate);
reaching_concatenated.event=events;
reaching_concatenated.urevent=urevents;
cd(namefolder);
save('reaching_concatenated.mat','reaching_concatenated');


%% Preprocessing on trials

% Trials visual inspection
disp('Trials visual inspection... ');
reaching_trials_visual_removal=reaching_concatenated;
eegplot(reaching_trials_visual_removal.data,'srate',sr_new,'events',...
    reaching_trials_visual_removal.event);
bad_trials=input('Bad trials (vector): ');
good_trials=setdiff(1:n_trials,bad_trials);
rm_data=[];
for nn=1:length(bad_trials)
    st=(bad_trials(nn)-1)*(n_frames+preframes)+1;
    rm_data=[rm_data st:st+n_frames+preframes-1];
end
reaching_trials_visual_removal.data(:,rm_data)=[];
reaching_trials_visual_removal.pnts=size(...
    reaching_trials_visual_removal.data,2);
reaching_trials_visual_removal.xmax=(...
    reaching_trials_visual_removal.pnts-1)/sr_new;
reaching_trials_visual_removal.times=linspace(...
    reaching_trials_visual_removal.xmin,...
    reaching_trials_visual_removal.xmax,...
    size(reaching_trials_visual_removal.data,2));
reaching_trials_visual_removal.event=...
    reaching_trials_visual_removal.event(good_trials);
reaching_trials_visual_removal.urevent=...
    reaching_trials_visual_removal.urevent(good_trials);
for nn=1:length(good_trials)
    reaching_trials_visual_removal.event(nn).latency=(nn-1)*(n_frames+...
        preframes)+1;
    reaching_trials_visual_removal.event(nn).init_index=nn;
    reaching_trials_visual_removal.event(nn).urevent=nn;
    reaching_trials_visual_removal.urevent(nn).latency=...
        reaching_trials_visual_removal.event(nn).latency;
    reaching_trials_visual_removal.urevent(nn).init_index=...
        reaching_trials_visual_removal.event(nn).init_index;
end
save('reaching_trials_visual_removal.mat','reaching_trials_visual_removal');

% Independent Component Analisys (ICA): Infomax
disp('Running Infomax ICA...');
tmpdata=reaching_trials_visual_removal.data-repmat(mean(...
    reaching_trials_visual_removal.data,2),[1 ...
    size(reaching_trials_visual_removal.data,2)]);
n_comp=min([rank(tmpdata),floor(sqrt((size(...
    reaching_trials_visual_removal.data,2)/rule_ICA)))]);
reaching_ICA2=pop_runica(reaching_trials_visual_removal,'icatype',...
    'runica','extended',1,'pca',n_comp,'stop',1e-7);
reaching_ICA2.setname=['Sub' subject '_reaching_ICA2'];
save('reaching_ICA2.mat','reaching_ICA2');

% Condition number revision
cond_num=cond(reaching_ICA2.icaweights);
disp(['Cond. number -> ' num2str(cond_num)]);

% Visual inspection and artifacts removal
disp('Independent Components Inspection...');
reaching_IC_visual_removal2=reaching_ICA2;
n_comp=size(reaching_IC_visual_removal2.icaweights,1);
pop_topoplot(reaching_IC_visual_removal2,0,1:n_comp,...
    'Independent Components');
eegplot(eeg_getdatact(reaching_IC_visual_removal2,'component',1:n_comp,...
    'reshape','2d'));
bad_IC=input('Bad componets (vector): ');
if(~isempty(bad_IC))
    projection=reaching_IC_visual_removal2.icawinv(:,bad_IC)*...
        eeg_getdatact(reaching_IC_visual_removal2,'component',bad_IC,...
        'reshape','2d');
    reaching_IC_visual_removal2.data=reaching_IC_visual_removal2.data-...
        projection;
    good_IC=setdiff(1:size(reaching_IC_visual_removal2.icaweights,1),...
        bad_IC);
    reaching_IC_visual_removal2.icawinv=...
        reaching_IC_visual_removal2.icawinv(:,good_IC);
    reaching_IC_visual_removal2.icaweights=...
        reaching_IC_visual_removal2.icaweights(good_IC,:);
end
reaching_IC_visual_removal2.reject=[];
reaching_IC_visual_removal2.setname=['Sub' subject ...
    '_reaching_IC_visual_removal2'];
save('reaching_IC_visual_removal2.mat','reaching_IC_visual_removal2');

% Dipole locations and bad ICs rejection
disp('Locating dipoles...');
reaching_dipole=reaching_IC_visual_removal2;
[~,c_transformation]=coregister(reaching_dipole.chanlocs,dip_path_chan,...
    'warp','auto','manual','off');
block=pop_dipfit_settings(reaching_dipole,'hdmfile',dip_path_hdm,'mrifile',...
    dip_path_mri,'chanfile',dip_path_chan,'coordformat',dip_coordformat,...
    'coord_transform',c_transformation);
block=pop_multifit(block,[],'rmout','on','threshold',dip_threshold);
dips=[];
dips=[dips block.dipfit.model(1,:).rv];
bad_IC=find(dips==1);
disp('Removing bad ICs...');
if(~isempty(bad_IC))
    projection=reaching_dipole.icawinv(:,bad_IC)*eeg_getdatact(...
        reaching_dipole,'component',bad_IC,'reshape','2d');
    reaching_dipole.data=reaching_dipole.data-projection;
    good_IC=setdiff(1:size(reaching_dipole.icaweights,1),bad_IC);
    reaching_dipole.icawinv=reaching_dipole.icawinv(:,good_IC);
    reaching_dipole.icaweights=reaching_dipole.icaweights(good_IC,:);
end
reaching_dipole.reject=[];
reaching_dipole.setname=['Sub' subject '_reaching_dipole'];
save('reaching_dipole.mat','reaching_dipole');

% Visual inspection and artifacts removal
disp('Independent Components Inspection...');
reaching_IC_visual_removal3=reaching_dipole;
n_comp=size(reaching_IC_visual_removal3.icaweights,1);
pop_topoplot(reaching_IC_visual_removal3,0,1:n_comp,...
    'Independent Components');
eegplot(eeg_getdatact(reaching_IC_visual_removal3,'component',1:n_comp,...
    'reshape','2d'));
bad_IC=input('Bad componets (vector): ');
if(~isempty(bad_IC))
    projection=reaching_IC_visual_removal3.icawinv(:,bad_IC)*...
        eeg_getdatact(reaching_IC_visual_removal3,'component',bad_IC,...
        'reshape','2d');
    reaching_IC_visual_removal3.data=reaching_IC_visual_removal3.data-...
        projection;
    good_IC=setdiff(1:size(reaching_IC_visual_removal3.icaweights,1),...
        bad_IC);
    reaching_IC_visual_removal3.icawinv=...
        reaching_IC_visual_removal3.icawinv(:,good_IC);
    reaching_IC_visual_removal3.icaweights=...
        reaching_IC_visual_removal3.icaweights(good_IC,:);
end
reaching_IC_visual_removal3.reject=[];
reaching_IC_visual_removal3.setname=['Sub' subject ...
    '_reaching_IC_visual_removal3'];
save('reaching_IC_visual_removal3.mat','reaching_IC_visual_removal3');

% Epoching and prestimilus baseline correction
disp('Epoching...');
reaching_epoched=reaching_IC_visual_removal3;
reaching_epoched=pop_epoch(reaching_epoched,{},[0 (preframes+n_frames-1)...
    /sr_new]);
n_epochs=size(reaching_epoched.data,3);
for nn=1:n_epochs
    [dataout,~]=rmbase(reaching_epoched.data(:,:,nn),[],1:preframes);
    reaching_epoched.data(:,:,nn)=dataout;
end
reaching_epoched.setname=['Sub' subject '_reaching_epoched'];
save('reaching_epoched.mat','reaching_epoched');

% Epochs visual inspection
disp('Epochs visual inspection...');
n_epochs=size(reaching_epoched.data,3);
eegplot(reaching_epoched.data,'srate',sr_new,'events',...
    reaching_epoched.event);
bad_epochs=input('Bad epochs (vector): ');
b_epochs=zeros(1,n_epochs);
b_epochs(bad_epochs)=1;
reaching_epochs_visual_removal=pop_rejepoch(reaching_epoched,b_epochs,0);
reaching_epochs_visual_removal.setname=['Sub' subject ...
    '_reaching_epochs_visual_removal'];
save('reaching_epochs_visual_removal.mat','reaching_epochs_visual_removal');

% Epochs averaging
disp('Averaging epochs...');
reaching_final=reaching_epochs_visual_removal;
n_epochs=size(reaching_final.data,3);
data_av=zeros(size(reaching_final.data,1),size(reaching_final.data,2));
for nn=1:n_epochs
    data_av=data_av+reaching_final.data(:,:,nn);
end
data_av=data_av/n_epochs;
reaching_final.data=data_av;
%reaching_final.event=struct([]);
%reaching_final.urevent=struct([]);
%reaching_final.epoch=struct([]);
reaching_final.setname=['Sub' subject '_reaching_final'];
save('reaching_final.mat','reaching_final');

% GFP
cd(c_dir);
reaching_GFP=struct('subject',subject,'GFP',GFP(reaching_final.data));
cd(namefolder);
save('reaching_GFP.mat','reaching_GFP');


%% Exporting

% CARTOOL
pop_writeeeg(reaching_final,['subject_' subject '_reaching.bdf'],'TYPE','BDF');
output=fopen(['subject_' subject '_reaching_electrodes.xyz'],'w');
len=length(reaching_final.chanlocs);
fprintf(output,'%u\t%u',len,1);
for nn=1:len
    fprintf(output,'\n%f\t%f\t%f\t%s',-reaching_final.chanlocs(nn).X,...
        -reaching_final.chanlocs(nn).Y,reaching_final.chanlocs(nn).Z,...
        reaching_final.chanlocs(nn).labels);
end
fclose(output);

% LORETA

cd(c_dir);