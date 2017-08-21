function preprocessing_baseline_NEW

c_dir=cd;


%% Future input arguments and global variables

w_dir1='E:\Data\NewData_EEG_EMG_grasping\New subjects\EEG';
w_dir2='E:\Data\NewData_EEG_EMG_grasping\New subjects\Events\EEG';
ch_file=[c_dir '\channels_location(64).ced'];
dip_path_hdm='C:\Program Files (x86)\MATLAB\R2013a\toolbox\eeglab12_0_2_5b\plugins\dipfit2.2\standard_BESA\standard_BESA.mat';
dip_path_mri='C:\Program Files (x86)\MATLAB\R2013a\toolbox\eeglab12_0_2_5b\plugins\dipfit2.2\standard_BESA\avg152t1.mat';
dip_path_chan='C:\Program Files (x86)\MATLAB\R2013a\toolbox\eeglab12_0_2_5b\plugins\dipfit2.2\standard_BESA\standard-10-5-cap385.elp';
dip_coordformat='Spherical';
% dip_path_hdm='E:\Jesus\eeglab11_0_5_4b\plugins\dipfit2.2\standard_BEM\standard_vol.mat';
% dip_path_mri='E:\Jesus\eeglab11_0_5_4b\plugins\dipfit2.2\standard_BEM\standard_mri.mat';
% dip_path_chan='E:\Jesus\eeglab11_0_5_4b\plugins\dipfit2.2\standard_BEM\elec\standard_1005.elc';
% dip_coordformat='MNI';
dip_threshold=15;
sr=2048;
sr_new=128;
low_cf=1;
high_cf=20;
filt_order=8;
rule_ICA=20;
subject='5';
blocks={'resting1','resting2'};
n_blocks=length(blocks);
namefolder=['Baseline preprocessing subject NEW' subject];
mkdir(namefolder);

%% Preprocessing

% Data loading
cd(['E:\Data\NewData_EEG_EMG_grasping\New subjects\Subject ', subject, '\Baseline preprocessing subject ', subject])
load('base_original.mat');

% Extra channels removing and channel locations reading
disp('Removing channels from 65 to 88...');
base_ch_removed1=base_original;
base_ch_removed1.data=base_ch_removed1.data(1:64,:);
base_ch_removed1.nbchan=64;
base_ch_removed1.chanlocs=readlocs(ch_file);
base_ch_removed1.setname=['Sub' subject '_resting_ch_removed1'];
save('base_ch_removed1.mat','base_ch_removed1');
clear base_original

% Data filtering
disp('Filtering data...');
base_filtered=base_ch_removed1;
base_filtered.data = double(base_filtered.data);
base_filtered.data=base_filtered.data-repmat(mean(base_filtered.data,2),...
    [1 size(base_filtered.data,2)]);
h1=fdesign.highpass('N,F3dB',filt_order,low_cf,sr_new);
d1=design(h1,'Butter');
h2=fdesign.lowpass('N,F3dB',filt_order,high_cf,sr_new);
d2=design(h2,'Butter');
freqz(d1)
freqz(d2)
for nn=1:size(base_filtered.data,1)
    base_filtered.data(nn,:)=filtfilt(d1.sosMatrix,d1.ScaleValues,...
        base_filtered.data(nn,:));
    base_filtered.data(nn,:)=filtfilt(d2.sosMatrix,d2.ScaleValues,...
        base_filtered.data(nn,:));
end
base_filtered.setname=['Sub' subject '_resting_filtered'];
save('base_filtered.mat','base_filtered','-v7.3');
clear base_ch_removed1


% Data resampling
disp('Resampling data...');
base_resampled=pop_resample(base_filtered,sr_new);        
base_resampled.setname=['Sub' subject '_resting_resampled'];          
save('base_resampled.mat','base_resampled','-v7.3');
clear base_filtered

% Channels visual removal
disp('Channels visual removal...');
base_ch_removed2=base_resampled;
eegplot(base_ch_removed2.data,'srate',sr_new);
bad_chans=input('Bad channels (vector): ');
good_chans=setdiff(1:64,bad_chans);
n_channels=length(good_chans);
base_ch_removed2.data=base_ch_removed2.data(good_chans,:);
base_ch_removed2.nbchan=n_channels;
base_ch_removed2.chanlocs=base_ch_removed2.chanlocs(good_chans);
base_ch_removed2.setname=['Sub' subject '_resting_ch_removed2'];
save('base_ch_removed2.mat','base_ch_removed2');

% Common average re-referencing
disp('Re-referencing data...');
base_rereferenced=base_ch_removed2;
base_rereferenced.data=base_rereferenced.data-repmat(mean...
    (base_rereferenced.data,1),[size(base_rereferenced.data,1) 1]);
base_rereferenced.setname=['Sub' subject '_resting_rereferenced'];
save('base_rereferenced.mat','base_rereferenced');

% Data scrolling
disp('Scrolling data...');
base_scrolled=base_rereferenced;
eegplot(base_scrolled.data,'srate',sr_new,'command','close');
pause;
rm_data=[];
for nn=1:size(TMPREJ,1)
    rm_data=[rm_data round(TMPREJ(nn,1)):round(TMPREJ(nn,2))];
end
base_scrolled.data(:,rm_data)=[];
base_scrolled.pnts=size(base_scrolled.data,2);
base_scrolled.xmax=(base_scrolled.pnts-1)/sr_new;
base_scrolled.times=linspace(base_scrolled.xmin,base_scrolled.xmax,...
    size(base_scrolled.data,2));
base_scrolled.setname=['Sub' subject '_resting_scrolled'];
save('base_scrolled.mat','base_scrolled');
save('rm_data1', 'rm_data')

% % Data scrolling
% disp('Scrolling data...');
% base_scrolled=base_rereferenced;
% eegplot(base_scrolled.data,'srate',sr_new,'command','close');
% pause;
% rm_data=[];
% TMPREJ= [];
% for nn=1:size(TMPREJ,1)
%     rm_data=[rm_data round(TMPREJ(nn,1)):round(TMPREJ(nn,2))];
% end
% data_rm = [];
% data_rm = base_scrolled.data(:,rm_data);
% base_scrolled.data=data_rm;
% base_scrolled.pnts=size(base_scrolled.data,2);
% base_scrolled.xmax=(base_scrolled.pnts-1)/sr_new;
% base_scrolled.times=linspace(base_scrolled.xmin,base_scrolled.xmax,...
%     size(base_scrolled.data,2));
% base_scrolled.setname=['Sub' subject '_resting_scrolled'];
% eegplot(base_scrolled.data,'srate',sr_new,'command','close');
% 
% save('base_scrolled.mat','base_scrolled');

% Independent Component Analisys (ICA): Infomax
disp('Running Infomax ICA...');
tmpdata=base_scrolled.data-repmat(mean(base_scrolled.data,2),...
    [1 size(base_scrolled.data,2)]);
n_comp=min([rank(tmpdata),floor(sqrt((size(base_scrolled.data,2)/...
    rule_ICA)))]);
base_ICA=pop_runica(base_scrolled,'icatype','runica','extended',1,'pca',...
    n_comp,'stop',1e-7);
base_ICA.setname=['Sub' subject '_resting_ICA'];
save('base_ICA.mat','base_ICA');

% Condition number revision
cond_num=cond(base_ICA.icaweights);
disp(['Cond. number -> ' num2str(cond_num)]);

% Visual inspection and artifacts removal
disp('Independent Components Inspection...');
base_IC_visual_removal1=base_ICA;
pop_topoplot(base_IC_visual_removal1,0,1:n_comp,'Independent Components');
eegplot(eeg_getdatact(base_IC_visual_removal1,'component',1:n_comp,...
    'reshape','2d'));
bad_IC=input('Bad componets (vector): ');
if(~isempty(bad_IC))
    projection=base_IC_visual_removal1.icawinv(:,bad_IC)*...
        eeg_getdatact(base_IC_visual_removal1,'component',bad_IC,...
        'reshape','2d');
    base_IC_visual_removal1.data=base_IC_visual_removal1.data-projection;
    good_IC=setdiff(1:size(base_IC_visual_removal1.icaweights,1),bad_IC);
    base_IC_visual_removal1.icawinv=...
        base_IC_visual_removal1.icawinv(:,good_IC);
    base_IC_visual_removal1.icaweights=...
        base_IC_visual_removal1.icaweights(good_IC,:);
end
base_IC_visual_removal1.reject=[];
base_IC_visual_removal1.setname=['Sub' subject '_resting_IC_visual_removal1'];
save('base_IC_visual_removal1.mat','base_IC_visual_removal1');

% Dipole locations and bad ICs rejection
disp('Locating dipoles...');
base_dipole=base_IC_visual_removal1;
[~,c_transformation]=coregister(base_dipole.chanlocs,dip_path_chan,...
    'warp','auto','manual','off');
block=pop_dipfit_settings(base_dipole,'hdmfile',dip_path_hdm,'mrifile',...
    dip_path_mri,'chanfile',dip_path_chan,'coordformat',dip_coordformat,...
    'coord_transform',c_transformation);
block=pop_multifit(block,[],'rmout','on','threshold',dip_threshold);
dips=[];
dips=[dips block.dipfit.model(1,:).rv];
bad_IC=find(dips==1);
disp('Removing bad ICs...');
if(~isempty(bad_IC))
    projection=base_dipole.icawinv(:,bad_IC)*eeg_getdatact(...
        base_dipole,'component',bad_IC,'reshape','2d');
    base_dipole.data=base_dipole.data-projection;
    good_IC=setdiff(1:size(base_dipole.icaweights,1),bad_IC);
    base_dipole.icawinv=base_dipole.icawinv(:,good_IC);
    base_dipole.icaweights=base_dipole.icaweights(good_IC,:);
end
base_dipole.reject=[];
base_dipole.setname=['Sub' subject '_resting_dipole'];
save('base_dipole.mat','base_dipole');

% Visual inspection and artifacts removal
disp('Independent Components Inspection...');
base_IC_visual_removal2=base_dipole;
n_comp=size(base_IC_visual_removal2.icaweights,1);
pop_topoplot(base_IC_visual_removal2,0,1:n_comp,'Independent Components');
eegplot(eeg_getdatact(base_IC_visual_removal2,'component',1:n_comp,...
    'reshape','2d'));
bad_IC=input('Bad componets (vector): ');
if(~isempty(bad_IC))
    projection=base_IC_visual_removal2.icawinv(:,bad_IC)*...
        eeg_getdatact(base_IC_visual_removal2,'component',bad_IC,...
        'reshape','2d');
    base_IC_visual_removal2.data=base_IC_visual_removal2.data-projection;
    good_IC=setdiff(1:size(base_IC_visual_removal2.icaweights,1),bad_IC);
    base_IC_visual_removal2.icawinv=...
        base_IC_visual_removal2.icawinv(:,good_IC);
    base_IC_visual_removal2.icaweights=...
        base_IC_visual_removal2.icaweights(good_IC,:);
end
base_IC_visual_removal2.reject=[];
base_IC_visual_removal2.setname=['Sub' subject '_resting_IC_visual_removal2'];
save('base_IC_visual_removal2.mat','base_IC_visual_removal2');

% Baseline correction (zero mean)
disp('Removing mean...');
base_zero=base_IC_visual_removal2;
base_zero.data=base_zero.data-repmat(mean(base_zero.data,2),...
    [1 size(base_zero.data,2)]);
base_zero.setname=['Sub' subject '_resting_final'];
save('base_zero.mat','base_zero');

% Region of interest selection
TMPREJ = [];
disp('Scrolling data...');
base_final=base_zero;
eegplot(base_final.data,'srate',sr_new,'command','close');
pause;
rm_data=[];
for nn=1:size(TMPREJ,1)
    rm_data=[rm_data round(TMPREJ(nn,1)):round(TMPREJ(nn,2))];
end
base_final.data(:,rm_data)=[];
base_final.pnts=size(base_final.data,2);
base_final.xmax=(base_final.pnts-1)/sr_new;
base_final.times=linspace(base_final.xmin,base_final.xmax,...
    size(base_final.data,2));
base_final.setname=['Sub' subject '_resting_scrolled'];
save('base_final.mat','base_final');
save('rm_data2', 'rm_data')

% GFP
base_GFP=struct('subject',subject,'GFP',GFP(base_final.data));
save('base_GFP.mat','base_GFP');

base_cartool = base_final;
base_cartool.data = base_final.data(:,end-60*sr_new:end-1);
base_cartool.pnts=size(base_cartool.data,2);
base_cartool.xmax=(base_cartool.pnts-1)/sr_new;
base_cartool.times=linspace(base_cartool.xmin,base_cartool.xmax,...
    size(base_cartool.data,2));
base_cartool.setname=['Sub' subject '_resting_scrolled'];
save('base_cartool.mat','base_cartool');


%% Exporting

% CARTOOL
pop_writeeeg(base_cartool,['subject_' subject '_resting.bdf'],'TYPE','BDF');
output=fopen(['subject_' subject '_resting_electrodes.xyz'],'w');
len=length(base_cartool.chanlocs);
fprintf(output,'%u\t%u',len,1);
for nn=1:len
    fprintf(output,'\n%f\t%f\t%f\t%s',-base_final.chanlocs(nn).X,...
        -base_final.chanlocs(nn).Y,base_final.chanlocs(nn).Z,...
        base_final.chanlocs(nn).labels);
end
fclose(output);

% LORETA

cd(c_dir);