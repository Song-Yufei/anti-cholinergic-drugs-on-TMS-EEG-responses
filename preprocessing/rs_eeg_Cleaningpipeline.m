% This script illustrates the steps of resting eeg data preprocessing for
% an individual subject. Some scripts were adpated from Nigel Rogasch
% https://github.com/nigelrogasch/DXM_TMS-EEG_paper/tree/master/pre_processing/resting_state_pipeline

% fieldtrip(fieldtrip-20230613);
% eeglab(v2023.0): TESA(v1.2.0),FastICA_25 should be installed in the plugins;
% requiredFunction folder;
% original_chanlocs.mat
%% Step 1: Import data to EEGLAB
% load data
 eeglab;
 EEG = pop_loadset('filepath','\example_dataset_before_cleaning\',...
'filename', 'CHODME_031_1_rsEEG_PRE.set');

[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
 EEG = eeg_checkset(EEG);
 eeglab redraw     
%% Step 2: Remove unused channels 
EEG = pop_select( EEG, 'nochannel',{'APBr', 'FDIr','VEOG','HEOG'});
%% Step 3: Downsample data (5000 Hz to 1000 Hz)
EEG = pop_resample( EEG, 1000);
%% Step 4: Filter 1~100 Hz
clearvars beffEEG affEEG
beffEEG = EEG.data;
affEEG  = ft_preproc_bandpassfilter(beffEEG,EEG.srate,[1,100],3,'but','twopass','hamming'); %1-100Hz was used in the first version
% order=3, https://www.fieldtriptoolbox.org/example/irasa/                   
EEG.data = [];
EEG.data = affEEG;
%% Step 5: Epoch into 2.5s segments 
% Generate triggers every 2.5 seconds (exclude first 10s and last 10s)
% 5s is to ensure 1Hz can have a 5 cycle long data;
% 2.5s is to ensure 2Hz can have a 5 cycle long data
s10 = EEG.srate*10;% Find 10 s from start
e10 = EEG.times(1,end)-EEG.srate*10; % Find 10 s from end
eventTimes = s10:EEG.srate*2.5:e10; % Calculate vector with 2.5 s intervals

% find out pre or pos and add event marker
if ~isempty( strfind(EEG.setname,'PRE'))
    trigName = '1';
elseif ~isempty( strfind(EEG.setname,'POS') )                      
    trigName = '2';
end 

idx = length(EEG.event);
for ev = 1:length(eventTimes)
    EEG.event(idx+ev).latency = eventTimes(ev);
    EEG.event(idx+ev).duration = NaN;
    EEG.event(idx+ev).channel = 0;
    EEG.event(idx+ev).type = trigName;
    EEG.event(idx+ev).urevent = idx+1;
    EEG.urevent(idx+ev).latency = eventTimes(ev);
    EEG.urevent(idx+ev).duration = NaN;
    EEG.urevent(idx+ev).channel = 0;
    EEG.urevent(idx+ev).type = trigName;
end

[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET);% creates a new dataset in variable ALLEEG.
EEG = eeg_checkset(EEG);
eeglab redraw  

% Epoch the data
EEG = pop_epoch( EEG, {trigName}, [0  2.5], 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset(EEG);
eeglab redraw;    

% Repeat Step1-5 for 'CHODME_031_1_rsEEG_POS.set'
% Concatenate Pre- and Pos- data
EEG = pop_mergeset(ALLEEG, 1:2, 0); 
%% Step 6: Remove heavily contaminated trial and channels (Visual inspection)
% visualize
pop_eegplot( EEG, 1, 1, 1);
% ---insert channel labels to remove
EEG.Badchan = {'P5','P6','AF7'}; % these three channels are removed from all subjects' data due a technical error
EEG = pop_select( EEG, 'nochannel',EEG.Badchan);

% ####Bad Trials####
% visualize
EEG = pop_jointprob(EEG,1,1:size(EEG.data,1) ,5,3,0,0);
pop_rejmenu(EEG,1);
pause_script = input('Highlight bad trials, update marks and then press enter');
% in the GUI, select scroll data in the top right corner. 
% hit update marks and close the main trial rejection window. 
% hit enter in the command window.
EEG.BadTr = unique([find(EEG.reject.rejjp==1) find(EEG.reject.rejmanual==1)]);
EEG = pop_rejepoch( EEG, EEG.BadTr ,0); 

%% Step 7: ICA to remove ocular artifact
EEG = pop_tesa_removedata( EEG, EEG.artCut);
EEG = runICA_TEP_on_EEG(EEG, []); % use defaut 35 ICs
% plot topographies and write down the indx for eye blink and movement.
figure; 
for i=1:35
    subplot(5,7,i)
    topoplot(EEG.icawinv(:,i), EEG.chanlocs);
end

EEG = pop_iclabel(EEG, 'default');
EEG.ICRemoved = [1 2]; % Insert number of ocular components 
EEG = pop_subcomp( EEG, EEG.ICRemoved, 0);
%% Step 8: Interpolate channels and re-reference
EEG = pop_interp( EEG,original_chanlocs, 'spherical');
EEG = pop_reref( EEG, []);

% save
