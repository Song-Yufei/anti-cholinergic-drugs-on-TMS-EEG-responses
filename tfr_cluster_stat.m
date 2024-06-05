% this script:
% 1) load individual tms-eeg data after preprocessing, and do TFR calculations for TIOs
% and SIOs
% Note: in SMAData structure, SHAM_...from each subject is the dataset for computing SIOs
% 2) run cluster based permuation F/t tests on tfrs
% 3) compute effect size
%% Prepare data
clear, clc
subjies={'002','003', '004', '007', '009', '010', '011', '012',...
    '013', '014', '015', '017', '018', '019', '020', '021', '022','024','025','027','028','030','031','032'};

TARG={'SMA'};
DRUGS={'Place', 'Scop', 'Bipe'};
[~,sessdrugs]=xlsread('Z:\Experimental Data\2020-08 CHODME (processed)\YS_Sisyphus\SUJsessionInfo2');

    for ii=1:length(subjies)
        i=subjies{ii};
        da_row=ismember(cellfun(@(x) x(5:end), sessdrugs(:,1), 'un', 0), {i});
        for sess=1:3 
            da_cell=find(ismember(cellfun(@(x) x(4:end), sessdrugs(da_row,2:4), 'un', 0), DRUGS(sess)));        
            SESSION_NUMBER=cell2mat( cellfun(@(x) x(2), sessdrugs(da_row,da_cell+1), 'un', 0) );
%%%% Load individual TMS-EEG data after preprocessing     
            load (['...\CHODME_' i '_' SESSION_NUMBER '_' TARG '_TEPs(PREPOS)_cleaned_ft'],...
                'data_SIR_FT') 
            
            trialinfo = data_SIR_FT.trialinfo;     
            cfg = [];
            cfg.reref         = 'yes';
            cfg.refchannel    = {'all'};
            cfg.refmethod     = 'avg';
            cfg.baselinewindow  = [-0.8 -0.1]; 
            EEG_interp = ft_preprocessing(cfg,data_SIR_FT);

            % sort trials
            cfg=[];
            cfg.keeptrials = 'yes';
            cfg.trials = find(ismember(trialinfo, 2)); % pre-drug, sham tms
            SOUND_SHAM_PRE = ft_timelockanalysis(cfg, EEG_interp); 
            cfg.trials = find(ismember(trialinfo, 1)); % pre-drug, active tms
            SOUND_REAL_PRE = ft_timelockanalysis(cfg, EEG_interp);
            cfg.trials = find(ismember(trialinfo, 3)); % pos-drug, active tms
            SOUND_REAL_POS = ft_timelockanalysis(cfg, EEG_interp);
            cfg.trials = find(ismember(trialinfo, 4)); % pos-drug, sham tms
            SOUND_SHAM_POS = ft_timelockanalysis(cfg, EEG_interp);

            % evoked potential 
            cfg=[];
            cfg.keeptrials = 'no';
            cfg.trials = find(ismember(trialinfo, 2)); 
            SOUND_SHAM_PRE_eatme = ft_timelockanalysis(cfg, EEG_interp); 
            cfg.trials = find(ismember(trialinfo, 1)); 
            SOUND_REAL_PRE_eatme = ft_timelockanalysis(cfg, EEG_interp);
            cfg.trials = find(ismember(trialinfo, 3)); 
            SOUND_REAL_POS_eatme = ft_timelockanalysis(cfg, EEG_interp);
            cfg.trials = find(ismember(trialinfo, 4)); 
            SOUND_SHAM_POS_eatme = ft_timelockanalysis(cfg, EEG_interp);
            
            % put into ft data structure
            SOUND_SHAM_PRE_eatme.trial(1,:,:) = SOUND_SHAM_PRE_eatme.avg;
            SOUND_REAL_PRE_eatme.trial(1,:,:) = SOUND_REAL_PRE_eatme.avg;
            SOUND_REAL_POS_eatme.trial(1,:,:) = SOUND_REAL_POS_eatme.avg;
            SOUND_SHAM_POS_eatme.trial(1,:,:) = SOUND_SHAM_POS_eatme.avg;

            SOUND_SHAM_PRE_eatme.dimord = SOUND_SHAM_PRE_eatme.dimord;
            SOUND_REAL_PRE_eatme.dimord = SOUND_REAL_PRE_eatme.dimord;
            SOUND_REAL_POS_eatme.dimord = SOUND_REAL_POS_eatme.dimord;
            SOUND_SHAM_POS_eatme.dimord = SOUND_SHAM_POS_eatme.dimord;

            clearvars EEG_interp 

%%%% Wavelet Convolution
               freq = 3:1:60; % only use 4-40Hz in later analysis
                for fr = 1:length(freq) 
                      ncycles = 9/73*(freq(fr))+156/73; % freq of interest 4:40Hz has a cycle of around 3-7 cycles
                      cfg = [];
                      cfg.method     = 'wavelet';
                      cfg.width      = ncycles; 
                      cfg.output     = 'pow';
                      cfg.foi        = freq(fr);
                      cfg.toi        = SOUND_SHAM_PRE.time(1):0.01:SOUND_SHAM_PRE.time(end); 

                        cfg.keeptrials = 'yes';
                        induced_SOUND_SHAM_PRE(fr) = ft_freqanalysis(cfg, SOUND_SHAM_PRE);
                        induced_SOUND_REAL_PRE(fr) = ft_freqanalysis(cfg, SOUND_REAL_PRE);
                        induced_SOUND_REAL_POS(fr) = ft_freqanalysis(cfg, SOUND_REAL_POS);
                        induced_SOUND_SHAM_POS(fr) = ft_freqanalysis(cfg, SOUND_SHAM_POS);

                        cfg.keeptrials = 'no';
                        evoked_SOUND_SHAM_PRE(fr) = ft_freqanalysis(cfg, SOUND_SHAM_PRE_eatme);
                        evoked_SOUND_REAL_PRE(fr) = ft_freqanalysis(cfg, SOUND_REAL_PRE_eatme);
                        evoked_SOUND_REAL_POS(fr) = ft_freqanalysis(cfg, SOUND_REAL_POS_eatme);
                        evoked_SOUND_SHAM_POS(fr) = ft_freqanalysis(cfg, SOUND_SHAM_POS_eatme);
                end

            clearvars SOUND_SHAM_PRE SOUND_REAL_PRE SOUND_REAL_POS SOUND_SHAM_POS SOUND_SHAM_PRE_eatme SOUND_REAL_PRE_eatme...
                SOUND_REAL_POS_eatme  SOUND_SHAM_POS_eatme

        % Total (Evoked+Induced)
         IND_SOUND_SHAM_PRE = induced_SOUND_SHAM_PRE(1);
         IND_SOUND_SHAM_PRE.powspctrm = cat(3,induced_SOUND_SHAM_PRE.powspctrm);    
         IND_SOUND_SHAM_PRE.freq = cat(2,induced_SOUND_SHAM_PRE.freq); 

         IND_SOUND_REAL_PRE = induced_SOUND_REAL_PRE(1);
         IND_SOUND_REAL_PRE.powspctrm = cat(3,induced_SOUND_REAL_PRE.powspctrm);    
         IND_SOUND_REAL_PRE.freq = cat(2,induced_SOUND_REAL_PRE.freq);   

         IND_SOUND_SHAM_POS = induced_SOUND_SHAM_POS(1);
         IND_SOUND_SHAM_POS.powspctrm = cat(3,induced_SOUND_SHAM_POS.powspctrm);    
         IND_SOUND_SHAM_POS.freq = cat(2,induced_SOUND_SHAM_POS.freq);   

         IND_SOUND_REAL_POS = induced_SOUND_REAL_POS(1);
         IND_SOUND_REAL_POS.powspctrm = cat(3,induced_SOUND_REAL_POS.powspctrm);    
         IND_SOUND_REAL_POS.freq = cat(2,induced_SOUND_REAL_POS.freq);  
         
         %  Evoked  
         EV_SOUND_SHAM_PRE = evoked_SOUND_SHAM_PRE(1);
         EV_SOUND_SHAM_PRE.powspctrm=[];
         EV_SOUND_SHAM_PRE.powspctrm(1,:,:,:) = cat(2,evoked_SOUND_SHAM_PRE.powspctrm);
         EV_SOUND_SHAM_PRE.powspctrm = repmat(EV_SOUND_SHAM_PRE.powspctrm,[size(IND_SOUND_SHAM_PRE.powspctrm,1) 1 1 1 ]);

         EV_SOUND_REAL_PRE = evoked_SOUND_REAL_PRE(1);
         EV_SOUND_REAL_PRE.powspctrm=[];
         EV_SOUND_REAL_PRE.powspctrm(1,:,:,:) = cat(2,evoked_SOUND_REAL_PRE.powspctrm);
         EV_SOUND_REAL_PRE.powspctrm = repmat(EV_SOUND_REAL_PRE.powspctrm,[size(IND_SOUND_REAL_PRE.powspctrm,1) 1 1 1 ]);

         EV_SOUND_SHAM_POS = evoked_SOUND_SHAM_POS(1);
         EV_SOUND_SHAM_POS.powspctrm=[];
         EV_SOUND_SHAM_POS.powspctrm(1,:,:,:) = cat(2,evoked_SOUND_SHAM_POS.powspctrm);
         EV_SOUND_SHAM_POS.powspctrm = repmat(EV_SOUND_SHAM_POS.powspctrm,[size(IND_SOUND_SHAM_POS.powspctrm,1) 1 1 1 ]);

         EV_SOUND_REAL_POS = evoked_SOUND_REAL_POS(1);
         EV_SOUND_REAL_POS.powspctrm=[];
         EV_SOUND_REAL_POS.powspctrm(1,:,:,:) = cat(2,evoked_SOUND_REAL_POS.powspctrm);
         EV_SOUND_REAL_POS.powspctrm = repmat(EV_SOUND_REAL_POS.powspctrm,[size(IND_SOUND_REAL_POS.powspctrm,1) 1 1 1 ]);

        % Induced
        IND_SOUND_SHAM_PRE.powspctrm = IND_SOUND_SHAM_PRE.powspctrm-EV_SOUND_SHAM_PRE.powspctrm;
        IND_SOUND_REAL_PRE.powspctrm = IND_SOUND_REAL_PRE.powspctrm-EV_SOUND_REAL_PRE.powspctrm;
        IND_SOUND_SHAM_POS.powspctrm = IND_SOUND_SHAM_POS.powspctrm-EV_SOUND_SHAM_POS.powspctrm;
        IND_SOUND_REAL_POS.powspctrm = IND_SOUND_REAL_POS.powspctrm-EV_SOUND_REAL_POS.powspctrm;

     clearvars  induced_SOUND_SHAM_PRE induced_SOUND_REAL_PRE induced_SOUND_SHAM_POS induced_SOUND_REAL_POS...
         evoked_SOUND_SHAM_PRE evoked_SOUND_REAL_PRE evoked_SOUND_SHAM_POS evoked_SOUND_REAL_POS...
         EV_SOUND_SHAM_PRE  EV_SOUND_REAL_PRE EV_SOUND_SHAM_POS EV_SOUND_REAL_POS

%%%% Power Normalize with respect to the whole epoch (contain pre- and post-stimulus time; trial by trial)  
        ERSPfullTBz_4D = IND_SOUND_SHAM_PRE;
        avgPerFreq = nanmean(IND_SOUND_SHAM_PRE.powspctrm,4);% average over time bins
        avgPerFreqMat = repmat(avgPerFreq,[1 1 1 size(IND_SOUND_SHAM_PRE.powspctrm,4)]); % set average value for all time bins
        sdPerFreq = nanstd(IND_SOUND_SHAM_PRE.powspctrm,0,4); % sd over time bins
        sdPerFreqMat =  repmat(sdPerFreq,[1 1 1 size(IND_SOUND_SHAM_PRE.powspctrm,4)]); % set sd value for all time bins
        ERSPfullTBz_4D.powspctrm = (IND_SOUND_SHAM_PRE.powspctrm - avgPerFreqMat) ./ sdPerFreqMat; % z-transform per frequency relative to full epoch
        ERSPfullTBz_SHAM_PRE = ft_freqdescriptives([],ERSPfullTBz_4D); % average ERSPfullTBz over trials (4D-->3D)
        clear ERSPfullTBz_4D

        ERSPfullTBz_4D = IND_SOUND_SHAM_POS;
        avgPerFreq = nanmean(IND_SOUND_SHAM_POS.powspctrm,4);
        avgPerFreqMat = repmat(avgPerFreq,[1 1 1 size(IND_SOUND_SHAM_POS.powspctrm,4)]); 
        sdPerFreq = nanstd(IND_SOUND_SHAM_POS.powspctrm,0,4); 
        sdPerFreqMat =  repmat(sdPerFreq,[1 1 1 size(IND_SOUND_SHAM_POS.powspctrm,4)]); 
        ERSPfullTBz_4D.powspctrm = (IND_SOUND_SHAM_POS.powspctrm - avgPerFreqMat) ./ sdPerFreqMat; 
        ERSPfullTBz_SHAM_POS = ft_freqdescriptives([],ERSPfullTBz_4D); 
        clear ERSPfullTBz_4D   
 
        ERSPfullTBz_4D = IND_SOUND_REAL_PRE;
        avgPerFreq = nanmean(IND_SOUND_REAL_PRE.powspctrm,4);
        avgPerFreqMat = repmat(avgPerFreq,[1 1 1 size(IND_SOUND_REAL_PRE.powspctrm,4)]); 
        sdPerFreq = nanstd(IND_SOUND_REAL_PRE.powspctrm,0,4); 
        sdPerFreqMat =  repmat(sdPerFreq,[1 1 1 size(IND_SOUND_REAL_PRE.powspctrm,4)]); 
        ERSPfullTBz_4D.powspctrm = (IND_SOUND_REAL_PRE.powspctrm - avgPerFreqMat) ./ sdPerFreqMat; 
        ERSPfullTBz_REAL_PRE = ft_freqdescriptives([],ERSPfullTBz_4D); 
        clear ERSPfullTBz_4D
      
        ERSPfullTBz_4D = IND_SOUND_REAL_POS;
        avgPerFreq = nanmean(IND_SOUND_REAL_POS.powspctrm,4);
        avgPerFreqMat = repmat(avgPerFreq,[1 1 1 size(IND_SOUND_REAL_POS.powspctrm,4)]); 
        sdPerFreq = nanstd(IND_SOUND_REAL_POS.powspctrm,0,4); 
        sdPerFreqMat =  repmat(sdPerFreq,[1 1 1 size(IND_SOUND_REAL_POS.powspctrm,4)]); 
        ERSPfullTBz_4D.powspctrm = (IND_SOUND_REAL_POS.powspctrm - avgPerFreqMat) ./ sdPerFreqMat; 
        ERSPfullTBz_REAL_POS = ft_freqdescriptives([],ERSPfullTBz_4D);
        clear ERSPfullTBz_4D   
                  
%%%% Baseline correction with respect to [-0.5 -0.2] pre-stimulus
        cfg = [];
        cfg.baselinetype = 'absolute';
        cfg.baseline = [-0.5 -0.2]; % from 500 to 100 ms pre-TMS
        INDUCED_SHAM_PRE = ft_freqbaseline(cfg, ERSPfullTBz_SHAM_PRE);
        clear ERSPfullTBz_SHAM_PRE
        INDUCED_SHAM_POS = ft_freqbaseline(cfg, ERSPfullTBz_SHAM_POS);
        clear ERSPfullTBz_SHAM_POS
        INDUCED_REAL_PRE = ft_freqbaseline(cfg, ERSPfullTBz_REAL_PRE);
        clear ERSPfullTBz_REAL_PRE
        INDUCED_REAL_POS = ft_freqbaseline(cfg, ERSPfullTBz_REAL_POS);
        clear ERSPfullTBz_REAL_POS

           %%%%%%%%%%%%%%%%%  
            if sess==1 % PLA
                INDUCED_SHAM_PLA_PRE{ii}=INDUCED_SHAM_PRE; % Placebo, PRE, Sham
                INDUCED_SHAM_PLA_POS{ii}=INDUCED_SHAM_POS;% Placebo, PRE, Active
                INDUCED_REAL_PLA_PRE{ii}=INDUCED_REAL_PRE; % Placebo, POS, Sham
                INDUCED_REAL_PLA_POS{ii}=INDUCED_REAL_POS; % Placebo, POS, Active

            elseif sess==2 % SCO
                INDUCED_SHAM_SCO_PRE{ii}=INDUCED_SHAM_PRE;  % Scopolamine, PRE, Sham
                INDUCED_SHAM_SCO_POS{ii}=INDUCED_SHAM_POS;% Scopolamine, PRE, Active
                INDUCED_REAL_SCO_PRE{ii}=INDUCED_REAL_PRE; % Scopolamine, POS, Sham
                INDUCED_REAL_SCO_POS{ii}=INDUCED_REAL_POS;  % Scopolamine, POS, Active

            elseif sess==3 % BIP
                INDUCED_SHAM_BIP_PRE{ii}=INDUCED_SHAM_PRE;% Biperiden, PRE, Sham
                INDUCED_SHAM_BIP_POS{ii}=INDUCED_SHAM_POS; % Biperiden, PRE, Active
                INDUCED_REAL_BIP_PRE{ii}=INDUCED_REAL_PRE; % Biperiden, POS, Sham
                INDUCED_REAL_BIP_POS{ii}=INDUCED_REAL_POS; % Biperiden, POS, Active 
            end
        end % end of session
    end % end of subjects
  
%%%% Save data        
data_TF = struct();
% Set field names based on the current brain region
field_names = {'SHAM_PLA_PRE', 'ACTIVE_PLA_PRE', 'SHAM_PLA_POS', 'ACTIVE_PLA_POS', ...
               'SHAM_SCO_PRE', 'ACTIVE_SCO_PRE', 'SHAM_SCO_POS', 'ACTIVE_SCO_POS', ...
               'SHAM_BIP_PRE', 'ACTIVE_BIP_PRE', 'SHAM_BIP_POS', 'ACTIVE_BIP_POS'};

% Populate the structure with data
for i = 1:numel(field_names)
    data_TF.(field_names{i}) = eval(['INDUCED_' TARG '_' field_names{i}]);
end
% Save the structure to a file
filepath = 'Z:\Experimental Data\2020-08 CHODME (processed)\YS_Sisyphus\GroupData\TFR_v5\';
filename = [TARG '_TFRs_PREPOS'];
%save([filepath filename], ['data_' regions{Ti}], '-v7.3');
%% Compute condition specific TFR average over trials
clc, clear
% Load data
load('...\SMA_TFRs_PREPOS.mat')

Data_TF = SMAData_TF; 

INDUCED_SHAM_PLA_PRE = Data_TF.SHAM_PLA_PRE; 
INDUCED_REAL_PLA_PRE = Data_TF.ACTIVE_PLA_PRE;
INDUCED_SHAM_PLA_POS = Data_TF.SHAM_PLA_POS;
INDUCED_REAL_PLA_POS = Data_TF.ACTIVE_PLA_POS;

INDUCED_SHAM_SCO_PRE = Data_TF.SHAM_SCO_PRE;
INDUCED_REAL_SCO_PRE = Data_TF.ACTIVE_SCO_PRE;
INDUCED_SHAM_SCO_POS = Data_TF.SHAM_SCO_POS;
INDUCED_REAL_SCO_POS = Data_TF.ACTIVE_SCO_POS;

INDUCED_SHAM_BIP_PRE = Data_TF.SHAM_BIP_PRE;
INDUCED_REAL_BIP_PRE = Data_TF.ACTIVE_BIP_PRE;
INDUCED_SHAM_BIP_POS = Data_TF.SHAM_BIP_POS;
INDUCED_REAL_BIP_POS = Data_TF.ACTIVE_BIP_POS;

cfg = [];
cfg.channel   = 'all';
cfg.keepindividual   = 'yes';
cfg.parameter = 'powspctrm';
% Placebo
GA_TFRs_SHAM_PLA_PRE = ft_freqgrandaverage(cfg,INDUCED_SHAM_PLA_PRE{:});  
GA_TFRs_SHAM_PLA_POS = ft_freqgrandaverage(cfg,INDUCED_SHAM_PLA_POS{:});  
GA_TFRs_ACTIVE_PLA_PRE = ft_freqgrandaverage(cfg,INDUCED_REAL_PLA_PRE{:});  
GA_TFRs_ACTIVE_PLA_POS = ft_freqgrandaverage(cfg,INDUCED_REAL_PLA_POS{:});  

% Scopolamine
GA_TFRs_SHAM_SCO_PRE = ft_freqgrandaverage(cfg,INDUCED_SHAM_SCO_PRE{:});  
GA_TFRs_SHAM_SCO_POS = ft_freqgrandaverage(cfg,INDUCED_SHAM_SCO_POS{:});  
GA_TFRs_ACTIVE_SCO_PRE = ft_freqgrandaverage(cfg,INDUCED_REAL_SCO_PRE{:});  
GA_TFRs_ACTIVE_SCO_POS = ft_freqgrandaverage(cfg,INDUCED_REAL_SCO_POS{:});  

% Biperiden
GA_TFRs_SHAM_BIP_PRE = ft_freqgrandaverage(cfg,INDUCED_SHAM_BIP_PRE{:});  
GA_TFRs_SHAM_BIP_POS = ft_freqgrandaverage(cfg,INDUCED_SHAM_BIP_POS{:});  
GA_TFRs_ACTIVE_BIP_PRE = ft_freqgrandaverage(cfg,INDUCED_REAL_BIP_PRE{:});  
GA_TFRs_ACTIVE_BIP_POS = ft_freqgrandaverage(cfg,INDUCED_REAL_BIP_POS{:});  
   
% cleaned
cfg=[];
cfg.parameter = 'powspctrm';
cfg.operation =  'subtract';           
GA_TFRs_clean_PLA_PRE = ft_math(cfg, GA_TFRs_ACTIVE_PLA_PRE, GA_TFRs_SHAM_PLA_PRE);            
GA_TFRs_clean_PLA_POS = ft_math(cfg, GA_TFRs_ACTIVE_PLA_POS, GA_TFRs_SHAM_PLA_POS);
GA_TFRs_clean_SCO_PRE = ft_math(cfg, GA_TFRs_ACTIVE_SCO_PRE, GA_TFRs_SHAM_SCO_PRE);
GA_TFRs_clean_SCO_POS = ft_math(cfg, GA_TFRs_ACTIVE_SCO_POS, GA_TFRs_SHAM_SCO_POS);
GA_TFRs_clean_BIP_PRE = ft_math(cfg, GA_TFRs_ACTIVE_BIP_PRE, GA_TFRs_SHAM_BIP_PRE);
GA_TFRs_clean_BIP_POS = ft_math(cfg, GA_TFRs_ACTIVE_BIP_POS, GA_TFRs_SHAM_BIP_POS);

% If the labels are from original NeurOne eeg structure, update
vars = who;
matchingVars = vars(startsWith(vars, 'GA'));
for i = 1:length(matchingVars)
    tempt = [];
    tempt = eval(matchingVars{i});
    if ismember('Afz', tempt.label)
        tempt.label(strcmpi('Afz', tempt.label)) = {'AFz'};
    end

    assignin('base', matchingVars{i}, tempt);
end

% Average PRE across three sessions
% make the basic fd structure
GA_TFRs_SHAM_ALL_PRE = GA_TFRs_SHAM_PLA_PRE;
GA_TFRs_ACTIVE_ALL_PRE = GA_TFRs_ACTIVE_PLA_PRE;
GA_TFRs_clean_ALL_PRE = GA_TFRs_clean_PLA_PRE;
for sj = 1:24
   GA_TFRs_SHAM_ALL_PRE.powspctrm(sj,:,:,:) =  mean( cat(1,GA_TFRs_SHAM_PLA_PRE.powspctrm(sj,:,:,:), GA_TFRs_SHAM_SCO_PRE.powspctrm(sj,:,:,:),GA_TFRs_SHAM_BIP_PRE.powspctrm(sj,:,:,:)), 1);
   GA_TFRs_ACTIVE_ALL_PRE.powspctrm(sj,:,:,:) = mean( cat(1,GA_TFRs_ACTIVE_PLA_PRE.powspctrm(sj,:,:,:), GA_TFRs_ACTIVE_SCO_PRE.powspctrm(sj,:,:,:), GA_TFRs_ACTIVE_BIP_PRE.powspctrm(sj,:,:,:)) ,1);
   GA_TFRs_clean_ALL_PRE.powspctrm(sj,:,:,:) = mean( cat(1,GA_TFRs_clean_PLA_PRE.powspctrm(sj,:,:,:), GA_TFRs_clean_SCO_PRE.powspctrm(sj,:,:,:),GA_TFRs_clean_BIP_PRE.powspctrm(sj,:,:,:)), 1);
end

%% Cluster F-test
clear COND1 COND2 COND3
load('...\neighbour_layout64_custom.mat')
clearvars COND1 COND2 COND3
% 
cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'powspctrm';
COND1 = ft_math(cfg, GA_TFRs_clean_PLA_POS, GA_TFRs_clean_PLA_PRE);
% COND1 = ft_math(cfg, GA_TFRs_SHAM_PLA_POS, GA_TFRs_SHAM_PLA_PRE); % for SIOs
COND2 = ft_math(cfg, GA_TFRs_clean_BIP_POS, GA_TFRs_clean_BIP_PRE);
COND3 = ft_math(cfg, GA_TFRs_clean_SCO_POS, GA_TFRs_clean_SCO_PRE);

cfg = [];
cfg.latency = [0.05 0.65]; %the whole post stimulation time window
cfg.avgovertime = 'no';
cfg.channel     = {'all'};
cfg.avgoverchan = 'no'; 
cfg.frequency= [4 40];
cfg.avgoverfreq = 'no';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesFunivariate';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05; 
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2; 
cfg.neighbours = neighbour_layout.custneighbor.neighbours;  
cfg.tail = 1; % For a F-statistic, it only make sense to calculate the right tail 
cfg.clustertail = 1;
cfg.correcttail = 'no'; 
cfg.alpha = 0.05;
cfg.numrandomization =2000;
cfg.parameter        = 'powspctrm';

subj = size(COND1.powspctrm,1); %enter number of participants
design = zeros(2,3*subj);
%%% Cond 1
for i = 1:subj
  design(1,i) = i;
end
%%% Cond 2
for i = 1:subj
  design(1,subj+i) = i;
end
%%% Cond 3
for i = 1:subj
  design(1,2*subj+i) = i;
end
design(2,1:subj)        = 1; 
design(2,subj+1:2*subj) = 2;
design(2,2*subj+1:3*subj)  = 3;
cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;
stat = ft_freqstatistics(cfg, COND1, COND2, COND3);
%% Follow-up analysis: extract cluster content (from cluster-F test), in freq and time and channel
% load stat 
stat_inuse = ; % stat results from cluster F-test

% positive
posClust=[];
if isfield(stat_inuse,'posclusters')
    if ~isempty(stat_inuse.posclusters)
        for clusx =1: length(stat_inuse.posclusters)
            if stat_inuse.posclusters(clusx).prob<0.05 % 
                clear clusterMask
                clusterMask = stat_inuse.posclusterslabelmat == clusx;
                posClust.selCh{clusx} = stat_inuse.label(any(clusterMask, [2, 3]));
                posClust.minT{clusx}  = min(stat_inuse.time(any(clusterMask, [1, 2])));
                posClust.maxT{clusx}  = max(stat_inuse.time(any(clusterMask, [1, 2])));
                posClust.minF{clusx}  = min(stat_inuse.freq(any(clusterMask, [1, 3])));
                posClust.maxF{clusx}  = max(stat_inuse.freq(any(clusterMask, [1, 3])));
            end
        end
    end
end
posClust

%%%% Cluster content 
clearvars Channels_in Freqs_in Times_in
Channels_in = stat_inuse.label(any(stat_inuse.posclusterslabelmat ==1, [2, 3])); 
Freqs_in= [min(stat_inuse.freq(any(stat_inuse.posclusterslabelmat ==1, [1, 3]))), ...
    max(stat_inuse.freq(any(stat_inuse.posclusterslabelmat ==1, [1, 3])))];
Times_in = [min(stat_inuse.time(any(stat_inuse.posclusterslabelmat ==1, [1, 2]))), ...
    max(stat_inuse.time(any(stat_inuse.posclusterslabelmat ==1, [1, 2])))];

%% Follow-up analysis: pairwise cluster t-test
load('...\neighbour_layout64_custom.mat')

% contrastsI example:
clearvars COND1 COND2
COND1 = GA_TFRs_clean_BIP_PRE;
COND2 = GA_TFRs_clean_BIP_POS;

cfg = [];
cfg.latency = Times_in; 
cfg.avgovertime = 'no';
cfg.channel     = Channels_in;
cfg.avgoverchan = 'no'; 
cfg.frequency= Freqs_in;
cfg.avgoverfreq = 'no';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;  
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2; 
cfg.neighbours = neighbour_layout.custneighbor.neighbours;  
cfg.tail = 0; 
cfg.clustertail = 0;
cfg.alpha = 0.05/3; % correct for number of comparisions 3
cfg.correcttail = 'alpha'; 
cfg.numrandomization =2000;
cfg.parameter        = 'powspctrm';

subj = size(COND1.powspctrm,1);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;
cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

stat =[];
stat = ft_freqstatistics(cfg, COND2, COND1);
%% Effect size following cluster t-test
clearvars COND1 COND2 COND3 COND4
%%%% Step1: compute data contrasts %%%%
cond = {'PRE','POS','POSminusPRE','SCOminusPLA','BIPminusPLA','SCOminusBIP'};
condName = {'PRE','POS','POS-PRE'};

cfg=[];
cfg.operation =  'subtract';
cfg.parameter =  'powspctrm';
COND1.(cond{3}) =ft_math(cfg, GA_TFRs_SHAM_PLA_POS,GA_TFRs_SHAM_PLA_PRE);
COND1.(cond{2}) = GA_TFRs_SHAM_PLA_POS;
COND1.(cond{1}) = GA_TFRs_SHAM_PLA_PRE;

COND2.(cond{3}) =ft_math(cfg, GA_TFRs_SHAM_BIP_POS,GA_TFRs_SHAM_BIP_PRE);
COND2.(cond{2}) = GA_TFRs_SHAM_BIP_POS;
COND2.(cond{1}) = GA_TFRs_SHAM_BIP_PRE;

COND3.(cond{3}) =ft_math(cfg, GA_TFRs_SHAM_SCO_POS,GA_TFRs_SHAM_SCO_PRE);
COND3.(cond{2}) = GA_TFRs_SHAM_SCO_POS;
COND3.(cond{1}) = GA_TFRs_SHAM_SCO_PRE;

COND4.(cond{4}) =ft_math(cfg, COND3.(cond{3}),COND1.(cond{3})); % Time:sco-pla
COND4.(cond{5}) =ft_math(cfg, COND2.(cond{3}),COND1.(cond{3})); % Time:bip-pla
COND4.(cond{6}) =ft_math(cfg, COND3.(cond{3}),COND2.(cond{3})); % Time:sco-pla

subj = 24;
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

%%%% Step2: load stat and extract cluster content in space and time %%%%
clearvars Freqs_in Channels_in Times_in stat_inuse Sigmask statSelectTOIFOIROI
% load stat result from followup t-tests 
stat_inuse = ; % stat from t-tests
% Two Info about LOADED stat: 
% 1) Know pos or neg; 2) Know crtical alpha level 
clusterGender = 'Pos'; % 'Pos' or 'Neg'
alpha_lim = 0.025/3;

if isequal(clusterGender, 'Pos')
        % positive
        posClust=[];
        if isfield(stat_inuse,'posclusters')
            if ~isempty(stat_inuse.posclusters)
                for clusx =1: length(stat_inuse.posclusters)
                    if stat_inuse.posclusters(clusx).prob<alpha_lim % 
                        clear clusterMask
                        % Channel x Frequency x Time
                        clusterMask = stat_inuse.posclusterslabelmat == clusx;
                        posClust.selCh{clusx} = stat_inuse.label(any(clusterMask, [2, 3]));
                        posClust.minT{clusx}  = min(stat_inuse.time(any(clusterMask, [1, 2])));
                        posClust.maxT{clusx}  = max(stat_inuse.time(any(clusterMask, [1, 2])));
                        posClust.minF{clusx}  = min(stat_inuse.freq(any(clusterMask, [1, 3])));
                        posClust.maxF{clusx}  = max(stat_inuse.freq(any(clusterMask, [1, 3])));
                    end
                end
            end
        end
        posClust
        
        %%%% Cluster content 
        clearvars Channels_in Freqs_in Times_in Sigmask
        Channels_in = stat_inuse.label(any(stat_inuse.posclusterslabelmat ==1, [2, 3])); 
        Freqs_in= [min(stat_inuse.freq(any(stat_inuse.posclusterslabelmat ==1, [1, 3]))), ...
            max(stat_inuse.freq(any(stat_inuse.posclusterslabelmat ==1, [1, 3])))];
        Times_in = [min(stat_inuse.time(any(stat_inuse.posclusterslabelmat ==1, [1, 2]))), ...
            max(stat_inuse.time(any(stat_inuse.posclusterslabelmat ==1, [1, 2])))];
        Sigmask = (stat_inuse.posclusterslabelmat==1);
else
        % negative
        negClust=[];
        if isfield(stat_inuse,'negclusters')
            if ~isempty(stat_inuse.negclusters)
                for clusx =1: length(stat_inuse.negclusters)
                    if stat_inuse.negclusters(clusx).prob<alpha_lim % 
                        clear clusterMask
                        % Channel x Frequency x Time
                        clusterMask = stat_inuse.negclusterslabelmat == clusx;
                        negClust.selCh{clusx} = stat_inuse.label(any(clusterMask, [2, 3]));
                        negClust.minT{clusx}  = min(stat_inuse.time(any(clusterMask, [1, 2])));
                        negClust.maxT{clusx}  = max(stat_inuse.time(any(clusterMask, [1, 2])));
                        negClust.minF{clusx}  = min(stat_inuse.freq(any(clusterMask, [1, 3])));
                        negClust.maxF{clusx}  = max(stat_inuse.freq(any(clusterMask, [1, 3])));
                    end
                end
            end
        end
        negClust
        clearvars Channels_in Freqs_in Times_in Sigmask
        Channels_in = stat_inuse.label(any(stat_inuse.negclusterslabelmat ==1, [2, 3])); 
        Freqs_in= [min(stat_inuse.freq(any(stat_inuse.negclusterslabelmat ==1, [1, 3]))), ...
            max(stat_inuse.freq(any(stat_inuse.negclusterslabelmat ==1, [1, 3])))];
        Times_in = [min(stat_inuse.time(any(stat_inuse.negclusterslabelmat ==1, [1, 2]))), ...
            max(stat_inuse.time(any(stat_inuse.negclusterslabelmat ==1, [1, 2])))];
        Sigmask = (stat_inuse.negclusterslabelmat==1);
end

clearvars -except Freqs_in Channels_in Times_in stat_inuse Sigmask...
    COND1 COND2 COND3 COND4 design cond condName clusterGender

%%%% Step3: Select data %%%%
% 1) keep data structure's time and channels content the same as stat's time and channels content)
% 2) Know which contrast corresponds to the LOADED stat in the begining (Here I listed all possible contrasts but only report the one intended)
contrast_cond = {'pla_posMinuspla_pre', 'sco_posMinussco_pre', 'bip_posMinusbip_pre'...
    'sco_contrastMinus_pla_contrast', 'bip_contrastMinus_pla_contrast', 'sco_contrastMinus_bip_contrast'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option1: Average Over Circumscribed rectangle (lower bound)
clearvars pla_posMinuspla_pre sco_posMinussco_pre bip_posMinusbip_pre ...
    sco_contrastMinus_pla_contrast bip_contrastMinus_pla_contrast sco_contrastMinus_bip_contrast
clearvars pla_pre pla_pos sco_pre sco_pos bip_pre bip_pos...
    pla_contrast sco_contrast bip_contrast
cfg = [];
cfg.frequency   = Freqs_in;
cfg.avgoverfreq = 'yes';
cfg.channel     = Channels_in;   
cfg.avgoverchan = 'yes';
cfg.latency     = Times_in; 
cfg.avgovertime = 'yes';
cfg.method      = 'analytic';
cfg.statistic   = 'cohensd';
cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;
pla_posMinuspla_pre = ft_freqstatistics(cfg, COND1.(cond{2}),COND1.(cond{1})); %pla
sco_posMinussco_pre = ft_freqstatistics(cfg, COND3.(cond{2}),COND3.(cond{1})); %sco
bip_posMinusbip_pre = ft_freqstatistics(cfg, COND2.(cond{2}),COND2.(cond{1})); %bip
sco_contrastMinus_pla_contrast = ft_freqstatistics(cfg, COND3.(cond{3}), COND1.(cond{3})); % sco-pla
bip_contrastMinus_pla_contrast = ft_freqstatistics(cfg, COND2.(cond{3}), COND1.(cond{3})); % bip-pla
sco_contrastMinus_bip_contrast = ft_freqstatistics(cfg, COND3.(cond{3}), COND2.(cond{3})); % sco-bip

if isequal(clusterGender, 'Pos')
    Pos.cohensd_rect = struct();
    for i =1:6    
        Pos.cohensd_rect.(contrast_cond{i}) = eval([contrast_cond{i}, '.cohensd']);
    end
else 
    Neg.cohensd_rect = struct();
    for i =1:6    
        Neg.cohensd_rect.(contrast_cond{i}) = eval([contrast_cond{i}, '.cohensd']);
    end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option2: Average Over the cluster (all time-channel pairs in the largest cluster)
% first make the same selection as used in the stat !!!!
clearvars pla_posMinuspla_pre sco_posMinussco_pre bip_posMinusbip_pre ...
    sco_contrastMinus_pla_contrast bip_contrastMinus_pla_contrast sco_contrastMinus_bip_contrast
clearvars pla_pre pla_pos sco_pre sco_pos bip_pre bip_pos...
    pla_contrast sco_contrast bip_contrast
cfg = [];
cfg.frequency   = [min(stat_inuse.freq) max(stat_inuse.freq)];
cfg.avgoverfreq = 'no';
cfg.channel     = stat_inuse.label;   
cfg.avgoverchan = 'no';
cfg.latency     = [min(stat_inuse.time) max(stat_inuse.time)];
cfg.avgoverfreq = 'no';
pla_pre =  ft_selectdata(cfg, COND1.(cond{1}));
pla_pos =  ft_selectdata(cfg, COND1.(cond{2}));
sco_pre =  ft_selectdata(cfg, COND3.(cond{1})); 
sco_pos =  ft_selectdata(cfg, COND3.(cond{2})); 
bip_pre =  ft_selectdata(cfg, COND2.(cond{1})); 
bip_pos =  ft_selectdata(cfg, COND2.(cond{2})); 
pla_contrast = ft_selectdata(cfg, COND1.(cond{3}));
sco_contrast = ft_selectdata(cfg, COND3.(cond{3}));
bip_contrast = ft_selectdata(cfg, COND2.(cond{3}));

subj = 24;
for ii = 1:subj
  % construct a 4-dimensional Boolean array to select the data from this participant
  sel4d = false(size(pla_pre.powspctrm));
  sel4d(ii,:,:,:) = Sigmask; 
  tmp = pla_pre.powspctrm(ii, sel4d(ii,:,:,:));
  pla_pre.avg(ii,1) = mean(tmp);

  % use the same logic for all the other dataset 
  pla_pos.avg(ii,1) = mean(pla_pos.powspctrm(ii, sel4d(ii,:,:,:)) );
  sco_pre.avg(ii,1) = mean(sco_pre.powspctrm(ii, sel4d(ii,:,:,:)) );
  sco_pos.avg(ii,1) = mean(sco_pos.powspctrm(ii, sel4d(ii,:,:,:)) );
  bip_pre.avg(ii,1) = mean(bip_pre.powspctrm(ii, sel4d(ii,:,:,:)) );
  bip_pos.avg(ii,1) = mean(bip_pos.powspctrm(ii, sel4d(ii,:,:,:)) );
  pla_contrast.avg(ii,1) = mean(pla_contrast.powspctrm(ii, sel4d(ii,:,:,:)) );
  sco_contrast.avg(ii,1) = mean(sco_contrast.powspctrm(ii, sel4d(ii,:,:,:)) );
  bip_contrast.avg(ii,1) = mean(bip_contrast.powspctrm(ii, sel4d(ii,:,:,:)) );   

  clear tem  sel4d 
end

pla_posMinuspla_pre = mean(pla_pos.avg-pla_pre.avg)./std((pla_pos.avg-pla_pre.avg), 0,1); %pla
sco_posMinussco_pre = mean(sco_pos.avg-sco_pre.avg)./std((sco_pos.avg-sco_pre.avg), 0,1); %sco
bip_posMinusbip_pre = mean(bip_pos.avg-bip_pre.avg)./std((bip_pos.avg-bip_pre.avg), 0,1); %bip
sco_contrastMinus_pla_contrast = mean(sco_contrast.avg-pla_contrast.avg)./std((sco_contrast.avg-pla_contrast.avg), 0,1); % sco-pla
bip_contrastMinus_pla_contrast = mean(bip_contrast.avg-pla_contrast.avg)./std((bip_contrast.avg-pla_contrast.avg), 0,1); % bip-pla
sco_contrastMinus_bip_contrast = mean(sco_contrast.avg-bip_contrast.avg)./std((sco_contrast.avg-bip_contrast.avg), 0,1); % sco-bip

if isequal(clusterGender, 'Pos')
        Pos.cohensd_jagged = struct();
        for i =1:6    
            Pos.cohensd_jagged.(contrast_cond{i}) = eval([contrast_cond{i}]);
        end
elseif isequal(clusterGender, 'Neg')

        Neg.cohensd_jagged = struct();
        for i =1:6    
            Neg.cohensd_jagged.(contrast_cond{i}) = eval([contrast_cond{i}]);
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option3: Maximum effect size (upper bound)
% Note: this method seems to be the upper bound, but Not Necessarily as
% cohen's d also dependens on the standad deviation across the subjects. 
clearvars pla_posMinuspla_pre sco_posMinussco_pre bip_posMinusbip_pre ...
    sco_contrastMinus_pla_contrast bip_contrastMinus_pla_contrast sco_contrastMinus_bip_contrast
clearvars pla_pre pla_pos sco_pre sco_pos bip_pre bip_pos...
    pla_contrast sco_contrast bip_contrast
cfg = [];
cfg.frequency   = [min(stat_inuse.freq) max(stat_inuse.freq)];
cfg.avgoverfreq = 'no';
cfg.channel     = stat_inuse.label;   
cfg.avgoverchan = 'no';
cfg.latency     = [min(stat_inuse.time) max(stat_inuse.time)];
cfg.avgoverfreq = 'no';
pla_pre =  ft_selectdata(cfg, COND1.(cond{1}));
pla_pos =  ft_selectdata(cfg, COND1.(cond{2}));
sco_pre =  ft_selectdata(cfg, COND3.(cond{1})); 
sco_pos =  ft_selectdata(cfg, COND3.(cond{2})); 
bip_pre =  ft_selectdata(cfg, COND2.(cond{1})); 
bip_pos =  ft_selectdata(cfg, COND2.(cond{2})); 
pla_contrast = ft_selectdata(cfg, COND1.(cond{3}));
sco_contrast = ft_selectdata(cfg, COND3.(cond{3}));
bip_contrast = ft_selectdata(cfg, COND2.(cond{3}));

subj = 24;
for ii = 1:subj
  % construct a 4-dimensional Boolean array to select the data from this participant
  clear  sel4d 
  sel4d = false(size(pla_pre.powspctrm));
  sel4d(ii,:,:,:) = Sigmask;  

  % select data in the cluster for this participant, represent it as a vector
  clearvars a b
  a = pla_pre.powspctrm(ii, sel4d(ii,:,:,:));
  b = pla_pos.powspctrm(ii, sel4d(ii,:,:,:));
  pla_posMinuspla_pre(ii,:) = b - a;

  clearvars a b 
  a = sco_pre.powspctrm(ii, sel4d(ii,:,:,:));
  b = sco_pos.powspctrm(ii, sel4d(ii,:,:,:));
  sco_posMinussco_pre(ii,:) = b-a;

  clearvars a b 
  a = bip_pre.powspctrm(ii, sel4d(ii,:,:,:));
  b = bip_pos.powspctrm(ii, sel4d(ii,:,:,:));
  bip_posMinusbip_pre(ii,:) = b-a;

  clearvars a b 
  a = pla_contrast.powspctrm(ii, sel4d(ii,:,:,:));
  b = sco_contrast.powspctrm(ii, sel4d(ii,:,:,:));
  sco_contrastMinus_pla_contrast(ii,:) = b-a;

  clearvars a b 
  a = pla_contrast.powspctrm(ii, sel4d(ii,:,:,:));
  b = bip_contrast.powspctrm(ii, sel4d(ii,:,:,:));
  bip_contrastMinus_pla_contrast(ii,:) = b-a;

  clearvars a b 
  a = bip_contrast.powspctrm(ii, sel4d(ii,:,:,:));
  b = sco_contrast.powspctrm(ii, sel4d(ii,:,:,:));
  sco_contrastMinus_bip_contrast(ii,:) = b-a;
end

for i = 1:size(pla_posMinuspla_pre,2)
 cohensd.pla_posMinuspla_pre(i) =( mean(pla_posMinuspla_pre(:,i)))./(std(pla_posMinuspla_pre(:,i),0,1));
 cohensd.sco_posMinussco_pre(i) =( mean(sco_posMinussco_pre(:,i)))./(std(sco_posMinussco_pre(:,i),0,1));
 cohensd.bip_posMinusbip_pre(i) =( mean(bip_posMinusbip_pre(:,i)))./(std(bip_posMinusbip_pre(:,i),0,1));
 cohensd.sco_contrastMinus_pla_contrast(i) =( mean(sco_contrastMinus_pla_contrast(:,i)))./...
     (std(sco_contrastMinus_pla_contrast(:,i),0,1));
 cohensd.bip_contrastMinus_pla_contrast(i) =( mean(bip_contrastMinus_pla_contrast(:,i)))./...
     (std(bip_contrastMinus_pla_contrast(:,i),0,1));
 cohensd.sco_contrastMinus_bip_contrast(i) =( mean(sco_contrastMinus_bip_contrast(:,i)))./...
     (std(sco_contrastMinus_bip_contrast(:,i),0,1));
end

[irow,jcol,kdept] = ind2sub(size(Sigmask), find(Sigmask(:) == 1)); %extract the coordinates of certain elements in a logical matrix. x = channel; y = freq; z =time (depth)
if isequal(clusterGender, 'Pos')
    Pos.cohensd_max = struct();
    Pos.cohensd_max_idx = struct();
    for i =1:6    
        [Pos.cohensd_max.(contrast_cond{i}), idx(i)]...
            = max(abs(cohensd.(contrast_cond{i})));
    end
    Pos.cohensd_max_idx = idx;
    % Determine maximum effect size and at which channel and time point Cohen's d is maximal
    % returns the row and column subscripts of each nonzero element in array X
    Pos.cohensd_max_irow = irow;
    Pos.cohensd_max_jcol = jcol;
    Pos.cohensd_max_kdept = kdept;
elseif isequal(clusterGender, 'Neg')
    Neg.cohensd_max = struct();
    Neg.cohensd_max_idx = struct();
    for i =1:6    
        [Neg.cohensd_max.(contrast_cond{i}), idx(i)]...
            = max(abs(cohensd.(contrast_cond{i})));
    end
    Neg.cohensd_max_idx = idx;
    Neg.cohensd_max_irow = irow;
    Neg.cohensd_max_jcol = jcol;
    Neg.cohensd_max_kdept = kdept;
end

% Display results
idx = 22731; %22731; 28755
fprintf('\n')
disp('~~~~~')
disp(['Maximum effect size is at channel: ' stat_inuse.label{irow(idx)} ' at frequency: ' num2str(stat_inuse.freq(jcol(idx))) ...
    ' Hz and at time ' num2str(stat_inuse.time(kdept(idx))) ' sec']);
disp('~~~~~')
fprintf('\n')