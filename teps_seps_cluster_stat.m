% this script:
% 1) load individual tms-eeg data after preprocessing
% Note: in SMAData structure, SHAM_...from each subject is the dataset for computing SEPs
% 2) run cluster based permuation F/t tests on teps and seps
% 3) compute effect size
%% Prepare data and save group data
clear;
clc;
% Define subject IDs and drug conditions
subjects = {'002','003', '004', '007', '009', '010', '011', '012', ...
    '013', '014', '015', '017', '018', '019', '020', '021', '022', '024', ...
    '025', '027', '028', '030', '031', '032'};
TARG = 'SMA';
DRUGS = {'Place', 'Scop', 'Bipe'};
SESSION_LABELS = {'PLA', 'SCO', 'BIP'};

% Read session and drug information
[~, sessdrugs] = xlsread('...\SUJsessionInfo2');

% Initialize storage for TEPs data
TEPs_SHAM_PRE = cell(1, length(subjects));
TEPs_ACTIVE_PRE = cell(1, length(subjects));
TEPs_SHAM_POS = cell(1, length(subjects));
TEPs_ACTIVE_POS = cell(1, length(subjects));

% Loop through each subject
for ii = 1:length(subjects)
    subj = subjects{ii};
    subj_rows = ismember(cellfun(@(x) x(5:end), sessdrugs(:, 1), 'UniformOutput', false), subj);

    for sess = 1:3
        % Find the appropriate session cell
        session_cell = find(ismember(cellfun(@(x) x(4:end), sessdrugs(subj_rows, 2:4), 'UniformOutput', false), DRUGS(sess)));
        session_number = cell2mat(cellfun(@(x) x(2), sessdrugs(subj_rows, session_cell + 1), 'UniformOutput', false));

        % Load individual tms-eeg data after preprocessing
        [SOUND_SHAM_PRE, SOUND_REAL_PRE, SOUND_SHAM_POS, SOUND_REAL_POS] = load_tmseeg(subj, session_number, TARG);

        % Store results based on session
        TEPs_SHAM_PRE{ii, sess} = SOUND_SHAM_PRE;
        TEPs_ACTIVE_PRE{ii, sess} = SOUND_REAL_PRE;
        TEPs_SHAM_POS{ii, sess} = SOUND_SHAM_POS;
        TEPs_ACTIVE_POS{ii, sess} = SOUND_REAL_POS;
    end
end

% Combine and save data
SMAData = struct();
for sess = 1:3
    SMAData.(['SHAM_' SESSION_LABELS{sess} '_PRE']) = TEPs_SHAM_PRE(:, sess);
    SMAData.(['ACTIVE_' SESSION_LABELS{sess} '_PRE']) = TEPs_ACTIVE_PRE(:, sess);
    SMAData.(['SHAM_' SESSION_LABELS{sess} '_POS']) = TEPs_SHAM_POS(:, sess);
    SMAData.(['ACTIVE_' SESSION_LABELS{sess} '_POS']) = TEPs_ACTIVE_POS(:, sess);
end

filepath = '...\GroupData\';
filename = 'SMA_TEPs_PREPOS'; 
% save([filepath filename], 'SMAData', '-v7.3');
%% Computes condition specific TEP/SEP average over trials
% Define the configuration for grand averaging
clear
% Load group data
load('...\SMA_TEPs_PREPOS.mat')
Data = SMAData; 
TEPs_SHAM_PLA_PRE = Data.SHAM_PLA_PRE;
TEPs_ACTIVE_PLA_PRE= Data.ACTIVE_PLA_PRE;
TEPs_SHAM_PLA_POS = Data.SHAM_PLA_POS;
TEPs_ACTIVE_PLA_POS = Data.ACTIVE_PLA_POS;

TEPs_SHAM_SCO_PRE = Data.SHAM_SCO_PRE;
TEPs_ACTIVE_SCO_PRE = Data.ACTIVE_SCO_PRE;
TEPs_SHAM_SCO_POS = Data.SHAM_SCO_POS;
TEPs_ACTIVE_SCO_POS = Data.ACTIVE_SCO_POS;

TEPs_SHAM_BIP_PRE = Data.SHAM_BIP_PRE;
TEPs_ACTIVE_BIP_PRE = Data.ACTIVE_BIP_PRE;
TEPs_SHAM_BIP_POS = Data.SHAM_BIP_POS;
TEPs_ACTIVE_BIP_POS = Data.ACTIVE_BIP_POS;

% Define the configuration for grand averaging
% $$$ SHAM and ACTIVE TEPs 
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
% Placebo
GA_TEPs_SHAM_PLA_PRE=ft_timelockgrandaverage(cfg,TEPs_SHAM_PLA_PRE{:}); 
GA_TEPs_ACTIVE_PLA_PRE=ft_timelockgrandaverage(cfg,TEPs_ACTIVE_PLA_PRE{:}); 
GA_TEPs_SHAM_PLA_POS=ft_timelockgrandaverage(cfg,TEPs_SHAM_PLA_POS{:});
GA_TEPs_ACTIVE_PLA_POS=ft_timelockgrandaverage(cfg,TEPs_ACTIVE_PLA_POS{:});

% Scopolamine
GA_TEPs_SHAM_SCO_PRE=ft_timelockgrandaverage(cfg,TEPs_SHAM_SCO_PRE{:});
GA_TEPs_ACTIVE_SCO_PRE=ft_timelockgrandaverage(cfg,TEPs_ACTIVE_SCO_PRE{:});
GA_TEPs_SHAM_SCO_POS=ft_timelockgrandaverage(cfg,TEPs_SHAM_SCO_POS{:}); 
GA_TEPs_ACTIVE_SCO_POS=ft_timelockgrandaverage(cfg,TEPs_ACTIVE_SCO_POS{:}); 

% Biperiden
GA_TEPs_SHAM_BIP_PRE=ft_timelockgrandaverage(cfg,TEPs_SHAM_BIP_PRE{:}); 
GA_TEPs_ACTIVE_BIP_PRE=ft_timelockgrandaverage(cfg,TEPs_ACTIVE_BIP_PRE{:});
GA_TEPs_SHAM_BIP_POS=ft_timelockgrandaverage(cfg,TEPs_SHAM_BIP_POS{:}); 
GA_TEPs_ACTIVE_BIP_POS=ft_timelockgrandaverage(cfg,TEPs_ACTIVE_BIP_POS{:});  
   
% Cleaned TEPs (Subtraction)
cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'individual';
% Placebo 
GA_TEPs_clean_PLA_PRE=ft_math(cfg, GA_TEPs_ACTIVE_PLA_PRE, GA_TEPs_SHAM_PLA_PRE);
GA_TEPs_clean_PLA_POS=ft_math(cfg, GA_TEPs_ACTIVE_PLA_POS, GA_TEPs_SHAM_PLA_POS);
% Scopolamine 
GA_TEPs_clean_SCO_PRE=ft_math(cfg, GA_TEPs_ACTIVE_SCO_PRE, GA_TEPs_SHAM_SCO_PRE);
GA_TEPs_clean_SCO_POS=ft_math(cfg, GA_TEPs_ACTIVE_SCO_POS, GA_TEPs_SHAM_SCO_POS);
% Biperiden 
GA_TEPs_clean_BIP_PRE=ft_math(cfg, GA_TEPs_ACTIVE_BIP_PRE, GA_TEPs_SHAM_BIP_PRE);
GA_TEPs_clean_BIP_POS=ft_math(cfg, GA_TEPs_ACTIVE_BIP_POS, GA_TEPs_SHAM_BIP_POS);

% If the labels are from original NeurOne eeg structure, update
vars = who;
matchingVars = vars(startsWith(vars, 'GA'));
for i = 1:length(matchingVars)
    tempt = [];
    tempt = eval(matchingVars{i});
    if ismember('Afz', tempt.label)
        tempt.label(strcmpi('Afz', tempt.label)) = {'AFz'};
    end

   if ismember('Afz', tempt.elec.label)
    tempt.elec.label(strcmpi('Afz',  tempt.elec.label)) = {'AFz'}; % Afz is not recognized
   end 

    assignin('base', matchingVars{i}, tempt);
end
%% Cluster F test
load('...\neighbour_layout64_custom.mat')
clearvars COND1 COND2 COND3

cfg=[];
cfg.operation =  'subtract';
cfg.parameter =  'individual';
COND1 =ft_math(cfg, GA_TEPs_clean_PLA_POS,GA_TEPs_clean_PLA_PRE);
% COND1 =ft_math(cfg, GA_TEPs_SHAM_PLA_POS,GA_TEPs_SHAM_PLA_PRE); % for SEPs
COND2 =ft_math(cfg, GA_TEPs_clean_BIP_POS, GA_TEPs_clean_BIP_PRE);
COND3 =ft_math(cfg, GA_TEPs_clean_SCO_POS, GA_TEPs_clean_SCO_PRE);

subj=size(COND1.individual,1);
target = 'sma';  

stat = [];
cfg = [];
cfg.latency = [0.02  0.3];
cfg.avgovertime = 'no';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesFunivariate';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.025; % for short-lived, focal effect: 0.025
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbour_layout.custneighbor.neighbours; 
cfg.tail = 1; % for a F-statistic, only right tail 
cfg.correcttail = 'no';
cfg.clustertail = 1;
cfg.alpha = 0.05;
cfg.numrandomization =10000;

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
stat = ft_timelockstatistics(cfg, COND1, COND2, COND3);

%% Follow-up analysis: get Cluster content (from cluster F-test), extent in space and time
% load stat results from cluster F test
stat = ; % stat results from cluster F-test
% positive (only positive from F-test)
posClust=[];
if isfield(stat,'posclusters')
  if ~isempty(stat.posclusters)
  for clusx =1: length(stat.posclusters)
    if stat.posclusters(clusx).prob<0.05
    posClust.selCh{clusx} =stat.label(any(stat.posclusterslabelmat==clusx,2)); % any(A,2) tests elements of each row and returns a column vector of logical 1s and 0s.
    posClust.minT{clusx} = min(stat.time( any(stat.posclusterslabelmat==clusx,1) )); % any(A,1) tests elements of each column and returns a row vector of logical 1s and 0s. 
    posClust.maxT{clusx} = max(stat.time( any(stat.posclusterslabelmat==clusx,1) ));
    end
  end
  end
end 
posClust

%%%% Cluster content
Channels_in = stat.label(any(stat.posclusterslabelmat==1,2)); %==1 1st cluster/largest cluster
Times_in = [min(stat.time( any(stat.posclusterslabelmat==1,1) ))...
    max(stat.time( any(stat.posclusterslabelmat==1,1) ))]; % ==1 1st cluster/largest cluster
%% Follow-up analysis: pairwise Cluster t-test
load('...\neighbour_layout64_custom.mat')
clearvars COND1 COND2

% contrastsI example:
COND1=GA_TEPs_clean_SCO_PRE;
COND2=GA_TEPs_clean_SCO_POS;

% contrastsII example:
cfg=[];
cfg.operation =  'subtract';
cfg.parameter =  'individual';
COND1 =ft_math(cfg, GA_TEPs_clean_PLA_POS,GA_TEPs_clean_PLA_PRE);
COND2 =ft_math(cfg, GA_TEPs_clean_SCO_POS, GA_TEPs_clean_SCO_PRE);

subj=size(COND1.individual,1);

    stat = [];
    cfg = [];
    cfg.latency = Times_in; 
    cfg.avgovertime = 'no'; 
    cfg.channel = Channels_in; 
    cfg.avgoverchan = 'no';
    cfg.method = 'montecarlo'; 
    cfg.statistic = 'ft_statfun_depsamplesT';
    cfg.correctm = 'cluster';
    cfg.clusteralpha = 0.025;% short-lived, focal effect
    cfg.clustertail = 0;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2; 
    cfg.neighbours = neighbour_layout.custneighbor.neighbours; 
    cfg.tail = 0; 
    cfg.correcttail = 'alpha';  
    cfg.alpha = 0.05/3; % corrected by the number of simple effects tests, eg., scopre vs scopos; bippre vs bippos; plapre vs plapos
    cfg.numrandomization = 10000;

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
    stat = ft_timelockstatistics(cfg, COND2, COND1);
%% Effect size following cluster t-test
clearvars COND1 COND2 COND3 COND4
%%%% Step1: sort data according experimental conditions %%%%
cond = {'PRE','POS','POSminusPRE','SCOminusPLA','BIPminusPLA','SCOminusBIP'};
condName = {'PRE','POS','POS-PRE'};

cfg=[];
cfg.operation =  'subtract';
cfg.parameter =  'individual';
COND1.(cond{3}) =ft_math(cfg, GA_TEPs_clean_PLA_POS,GA_TEPs_clean_PLA_PRE);
COND1.(cond{2}) = GA_TEPs_clean_PLA_POS;
COND1.(cond{1}) = GA_TEPs_clean_PLA_PRE;

cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'individual';
COND2.(cond{3}) =ft_math(cfg, GA_TEPs_clean_BIP_POS, GA_TEPs_clean_BIP_PRE);
COND2.(cond{2}) = GA_TEPs_clean_BIP_POS;
COND2.(cond{1}) = GA_TEPs_clean_BIP_PRE;

cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'individual';
COND3.(cond{3}) =ft_math(cfg, GA_TEPs_clean_SCO_POS, GA_TEPs_clean_SCO_PRE);
COND3.(cond{2}) = GA_TEPs_clean_SCO_POS;
COND3.(cond{1}) = GA_TEPs_clean_SCO_PRE;

cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'individual';
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
clearvars Times_in Channels_in stat Sigmask

% load stat result from followup t-tests 
stat = ; % stat results from t-tests
% Two Info about LOADED stat: 
% 1) Know pos or neg; 2) Know crtical alpha level 
clusterGender = 'Pos'; % 'Pos' or 'Neg'
alpha_lim = 0.025/3;

if isequal(clusterGender, 'Pos')
    posClust=[];
    if isfield(stat,'posclusters')
      if ~isempty(stat.posclusters)
      for clusx =1: length(stat.posclusters)
        if stat.posclusters(clusx).prob<alpha_lim
        posClust.selCh{clusx} =stat.label(any(stat.posclusterslabelmat==clusx,2)); 
        posClust.minT{clusx} = min(stat.time( any(stat.posclusterslabelmat==clusx,1) )); 
        posClust.maxT{clusx} = max(stat.time( any(stat.posclusterslabelmat==clusx,1) ));
        end
      end
      end
    end 
    posClust

    Channels_in = stat.label(any(stat.posclusterslabelmat==1,2));
    Times_in = [min(stat.time( any(stat.posclusterslabelmat==1,1) ))...
        max(stat.time( any(stat.posclusterslabelmat==1,1) ))]; % ==1 1st cluster/largest cluster
 Sigmask = (stat.posclusterslabelmat==1);

elseif isequal(clusterGender, 'Neg')
    
    negClust=[];
    if isfield(stat,'negclusters')
        if ~isempty(stat.negclusters)
        for clusx =1: length(stat.negclusters)
            if stat.negclusters(clusx).prob<alpha_lim
            negClust.selCh{clusx} =stat.label(any(stat.negclusterslabelmat==clusx,2)); 
            negClust.minT{clusx} = min(stat.time( any(stat.negclusterslabelmat==clusx,1) ));
            negClust.maxT{clusx} = max(stat.time( any(stat.negclusterslabelmat==clusx,1) ));
            end
        end
        end
    end
    negClust
    
    Channels_in = stat.label(any(stat.negclusterslabelmat==1,2));
    Times_in = [min(stat.time( any(stat.negclusterslabelmat==1,1) ))...
        max(stat.time( any(stat.negclusterslabelmat==1,1) ))]; % ==1 1st cluster/largest cluster
    Sigmask = (stat.negclusterslabelmat==1);
end 
clearvars -except Times_in Channels_in stat Sigmask...
    COND1 COND2 COND3 COND4 design cond condName clusterGender

%%%% Step3: Select data %%%%
% 1) keep data structure's time and channels content the same as stat's time and channels content)
% 2) Know which contrast corresponds to the LOADED stat in the begining (Here I listed all possible contrasts,but only report the one intended)
contrast_cond = {'pla_posMinuspla_pre', 'sco_posMinussco_pre', 'bip_posMinusbip_pre'...
    'sco_contrastMinus_pla_contrast', 'bip_contrastMinus_pla_contrast', 'sco_contrastMinus_bip_contrast'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option1: Average Over Circumscribed rectangle (lower bound)
clearvars pla_posMinuspla_pre sco_posMinussco_pre bip_posMinusbip_pre ...
    sco_contrastMinus_pla_contrast bip_contrastMinus_pla_contrast sco_contrastMinus_bip_contrast
clearvars pla_pre pla_pos sco_pre sco_pos bip_pre bip_pos...
    pla_contrast sco_contrast bip_contrast
cfg = [];
cfg.latency   = Times_in;
cfg.avgovertime = 'yes';
cfg.channel     = Channels_in;   
cfg.avgoverchan = 'yes';
cfg.method      = 'analytic';
cfg.statistic   = 'cohensd';
cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;
pla_posMinuspla_pre = ft_timelockstatistics(cfg, COND1.(cond{2}),COND1.(cond{1})); %pla
sco_posMinussco_pre = ft_timelockstatistics(cfg, COND3.(cond{2}),COND3.(cond{1})); %sco
bip_posMinusbip_pre = ft_timelockstatistics(cfg, COND2.(cond{2}),COND2.(cond{1})); %bip
sco_contrastMinus_pla_contrast = ft_timelockstatistics(cfg, COND3.(cond{3}), COND1.(cond{3})); % sco-pla
bip_contrastMinus_pla_contrast = ft_timelockstatistics(cfg, COND2.(cond{3}), COND1.(cond{3})); % bip-pla
sco_contrastMinus_bip_contrast = ft_timelockstatistics(cfg, COND3.(cond{3}), COND2.(cond{3})); % sco-bip

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
cfg.latency   = [min(stat.time) max(stat.time)];
cfg.avgovertime = 'no';
cfg.channel     = stat.label;   
cfg.avgoverchan = 'no';
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
  % construct a 3-dimensional Boolean array to select the data from this participant
  sel3d = false(size(pla_pre.individual));
  sel3d(ii,:,:) = Sigmask; 
  tmp = pla_pre.individual(ii, sel3d(ii, :, :));
  pla_pre.avg(ii,1) = mean(tmp);

  % use the same logic for all the other dataset 
  pla_pos.avg(ii,1) = mean(pla_pos.individual(ii, sel3d(ii, :, :)) );
  sco_pre.avg(ii,1) = mean(sco_pre.individual(ii, sel3d(ii, :, :)) );
  sco_pos.avg(ii,1) = mean(sco_pos.individual(ii, sel3d(ii, :, :)) );
  bip_pre.avg(ii,1) = mean(bip_pre.individual(ii, sel3d(ii, :, :)) );
  bip_pos.avg(ii,1) = mean(bip_pos.individual(ii, sel3d(ii, :, :)) );
  pla_contrast.avg(ii,1) = mean(pla_contrast.individual(ii, sel3d(ii, :, :)) );
  sco_contrast.avg(ii,1) = mean(sco_contrast.individual(ii, sel3d(ii, :, :)) );
  bip_contrast.avg(ii,1) = mean(bip_contrast.individual(ii, sel3d(ii, :, :)) );   

  clear tem  sel3d 
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
cfg.latency   =  [min(stat.time) max(stat.time)];
cfg.avgovertime = 'no';
cfg.channel     = stat.label;   
cfg.avgoverchan = 'no';
pla_pre =  ft_selectdata(cfg, COND1.(cond{1}));
pla_pos =  ft_selectdata(cfg, COND1.(cond{2}));
sco_pre =  ft_selectdata(cfg, COND3.(cond{1})); % SCO
sco_pos =  ft_selectdata(cfg, COND3.(cond{2})); % SCO
bip_pre =  ft_selectdata(cfg, COND2.(cond{1})); % BIP 
bip_pos =  ft_selectdata(cfg, COND2.(cond{2})); % BIP 
pla_contrast = ft_selectdata(cfg, COND1.(cond{3}));
sco_contrast = ft_selectdata(cfg, COND3.(cond{3}));
bip_contrast = ft_selectdata(cfg, COND2.(cond{3}));

subj = 24;
for ii = 1:subj
  % construct a 3-dimensional Boolean array to select the data from this participant
  clear  sel3d 
  sel3d = false(size(pla_pre.individual));
  sel3d(ii,:,:) = Sigmask;  

  % select data in the cluster for this participant, represent it as a vector
  clearvars a b
  a = pla_pre.individual(ii, sel3d(ii, :, :));
  b = pla_pos.individual(ii, sel3d(ii, :, :));
  pla_posMinuspla_pre(ii,:) = b-a;

  clearvars a b 
  a = sco_pre.individual(ii, sel3d(ii, :, :));
  b = sco_pos.individual(ii, sel3d(ii, :, :));
  sco_posMinussco_pre(ii,:) = b-a;

  clearvars a b 
  a = bip_pre.individual(ii, sel3d(ii, :, :));
  b = bip_pos.individual(ii, sel3d(ii, :, :));
  bip_posMinusbip_pre(ii,:) = b-a;

  clearvars a b 
  a = pla_contrast.individual(ii, sel3d(ii, :, :));
  b = sco_contrast.individual(ii, sel3d(ii, :, :));
  sco_contrastMinus_pla_contrast(ii,:) = b-a;

  clearvars a b 
  a = pla_contrast.individual(ii, sel3d(ii, :, :));
  b = bip_contrast.individual(ii, sel3d(ii, :, :));
  bip_contrastMinus_pla_contrast(ii,:) = b-a;

  clearvars a b 
  a = bip_contrast.individual(ii, sel3d(ii, :, :));
  b = sco_contrast.individual(ii, sel3d(ii, :, :));
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

[irow,jcol] = ind2sub(size(Sigmask), find(Sigmask(:) == 1));
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
end
