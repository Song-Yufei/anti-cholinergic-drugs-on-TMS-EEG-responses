% Function to load data
function [SOUND_SHAM_PRE, SOUND_REAL_PRE, SOUND_SHAM_POS, SOUND_REAL_POS] = load_tmseeg(subj, session_number, TARG)
    data_file = ['...'\CHODME_' subj '_' session_number '_' TARG '_TEPs(PREPOS)_cleaned_ft'];
    load(data_file, 'data_SIR_FT');

    trialinfo = data_SIR_FT.trialinfo;

    % base line correction
    cfg = [];
    cfg.reref = 'yes';
    cfg.refchannel = {'all'};
    cfg.refmethod = 'avg';
    cfg.baselinewindow = [-0.8 -0.01];
    EEG_interp = ft_preprocessing(cfg, data_SIR_FT);

    % Low pass filter for all data
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 45;
    cfg.covariance = 'yes';
    cfg.covariancewindow = 'all';
    cfg.keeptrials = 'no';

    % sort data for each condition
    cfg.trials = find(ismember(trialinfo, 2)); % pre-drug, sham tms
    SOUND_SHAM_PRE = ft_timelockanalysis(cfg, EEG_interp);

    cfg.trials = find(ismember(trialinfo, 1)); % pre-drug, active tms
    SOUND_REAL_PRE = ft_timelockanalysis(cfg, EEG_interp);

    cfg.trials = find(ismember(trialinfo, 3)); % pos-drug, active tms
    SOUND_REAL_POS = ft_timelockanalysis(cfg, EEG_interp);

    cfg.trials = find(ismember(trialinfo, 4)); % pos-drug, sham tms
    SOUND_SHAM_POS = ft_timelockanalysis(cfg, EEG_interp);
end