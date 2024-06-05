% this script: 
% 1) plot grand average TFRs across all conditions
% Note: in SMAData structure, SHAM_...from each subject is the dataset for computing SIOs
% 2) illustrate topography of TFRs at selected (or average over) freq and time of interest
% 3) illustrate cluster F-test results; illustrate cluster t-test results and at selected channel and timing
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
   
% clean
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

%    if ismember('Afz', tempt.elec.label)
%     tempt.elec.label(strcmpi('Afz',  tempt.elec.label)) = {'AFz'}; % Afz is not recognized
%    end 

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

%% Plot grand average TFRs
clear COND1 COND2 COND3 COND4
%%%% Step1: compute data of interest %%%%
load('...\original_chanlocs.mat');
cond = {'PRE','POS','POSminusPRE','SCOminusPLA','BIPminusPLA','SCOminusBIP'};
condName = {'PRE','POS','POS-PRE'};

cfg=[];
cfg.operation =  'subtract';
cfg.parameter =  'powspctrm';
COND1.(cond{3}) =ft_math(cfg, GA_TFRs_SHAM_PLA_POS,GA_TFRs_SHAM_PLA_PRE);
COND1.(cond{2}) = GA_TFRs_SHAM_PLA_POS;
COND1.(cond{1}) = GA_TFRs_SHAM_PLA_PRE;

COND2.(cond{3}) =ft_math(cfg, GA_TFRs_SHAM_BIP_POS, GA_TFRs_SHAM_BIP_PRE);
COND2.(cond{2}) = GA_TFRs_SHAM_BIP_POS;
COND2.(cond{1}) = GA_TFRs_SHAM_BIP_PRE;

COND3.(cond{3}) =ft_math(cfg, GA_TFRs_SHAM_SCO_POS, GA_TFRs_SHAM_SCO_PRE);
COND3.(cond{2}) = GA_TFRs_SHAM_SCO_POS;
COND3.(cond{1}) = GA_TFRs_SHAM_SCO_PRE;

COND4.(cond{4}) =ft_math(cfg, COND3.(cond{3}),COND1.(cond{3})); % Time:sco-pla
COND4.(cond{5}) =ft_math(cfg, COND2.(cond{3}),COND1.(cond{3})); % Time:bip-pla
COND4.(cond{6}) =ft_math(cfg, COND3.(cond{3}),COND2.(cond{3})); % Time:sco-pla

%%%% Step2: layout %%%%
fig = figure('color','w');
set(gcf,'position',[0.1 0.02 1.4 0.6]*1000);
c = get(0, 'DefaultAxesColorOrder');
tiledlayout(fig,3,3,"TileSpacing","tight","TileIndexing","columnmajor")

chan2plot = 'FC1'; %'FC1';'P3';'Fz'
chpIndx = find(strcmpi(chan2plot, {original_chanlocs.labels}));
times = COND1.(cond{1}).time; 
freqs = COND1.(cond{1}).freq;
% Draw a black box to mask the time range -0.05s to 0.05s
timeRangeStart = -0.05;
timeRangeEnd = 0.05;

% 1st column: PLA
for i =1:3
ax(i) = nexttile;
clear meanData
% meanData = squeeze(mean(mean(COND1.(cond{i}).powspctrm(:,:,:,:),2),1)); 
meanData = squeeze(mean(mean(COND1.(cond{i}).powspctrm(:,chpIndx,:,:),2),1)); 
contourf(times,freqs,meanData,40,'linecolor','none');
ax(i).CLim = [-0.5 0.5];
ax(i).YLim = [4 40];
ax(i).YDir = 'normal';
ax(i).XLim = [-0.2 0.65];
ax(i).XLabel.String = 'Times (ms)';
ax(i).YLabel.String = 'Frequency (Hz)';
ax(i).YAxis.TickDirection = 'out';
ax(i).YAxis.TickLength = [0.005 0.02];
ax(i).FontSize =14;
ax(i).Box  = "off";
ax(i).LineWidth = 1;
ax(1).Title.String = 'PLA';
boxX = [timeRangeStart, timeRangeEnd, timeRangeEnd, timeRangeStart, timeRangeStart];
boxY = [ax(i).YLim(1), ax(i).YLim(1), ax(i).YLim(2), ax(i).YLim(2), ax(i).YLim(1)];
rectangle(ax(i), 'Position', [timeRangeStart, ax(1).YLim(1), timeRangeEnd - timeRangeStart, ax(1).YLim(2) - ax(1).YLim(1)], 'FaceColor', 'k', 'EdgeColor', 'none');
end

[map, ~, ~] = colorcet('D1A'); %'D1A' 'L20'
% map = brewermap([], '-RdYlBu');
colormap(map)
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Ticks = [min(ax(1).CLim);0;max(ax(1).CLim)];
    cb.LineWidth = 0.5;
    cb.FontSize = 14;
    title(cb,'z','fontsize',14);

% 2nd column: BIP
for i =1:3
ax(i+3) = nexttile;
clear meanData
% meanData = squeeze(mean(mean(COND2.(cond{i}).powspctrm(:,:,:,:),2),1)); 
meanData = squeeze(mean(mean(COND2.(cond{i}).powspctrm(:,chpIndx,:,:),2),1)); 
contourf(times,freqs,meanData,40,'linecolor','none');
ax(i+3).CLim = [-0.5 0.5];
ax(i+3).YLim = [4 40];
ax(i+3).YDir = 'normal';
ax(i+3).XLim = [-0.2 0.65];
ax(i+3).XLabel.String = 'Times (ms)';
ax(i+3).YLabel.String = 'Frequency (Hz)';
ax(i+3).YAxis.TickDirection = 'out';
ax(i+3).YAxis.TickLength = [0.005 0.02];
ax(i+3).FontSize =14;
ax(i+3).Box  = "off";
ax(i+3).LineWidth = 1;
ax(4).Title.String = 'BIP';
boxX = [timeRangeStart, timeRangeEnd, timeRangeEnd, timeRangeStart, timeRangeStart];
boxY = [ax(i+3).YLim(1), ax(i+3).YLim(1), ax(i+3).YLim(2), ax(i+3).YLim(2), ax(i+3).YLim(1)];
rectangle(ax(i+3), 'Position', [timeRangeStart, ax(1).YLim(1), timeRangeEnd - timeRangeStart, ax(1).YLim(2) - ax(1).YLim(1)], 'FaceColor', 'k', 'EdgeColor', 'none');
end

% 3rd column: SCO
for i =1:3
ax(i+6) = nexttile;
clear meanData
% meanData = squeeze(mean(mean(COND3.(cond{i}).powspctrm(:,:,:,:),2),1)); 
meanData = squeeze(mean(mean(COND3.(cond{i}).powspctrm(:,chpIndx,:,:),2),1)); 
contourf(times,freqs,meanData,40,'linecolor','none');
ax(i+6).CLim = [-0.5 0.5];
ax(i+6).YLim = [4 40];
ax(i+6).YDir = 'normal';
ax(i+6).XLim = [-0.2 0.65];
ax(i+6).XLabel.String = 'Times (ms)';
ax(i+6).YLabel.String = 'Frequency (Hz)';
ax(i+6).YAxis.TickDirection = 'out';
ax(i+6).YAxis.TickLength = [0.005 0.02];
ax(i+6).FontSize =14;
ax(i+6).Box  = "off";
ax(i+6).LineWidth = 1;
ax(7).Title.String = 'SCO';
boxX = [timeRangeStart, timeRangeEnd, timeRangeEnd, timeRangeStart, timeRangeStart];
boxY = [ax(i+6).YLim(1), ax(i+6).YLim(1), ax(i+6).YLim(2), ax(i+6).YLim(2), ax(i+6).YLim(1)];
rectangle(ax(i+6), 'Position', [timeRangeStart, ax(1).YLim(1), timeRangeEnd - timeRangeStart, ax(1).YLim(2) - ax(1).YLim(1)], 'FaceColor', 'k', 'EdgeColor', 'none');
end

%% Plot TFR in space (avg over selected frequency band and time window)
load('...\neighbour_layout64_custom.mat')
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

cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'powspctrm';
COND2.(cond{3}) =ft_math(cfg, GA_TFRs_SHAM_BIP_POS, GA_TFRs_SHAM_BIP_PRE);
COND2.(cond{2}) = GA_TFRs_SHAM_BIP_POS;
COND2.(cond{1}) = GA_TFRs_SHAM_BIP_PRE;

cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'powspctrm';
COND3.(cond{3}) =ft_math(cfg, GA_TFRs_SHAM_SCO_POS, GA_TFRs_SHAM_SCO_PRE);
COND3.(cond{2}) = GA_TFRs_SHAM_SCO_POS;
COND3.(cond{1}) = GA_TFRs_SHAM_SCO_PRE;

cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'powspctrm';
COND4.(cond{4}) =ft_math(cfg, COND3.(cond{3}),COND1.(cond{3})); % Time:sco-pla
COND4.(cond{5}) =ft_math(cfg, COND2.(cond{3}),COND1.(cond{3})); % Time:bip-pla
COND4.(cond{6}) =ft_math(cfg, COND3.(cond{3}),COND2.(cond{3})); % Time:sco-pla

%%%% Step2 plot %%%%
fig = figure('color','w');
set(gcf,'position',[0.1 0.02 0.6 0.6]*1000);
tiledlayout(fig,3,3,"TileSpacing","tight","TileIndexing","columnmajor")

% pre-assign some parameters
freq2plot = 19; % in Hz: 6 or [4,12]
time2plot = 0.44; % in sec: 0.44 or [0.25, 0.65]
plotZvalue = [-0.15,0.15];
clearvars plotCond plotStruc
plotStruc.time = 1;
plotStruc.dimord = 'chan';
plotStruc.label = COND1.(cond{1}).label;
fqIndx = dsearchn(COND1.(cond{1}).freq', freq2plot'); 
tpIndx = dsearchn(COND1.(cond{1}).time', time2plot'); 

% 1st column
clearvars plotCond ax
plotCond = COND2;
for fi =1:3
clear meanData plotStruc.avg
meanData = squeeze(mean(plotCond.(cond{fi}).powspctrm,1)); 
if size(freq2plot,2)==1 && size(time2plot,2)==1
    plotStruc.avg = mean(meanData(:,fqIndx(1), tpIndx(1)),[2 3]);% specific freq and time
else
    plotStruc.avg = mean(meanData(:,fqIndx(1):fqIndx(2), tpIndx(1):tpIndx(2)),[2 3]); %freq band and time window
end
ax(fi) = nexttile;
    cfg = [];
    cfg.layout = neighbour_layout.layout;
    cfg.comment = 'no';
    cfg.interactive = 'no';
    cfg.markersymbol = '.';
    cfg.zlim = plotZvalue;
    cfg.figure      = 'gca';
    ax(1).Title.String = ['PLA ' num2str(freq2plot) 'Hz' ];
    ax(1).FontSize =14;
    ft_topoplotER(cfg,plotStruc);
end       


% 2nd column
clearvars plotCond ax
plotCond = COND2;
for fi =1:3
clear meanData plotStruc.avg
meanData = squeeze(mean(plotCond.(cond{fi}).powspctrm,1)); 
if size(freq2plot,2)==1 && size(time2plot,2)==1
    plotStruc.avg = mean(meanData(:,fqIndx(1), tpIndx(1)),[2 3]);% specific freq and time
else
    plotStruc.avg = mean(meanData(:,fqIndx(1):fqIndx(2), tpIndx(1):tpIndx(2)),[2 3]);%freq band and time window
end
ax(fi) = nexttile;
    cfg = [];
    cfg.layout = neighbour_layout.layout;
    cfg.comment = 'no';
    cfg.interactive = 'no';
    cfg.markersymbol = '.';
    cfg.zlim = plotZvalue;
    cfg.figure      = 'gca';
    ax(1).Title.String = ['BIP ' num2str(freq2plot) 'Hz' ];
    ax(1).FontSize =14;
    ft_topoplotER(cfg,plotStruc);
end 

% 3rd column
clearvars plotCond ax
plotCond = COND3;
for fi =1:3
clear meanData plotStruc.avg
meanData = squeeze(mean(plotCond.(cond{fi}).powspctrm,1)); 
if size(freq2plot,2)==1 && size(time2plot,2)==1
    plotStruc.avg = mean(meanData(:,fqIndx(1), tpIndx(1)),[2 3]);% specific freq and time
else
    plotStruc.avg = mean(meanData(:,fqIndx(1):fqIndx(2), tpIndx(1):tpIndx(2)),[2 3]);%freq band and time window
end
ax(fi) = nexttile;
    cfg = [];
    cfg.layout = neighbour_layout.layout;
    cfg.comment = 'no';
    cfg.interactive = 'no';
    cfg.markersymbol = '.';
    cfg.zlim = plotZvalue;
    cfg.figure      = 'gca';
    ax(1).Title.String = ['SCO ' num2str(freq2plot) 'Hz' ];
    ax(1).FontSize =14;
    ft_topoplotER(cfg,plotStruc);
end       

    map = brewermap([], '-RdYlBu');
    colormap(map)
    cb = colorbar;
    cb.Layout.Tile = 'south';
    cb.Ticks = [min(plotZvalue);0;max(plotZvalue)];
    cb.LineWidth = 0.5;
    cb.FontSize = 14;
    title(cb,'z','fontsize',14);

%% Plot F/t values from cluster F-test/t-tests
clear 
load('...\neighbour_layout64_custom.mat')
%%%% Step1 %%%%
% One Info from LOADED stat:
%  which value ( 'F-value' or 't-value')
WhichValue = 't-value'; 

stat_select = ; % stat results from cluster t-tests
stat_brod = ; % stat results from cluster F-tests

if isequal(WhichValue, 'F-value')
    stat_inuse = stat_brod;
    alpha_lim = 0.05;
else
    stat_inuse = stat_select;
    alpha_lim = 0.025/3; % neg or pos, followed by 3 pairwise tests
end

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

%%%% Step2 F/t values in channel-time_frquency axis %%%%
% layout
fig = figure('color','w');
set(gcf,'position',[0.1 0.02 1.2 0.6]*1000);
tiledlayout(fig,1,2,"TileSpacing","tight")

clearvars plot_clus plot_stat
plot_stat = stat_inuse;
plot_clus = zeros(size(stat_inuse.prob));
if isempty(posClust) && ~isempty(negClust)
     plot_clus(stat_inuse.negclusterslabelmat==1) =  -1; % negative cluster

elseif isempty(negClust) && ~isempty(posClust)
     plot_clus(stat_inuse.posclusterslabelmat==1) =  1; % positive cluster

elseif ~isempty(negClust) && ~isempty(posClust)
     plot_clus(stat_inuse.posclusterslabelmat==1) =  1; % positive cluster
     plot_clus(stat_inuse.negclusterslabelmat==1) =  -1; % negative cluster
end
plot_stat.stat = plot_clus;

if ~isequal(WhichValue,  'F-value') % t-value       
    % plot statistc value
        ax(1) = nexttile;
        cfg = [];
        cfg.marker  = 'on';
        cfg.layout  = neighbour_layout.layout;
        cfg.parameter  = 'stat';
        cfg.interactive ='yes';
        cfg.colormap = 'coolwarm'; % check ft_colormap
        cfg.zlim   = [-4,4];
        cfg.box  = 'yes';
        cfg.showscale    = 'no' ;
        cfg.comment   ='no';
        cfg.showcomment  = 'no' ;
        cfg.figure = gca;
        ft_multiplotTFR(cfg, stat_inuse);
            cb = colorbar;
            cb.Ticks = [-4;0;4];
            cb.LineWidth = 0.5;
            cb.FontSize = 14;
            cb.Title.String = 't-values';    
            cb.Title.FontSize = 14;
            cb.Location ="northoutside";
    % plot cluster extent
        ax(2) = nexttile;
        cfg = [];
        cfg.marker  = 'on';
        cfg.layout  = neighbour_layout.layout;
        cfg.parameter  = 'stat';
        cfg.interactive ='yes';
        cfg.maskparameter = 'mask';  
        cfg.colormap = 'coolwarm';
        cfg.zlim   = [-1,1];
        cfg.maskstyle = 'opacity';
        cfg.box  = 'yes';
        cfg.showscale    = 'no' ;
        cfg.comment   ='no';
        cfg.showcomment  = 'no' ;
        cfg.figure = gca;
        ft_multiplotTFR(cfg, plot_stat);
else % F-value  
    % plot statistc value
        ax(1) = nexttile;
        cfg = [];
        cfg.marker  = 'on';
        cfg.layout  = neighbour_layout.layout;
        cfg.parameter  = 'stat';
        cfg.interactive ='yes';
        cfg.colormap = 'OrRd';
        cfg.zlim   = [0,12];
        cfg.box  = 'yes';
        cfg.showscale    = 'no' ;
        cfg.comment   ='no';
        cfg.showcomment  = 'no' ;
        cfg.figure = gca;
        ft_multiplotTFR(cfg, stat_inuse);
            cb = colorbar;
            cb.Ticks = [0;6;12];
            cb.LineWidth = 0.5;
            cb.FontSize = 14;
            cb.Title.String = 'F-values';    
            cb.Title.FontSize = 14;
            cb.Location ="northoutside";
    % plot cluster extent
        ax(2) = nexttile;
        cfg = [];
        cfg.marker  = 'on';
        cfg.layout  = neighbour_layout.layout;
        cfg.parameter  = 'stat';
        cfg.interactive ='yes';
        cfg.maskparameter = 'mask';  
        cfg.colormap = 'OrRd';
        cfg.zlim   = [0,1];
        cfg.maskstyle = 'opacity';
        cfg.box  = 'yes';
        cfg.showscale    = 'no' ;
        cfg.comment   ='no';
        cfg.showcomment  = 'no' ;
        cfg.figure = gca;
        ft_multiplotTFR(cfg, plot_stat);
end

%%%% Step3 F/t values in space (topograhy) %%%%
clear valMap
if ~isequal(WhichValue,  'F-value') % t-value
            if size(stat_inuse.prob,2) >1 % stat: avgoverfreq = no
                % Initialize matrix valMap with 0
                valMap = zeros(size(stat_brod.prob));                
                % Create a mapping between Channels, Freqs and Times in stat_select.stat and stat_brod.stat
                [~, channel_indices_select_in_brod] = ismember(stat_inuse.label, stat_brod.label);
                [~, freq_indices_select_in_brod] = ismember(stat_inuse.freq, stat_brod.freq);
                [~, time_indices_select_in_brod] = ismember(stat_inuse.time, stat_brod.time);      
                for i = 1:size(stat_inuse.label,1)
                    for j = 1:size(stat_inuse.freq,2)
                        for k = 1: size(stat_inuse.time,2)
                        channel_idx = channel_indices_select_in_brod(i);
                        freq_idx = freq_indices_select_in_brod(j);
                        time_idx = time_indices_select_in_brod(k);
                        valMap(channel_idx, freq_idx,time_idx) = stat_inuse.stat(i,j,k); % channel xfreqx time
                        end
                    end
                end

            else % stat: avgoverfreq = yes
                valMap = zeros(size(stat_brod.prob,1),1);
                [~, channel_idx] = ismember(stat_inuse.label, stat_brod.label);
                valMap(channel_idx) =stat_inuse.stat;
            end

else  % stat: F test
       valMap = stat_inuse.stat;

end 

% layout
fig = figure('color','w');
set(gcf,'position',[0.1 0.01 0.3 0.3]*1000);
c = get(0, 'DefaultAxesColorOrder');
tiledlayout(fig,1,1)

load('...\neighbour_layout64_custom.mat')

if size(stat_inuse.prob,2) >1 % stat: avgoverfreq = no
        freq2plot = 6; %in Hz
        time2plot = 0.44; %in sec
        plotStruc = [];
        plotStruc.time = 1;
        plotStruc.dimord = 'chan';
        plotStruc.label = stat_brod.label;
        [~,fqIndx_in_select] = min(abs(stat_inuse.freq-freq2plot));
        [~,fqIndx_in_brod] = min(abs(stat_brod.freq-freq2plot));
        [~,tpIndx_in_select] = min(abs(stat_inuse.time-time2plot));
        [~,tpIndx_in_brod] = min(abs(stat_brod.time-time2plot));
        plotStruc.avg = valMap(:,fqIndx_in_brod,tpIndx_in_brod);

else
        plotStruc = [];
        plotStruc.time = 1;
        plotStruc.dimord = 'chan';
        plotStruc.label = stat_brod.label;
        fqIndx_in_select = 1;
        plotStruc.avg = valMap;
end

if ~isequal(WhichValue,  'F-value')
    cfg = [];
    cfg.layout = neighbour_layout.layout;
    cfg.comment = 'no';
    cfg.interactive = 'no';
    cfg.markersymbol = '.';
    cfg.zlim = [-4,4];
    cfg.figure      = 'gca';
    cfg.interpolation      = 'nearest'; %'linear','cubic','nearest','v4' 
    cfg.style = 'straight';

    if isfield(stat_inuse,'posclusters')
        if ~isempty(stat_inuse.posclusters)
            if stat_inuse.posclusters(1).prob < alpha_lim % 
                cfg.highlight = 'on';
                cfg.highlightchannel   =  stat_inuse.label(stat_inuse.posclusterslabelmat(:,fqIndx_in_select,tpIndx_in_select)==1); % ==1 the 1st cluster
                cfg.highlightsymbol    = '.';
                cfg.highlightsize      = 16;
                cfg.highlightcolor   = [0 0 0];
            end
        end
    end
    if isfield(stat_inuse,'negclusters')
        if ~isempty(stat_inuse.negclusters)
            if stat_inuse.negclusters(1).prob < alpha_lim 
               
                if isfield(cfg,'highlightchannel')
                    if ~isempty(cfg.highlightchannel) % positive + negative channels highlighted
                        cfg.highlight = {'on','on'};
                        cfg.highlightchannel   =  {cfg.highlightchannel,stat_inuse.label(stat_inuse.negclusterslabelmat(:,fqIndx_in_select,tpIndx_in_select)==1)};% ==1 the 1st cluster
                        cfg.highlightsymbol    = {'.','.'};
                        cfg.highlightsize      = {16,16};
                        cfg.highlightcolor   = {[0 0 0],[0.99 0.99 0.99]};
                    end
                else
                    cfg.highlight = 'on';            % only neative channels highlighted
                    cfg.highlightchannel   =  stat_inuse.label(stat_inuse.negclusterslabelmat(:,fqIndx_in_select,tpIndx_in_select)==1);% ==1 the 1st cluster
                    cfg.highlightsymbol    = '.';
                    cfg.highlightsize      = 16;
                    cfg.highlightcolor   = [ 0.99 0.99 0.99];
                end
            end
        end
    end

    ax = nexttile;
    ft_topoplotER(cfg,plotStruc);
    ax.Title.String = [num2str(freq2plot) 'Hz ' num2str(time2plot) 'sec'];
    ax.FontSize =14;
%     [map, ~, ~] = colorcet('D1A');
    map = brewermap([], '-RdBu');
    colormap(map)
    cb = colorbar;
    cb.Ticks = [-4;0;4];
    cb.LineWidth = 0.5;
    cb.Title.String = 't-values';
    cb.Title.FontSize = 14;  
    cb.Layout.Tile = 'south';

else  % F-value

    cfg = [];
    cfg.layout = neighbour_layout.layout;
    cfg.comment = 'no';
    cfg.interactive = 'no';
    cfg.markersymbol = '.';
    cfg.zlim = [0,12];
    cfg.figure      = 'gca';
    cfg.interpolation      = 'nearest'; %'linear','cubic','nearest','v4' 
    cfg.style = 'straight';

    if isfield(stat_inuse,'posclusters')
        if ~isempty(stat_inuse.posclusters)
            if stat_inuse.posclusters(1).prob < alpha_lim 
                cfg.highlight = 'on';
                cfg.highlightchannel   =  stat_inuse.label(stat_inuse.posclusterslabelmat(:,fqIndx_in_select,tpIndx_in_select)==1); % ==1 the 1st cluster; ==2 2nd cluster
                cfg.highlightsymbol    = '.';
                cfg.highlightsize      = 16;
                cfg.highlightcolor   = [0 0 0];
            end
        end
    end
    % F-value only has positive
    ax = nexttile;
    ft_topoplotER(cfg,plotStruc);
    ax.Title.String = [num2str(freq2plot) 'Hz'];
    ax.FontSize =14;

    map = brewermap([], 'OrRd');%'-RdYlBu'
    colormap(map)
    cb = colorbar;
    cb.Ticks = [0;6;12];
    cb.LineWidth = 0.5;
    cb.FontSize = 14;
    cb.Title.String = 'F-values';    
    cb.Title.FontSize = 14;
end

