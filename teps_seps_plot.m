% this script: 
% 1) plot multi-channel TEPs/SEPs for all conditions
% Note: in SMAData structure, SHAM_...from each subject is the dataset for computing SEPs
% 2) plot TEPs/SEPs topography
% 2) illustrate cluster F-test results; illustrate cluster t-test results
%% Compute condition specific TEP average over trials
clc, clear
% Load data
load('..\SMA_TEPs_PREPOS.mat')
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
%% Plot grand average TEPs in time 
%%%% Step1: compute data of interest %%%%
load('...\original_chanlocs.mat');
cond = {'PRE','POS','POSminusPRE','SCOminusPLA','BIPminusPLA','SCOminusBIP'};
condName = {'PRE','POS','POS-PRE'};

cfg=[];
cfg.operation =  'subtract';
cfg.parameter =  'individual';
COND1.(cond{3}) =ft_math(cfg, GA_TEPs_clean_PLA_POS,GA_TEPs_clean_PLA_PRE);
%COND1.(cond{3}) =ft_math(cfg, GA_TEPs_SHAM_PLA_POS,GA_TEPs_SHAM_PLA_PRE); % for SEPs
COND1.(cond{2}) = GA_TEPs_clean_PLA_POS;
COND1.(cond{1}) = GA_TEPs_clean_PLA_PRE;

COND2.(cond{3}) =ft_math(cfg, GA_TEPs_clean_BIP_POS, GA_TEPs_clean_BIP_PRE);
COND2.(cond{2}) = GA_TEPs_clean_BIP_POS;
COND2.(cond{1}) = GA_TEPs_clean_BIP_PRE;

COND3.(cond{3}) =ft_math(cfg, GA_TEPs_clean_SCO_POS, GA_TEPs_clean_SCO_PRE);
COND3.(cond{2}) = GA_TEPs_clean_SCO_POS;
COND3.(cond{1}) = GA_TEPs_clean_SCO_PRE;

COND4.(cond{4}) =ft_math(cfg, COND3.(cond{3}),COND1.(cond{3})); % Time:sco-pla
COND4.(cond{5}) =ft_math(cfg, COND2.(cond{3}),COND1.(cond{3})); % Time:bip-pla
COND4.(cond{6}) =ft_math(cfg, COND3.(cond{3}),COND2.(cond{3})); % Time:sco-pla

%%%% Step2: layout %%%%
fig = figure('color','w');
set(gcf,'position',[0.1 0.02 0.4 0.6]*1000);
c = get(0, 'DefaultAxesColorOrder');
tiledlayout(fig,3,1,"TileSpacing","tight")

chan2plot = 'Fz'; %'FC1';'P3';'Fz'
chpIndx = find(strcmpi(chan2plot, {original_chanlocs.labels}));
times = COND1.(cond{3}).time; 
time2plot = [45,61,110,182]./1000;%!!!!!!!!!!!!!! mpfc
%  time2plot = [41,62,110,182]./1000;%!!!!!!!!!!!!!! sma
% time2plot = [36,68,109,182]./1000;%!!!!!!!!!!!!!! ag
[~,t1] = min(abs(times*1000 - -500));
[~,t2] = min(abs(times*1000 - -5));
[~,t3] = min(abs(times*1000 - 20));
[~,t4] = min(abs(times*1000 - 500));

% 1st row
ax(1) = nexttile;
clear meanData
meanData = squeeze(mean(COND3.(cond{1}).individual,1)); 
plot(times(t1:t2),meanData(:,(t1:t2)),'Linewidth',1.5,'Color',[120,120,120]./255);
hold on
plot(times(t1:t2),meanData(chpIndx,t1:t2),'Linewidth',1.5,'Color',[53,94,182]./255); 
plot(times(t3:t4),meanData(:,(t3:t4)),'Linewidth',1.5,'Color',[120,120,120]./255);
plot(times(t3:t4),meanData(chpIndx,t3:t4),'Linewidth',1.5,'Color',[53,94,182]./255); 
ax(1).XLabel.String = 'Time (s)';
ax(1).YLabel.String = 'Amplitude (\muV)';
ax(1).XLim = [-50,300]./1000;
ax(1).YLim = [-5,5];
ax(1).YTick = -5:2.5:5;
ax(1).YAxis.TickDirection = 'out';
ax(1).YAxis.TickLength = [0.005 0.02];
ax(1).FontSize =14;
ax(1).Box  = "off";
ax(1).LineWidth = 1;
% title(ax(1),'PRE')
% plot dashline
plot([0,0],get(gca,'ylim'),'-k','linewidth',1.5);   
for dl =1: length(time2plot)
 hold on
 plot([time2plot(dl),time2plot(dl)],get(gca,'ylim'),'k--','linewidth',0.9,'Color',[144,144,144]./255);
end 

% 2nd row
ax(2) = nexttile;
clear meanData
meanData = squeeze(mean(COND3.(cond{2}).individual,1)); 
plot(times(t1:t2),meanData(:,(t1:t2)),'Linewidth',1.5,'Color',[120,120,120]./255);
hold on
plot(times(t1:t2),meanData(chpIndx,t1:t2),'Linewidth',1.5,'Color',[53,94,182]./255); 
plot(times(t3:t4),meanData(:,(t3:t4)),'Linewidth',1.5,'Color',[120,120,120]./255);
plot(times(t3:t4),meanData(chpIndx,t3:t4),'Linewidth',1.5,'Color',[53,94,182]./255); 
ax(2).XLabel.String = 'Time (s)';
ax(2).YLabel.String = 'Amplitude (\muV)';
ax(2).XLim = [-50,300]./1000;
ax(2).YLim = [-5,5];
ax(2).YTick = -5:2.5:5;
ax(2).YAxis.TickDirection = 'out';
ax(2).YAxis.TickLength = [0.005 0.02];
ax(2).FontSize =14;
ax(2).Box  = "off";
ax(2).LineWidth = 1;
% title(ax(2),'POS')
% plot dashline
plot([0,0],get(gca,'ylim'),'-k','linewidth',1.5);   
for dl =1: length(time2plot)
 hold on
 plot([time2plot(dl),time2plot(dl)],get(gca,'ylim'),'k--','linewidth',0.9,'Color',[144,144,144]./255);
end 

% 3rd row
ax(3) = nexttile;
clear meanData
meanData = squeeze(mean(COND3.(cond{3}).individual,1)); 
plot(times(t1:t2),meanData(:,(t1:t2)),'Linewidth',1.5,'Color',[120,120,120]./255);
hold on
plot(times(t1:t2),meanData(chpIndx,t1:t2),'Linewidth',1.5,'Color',[53,94,182]./255); 
plot(times(t3:t4),meanData(:,(t3:t4)),'Linewidth',1.5,'Color',[120,120,120]./255);
plot(times(t3:t4),meanData(chpIndx,t3:t4),'Linewidth',1.5,'Color',[53,94,182]./255); 
ax(3).XLabel.String = 'Time (s)';
ax(3).YLabel.String = 'Amplitude (\muV)';
ax(3).XLim = [-50,300]./1000;
ax(3).YLim = [-5,5];
ax(3).YTick = -5:2.5:5;
ax(3).YAxis.TickDirection = 'out';
ax(3).YAxis.TickLength = [0.005 0.02];
ax(3).FontSize =14;
ax(3).Box  = "off";
ax(3).LineWidth = 1;
% title(ax(3),'POS-PRE')
% plot dashline
plot([0,0],get(gca,'ylim'),'-k','linewidth',1.5);   
for dl =1: length(time2plot)
 hold on
 plot([time2plot(dl),time2plot(dl)],get(gca,'ylim'),'k--','linewidth',0.9,'Color',[144,144,144]./255);
end 
%% Plot grand average TEPs in space
%%%% Step1: compute data of interest %%%%
load('...\neighbour_layout64_custom.mat')
cond = {'PRE','POS','POSminusPRE','SCOminusPLA','BIPminusPLA','SCOminusBIP'};
condName = {'PRE','POS','POS-PRE'};

cfg=[];
cfg.operation =  'subtract';
cfg.parameter =  'individual';
COND1.(cond{3}) =ft_math(cfg, GA_TEPs_clean_PLA_POS,GA_TEPs_clean_PLA_PRE);
COND1.(cond{2}) = GA_TEPs_clean_PLA_POS;
COND1.(cond{1}) = GA_TEPs_clean_PLA_PRE;

COND2.(cond{3}) =ft_math(cfg, GA_TEPs_clean_BIP_POS, GA_TEPs_clean_BIP_PRE);
COND2.(cond{2}) = GA_TEPs_clean_BIP_POS;
COND2.(cond{1}) = GA_TEPs_clean_BIP_PRE;

COND3.(cond{3}) =ft_math(cfg, GA_TEPs_clean_SCO_POS, GA_TEPs_clean_SCO_PRE);
COND3.(cond{2}) = GA_TEPs_clean_SCO_POS;
COND3.(cond{1}) = GA_TEPs_clean_SCO_PRE;

COND4.(cond{4}) =ft_math(cfg, COND3.(cond{3}),COND1.(cond{3})); % Time:sco-pla
COND4.(cond{5}) =ft_math(cfg, COND2.(cond{3}),COND1.(cond{3})); % Time:bip-pla
COND4.(cond{6}) =ft_math(cfg, COND3.(cond{3}),COND2.(cond{3})); % Time:sco-pla

%%%% Step2: layout %%%%
fig = figure('color','w');
set(gcf,'position',[0.1 0.02 0.3 0.5]*1000);
c = get(0, 'DefaultAxesColorOrder');
tiledlayout(fig,3,4,"TileSpacing","none","TileIndexing","rowmajor")

time2plot = [45,61,110,182]./1000;%!!!!!!!!!!!!!! mpfc
%  time2plot = [41,62,110,182]./1000;%!!!!!!!!!!!!!! sma
% time2plot = [36,68,109,182]./1000;%!!!!!!!!!!!!!! ag
tpIndx = dsearchn(COND3.(cond{1}).time', time2plot'); 

% 1st row--- PRE
clear meanData plotStruc
meanData = squeeze(mean(COND3.(cond{1}).individual,1)); 
plotStruc.time = 1;
plotStruc.dimord = 'chan';
plotStruc.label = COND3.(cond{1}).label;
for ti =1:4
clear plotStruc.avg
plotStruc.avg = meanData(:,tpIndx(ti));
ax(ti) = nexttile;
    cfg = [];
    cfg.layout = neighbour_layout.layout;
    cfg.comment = 'no';
    cfg.interactive = 'no';
    cfg.marker = 'off';
    cfg.zlim = [-2,2];
    cfg.figure      = 'gca';
    cfg.style = 'straight';
    ft_topoplotER(cfg,plotStruc);  
%     ax(ti).Title.String = [num2str(time2plot(ti)) 's' ];
    ax(ti).FontSize =14;
end
    map = brewermap([], '-RdYlBu');
    colormap(map)
       
% 2nd row --- POS
clear meanData plotStruc
meanData = squeeze(mean(COND3.(cond{2}).individual,1)); 
plotStruc.time = 1;
plotStruc.dimord = 'chan';
plotStruc.label = COND2.(cond{1}).label;
for ti =1:4
clear plotStruc.avg
plotStruc.avg = meanData(:,tpIndx(ti));
ax(ti+4) = nexttile;
    cfg = [];
    cfg.layout = neighbour_layout.layout;
    cfg.comment = 'no';
    cfg.interactive = 'no';
    cfg.marker = 'off';
    cfg.zlim = [-2,2];
    cfg.figure      = 'gca';
    cfg.style = 'straight';
    ft_topoplotER(cfg,plotStruc);  
    ax(ti+4).FontSize =14;
end
    map = brewermap([], '-RdYlBu');
    colormap(map)

% 3rd row --- POS-PRE
clear meanData plotStruc
meanData = squeeze(mean(COND3.(cond{3}).individual,1)); 
plotStruc.time = 1;
plotStruc.dimord = 'chan';
plotStruc.label = COND3.(cond{1}).label;
for ti =1:4
clear plotStruc.avg
plotStruc.avg = meanData(:,tpIndx(ti));
ax(ti+8) = nexttile;
    cfg = [];
    cfg.layout = neighbour_layout.layout;
    cfg.comment = 'no';
    cfg.interactive = 'no';
    cfg.marker = 'off';
    cfg.zlim = [-2,2];
    cfg.figure      = 'gca';
    cfg.style = 'straight';
    ft_topoplotER(cfg,plotStruc);  
    ax(ti+8).FontSize =14;
end
    map = brewermap([], '-RdYlBu');
    colormap(map)
%% Plot TEPs in space (at selected time point)
%%%% Step1: compute data of interest %%%%
load('...\neighbour_layout64_custom.mat')
cond = {'PRE','POS','POSminusPRE','SCOminusPLA','BIPminusPLA','SCOminusBIP'};
condName = {'PRE','POS','POS-PRE'};

cfg=[];
cfg.operation =  'subtract';
cfg.parameter =  'individual';
COND1.(cond{3}) =ft_math(cfg, GA_TEPs_clean_PLA_POS,GA_TEPs_clean_PLA_PRE);
COND1.(cond{2}) = GA_TEPs_clean_PLA_POS;
COND1.(cond{1}) = GA_TEPs_clean_PLA_PRE;

COND2.(cond{3}) =ft_math(cfg, GA_TEPs_clean_BIP_POS, GA_TEPs_clean_BIP_PRE);
COND2.(cond{2}) = GA_TEPs_clean_BIP_POS;
COND2.(cond{1}) = GA_TEPs_clean_BIP_PRE;

COND3.(cond{3}) =ft_math(cfg, GA_TEPs_clean_SCO_POS, GA_TEPs_clean_SCO_PRE);
COND3.(cond{2}) = GA_TEPs_clean_SCO_POS;
COND3.(cond{1}) = GA_TEPs_clean_SCO_PRE;

COND4.(cond{4}) =ft_math(cfg, COND3.(cond{3}),COND1.(cond{3})); % Time:sco-pla
COND4.(cond{5}) =ft_math(cfg, COND2.(cond{3}),COND1.(cond{3})); % Time:bip-pla
COND4.(cond{6}) =ft_math(cfg, COND3.(cond{3}),COND2.(cond{3})); % Time:sco-pla

%%%% Step2 plot %%%%
fig = figure('color','w');
set(gcf,'position',[0.1 0.02 0.3 0.8]*1000);
c = get(0, 'DefaultAxesColorOrder');
tiledlayout(fig,3,1)

time2plot = 0.047;
% 1st row--- PRE
clear meanData plotStruc
meanData = squeeze(mean(COND3.(cond{1}).individual,1)); 
plotStruc.time = 1;
plotStruc.dimord = 'chan';
plotStruc.label = COND3.(cond{1}).label;
tpIndx = dsearchn(COND3.(cond{1}).time', time2plot); 
plotStruc.avg = meanData(:,tpIndx);
ax(1) = nexttile;
    cfg = [];
    cfg.layout = neighbour_layout.layout;
    cfg.comment = 'no';
    cfg.interactive = 'no';
    cfg.markersymbol = '.';
    cfg.zlim = [-2,2];
    cfg.figure      = 'gca';
    cfg.style = 'straight';
    ax(1).Title.String = [num2str(time2plot) 's' ];
    ax(1).FontSize =14;
    ft_topoplotER(cfg,plotStruc);
    map = brewermap([], '-RdYlBu');
    colormap(map)
    cb = colorbar;
    cb.Layout.Tile = 'south';
    cb.Ticks = [-2;0;2];
    cb.LineWidth = 0.5;
    cb.FontSize = 14;
    title(cb,'\muV','fontsize',14);
       
% 2nd row --- POS
clear meanData plotStruc
meanData = squeeze(mean(COND3.(cond{2}).individual,1)); 
plotStruc.time = 1;
plotStruc.dimord = 'chan';
plotStruc.label = COND1.(cond{1}).label;
tpIndx = dsearchn(COND1.(cond{1}).time', time2plot); 
plotStruc.avg = meanData(:,tpIndx);
ax(2) = nexttile;
    cfg = [];
    cfg.layout = neighbour_layout.layout;
    cfg.comment = 'no';
    cfg.interactive = 'no';
    cfg.markersymbol = '.';
    cfg.zlim = [-2,2];
    cfg.figure      = 'gca';
    cfg.style = 'straight';
    ft_topoplotER(cfg,plotStruc);
    map = brewermap([], '-RdYlBu');
    colormap(map)

% 3rd row --- POS-PRE
clear meanData plotStruc
meanData = squeeze(mean(COND3.(cond{3}).individual,1)); 
plotStruc.time = 1;
plotStruc.dimord = 'chan';
plotStruc.label = COND1.(cond{1}).label;
tpIndx = dsearchn(COND1.(cond{1}).time', time2plot); 
plotStruc.avg = meanData(:,tpIndx);
ax(3) = nexttile;
    cfg = [];
    cfg.layout = neighbour_layout.layout;
    cfg.comment = 'no';
    cfg.interactive = 'no';
    cfg.markersymbol = '.';
    cfg.zlim = [-2,2];
    cfg.figure      = 'gca';
    cfg.style = 'straight';
    ft_topoplotER(cfg,plotStruc);
    map = brewermap([], '-RdYlBu');
    colormap(map)

%% Plot F/t values from cluster F-test/t-tests
clear 
%%%% Step1 %%%%
stat_select = ; % stat results from cluster t-tests
stat_brod = ; % stat results from cluster F-tests
% One Info from LOADED stat:
%  which value ( 'F-value' or 't-value')
WhichValue = 't-value'; 

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
    if stat_inuse.posclusters(clusx).prob<alpha_lim
    posClust.selCh{clusx} =stat_inuse.label(any(stat_inuse.posclusterslabelmat==clusx,2)); 
    posClust.minT{clusx} = min(stat_inuse.time( any(stat_inuse.posclusterslabelmat==clusx,1) )); 
    posClust.maxT{clusx} = max(stat_inuse.time( any(stat_inuse.posclusterslabelmat==clusx,1) ));
    end
  end
  end
end 
posClust
% Cluster content 
Channels_in_p =stat_inuse.label(any(stat_inuse.posclusterslabelmat==1,2));
Times_in_p = [min(stat_inuse.time( any(stat_inuse.posclusterslabelmat==1,1) ))...
    max(stat_inuse.time( any(stat_inuse.posclusterslabelmat==1,1) ))]; % ==1 1st cluster/largest cluster

% negative
negClust=[];
if isfield(stat_inuse,'negclusters')
    if ~isempty(stat_inuse.negclusters)
    for clusx =1: length(stat_inuse.negclusters)
        if stat_inuse.negclusters(clusx).prob<alpha_lim
        negClust.selCh{clusx} =stat_inuse.label(any(stat_inuse.negclusterslabelmat==clusx,2)); % 2 dimension
        negClust.minT{clusx} = min(stat_inuse.time( any(stat_inuse.negclusterslabelmat==clusx,1) ));
        negClust.maxT{clusx} = max(stat_inuse.time( any(stat_inuse.negclusterslabelmat==clusx,1) ));
        end
    end
    end
end
negClust
% Cluster content
Channels_in_n = stat_inuse.label(any(stat_inuse.negclusterslabelmat==1,2)); %==1 1st cluster/largest cluster
Times_in_n = [min(stat_inuse.time( any(stat_inuse.negclusterslabelmat==1,1) ))...
    max(stat_inuse.time( any(stat_inuse.negclusterslabelmat==1,1) ))]; % ==1 1st cluster/largest cluster

%%%% Step2 F/t values in channel-time axis %%%%
% layout
fig = figure('color','w');
set(gcf,'position',[0.1 0.02 0.8 0.3]*1000);
tiledlayout(fig,1,2,"TileSpacing","loose")

time2plot = 0.049; % in s
clearvars plot_clus Channels_in

plot_clus = zeros(size(stat_inuse.prob));
if isempty(posClust) && ~isempty(negClust)
     Channels_in = Channels_in_n;
     plot_clus(stat_inuse.negclusterslabelmat==1) =  -1; % negative cluster

elseif isempty(negClust) && ~isempty(posClust)
     Channels_in = Channels_in_p;
     plot_clus(stat_inuse.posclusterslabelmat==1) =  1; % positive cluster


elseif ~isempty(negClust) && ~isempty(posClust)
     Channels_in = union(Channels_in_n,Channels_in_p);
     plot_clus(stat_inuse.posclusterslabelmat==1) =  1; % positive cluster
     plot_clus(stat_inuse.negclusterslabelmat==1) =  -1; % negative cluster
end

if ~isequal(WhichValue,  'F-value') % t-value       
        % left
        ax(2) = nexttile;
        imagesc(stat_inuse.time, 1:size(stat_inuse.label,1),  stat_inuse.stat );
        hold on
        line([time2plot, time2plot], [0,size(stat_inuse.label,1)+1],'Color', 'black', 'LineStyle','--','LineWidth',1 );
        ax(2).XLabel.String = 'Time (s)';
        ax(2).YLabel.String = 'Channel number';
        ax(2).XLim = [0.02 0.3];
        ax(2).CLim =[-4 4];
        ax(2).XTick = [0.02, 0.1, 0.2, 0.3];
        ax(2).YTick = [1,11, 22];
        ax(2).YAxis.TickDirection = 'out';
        ax(2).YAxis.TickLength = [0.005 0.02];
        ax(2).FontSize =14;
        
        [map, ~, ~] = colorcet('D1A');
%         map = brewermap([], 'RdBu');%'-RdYlBu'
        colormap(map)
        cb = colorbar;
        cb.Ticks = [-4;0;4];
        cb.LineWidth = 0.5;
        cb.FontSize = 14;
        cb.Title.String = 't-values';    
        cb.Title.FontSize = 14;
        cb.Location ="northoutside";

        % right
        ax(1) = nexttile;
        imagesc(stat_inuse.time, 1:size(stat_inuse.label,1),plot_clus)
        ax(1).XLabel.String = 'Time (s)';
        ax(1).YLabel.String = 'Channel number';
        ax(1).XLim = [0.02 0.3];
        ax(1).XTick = [0.02, 0.1, 0.2, 0.3];
        ax(1).CLim =[-1 1];
        ax(1).YTick = [1,11,22];
        ax(1).YAxis.TickDirection = 'out';
        ax(1).YAxis.TickLength = [0.005 0.02];
       ax(1).FontSize =14;
else 
 
        % left
        ax(2) = nexttile;
        imagesc(stat_inuse.time, 1:size(stat_inuse.label,1),  stat_inuse.stat );
        hold on
        line([time2plot, time2plot], [0,size(stat_inuse.label,1)+1],'Color', 'black', 'LineStyle','--','LineWidth',1 );
        ax(2).XLabel.String = 'Time (s)';
        ax(2).YLabel.String = 'Channel number';
        ax(2).XLim = [0.02 0.3];
        ax(2).XTick = [0.02, 0.1, 0.2, 0.3];
        ax(2).CLim =[0 12];
        ax(2).YTick = [1, 32, 64];
        ax(2).YAxis.TickDirection = 'out';
        ax(2).YAxis.TickLength = [0.005 0.02];
        ax(2).FontSize =14;
        map = brewermap([], 'OrRd');%'-RdYlBu'
        colormap(map)
            cb = colorbar;
            cb.Ticks = [0;6;12];
            cb.LineWidth = 0.5;
            cb.FontSize = 14;
            cb.Title.String = 'F-values';    
            cb.Title.FontSize = 14;
            cb.Location ="northoutside";

         % right
        ax(1) = nexttile;
        imagesc(stat_inuse.time, 1:size(stat_inuse.label,1),plot_clus)
        ax(1).XLabel.String = 'Time (s)';
        ax(1).YLabel.String = 'Channel number';
        ax(1).XLim = [0.02 0.3];
        ax(1).XTick = [0.02, 0.1, 0.2, 0.3];
        ax(1).CLim =[0 1];
        ax(1).YTick = [1,32,64];
        ax(1).YAxis.TickDirection = 'out';
        ax(1).YAxis.TickLength = [0.005 0.02];
        ax(1).FontSize =14;
end

%%%% Step3 F/t values in space (topograhy) %%%%
clear valMap
if ~isequal(WhichValue,  'F-value') % t-value
        % Initialize matrix valMap(C) with 0
        valMap = zeros(size(stat_brod.prob));
        
        % Create a mapping between Channels and Times in stat_select.stat and stat_brod.stat
        [~, channel_indices_select_in_brod] = ismember(stat_inuse.label, stat_brod.label);
        [~, time_indices_select_in_brod] = ismember(stat_inuse.time, stat_brod.time);
        
        % Map the t-values from A to their corresponding positions in C
        for i = 1:size(stat_inuse.label,1)
            for j = 1:size(stat_inuse.time,2)
                channel_idx = channel_indices_select_in_brod(i);
                time_idx = time_indices_select_in_brod(j);
                valMap(channel_idx, time_idx) = stat_inuse.stat(i, j); % channel x time
            end
        end
else 
       valMap = stat_inuse.stat;
end 

% layout
fig = figure('color','w');
set(gcf,'position',[0.1 0.02 0.3 0.3]*1000);
c = get(0, 'DefaultAxesColorOrder');
tiledlayout(fig,1,1)

load('...\neighbour_layout64_custom.mat')
time2plot = 0.049; % in s

plotStruc = [];
plotStruc.time = 1;
plotStruc.dimord = 'chan';
plotStruc.label = stat_brod.label;
[~,tpIndx_in_select] = min(abs(stat_inuse.time-time2plot));
[~,tpIndx_in_brod] = min(abs(stat_brod.time-time2plot));
plotStruc.avg = valMap(:,tpIndx_in_brod);

if ~isequal(WhichValue,  'F-value')
    cfg = [];
    cfg.layout = neighbour_layout.layout;
    cfg.comment = 'no';
    cfg.interactive = 'no';
    cfg.markersymbol = '.';
    cfg.zlim = [-4,4];
    cfg.figure      = 'gca';
    cfg.interpolation      = 'nearest'; 
    cfg.style = 'straight';

    if isfield(stat_inuse,'posclusters')
        if ~isempty(stat_inuse.posclusters)
            if stat_inuse.posclusters(1).prob < alpha_lim 
                cfg.highlight = 'on';
                cfg.highlightchannel   =  stat_inuse.label(stat_inuse.posclusterslabelmat(:,tpIndx_in_select)==1); % ==1 the 1st cluster
                cfg.highlightsymbol    = '.';
                cfg.highlightsize      = 16;
                cfg.highlightcolor   = [0 0 0];
            end
        end
    end
    if isfield(stat_inuse,'negclusters')
        if ~isempty(stat_inuse.negclusters)
            if stat_inuse.negclusters(1).prob < alpha_lim 
                if isfield(cfg,'highlightchannel') && length(cfg.highlightchannel) ~= 0

                    if ~isempty(cfg.highlightchannel) % positive + negative channels highlighted
                        cfg.highlight = {'on','on'};
                        cfg.highlightchannel   =  {cfg.highlightchannel,stat_inuse.label(stat_inuse.negclusterslabelmat(:,tpIndx_in_select)==1)};% ==1 the 1st cluster
                        cfg.highlightsymbol    = {'.','.'};
                        cfg.highlightsize      = {16,16};
                        cfg.highlightcolor   = {[0 0 0],[0.99 0.99 0.99]};
                    end
                else
                    cfg.highlight = 'on';            % only neative channels highlighted
                    cfg.highlightchannel   =  stat_inuse.label(stat_inuse.negclusterslabelmat(:,tpIndx_in_select)==1);% ==1 the 1st cluster
                    cfg.highlightsymbol    = '.';
                    cfg.highlightsize      = 16;
                    cfg.highlightcolor   = [ 0.99 0.99 0.99];
                end
            end
        end
    end

    ax = nexttile;
    ft_topoplotER(cfg,plotStruc);
    ax.Title.String = [num2str(time2plot) 's'];
    ax.FontSize =14;

%     [map, ~, ~] = colorcet('D1A');
    map = brewermap([], '-RdBu');
    colormap(map)
    cb = colorbar;
    cb.Ticks = [-4;0;4];
    cb.LineWidth = 0.5;
    cb.Title.String = 't-values';
    cb.Title.FontSize = 14;
else 

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
            if stat_inuse.posclusters(1).prob < alpha_lim % ONLY 1 cluster may have the chance to be less than 0.025
                cfg.highlight = 'on';
                cfg.highlightchannel   =  stat_inuse.label(stat_inuse.posclusterslabelmat(:,tpIndx_in_select)==1); % ==1 the 1st cluster
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
                        cfg.highlightchannel   =  {cfg.highlightchannel,stat_inuse.label(stat_inuse.negclusterslabelmat(:,tpIndx_in_select)==1)};% ==1 the 1st cluster
                        cfg.highlightsymbol    = {'.','.'};
                        cfg.highlightsize      = {16,16};
                        cfg.highlightcolor   = {[0 0 0],[0.99 0.99 0.99]};
                    end
                else
                    cfg.highlight = 'on';            % only neative channels highlighted
                    cfg.highlightchannel   =  stat_inuse.label(stat_inuse.negclusterslabelmat(:,tpIndx_in_select)==1);% ==1 the 1st cluster
                    cfg.highlightsymbol    = '.';
                    cfg.highlightsize      = 16;
                    cfg.highlightcolor   = [ 0.99 0.99 0.99];
                end
            end
        end
    end

    ax = nexttile;
    ft_topoplotER(cfg,plotStruc);
    ax.Title.String = [num2str(time2plot) 's'];
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