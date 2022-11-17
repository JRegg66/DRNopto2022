%% Figure 3
% plotting script
% set default plot settings
set(groot,...
    'DefaultFigureRenderer','painters',...
    'DefaultFigureColor','w',...
    'DefaultAxesLineWidth',.5,...
    'DefaultAxesXColor','k',...
    'DefaultAxesYColor','k',...
    'DefaultAxesFontUnits','points',...
    'DefaultAxesFontSize',7,...
    'DefaultAxesFontName','Arial',...
    'DefaultAxesBox','off',...
    'DefaultAxesTickDir','in',...
    'DefaultAxesTickDirMode','manual',...
    'DefaultLineLineWidth',1,...
    'DefaultHistogramLineWidth',1,...
    'DefaultTextFontUnits','Points',...
    'DefaultTextFontSize',7,...
    'DefaultTextFontName','Arial');

%% load data
clear
PlotColours
% load WholeLGN_Data_Nov2022.mat

%% pull out data

clearvars -except data ColourMap
fr = 15.5; QI_thresh = .15;
QIbar = []; QIff = []; SBC = []; Allff = []; QI = []; Base = []; BigRespAll = []; XY = []; FOV = []; P = []; Ridx = []; Sidx = [];
BaseDR = []; BaseOptoDR = []; BigRespAllDR = []; OOidx = []; VertPref = []; HorzPref = []; OnOff =  []; RetResp = [];

for i = 1:6
    boutons = size(data(i).OptoMean,2);
    OptoQI = nan(boutons,length(data(i).OptoTraces));
    OptoMean = nan(204,boutons,3);
    for k = 1:length(data(i).OptoTraces)
        OptoTraces{k} = data(i).OptoTraces{k}(:,data(i).Arousal.Opto.FractClosedEye{k} < .15,:);
        for b = 1:boutons
            OptoTraces{k}(:,:,b) = smoothdata(OptoTraces{k}(:,:,b),1,'movmean',7);
        end
        OptoMean(:,:,k) = squeeze(nanmean(OptoTraces{k},2));
        meanresp = squeeze(nanmean(OptoTraces{k}(round(3*fr):end,:,:),2));
        varmeanresp = squeeze(var(meanresp,[],1));
        varresp = squeeze(var(OptoTraces{k}(round(3*fr):end,:,:),[],1,'omitnan'));
        meanvarresp = squeeze(nanmean(varresp,1));
        OptoQI(:,k) = varmeanresp./meanvarresp;
        clear meanresp varmeanresp varresp meanvarresp
    end
    % 20th and 80th percentile for normalizing
    perc = [20,80]; TotResp = [];
    xstart = [4,5.1,6,7.1,8,9.1,10,11.1]; xend = [5,6,7,8,9,10,11,12];
    for p = 1:length(perc)
        TotResp = cat(1,TotResp,(prctile(OptoMean(1:round(2*fr),:,1),perc(p))));
        for k = 1 % if you want to include opto trials, make this 1:2
            for j = 1:length(xstart) % ctrl lum
                TotResp = cat(1,TotResp,(prctile(OptoMean(round(xstart(j)*fr):round(xend(j)*fr),:,k),perc(p))));
            end
        end
        TotResp = cat(1,TotResp,permute(squeeze(prctile(data(i).RetMean(16:31,:,:),perc(p))),[2,1]),...
            permute(squeeze(prctile(data(i).RetMean(32:47,:,:),perc(p))),[2,1]));
    end
    bot = min(TotResp,[],1); TotResp = TotResp - min(TotResp,[],1);
    top = max(TotResp,[],1); TotResp = TotResp./max(TotResp,[],1);
    BigRespAll = cat(2,BigRespAll,cat(1,reshape(permute(OptoMean,[1,3,2]),204*3,boutons),...
        reshape(permute(data(i).RetMean(1:65,:,:),[1,3,2]),65*18,boutons)));
    BigRespAllDR = cat(2,BigRespAllDR,(cat(1,reshape(permute(OptoMean,[1,3,2]),204*3,boutons),...
        reshape(permute(data(i).RetMean(1:65,:,:),[1,3,2]),65*18,boutons))-bot)./top);
    Base = cat(2,Base,squeeze(mean(mean(OptoMean(round(4*fr):round(5*fr),:,[2:3]),1),3)));
    BaseDR = cat(2,BaseDR,(squeeze(mean(mean(OptoMean(1:round(2*fr),:,:),1),3))-bot)./top);
    BaseOptoDR = cat(2,BaseOptoDR,(squeeze(mean(mean(OptoMean(round(4*fr):round(5*fr),:,[2:3]),1),3))-bot)./top);
    QIbar = cat(1,QIbar,sum([max(data(i).Ret.QI(:,1:8),[],2) > QI_thresh,max(data(i).Ret.QI(:,9:16),[],2) > QI_thresh],2));
    QIff = cat(1,QIff,max(data(i).Ret.QI(:,17:18),[],2) > QI_thresh);
    Ret = squeeze(mean(data(i).RetMean(16:47,:,1:16),1)); 
    RetResp = cat(1,RetResp,Ret);
    HorzPref = cat(1,HorzPref,cat(1,data(i).Ret.Retinotopy.horzPref));
    VertPref = cat(1,VertPref,cat(1,data(i).Ret.Retinotopy.vertPref));
    OnOff = cat(1,OnOff,squeeze(mean(data(i).RetMean(16:47,:,17:18),1)));
    bad = data(i).Ret.QI(:,1:16) < QI_thresh;
    Ret(bad) = NaN;
    S = cat(2,sum(Ret(:,1:8) < 0 & ~isnan(Ret(:,1:8)),2)>0 & sum(Ret > 0 & ~isnan(Ret),2) == 0,...
        sum(Ret(:,9:16) < 0 & ~isnan(Ret(:,9:16)),2)>0 & sum(Ret > 0 & ~isnan(Ret),2) == 0);
    SBC = cat(1,SBC,sum(S,2)>=2);
    Allff = cat(2,Allff,data(i).Ret.FFpref');
    OOidx = cat(1,OOidx,data(i).Ret.OnOffPref);
    Ridx = cat(1,Ridx,data(i).Ret.RampIdx);
    Sidx = cat(1,Sidx,data(i).Ret.SustainIdx);
    QI = cat(1,QI,cat(2,OptoQI(:,1:2),...
        max(data(i).Ret.QI(:,1:8),[],2),max(data(i).Ret.QI(:,9:16),[],2),data(i).Ret.QI(:,17:18)));
    XY = cat(1,XY,data(i).ROIxyReg);
    FOV = cat(1,FOV, ones(boutons,1).*i);
    P = cat(2,P,data(i).posterior);
end

toss = ~((QIbar == 2 | QIff == 1)); % toss = toss | SBC; FOV2 = FOV;
SBC(toss) = []; Allff(toss) = []; QIbar(toss) = []; QIff(toss) = []; QI(toss,:) = []; BigRespAllDR(:,toss) = []; OOidx(QIff == 0) = NaN; OOidx(toss) = [];
Base(toss) = []; BigRespAll(:,toss) = []; XY(toss,:) = []; FOV(toss) = []; BaseDR(toss) = []; BaseOptoDR(toss) = [];
HorzPref(toss) = []; VertPref(toss) = []; OnOff(toss,:) = []; RetResp(toss,:) = []; P(toss) = []; Ridx(toss) = []; Sidx(toss) = [];
clearvars -except data SBC QIbar QIff QI Base Allff BigRespAll XY ColourMap ...
    FOV BaseOptoDR BaseDR BigRespAllDR OOidx OnOff HorzPref VertPref RetResp FOV2 toss P Ridx Sidx

% define the 4 categories
newFF = NaN(size(Allff)); 
newFF(QIff == 1 & QIbar < 2) = 1; % FF
newFF(QIff == 1 & QIbar == 2) = 2; % FF+Bar
newFF(QIff == 0 & QIbar == 2) = 3; % bar
newFF(SBC == 1) = 4; % SBC

titles = {'FF','FF+bar','bar','SBC'};
cmap = cat(1,ColourMap.FF,ColourMap.Both,ColourMap.Bar,ColourMap.SBC);

%% fig 3b. example traces

fr = 15.5;
FFexamples = find(newFF == 1 & BaseDR > .5 & (BaseOptoDR - BaseDR) < -.2 & QI(:,1)' > .25);
Bothexample = find(newFF == 2 & (BaseOptoDR - BaseDR) < -.1 & max(QI,[],2)' > .35 & ~SBC');
Barexample = find(newFF == 3 & BaseDR < .5 & (BaseOptoDR - BaseDR) < 0 & max(QI(:,3:4),[],2)' > .25);
SBCexample = find(newFF == 4 & (BaseOptoDR - BaseDR) < -.1 & max(QI,[],2)' > .35);
Ball = BigRespAll; Ball(205:408,:) = [];

Stim = NaN(size(Ball,1),1);
Lum = cat(1,zeros(ceil(5*fr),1),ones(round(2*fr),1)*-1,ones(round(2*fr),1),ones(round(2*fr),1)*-1,zeros(round(2*fr)+2,1));
Stim(1:204) = Lum; Stim(205:408) = 0; clear Lum
start = 1+408+65*16; Stim(start:start+16) = 0; Stim(start+16:start+47) = -1; Stim(start+48:start+65) = 0;
start = 66+408+65*16; Stim(start:start+16) = 0; Stim(start+16:start+47) = 1; Stim(start+48:start+65) = 0;
clear start

x = [204,65]; xpos = 0;
for j = 1:20
    if j <= 2
        xpos(j+1) = xpos(j)+x(1);
    else
        xpos(j+1) = xpos(j)+x(2);
    end
end
figure(5)
yyaxis left
patch([205+round(2*fr) 205+round(12*fr) 205+round(12*fr) 205+round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none'); hold on
for j = 1:length(xpos); 
    if j > 2 & j < 19
        patch([xpos(j)+(fr) xpos(j)+(3*fr) xpos(j)+(3*fr) xpos(j)+(fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
    end
end
plot(Stim,'-k','LineWidth',1);
box off; set(gca,'visible','off'); xlim([1 size(Ball,1)]); ylim([-1.1 1.1])
set(gcf,'units','centimeters','position',[10 10 7 1]); 
set(gca,'units','centimeters','position',[.2 .2 6.6 .8]); 

h = figure(1) % FF
i = FFexamples(1);
yyaxis left
patch([205+round(2*fr) 205+round(12*fr) 205+round(12*fr) 205+round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none'); hold on
plot(smooth(Ball(:,i),7.5),'-','Color',cmap(1,:),'LineWidth',1);
xlim([0 xpos(end)]); ylim([min(smooth(Ball(:,i),7.5)) max(smooth(Ball(:,i),7.5))]); box off
ylabel('\Delta F/F')
yyaxis right
for j = 1:length(xpos); 
    if j < 2
        patch([xpos(j)+(5*fr) xpos(j)+(11*fr) xpos(j)+(11*fr) xpos(j)+(5*fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
    elseif j > 2
        patch([xpos(j)+(fr) xpos(j)+(3*fr) xpos(j)+(3*fr) xpos(j)+(fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
    end
end

B = smooth(Ball(:,i),7.5); B = (B - min(B))./max(B - min(B));
plot(B,'-','Color',cmap(1,:),'LineWidth',1);
xlim([0 xpos(end)]); ylim([min(B) max(B)]); box off
ylabel('dynamic range')
set(gca,'xColor','none');
set(gcf,'units','centimeters','position',[10 12 7 2.8]); 

h = figure(2) % FF +Bar
i = Bothexample(1);
yyaxis left
patch([205+round(2*fr) 205+round(12*fr) 205+round(12*fr) 205+round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none'); hold on
plot(smooth(Ball(:,i),7.5),'-','Color',cmap(2,:),'LineWidth',1);
xlim([0 xpos(end)]); ylim([min(smooth(Ball(:,i),7.5)) max(smooth(Ball(:,i),7.5))]); box off
ylabel('\Delta F/F')
yyaxis right
for j = 1:length(xpos); 
    if j < 2
        patch([xpos(j)+(5*fr) xpos(j)+(11*fr) xpos(j)+(11*fr) xpos(j)+(5*fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
    elseif j > 2
        patch([xpos(j)+(fr) xpos(j)+(3*fr) xpos(j)+(3*fr) xpos(j)+(fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
    end
end
B = smooth(Ball(:,i),7.5); B = (B - min(B))./max(B - min(B));
plot(B,'-','Color',cmap(2,:),'LineWidth',1);
xlim([0 xpos(end)]); ylim([min(B) max(B)]); box off
ylabel('dynamic range')
set(gca,'xColor','none');
set(gcf,'units','centimeters','position',[10 15 7 2.8]); 

h = figure(3) % bar
i = Barexample(1);
yyaxis left
patch([205+round(2*fr) 205+round(12*fr) 205+round(12*fr) 205+round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none'); hold on
plot(smooth(Ball(:,i),7.5),'-','Color',cmap(3,:),'LineWidth',1);
xlim([0 xpos(end)]); ylim([min(smooth(Ball(:,i),7.5)) max(smooth(Ball(:,i),7.5))]); box off
ylabel('\Delta F/F')
yyaxis right
for j = 1:length(xpos); 
    if j < 2
        patch([xpos(j)+(5*fr) xpos(j)+(11*fr) xpos(j)+(11*fr) xpos(j)+(5*fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
    elseif j > 2
        patch([xpos(j)+(fr) xpos(j)+(3*fr) xpos(j)+(3*fr) xpos(j)+(fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
    end
end
B = smooth(Ball(:,i),7.5); B = (B - min(B))./max(B - min(B));
plot(B,'-','Color',cmap(3,:),'LineWidth',1);
line([1500 1655],[-.1 -.1]);
xlim([0 xpos(end)]); ylim([min(B) max(B)]); box off
ylabel('dynamic range')
set(gca,'xColor','none');
set(gcf,'units','centimeters','position',[10 18 7 2.8]); 

SBCexample = find(newFF == 4 & (BaseOptoDR - BaseDR) < -.1 & max(QI,[],2)' > .35);
h = figure(4); clf % SBC
i = SBCexample(10);
yyaxis left
patch([205+round(2*fr) 205+round(12*fr) 205+round(12*fr) 205+round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none'); hold on
plot(smooth(Ball(:,i),7.5),'-','Color',cmap(4,:),'LineWidth',1);
xlim([0 xpos(end)]); ylim([min(smooth(Ball(:,i),7.5)) max(smooth(Ball(:,i),7.5))]); box off
ylabel('\Delta F/F')
yyaxis right
for j = 1:length(xpos); 
    if j < 2
        patch([xpos(j)+(5*fr) xpos(j)+(11*fr) xpos(j)+(11*fr) xpos(j)+(5*fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
    elseif j > 2 
        patch([xpos(j)+(fr) xpos(j)+(3*fr) xpos(j)+(3*fr) xpos(j)+(fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
    end
end
B = smooth(Ball(:,i),7.5); B = (B - min(B))./max(B - min(B));
plot(B,'-','Color',cmap(4,:),'LineWidth',1);
line([1500 1655],[-.1 -.1]);
xlim([0 xpos(end)]); ylim([min(B) max(B)]); box off
ylabel('dynamic range')
set(gca,'xColor','none');
set(gcf,'units','centimeters','position',[10 18 7 2.8]); 

%% fig 3c. cdf of suppression for FF, both and bar

close all
bins = [-.2:.01:.05]; plotbins = [-.19:.01:.05];

for i = 1:max(newFF)
    figure(1); clf; h = histogram(Base(newFF == i),bins,'Normalization','cdf');
    [~,p(i)] = kstest(Base(newFF == i));
    figure(2); plot(plotbins,h.Values,'-','Color',cmap(i,:)); hold on
end
plot([0 0],[0 1],'--k');
box off; axis square
xlabel('\Delta F/F'); xlim([-.2 .05]); xticks([-.1 0])
ylabel('cumulative probability'); yticks([0 .5 1])

% for the inset
for i = 1:max(newFF)
    for j = 1:max(FOV)
        FOVbase(i,j) = nanmean(Base(newFF' == i & FOV == j));
    end
end
[p,~,stats] = kruskalwallis(FOVbase');
c = multcompare(stats);

figure(2); inset1 = axes('Position',[.25 .5 .2 .35]);
for j = 1:size(FOVbase,1)
    bar(j,nanmean(FOVbase(j,:),2),'FaceColor',cmap(j,:),'EdgeColor','none'); hold on
end
for j = 1:size(FOVbase,2)
    plot([1:size(FOVbase,1)],FOVbase(:,j),'-k','LineWidth',.5);
end
plot([0 size(FOVbase,1)+1],[0 0],'-k');
plot([1 3],[-.046 -.046],'-k'); text(1.9,-.049,'*');
box off; xticks([1:size(FOVbase,1)]); xticklabels(titles(1:size(FOVbase,1))); xtickangle(45)
set(gca, 'XAxisLocation', 'top'); set(gca, 'YAxisLocation', 'right');
xlim([0 size(FOVbase,1)+1]); ylabel('\Delta F/F')
ylim([-.05 0]); yticks([-.04 -.02 0]); yticklabels({'-.04', '-.02','0'})
set(gcf,'units','centimeters','position',[10 10 5 4]); 

%% LME
Suppression = Base;
Mouse = nominal(FOV);

% FF reference
Type = strsplit(num2str(newFF))';
Type(newFF == 1) = {'1FF'}; Type(newFF == 2) = {'2FF+Bar'}; Type(newFF == 3) = {'3Bar'}; Type(newFF == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,Mouse);
lmeFF = fitlme(Opto,'Suppression ~ 1 + Type + (1|Mouse)')

% FF reference
Type = strsplit(num2str(newFF))';
Type(newFF == 1) = {'2FF'}; Type(newFF == 2) = {'1FF+Bar'}; Type(newFF == 3) = {'3Bar'}; Type(newFF == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,Mouse);
lmeFFBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|Mouse)')

% FF reference
Type = strsplit(num2str(newFF))';
Type(newFF == 1) = {'3FF'}; Type(newFF == 2) = {'2FF+Bar'}; Type(newFF == 3) = {'1Bar'}; Type(newFF == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,Mouse);
lmeBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|Mouse)')

% FF reference
Type = strsplit(num2str(newFF))';
Type(newFF == 1) = {'4FF'}; Type(newFF == 2) = {'2FF+Bar'}; Type(newFF == 3) = {'3Bar'}; Type(newFF == 4) = {'1SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,Mouse);
lmeSBC = fitlme(Opto,'Suppression ~ Type + (1|Mouse)')

%% fig 3d. cdf + bar plot of baseline suppression in dynamic range 

close all
bins = [-.4:.05:.1]; plotbins = [-.375:.05:.1];
B = BaseOptoDR - BaseDR;
for i = 1:max(newFF)
    figure(1); h = histogram(B(newFF == i),bins,'Normalization','cdf');
    [~,p(i)] = kstest(B(newFF == i));
    figure(2); plot(plotbins,h.Values,'-','Color',cmap(i,:)); hold on
end
plot([0 0],[0 1],'--k');
axis square; box off; 
xlabel('opto-ctrl'); ylabel('cumulative probability')
xlim([-.4 .1]); xticks([-.2 0]); ylim([0 1]); yticks([0 .5 1])

% inset
B = BaseOptoDR - BaseDR;
clear FOVbase
for i = 1:max(newFF)
    for j = 1:max(FOV)
        FOVbase(i,j) = mean(B(newFF' == i & FOV == j));
    end
end
[~,~,stats] = kruskalwallis(FOVbase');
c = multcompare(stats);

figure(2); inset3 = axes('Position',[.25 .5 .2 .35]);
for j = 1:size(FOVbase,1)
    bar(j,nanmean(FOVbase(j,:),2),'FaceColor',cmap(j,:),'EdgeColor','none'); hold on
end
for j = 1:size(FOVbase,2)
    plot([1:size(FOVbase,1)],FOVbase(:,j),'-k','LineWidth',.5);
end
plot([0 size(FOVbase,1)+1],[0 0],'-k');
plot([1 3],[-.14 -.14],'-k'); text(1.6,-.15,'**');
plot([1 2],[-.16 -.16],'-k'); text(1.2,-.17,'*');
box off; xticks([1:size(FOVbase,1)]); xticklabels(titles(1:size(FOVbase,1))); xtickangle(45)
set(gca, 'XAxisLocation', 'top'); set(gca, 'YAxisLocation', 'right');
xlim([0 size(FOVbase,1)+1]); ylabel('opto-ctrl')
ylim([-.17 0]); yticks([-.1 0])
set(gcf,'units','centimeters','position',[10 10 5 4]); 

%% LME

Suppression = BaseOptoDR' - BaseDR';
Mouse = nominal(FOV);

% FF reference
Type = strsplit(num2str(newFF))';
Type(newFF == 1) = {'1FF'}; Type(newFF == 2) = {'2FF+Bar'}; Type(newFF == 3) = {'3Bar'}; Type(newFF == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,Mouse);
lmeFF = fitlme(Opto,'Suppression ~ 1 + Type + (1|Mouse)')

% FF+Bar reference
Type = strsplit(num2str(newFF))';
Type(newFF == 1) = {'2FF'}; Type(newFF == 2) = {'1FF+Bar'}; Type(newFF == 3) = {'3Bar'}; Type(newFF == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,Mouse);
lmeFFBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|Mouse)')

% Bar reference
Type = strsplit(num2str(newFF))';
Type(newFF == 1) = {'3FF'}; Type(newFF == 2) = {'2FF+Bar'}; Type(newFF == 3) = {'1Bar'}; Type(newFF == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,Mouse);
lmeBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|Mouse)')

% SBC reference
Type = strsplit(num2str(newFF))';
Type(newFF == 1) = {'4FF'}; Type(newFF == 2) = {'2FF+Bar'}; Type(newFF == 3) = {'3Bar'}; Type(newFF == 4) = {'1SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,Mouse);
lmeSBC = fitlme(Opto,'Suppression ~ Type + (1|Mouse)')

%% fig 3e. cdf plus bar plot of baseline level

close all
bins = [0:.1:1]; plotbins = [.05:.1:1];
for i = 1:max(newFF)
    figure(1); h = histogram(BaseDR(newFF == i),bins,'Normalization','cdf');
    figure(2); plot(plotbins,h.Values,'-','Color',cmap(i,:)); hold on
end
axis square; box off; xticks([0 .5 1]); yticks([0 .5 1]);
ylim([0 1]); xlim([0 1]); xlabel('baseline level'); ylabel('cumulative probability')

% inset
clear FOVbase
for i = 1:max(newFF)
    for j = 1:max(FOV)
        FOVbase(i,j) = nanmean(BaseDR(newFF' == i & FOV == j));
    end
end
[p,~,stats] = kruskalwallis(FOVbase');
c = multcompare(stats);

figure(2); inset2 = axes('Position',[.6 .4 .2 .35]);
for j = 1:size(FOVbase,1)
    bar(j,nanmean(FOVbase(j,:),2),'FaceColor',cmap(j,:),'EdgeColor','none'); hold on
end
for j = 1:size(FOVbase,2)
    plot([1:size(FOVbase,1)],FOVbase(:,j),'-k','LineWidth',.5);
end
plot([0 size(FOVbase,1)+1],[0 0],'-k');
plot([2 4],[.8 .8],'-k'); text(2.8,.83,'**');
plot([3 4],[.9 .9],'-k'); text(3,.93,'**');
box off; xticks([1:size(FOVbase,1)]); xticklabels(titles(1:size(FOVbase,1))); xtickangle(45)
xlim([0 size(FOVbase,1)+1]); ylabel('baseline level')
ylim([0 1]); yticks([0 .5 1])
set(gca, 'YAxisLocation', 'right');
set(gcf,'units','centimeters','position',[10 10 5 4]); 

%% LME
BaseLevel = BaseDR';
Mouse = nominal(FOV);

% FF reference
Type = strsplit(num2str(newFF))';
Type(newFF == 1) = {'1FF'}; Type(newFF == 2) = {'2FF+Bar'}; Type(newFF == 3) = {'3Bar'}; Type(newFF == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(BaseLevel,Type,Mouse);
lmeFF = fitlme(Opto,'BaseLevel ~ 1 + Type + (1|Mouse)')

% FF+Bar reference
Type = strsplit(num2str(newFF))';
Type(newFF == 1) = {'2FF'}; Type(newFF == 2) = {'1FF+Bar'}; Type(newFF == 3) = {'3Bar'}; Type(newFF == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(BaseLevel,Type,Mouse);
lmeFFBar = fitlme(Opto,'BaseLevel ~ 1 + Type + (1|Mouse)')

% Bar reference
Type = strsplit(num2str(newFF))';
Type(newFF == 1) = {'3FF'}; Type(newFF == 2) = {'2FF+Bar'}; Type(newFF == 3) = {'1Bar'}; Type(newFF == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(BaseLevel,Type,Mouse);
lmeBar = fitlme(Opto,'BaseLevel ~ 1 + Type + (1|Mouse)')

% SBC reference
Type = strsplit(num2str(newFF))';
Type(newFF == 1) = {'4FF'}; Type(newFF == 2) = {'2FF+Bar'}; Type(newFF == 3) = {'3Bar'}; Type(newFF == 4) = {'1SBC'}; 
Type = nominal(Type);
Opto = table(BaseLevel,Type,Mouse);
lmeSBC = fitlme(Opto,'BaseLevel ~ Type + (1|Mouse)')

%% fig 3f. scatter histogram

figure;
x = BaseDR;
y = BaseOptoDR - BaseDR;
b = corrcoef(x,y);
scatterhist(x,y,'Group',newFF,'Kernel','on','Location','NorthEast',...
    'Direction','out','Marker','.','MarkerSize',.2,'Color',cmap(1:max(newFF),:),...
    'LineStyle',{'-'},'LineWidth',1);
axis square; box off
legend(titles)
xlabel('baseline level'); ylabel('baseline opto - ctrl')
xticks([0 .5 1]); 
hold on; 
plot([0 1],[0 0],'--k');
plot([0 1],[0 -1],'--k');
xlim([0 1]); ylim([-.8 .2])
set(gcf,'units','centimeters','position',[10 10 5 4]); 

%% fig 3g. cdf + bar plot dynamic range of suppression POSTERIOR

bins = [-.4:.05:.1]; plotbins = [-.375:.05:.1];
B = BaseOptoDR - BaseDR;
for i = 1:max(newFF)
    figure(1); h = histogram(B(newFF == i & P == 1),bins,'Normalization','cdf');
    [~,p(i)] = kstest(B(newFF == i & P == 1));
    figure(2); plot(plotbins,h.Values,'-','Color',cmap(i,:)); hold on
end
plot([0 0],[0 1],'--k');
axis square; box off; 
xlabel('opto-ctrl'); ylabel('cumulative probability')
xlim([-.4 .1]); xticks([-.2 0]); ylim([0 1]); yticks([0 .5 1])

% inset
B = BaseOptoDR - BaseDR;
clear FOVbase
for i = 1:max(newFF)
    for j = 1:max(FOV)
        FOVbase(i,j) = mean(B(newFF' == i & FOV == j & P' == 1));
    end
end
[~,~,stats] = kruskalwallis(FOVbase');
c = multcompare(stats);

figure(2); inset3 = axes('Position',[.25 .5 .2 .35]);
for j = 1:size(FOVbase,1)
    bar(j,nanmean(FOVbase(j,:),2),'FaceColor',cmap(j,:),'EdgeColor','none'); hold on
end
for j = 1:size(FOVbase,2)
    plot([1:size(FOVbase,1)],FOVbase(:,j),'-k','LineWidth',.5);
end
plot([0 size(FOVbase,1)+1],[0 0],'-k');
box off; xticks([1:size(FOVbase,1)]); xticklabels(titles(1:size(FOVbase,1))); xtickangle(45)
set(gca, 'XAxisLocation', 'top'); set(gca, 'YAxisLocation', 'right');
xlim([0 size(FOVbase,1)+1]); ylabel('opto-ctrl')
ylim([-.1 0]); yticks([-.1 0])
set(gcf,'units','centimeters','position',[10 10 5 4]); 

%% LME

Suppression = BaseOptoDR(P' == 1)' - BaseDR(P' == 1)';
Mouse = nominal(FOV(P' == 1));
ff = newFF(P == 1);

% FF reference
Type = strsplit(num2str(ff))';
Type(ff == 1) = {'1FF'}; Type(ff == 2) = {'2FF+Bar'}; Type(ff == 3) = {'3Bar'}; Type(ff == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,Mouse);
lmeFF = fitlme(Opto,'Suppression ~ 1 + Type + (1|Mouse)')

% FF+Bar reference
Type = strsplit(num2str(ff))';
Type(ff == 1) = {'2FF'}; Type(ff == 2) = {'1FF+Bar'}; Type(ff == 3) = {'3Bar'}; Type(ff == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,Mouse);
lmeFFBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|Mouse)')

% Bar reference
Type = strsplit(num2str(ff))';
Type(ff == 1) = {'3FF'}; Type(ff == 2) = {'2FF+Bar'}; Type(ff == 3) = {'1Bar'}; Type(ff == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,Mouse);
lmeBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|Mouse)')

% SBC reference
Type = strsplit(num2str(ff))';
Type(ff == 1) = {'4FF'}; Type(ff == 2) = {'2FF+Bar'}; Type(ff == 3) = {'3Bar'}; Type(ff == 4) = {'1SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,Mouse);
lmeSBC = fitlme(Opto,'Suppression ~ Type + (1|Mouse)')

%% --------- SUPPLEMENTARY PANELS

%% supp fig 3A big heatmap

AllB = []; boutons = 0; Mid = [];
for i = 1:max(newFF)
    for j = 1:6
        I = newFF == i & FOV' == j;
        F = FOV(I); 
        thisB = cat(1,BigRespAll([1:204],I),...
            BigRespAll([409:612],I),...
            BigRespAll([613:end],I));
        [~,idx] = sort(Base(I));
        AllB = cat(2,AllB,thisB(:,idx));
        Mid = cat(1,Mid,F(idx));
    end
    I = newFF == i; boutons(i+1) = sum(I) + boutons(i);
end

figure(1);
imagesc(smoothdata(AllB,1)'); caxis([-.3 .3]); colormap('bluewhitered'); colorbar('northoutside')
hold on
for i = 2:length(boutons)
    plot([1 size(AllB,1)],[boutons(i) boutons(i)],'-k');
end
xticks([])
ylabel('boutons')
set(gcf,'units','centimeters','position',[5 5 15 15]); 
set(gca,'units','centimeters','position',[1.5 .2 13 13.5]); 

figure(2)
Type = ones(1,max(boutons));
for j = 2:4
    Type(boutons(j):boutons(j+1)) = j;
end
imagesc(Type'); colormap(cmap)
xticks([]); yticks([])
set(gcf,'units','centimeters','position',[20 5 .5 15]); 

figure(3)
imagesc(Mid); colormap('jet')
xticks([]); yticks([])
set(gcf,'units','centimeters','position',[25 5 .5 15]); 

%% supp fig 3c - sptial type quantification

% AP and ML distribution of FF vs Bar

clear B_AP
for i = 1:6
    AP = XY(FOV == i,1)'; FF = newFF(FOV == i);
    AP = discretize(AP,[0:60:550]);
    for j = 1:max(AP)
        boutons = sum(AP' == j & sum(FF' == [1,3],2)>0);% boutons = sum(AP == j);
        for k = [1,3]
            f = FF == k;
            B_AP(j,k,i) = sum(f(AP == j))./boutons;
        end
    end
end

% same for ML
clear B_ML
for i = 1:6
    ML = XY(FOV == i,1)'; FF = newFF(FOV == i);
    ML = discretize(ML,[0:60:400]);
    for j = 1:max(ML)
        boutons = sum(ML' == j & sum(FF' == [1,3],2)>0);% boutons = sum(ML == j);
        for k = [1,3]
            f = FF == k;
            B_ML(j,k,i) = sum(f(ML == j))./boutons;
        end
    end
end

figure(1); clf
for i = 1:6
    for k = [1,3]
        plot([30:60:520]*(5/3),B_AP(:,k,i),'Color',cmap(k,:)); hold on
    end
end
xlim([30 520]*(5/3))
box off; axis square
xlabel('dist. from anterior edge (\mum)')
ylabel('fraction of boutons')
set(gcf,'units','centimeters','position',[10 10 4 4]); 

figure(2);
for i = 1:6
    for k = [1,3]
        plot([30:60:370]*(5/3),B_ML(:,k,i),'Color',cmap(k,:)); hold on
    end
end
xlim([30 370]*(5/3))
box off; axis square
xlabel('dist. from medial edge (\mum)')
ylabel('fraction of boutons')
set(gcf,'units','centimeters','position',[10 10 4 4]); 

%% supp fig 3d -  Allen data

slices = [82:120];
x = [270,325]; y = [280, 325]; z = [6,9];

% vglut2 ipsi and contra transverse AVERAGE
clear Kmod Fmod F_stack K_stack
% make path to these files!
% F = 'foxp2_avg_rottransverse.tif';
% K = 'kcng4_avg_rottransverse.tif';
t = 1;
for i = slices
    F_stack(:,:,t) = rgb2gray(imread(F,i));
    K_stack(:,:,t) = rgb2gray(imread(K,i));
    t = t+1;
end

Kmod = mean(mat2gray(K_stack(x(1):x(2),y(1):y(2),z(1):z(2))),3)'; AP_allen = size(Kmod,2)*25; ML_allen = size(Kmod,1)*25; % 25um per pix
Fmod = mean(mat2gray(F_stack(x(1):x(2),y(1):y(2),z(1):z(2))),3)'; 
Gmod = zeros(size(Kmod));
RGBmod = cat(3,Gmod,Fmod,Kmod);
figure(2); imshow(RGBmod); hold on
scalebarLength = 10;
plot([2 2+scalebarLength],[size(RGBmod,1)-2 size(RGBmod,1)-2],'-w','Linewidth',1); text(5,size(RGBmod,1)-5,'250\mum','Color','w')
set(gcf,'units','centimeters','position',[10 10 4 4]); 
set(gca,'units','centimeters','position',[0 0 4 4]); 

% AP averaging
clear AP_*
AP_f = squeeze(nanmean(Fmod,1)); AP_k = squeeze(nanmean(Kmod,1)); 
bins = [0:4:size(AP_f,2)]; if max(bins) < size(AP_f,2); bins(end+1) = size(AP_f,2); end
plotbins = bins(1:end-1)+(diff(bins)/2);
AP = [1:size(AP_f,2)]; AP = discretize(AP,bins);
for j = 1:max(AP)
    AP_foxp2(j) = nanmean(AP_f(AP == j)); 
    AP_kcng(j) = nanmean(AP_k(AP == j)); 
end
AP_foxp2 = mat2gray(AP_foxp2);
AP_kcng = mat2gray(AP_kcng);

% ML averaging
clear ML_*
ML_f = squeeze(nanmean(Fmod,2)); ML_k = squeeze(nanmean(Kmod,2)); 
binsML = [0:4:size(ML_f,1)]; if max(bins) < size(ML_f,1); bins(end+1) = size(ML_f,1); end
plotbinsML = binsML(1:end-1)+(diff(binsML)/2);
ML = [1:size(ML_f,1)]; ML = discretize(ML,binsML);
for j = 1:max(ML)
    ML_foxp(j) = nanmean(ML_f(ML == j)); 
    ML_kcng(j) = nanmean(ML_k(ML == j)); 
end
ML_foxp = mat2gray(ML_foxp);
ML_kcng = mat2gray(ML_kcng);

figure(1); clf % AP average
plot(plotbins*25,AP_foxp2,'Color',ColourMap.foxp2); hold on
plot(plotbins*25,AP_kcng,'Color',ColourMap.alpha); hold on
xlim([min(plotbins) max(plotbins)]*25)
yticks([0 1])
box off; axis square
xlabel('dist. anterior edge (\mum)')
ylabel('norm. fluorescence (a.u.)')
set(gcf,'units','centimeters','position',[10 10 4 4]); 

figure(2);clf % ML average
plot(plotbinsML*25,ML_foxp,'Color',ColourMap.foxp2); hold on
plot(plotbinsML*25,ML_kcng,'Color',ColourMap.alpha); hold on
xlim([min(plotbinsML) max(plotbinsML)]*25)
yticks([0 1])
box off; axis square
xlabel('dist. medial edge (\mum)')
set(gcf,'units','centimeters','position',[10 10 4 4]); 

%% fig 3e and 3f1 scatter maps of retinotopy and on off

siz = .2;
keeps = (newFF == 3 | newFF == 2)' & HorzPref > 0 & HorzPref < 9;
figure(1);
scatter(XY(~keeps,1),XY(~keeps,2),.15,[.5 .5 .5],'filled','MarkerFaceAlpha',.2); hold on
H = discretize(HorzPref,[0:.5:9]); cmap = jet(max(H)-1);
for j = 1:length(cmap)
    scatter(XY(keeps & H == j,1),XY(keeps & H == j,2),siz,cmap(j,:),'filled','MarkerEdgeColor',[.7 .7 .7],'LineWidth',.05);
end
set(gca,'YDir','reverse'); daspect([1 1 1]); set(gca,'visible','off'); 
set(gcf,'units','centimeters','position',[10 10 4.5 3]); 

keeps = (newFF == 3 | newFF == 2)' & VertPref > 0 & VertPref < 9;
figure(2);
scatter(XY(~keeps,1),XY(~keeps,2),.15,[.5 .5 .5],'filled','MarkerFaceAlpha',.2); hold on
H = discretize(VertPref,[0:.5:9]); cmap = jet(max(H)-1);
for j = 1:length(cmap)
    scatter(XY(keeps & H == j,1),XY(keeps & H == j,2),siz,cmap(j,:),'filled','MarkerEdgeColor',[.7 .7 .7],'LineWidth',.05);
end
set(gca,'YDir','reverse'); daspect([1 1 1]); set(gca,'visible','off');
set(gcf,'units','centimeters','position',[20 10 4.5 3]); 

figure(6); % this is a colormap!!
imagesc([.5:.5:9]'); colormap(jet); ylabel('\DeltaF/F')
set(gca,'xtick',[])
set(gcf,'units','centimeters','position',[20 10 .1 2]); 

for i = 1:2
    figure(2+i);
    clear cmap
    edges = [-.2:.025:.2]; ref = edges(1:end-1)+.0125;
    bwr = [0 0 1; 1 1 1; 1 0 0];
    for j = 1:3
        cmap(:,j) = cat(2,interp1([1:2],bwr(1:2,j),[1:2/length(edges):2]),interp1([1:2],bwr(2:3,j),[1:2/length(edges):2]));
    end

    A = OnOff(:,i);
    A(A >= .2) = .2; A(A <= -.2) = -.2; 
    A = discretize(A,edges);
    for j = 1:max(A)
        scatter(XY(A == j,1),XY(A == j,2),siz,cmap(j+1,:),'filled'); hold on
    end

    set(gca,'YDir','reverse'); daspect([1 1 1]); set(gca,'visible','off');
    set(gcf,'units','centimeters','position',[10 10 4.5 3]); 
end

figure(5); % this is a colormap!!
imagesc(ref'); colormap(bluewhitered); ylabel('\DeltaF/F')
set(gca,'xtick',[])
yticks([1 length(ref)])
yticklabels({'-.2', '.2'})
set(gcf,'units','centimeters','position',[20 10 .1 2]); 

%% supp fig s3f2 what is fraction of ON and OFF types

OO = OOidx; OO(OO<-1) = -1; OO(OO>1) = 1;
bins = [-1,-.9:.3:.9,1]; plotbins = [-.95,-.75:.3:.9,.95];
for i = 1:6
    figure(1); h = histogram(OO(newFF == 1 & FOV' == i),bins,'Normalization','probability');
    figure(2); plot(plotbins,h.Values,'-k'); hold on
end
axis square; box off
xlabel('On-Off index')
ylabel('bouton fraction')
set(gcf,'units','centimeters','position',[10 10 4 4]); 

