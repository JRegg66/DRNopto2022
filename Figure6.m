%% DRN axon effect

% plotting script
% set default plot settings
set(groot,...
    'DefaultFigureRenderer','painters',...
    'DefaultFigureColor','w',...
    'DefaultAxesLineWidth',1.5,...
    'DefaultAxesXColor','k',...
    'DefaultAxesYColor','k',...
    'DefaultAxesFontUnits','points',...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesBox','off',...
    'DefaultAxesTickDir','in',...
    'DefaultAxesTickDirMode','manual',...
    'DefaultLineLineWidth',1,...
    'DefaultHistogramLineWidth',1.5,...
    'DefaultTextFontUnits','Points',...
    'DefaultTextFontSize',13,...
    'DefaultTextFontName','Helvetica');

%% ----- WHOLE LGN DATA
%% load data
% load WholeLGN_Data_Nov2022.mat
PlotColours

%% pull out data

PlotColours
clearvars -except data ColourMap
fr = 15.5; QI_thresh = .15;
QIbar = []; QIff = []; SBC = []; Allff = []; QI = []; Base = []; BigRespAll = []; XY = []; FOV = []; DRN = []; DRN2 = [];
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
    QI = cat(1,QI,cat(2,OptoQI(:,1:2),...
        max(data(i).Ret.QI(:,1:8),[],2),max(data(i).Ret.QI(:,9:16),[],2),data(i).Ret.QI(:,17:18)));
    XY = cat(1,XY,data(i).ROIxyReg);
    FOV = cat(1,FOV, ones(boutons,1).*i);
    DRN = cat(2,DRN,data(i).DRNintensity);
    DRN2 = cat(2,DRN2,data(i).DRNintensityROI);
end

toss = ~((QIbar == 2 | QIff == 1)); % toss = toss | SBC;
SBC(toss) = []; Allff(toss) = []; QIbar(toss) = []; QIff(toss) = []; QI(toss,:) = []; BigRespAllDR(:,toss) = []; OOidx(QIff == 0) = NaN; OOidx(toss) = [];
Base(toss) = []; BigRespAll(:,toss) = []; XY(toss,:) = []; FOV(toss) = []; BaseDR(toss) = []; BaseOptoDR(toss) = [];
HorzPref(toss) = []; VertPref(toss) = []; OnOff(toss,:) = []; RetResp(toss,:) = []; DRN(:,toss) = []; DRN2(:,toss) = [];
clearvars -except data SBC QIbar QIff QI Base Allff BigRespAll XY ColourMap ...
    FOV BaseOptoDR DRN  DRN2 BaseDR BigRespAllDR OOidx OnOff HorzPref VertPref RetResp
normBase = BaseOptoDR - BaseDR;

% define the 4 categories
newFF = NaN(size(Allff)); 
newFF(QIff == 1 & QIbar < 2) = 1; % FF
newFF(QIff == 1 & QIbar == 2) = 2; % both
newFF(QIff == 0 & QIbar == 2) = 3; % bar
newFF(SBC == 1) = 4; % SBC

titles = {'FF','FF+bar','bar','SBC'};
cmap = cat(1,ColourMap.FF,ColourMap.Both,ColourMap.Bar,ColourMap.SBC);
DRN_colour = [148 26 28]./255;

%% fig 6b look at density effects for all RGCs

i = 3; % 10um circle size
figure;
bins = [0:.1:.55];
F = unique(FOV);
for k = 1:length(F)
    for b = 2:length(bins)
        h(b-1,k) = mean(normBase(FOV == F(k) & DRN(i,:)' > bins(b-1) & DRN(i,:)' < bins(b)));
    end
    plot(bins(1:end-1)+.05,h(:,k),'-','Color',DRN_colour+.1,'LineWidth',.5); hold on
end
plot([0 1],[0 0],'--k')
plot(bins(1:end-1)+.05, nanmean(h,2),'Color',DRN_colour,'LineWidth',2);
clear l h
xlim([0 .5]); ylim([-.2 .05]);
yticks([-.2 -.1 0]); xticks([0 .25 .5])
xlabel('DRN axon density'); ylabel('\DeltaF/F')
axis square; box off
set(gcf,'units','centimeters','position',[10 10 4 4]); 

%% figure 6c per bins of DRN axon density FF vs Bar

i = 3; % 10um ring size
j = [1,3];

bins = [0:.1:.55];
F = unique(FOV);
clear h
figure(1); hold on
for k = 1:length(F)
    for b = 2:length(bins)
        h(b-1,k) = mean(normBase(newFF' == 1 & FOV == F(k) & DRN(i,:)' > bins(b-1) & DRN(i,:)' < bins(b)));
    end
    plot(bins(1:end-1)+.05,h(:,k),'-','Color',cmap(1,:),'LineWidth',.5);
end
plot([0 1],[0 0],'--r')
plot(bins(1:end-1)+.05, nanmean(h,2),'Color',cmap(1,:),'LineWidth',2);
clear l h

for k = 1:length(F)
    for b = 2:length(bins)
        h(b-1,k) = mean(normBase(newFF' == 3 & FOV == F(k) & DRN(i,:)' > bins(b-1) & DRN(i,:)' < bins(b)));
    end
    plot(bins(1:end-1)+.05,h(:,k),'-','Color',cmap(3,:),'LineWidth',.5);
end
plot([0 1],[0 0],'--k')
plot(bins(1:end-1)+.05, nanmean(h,2),'Color',cmap(3,:),'LineWidth',2);
clear l h
xlim([0 .5]); ylim([-.2 .05]); axis square
yticks([-.2 -.1 0]); xticks([0 .25 .5])
xlabel('DRN density'); ylabel('\DeltaF/F')
set(gcf,'units','centimeters','position',[10 10 4 4]); 

%% figure 6d difference in baseline suppression for FF and Bar for all boutons and DRN matched

c = hsv(length(F));
for j = 1:length(F)
    % FF opto
    VisIdx = newFF == 1 & FOV' == F(j); FFnum(j) = sum(VisIdx)/sum(FOV == F(j));
    FF(j) = mean(normBase(VisIdx));
    VisIdx = newFF == 1 & FOV' == F(j) & DRN(3,:) > .25 & DRN(3,:) < .4;  FFdrn(j) = sum(VisIdx)/sum(FOV == F(j));
    FFD(j) = mean(normBase(VisIdx));
    % bar opto
    VisIdx = newFF == 3 & FOV' == F(j); Bnum(j) = sum(VisIdx)/sum(FOV == F(j));
    Bar(j) = mean(normBase(VisIdx));
    VisIdx = newFF == 3 & FOV' == F(j) & DRN(3,:) > .25 & DRN(3,:) < .4;  Bardrn(j) = sum(VisIdx)/sum(FOV == F(j));
    BarD(j) = mean(normBase(VisIdx));
end

All = FF - Bar;
Matched = FFD - BarD;
figure(5);
scatter(All,Matched,5,DRN_colour,'filled'); hold on
plot([-.2 0],[-.2 0],'--k');
xlim([-.2 0]); ylim([-.2 0]); axis square
set(gca, 'XAxisLocation', 'top'); set(gca, 'YAxisLocation', 'right');
xlabel({'FF - Bar', 'all'}); ylabel({'FF - Bar', 'DRN density matched'}); 
set(gcf,'units','centimeters','position',[10 10 4 4]); 

%% ---- SUPPLEMENTARY PANELS

%% supp fig s6a DRN density coloured scatter plot

clear cmap2
edges = [0:.05:.5]; ref = edges(1:end-1)+.025;
cmap2 = [[0:.1:1]',zeros(11,1),zeros(11,1)];

A = DRN(3,:); A(A >= .5) = .5;
A = discretize(A,edges);
figure(1);
for j = 1:max(A)
    scatter(XY(A == j,1),XY(A == j,2),.5,cmap2(j+1,:),'filled'); hold on
end
    
figure(1); set(gca,'YDir','reverse'); daspect([1 1 1]); xlim([0 600]); ylim([0 400]); set(gca,'visible','off');
set(gcf,'units','centimeters','position',[10 10 4.5 3]); 
set(gca,'units','centimeters','position',[.1 .1 4.3 2.8]); 

% just the example FOV, figure s6a1
figure(2);
f = 1;
for j = 1:max(A)
    scatter(XY(A == j & FOV' == f,1),XY(A == j & FOV' == f,2),.5,cmap2(j+1,:),'filled'); hold on
end
    
figure(2); set(gca,'YDir','reverse'); daspect([1 1 1]); xlim([0 600]); ylim([0 400]); set(gca,'visible','off');
set(gcf,'units','centimeters','position',[15 10 4.5 3]); 
set(gca,'units','centimeters','position',[.1 .1 4.3 2.8]);

figure(3); % this is a colormap!!
ref = [0:.1:1]; ref = ref(1:end-1)+.05;
imagesc(flipud(ref')); colormap(cmap2); ylabel('DRN density')
set(gca,'xtick',[])
yticks([1 length(ref)])
yticklabels({'.5', '0'})
set(gcf,'units','centimeters','position',[20 10 .5 2]); 
set(gcf,'units','centimeters','position',[20 10 .1 2]); 

%% supp figure s6b DRN density along A-P axis and M-L axis

F = FOV; cc = hsv(6);
clear B_AP
for i = 1:max(F)
    AP = XY(F == i,1)'; B = DRN(3,F == i);
    AP = discretize(AP,[0:60:550]);
    for j = 1:max(AP)
        B_AP(j,i) = nanmean(B(AP == j));
    end
end

figure(1); clf
for i = 1:6
    plot([30:60:520]*(5/3),B_AP(:,i),'Color',cc(i,:)); hold on
end
xlim([30 520]*(5/3)); ylim([0 .35])
yticks([.1 .3])
box off
xlabel('\mum')
ylabel('DRN density')
set(gcf,'units','centimeters','position',[10 10 4 4]); 

% and M-L axis
clear B_ML
for i = 1:max(F)
    ML = XY(F == i,2)'; B = DRN(3,F == i);
    ML = discretize(ML,[0:60:400]);
    for j = 1:max(ML)
        B_ML(j,i) = nanmean(B(ML == j));
    end
end

figure(2); clf
ypos = .15;
for i = 1:6
    plot([30:60:370]*(5/3),B_ML(:,i),'Color',cc(i,:)); hold on
    plot([100 150],[ypos ypos],'-','Color',cc(i,:));
    text(170,ypos,strcat('M',num2str(i)));
    ypos = ypos-.015;
end
xlim([30 370]*(5/3)); ylim([0 .35])
yticks([.1 .3])
box off
xlabel('\mum')
ylabel('DRN density')
set(gcf,'units','centimeters','position',[10 10 4 4]); 

%% ----- ZOOMED IN DATASET

clear
% load Zoom_2p_Data_Nov2022
PlotColours

%% pull out data

clearvars -except data ColourMap
fr = 15.5; QI_thresh = .15; 
BigRespAll = []; BigStd = []; QIbar = []; QIff = []; Allff = []; SBC = []; Reference = []; QI = [];  FOV = [];
PupilTrials = []; Base = []; ArousalCoeff = []; BigRespTrials = []; PupilAll = []; ArousalCoefR= [];
HighLowPrct = []; HighLow = []; DRN = [];

for i = 1:length(data)
    boutons = size(data(i).OptoMean,2);
    FOV = cat(1,FOV,ones(boutons,1)*i);
    clear OptoTraces OptoMean
    OptoQI = nan(boutons,length(data(i).OptoTraces));
    OptoMean = nan(204,boutons,5); OptoStd = nan(204,boutons,5);
    OptoTraces = nan(204,40,boutons,5); Pupil = nan(40,5);
    for k = 1:length(data(i).OptoTraces)
        open = find(data(i).Arousal.Opto.FractClosedEye{k} < .15);
        OptoTraces(:,1:length(open),:,k) = data(i).OptoTraces{k}(:,open,:);
        Pupil(1:length(open),k) = data(i).Arousal.Opto.PupTrial{k}(open);
        for b = 1:boutons
            OptoTraces(:,:,b,k) = smoothdata(OptoTraces(:,:,b,k),1,'movmean',7);
        end
        OptoMean(:,:,k) = squeeze(nanmean(OptoTraces(:,:,:,k),2));
        OptoStd(:,:,k) = squeeze(nanstd(OptoTraces(:,:,:,k),[],2))./sqrt(size(OptoTraces(:,:,:,k),2));
        meanresp = squeeze(nanmean(OptoTraces(round(3*fr):end,:,:,k),2));
        varmeanresp = squeeze(var(meanresp,[],1));
        varresp = squeeze(var(OptoTraces(round(3*fr):end,:,:,k),[],1,'omitnan'));
        meanvarresp = squeeze(nanmean(varresp,1));
        OptoQI(:,k) = varmeanresp./meanvarresp;
        clear meanresp varmeanresp varresp meanvarresp
    end
    Base = cat(2,Base,squeeze(mean(mean(OptoMean(round(4*fr):round(5*fr),:,[2,4,5]),1),3)));
    Resp = (mean(OptoMean(1:round(2*fr),:,1),1));
    xstart = [4,5.1,6,7.1,8,9.1,10,11.1]; xend = [5,6,7,8,9,10,11,12];
    for k = 1:2
        for j = 1:length(xstart) % ctrl lum
            Resp = cat(1,Resp,(mean(OptoMean(round(xstart(j)*fr):round(xend(j)*fr),:,0+k),1)));
        end
        for j = 2:length(xstart)-1 % ctrl bar
            Resp = cat(1,Resp,(mean(OptoMean(round(xstart(j)*fr):round(xend(j)*fr),:,2+k),1)));
        end
    end
    % 20th and 80th percentile for normalizing
    perc = [20,80]; TotResp = [];
    xstart = [4,5.1,6,7.1,8,9.1,10,11.1]; xend = [5,6,7,8,9,10,11,12];
    for p = 1:length(perc)
        TotResp = cat(1,TotResp,(prctile(OptoMean(1:round(2*fr),:,1),perc(p))));
        for k = 1 % if you want to include opto trials, make this 1:2
            for j = 1:length(xstart) % ctrl lum
                TotResp = cat(1,TotResp,(prctile(OptoMean(round(xstart(j)*fr):round(xend(j)*fr),:,0+k),perc(p))));
            end
            for j = 2:length(xstart)-1 % ctrl bar
                TotResp = cat(1,TotResp,(prctile(OptoMean(round(xstart(j)*fr):round(xend(j)*fr),:,2+k),perc(p))));
            end
        end
        TotResp = cat(1,TotResp,permute(squeeze(prctile(data(i).RetMean(15:30,:,:),perc(p))),[2,1]),...
            permute(squeeze(prctile(data(i).RetMean(31:46,:,:),perc(p))),[2,1]),...
            permute(reshape(squeeze(prctile(data(i).OriMean(15:30,:,:,:),perc(p))),boutons,24),[2,1]),...
            permute(reshape(squeeze(prctile(data(i).OriMean(31:46,:,:,:),perc(p))),boutons,24),[2,1]));
    end
    TR = TotResp; minR = min(TotResp,[],1); TR = TR - minR; maxR = max(TR,[],1); clear TR
    % correlation coefficient for arousal modulation on ctrl 
    P = reshape(Pupil,200,1); P = repmat(P,1,boutons); Pt = nan(200,boutons); Pt(1:size(P,1),:) = P;
    PupilTrials = cat(2,PupilTrials,Pt);
    t = 1; AC = NaN(14,boutons); Pt = NaN(14,40,boutons); Rt = NaN(14,40,boutons); ACsig = NaN(14,boutons);
    HminL = NaN(14,boutons); HminL_prct = NaN(14,boutons);
    for k = [1,3]
        for j = 2:length(xstart)
            clear R
            R = squeeze(mean(OptoTraces(round(xstart(j)*fr):round(xend(j)*fr),:,:,k),1)); 
            R = R - minR; R = R./maxR;
            for b = 1:boutons
                [c,d] = corrcoef(R(:,b),Pupil(:,k));
                AC(t,b) = c(2,1);
                ACsig(t,b) = d(2,1);
                HminL(t,b) = nanmean(R(Pupil(:,k)>.5,b)) - nanmean(R(Pupil(:,k)<.5,b));
                HminL_prct(t,b) = nanmean(R(Pupil(:,k)>prctile(Pupil(:,k),50),b)) - nanmean(R(Pupil(:,k)<prctile(Pupil(:,k),50),b));
            end
            Pt(t,1:length(Pupil(:,k)),:) = repmat(Pupil(:,k),[1,boutons]);
            Rt(t,1:length(Pupil(:,k)),:) = R;
            t = t+1;
        end
    end
    AC(t-1,:) = []; Rt(t-1,:,:) = []; Pt(t-1,:,:) = []; ACsig(t-1,:) = []; HminL(t-1,:) = []; HminL_prct(t-1,:) = [];
    ArousalCoeff = cat(2,ArousalCoeff,AC); ArousalCoefR = cat(2,ArousalCoefR,ACsig); 
    HighLow = cat(2,HighLow,HminL); HighLowPrct = cat(2,HighLowPrct,HminL_prct);
    % response matrix
    BigStd = cat(2,BigStd,reshape(permute(OptoStd,[1,3,2]),204*5,boutons));
    BigResp = cat(1,reshape(permute(OptoMean,[1,3,2]),204*5,boutons),...
            reshape(permute(data(i).RetMean(1:62,:,:),[1,3,2]),62*18,boutons),...
            reshape(permute(data(i).OriMean(1:62,:,:,:),[1,3,4,2]),62*24,boutons));
    %normalize
    RespTrials = reshape(permute(OptoTraces,[1,4,2,3]),1020,40,boutons);
    RespTrials = RespTrials - permute(min(TotResp,[],1),[1,3,2]); 
    Resp = Resp - min(TotResp,[],1); TotResp = TotResp - min(TotResp,[],1);
    RespTrials = RespTrials./permute(max(TotResp,[],1),[1,3,2]);
    Resp = Resp./max(TotResp,[],1); TotResp = TotResp./max(TotResp,[],1);
    Reference = cat(2,Reference,Resp);
    BigRespAll = cat(2,BigRespAll,BigResp); clear TotResp
    BigRespTrials = cat(3,BigRespTrials,RespTrials); PupilAll = cat(3,PupilAll,repmat(Pupil,[1,1,boutons]));
    QIbar = cat(1,QIbar,sum([max(data(i).Ret.QI(:,1:8),[],2) >= QI_thresh,max(data(i).Ret.QI(:,9:16),[],2) >= QI_thresh],2));
    QIff = cat(1,QIff,max(data(i).Ret.QI(:,17:18),[],2) >= QI_thresh);
    Ret = cat(2,squeeze(mean(data(i).RetMean(16:47,:,1:16),1)),reshape(squeeze(mean(data(i).OriMean(16:47,:,:,:),1)),boutons,24));
    bad = cat(2,data(i).Ret.QI(:,1:16),reshape(data(i).Ori.QI,boutons,24)) < QI_thresh;
    Ret(bad) = NaN;
    S = cat(2,sum(Ret(:,1:8) < 0 & ~isnan(Ret(:,1:8)),2)>0 & sum(Ret > 0 & ~isnan(Ret),2) == 0,...
        sum(Ret(:,9:16) < 0 & ~isnan(Ret(:,9:16)),2)>0 & sum(Ret > 0 & ~isnan(Ret),2) == 0,...
        sum(Ret(:,17:end) < 0 & ~isnan(Ret(:,17:end)),2)>0 & sum(Ret > 0 & ~isnan(Ret),2) == 0);
    SBC = cat(1,SBC,sum(S,2)>=2);
    QI = cat(1,QI,cat(2,OptoQI(:,1:4),...
        max(data(i).Ret.QI(:,1:8),[],2),max(data(i).Ret.QI(:,9:16),[],2),data(i).Ret.QI(:,17:18),...
        squeeze(max(data(i).Ori.QI,[],2))));
    DRN = cat(2,DRN, data(i).DRNintensity_new);
end

optoQIff = max(QI(:,1:2),[],2); optoQIbar = max(QI(:,3:4),[],2);
toss = ~((QIbar == 2 | QIff == 1) & (optoQIff >= QI_thresh | optoQIbar >= QI_thresh));
SBC(toss) = []; BigRespAll(:,toss) = []; BigStd(:,toss) = []; Reference(:,toss) = []; QIbar(toss) = []; QIff(toss) = []; QI(toss,:) = []; 
FOV2 = FOV; FOV(toss) = [];
optoQIff(toss) = []; optoQIbar(toss) = []; Base(:,toss) = []; PupilTrials(:,toss) = []; 
ArousalCoeff(:,toss) = []; BigRespTrials(:,:,toss) = []; PupilAll(:,:,toss) = []; ArousalCoefR(:,toss) = [];
HighLowPrct(:,toss) = []; HighLow(:,toss) = []; DRN(:,toss) = [];
clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials Base ...
    dataRaw ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF DRN toss FOV2

% define the 4 categories
QI_thresh = .15;
newFF = NaN(size(SBC)); 
newFF(QIff == 1 & QIbar < 2 & optoQIff >= QI_thresh & optoQIbar < QI_thresh) = 1; % FF
newFF(QIff == 1 & QIbar == 2 & optoQIff >= QI_thresh & optoQIbar >= QI_thresh) = 2; % both
newFF(QIff == 0 & QIbar == 2 & optoQIff < QI_thresh & optoQIbar >= QI_thresh) = 3; % bar
newFF(SBC == 1) = 4; % SBC

M{1} = [1,2,16]; % JR220
M{2} = [3,4,15];% JR227
M{3} = [5];% JR228
M{4} = [6,14];% JR230
M{5} = [7,13];% JR218
M{6} = [8,10];% JR219
M{7} = [9,11];% JR221
M{8} = [12,17];% JR212
M{9} = [18,19,20];% JR236

titles = {'FF','both','bar','SBC'};
cmap = cat(1,ColourMap.FF,ColourMap.Both,ColourMap.Bar,ColourMap.SBC);
DRN_colour = [148 26 28]./255;
DRN_match = [.15 .3];

%% fig 6F, DRN axon intensity for the 4 categories
sizes = 10;
F = unique(FOV);
for i = 1; 
    clear b
    for j = 1:4
        for k = 1:length(F)
            b(k,j) = mean(DRN(i,newFF == j & FOV == F(k)));
        end
        figure(3); % subplot(1,4,i); 
        bar(j,nanmean(b(:,j)),'FaceColor',cmap(j,:)); hold on
        scatter(ones(1,size(b,1)).*j,b(:,j),'.k');
    end
    figure(3); % subplot(1,4,i);
    plot(b','-k');
    title(strcat('radius ',num2str(sizes(i)),'um'));
    figure(4); % subplot(1,4,i);
    scatter(b(:,1),b(:,3),5,DRN_colour,'filled'); hold on
    plot([0 .4],[0 .4],'--k')
    xlim([0 .4]); ylim([0 .4]); axis square
end

figure(3); ylabel('DRN axon intensity')
figure(4); xlabel('FF DRN axon density'); ylabel('Bar DRN axon density'); box off
xticks([0 .2 .4]); yticks([0 .2 .4]);
set(gcf,'units','centimeters','position',[10 10 4 4]); 

%% fig 6G-H, FF vs Bar suppression, all boutons and DRN density matched boutons
% Only use FOV with sufficient FF and Bar boutons...

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials Base ...
    dataRaw ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF DRN toss FOV2 newFF cmap DRN_colour DRN_match

ctrl = [2:15]; opto = [16:29]; % all steps of the response
DR = [.5 .7];

c = hsv(length(M)); t = 1;
for j = 1:length(M)
    F = FOV == M{j};
    for i = 1:size(F,2)
        % FF opto
        VisIdx = newFF == 1 & F(:,i); 
        FFnum(t) = sum(VisIdx)/sum(F(:,i)); 
        C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
        O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); O(isnan(O)) = [];
        FF(t,1) = mean(O);
        VisIdx = newFF == 1 & F(:,i) & DRN(1,:)' > DRN_match(1) & DRN(1,:)' < DRN_match(2); 
        FFdrn(t) = sum(VisIdx)/sum(F(:,i));
        C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
        O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); O(isnan(O)) = [];
        FFD(t,1) = mean(O);
        % bar opto
        VisIdx = newFF == 3 & F(:,i); 
        Bnum(t) = sum(VisIdx)/sum(F(:,i));
        C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
        O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); O(isnan(O)) = [];
        Bar(t,1) = mean(O);
        VisIdx = newFF == 3 & F(:,i) & DRN(1,:)' > DRN_match(1) & DRN(1,:)' < DRN_match(2); 
        Bardrn(t) = sum(VisIdx)/sum(F(:,i));
        C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
        O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); O(isnan(O)) = [];
        BarD(t,1) = mean(O);
        % toss FOV with too few FF or bar
        whichM(t) = j;
        if FFnum(t) < .05 | Bnum(t) < .05; Bar(t,:) = NaN; FF(t,:) = NaN; BarD(t,:) = NaN; FFD(t,:) = NaN; end
        if FFdrn(t) < .05 | Bardrn(t) < .05; BarD(t,:) = NaN; FFD(t,:) = NaN; end
        t = t+1;
    end
end

figure(4); clf % fig 6G
for j = 1:length(M)
    B_m(j,:) = nanmean(BarD(whichM == j,:),1);
    F_m(j,:) = nanmean(FFD(whichM == j,:),1);
end
Matched = F_m - B_m;
bar(1,nanmean(F_m(:,1)),'FaceColor',cmap(1,:),'EdgeColor','none'); hold on
bar(2,nanmean(B_m(:,1)),'FaceColor',cmap(3,:),'EdgeColor','none');
for i = 1:length(M)
    plot([1,2],cat(2,F_m(i,1),B_m(i,1)),'-k');
end
[h,p] = ttest2(F_m,B_m);
box off
plot([1 2],[.05 .05],'-k','LineWidth',.5); text(1.4,.06,'*');
xticks([1,2]); xticklabels({'FF','bar'})
ylim([-.3 .07]); yticks([-.2 0]); yticklabels({'-.2','0'})
xlim([.2 2.7])
set(gca,'YColor','k')
ylabel('R_{opto} - R_{ctrl}')
set(gcf,'units','centimeters','position',[15 10 2.5 4]); 

for j = 1:length(M)
    B_m(j,:) = nanmean(Bar(whichM == j,:),1);
    F_m(j,:) = nanmean(FF(whichM == j,:),1);
end
All = F_m - B_m;

figure(5); % fig 6H
scatter(All,Matched,5,DRN_colour,'filled'); hold on
plot([-.1 .05],[-.1 .05],'--k');
xlim([-.1 .05]); ylim([-.1 .05]); axis square
xlabel({'FF - Bar'}); ylabel({'FF - Bar'}); 
set(gca, 'XAxisLocation', 'top'); set(gca, 'YAxisLocation', 'right');
set(gcf,'units','centimeters','position',[10 10 4 4]); 

%% supp fig S1C and S1D, response vs suppression for FF and bar

ctrl = [2:15]; opto = [16:29]; % all steps of the response

% fig 4f cdf of suppression between .5 and .7 for FF, both, bar and SBC
DR = [.5 .7]; bins = [-.6:.05:.2]; plotbins = [-.55:.05:.2];
C = Reference(ctrl,:); O = Reference(opto,:)-Reference(ctrl,:);
O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); %O(isnan(O)) = [];

i = 1; % fig S1C1
figure(1); clf 
scatter(DRN(i,newFF == 1),O(newFF == 1),1,cmap(1,:),'filled'); hold on
plot([0 1],[0 0],'--k')
xlim([0 .5]); ylim([-.8 .3]);
xticks([0 .25 .5]); yticks([-.6 -.3 0])
xlabel('DRN density'); ylabel('R_{opto} - R_{ctrl}')
axis square; box off
set(gcf,'units','centimeters','position',[10 10 4 4]); 

figure(2); % fig S1C2
scatter(DRN(i,newFF == 3),O(newFF == 3),1,cmap(3,:),'filled'); hold on
plot([0 1],[0 0],'--k')
xlim([0 .5]); ylim([-.8 .3]);
xticks([0 .25 .5]); yticks([-.6 -.3 0])
xlabel('DRN density'); ylabel('R_{opto} - R_{ctrl}')
axis square; box off
set(gcf,'units','centimeters','position',[15 10 4 4]); 

% per mouse averaged lines, fig S1D
clear h
bins = [0:.1:.5];
figure(3); hold on
patch([DRN_match(1) DRN_match(2) DRN_match(2) DRN_match(1)],[-.4 -.4 .1 .1],[.7 .7 .7],'EdgeColor','none'); hold on
for k = 1:length(M)
    Fs = sum(FOV == M{k},2)>0;
    for b = 2:length(bins)
        h(b-1,k) = nanmean(O(newFF == 1 & Fs & DRN(i,:)' > bins(b-1) & DRN(i,:)' < bins(b)));
    end
    plot(bins(1:end-1)+.05,h(:,k),'-','Color',cmap(1,:),'LineWidth',.5);
end
plot(bins(1:end-1)+.05, nanmean(h,2),'Color',cmap(1,:),'LineWidth',2);
clear l h

for k = 1:length(M)
    Fs = sum(FOV == M{k},2)>0;
    for b = 2:length(bins)
        h(b-1,k) = nanmean(O(newFF == 3 & Fs & DRN(i,:)' > bins(b-1) & DRN(i,:)' < bins(b)));
    end
    plot(bins(1:end-1)+.05,h(:,k),'-','Color',cmap(3,:),'LineWidth',.5);
end
plot([0 1],[0 0],'--k')
plot(bins(1:end-1)+.05, nanmean(h,2),'Color',cmap(3,:),'LineWidth',2);
clear l h
xlim([0 .5]); ylim([-.4 .1]);
xticks([0 .25 .5]); yticks([-.4 -.2 0])
xlabel('DRN density'); ylabel('R_{opto} - R_{ctrl}')
axis square; box off
set(gcf,'units','centimeters','position',[20 10 4 4]); 