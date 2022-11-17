%% Figure 4

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
% load Zoom_2p_Data_Nov2022.mat
PlotColours

%% pull out data

clearvars -except data ColourMap
fr = 15.5; QI_thresh = .15; 
BigRespAll = []; QIbar = []; QIff = []; SBC = []; Reference = []; QI = [];  FOV = [];
PupilTrials = []; BaseTrials = []; ArousalCoeff = []; BigRespTrials = []; PupilAll = [];
ASI = []; DSI = []; ASIangle = []; DSIangle = []; SFpref = [];
BaseDR = []; BaseOptoDR = [];

for i = 10:31
    if ~isempty(data(i).Opto)
    boutons = size(data(i).OptoMean,2);
    FOV = cat(1,FOV,ones(boutons,1)*i);
    clear OptoTraces OptoMean
    OptoQI = nan(boutons,length(data(i).OptoTraces));
    OptoMean = nan(204,boutons,5); OptoStd = nan(204,boutons,5);
    OptoTraces = nan(204,40,boutons,5);
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
    R = reshape(permute(OptoTracesRaw,[1,3,2,4]),204,boutons,200);
    R  =squeeze(mean(R(round(fr):2*fr,:,:),1)); R = R'; Rt = nan(200,boutons); Rt(1:size(R,1),:) = R;
    P = reshape(Pupil,200,1); P = repmat(P,1,boutons); Pt = nan(200,boutons); Pt(1:size(P,1),:) = P;
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
    ArousalCoeff = cat(2,ArousalCoeff,AC); 
    % response matrix
    BigResp = cat(1,reshape(permute(OptoMean,[1,3,2]),204*5,boutons),...
            reshape(permute(data(i).RetMean(1:62,:,:),[1,3,2]),62*18,boutons),...
            reshape(permute(data(i).OriMean(1:62,:,:,:),[1,3,4,2]),62*24,boutons));
    %normalize
    BigRespAll = cat(2,BigRespAll,BigResp);
    RespTrials = reshape(permute(OptoTraces,[1,4,2,3]),1020,40,boutons);
    BigResp = BigResp - min(TotResp,[],1); 
    RespTrials = RespTrials - permute(min(TotResp,[],1),[1,3,2]); 
    Resp = Resp - min(TotResp,[],1); bot = min(TotResp,[],1); TotResp = TotResp - min(TotResp,[],1); 
    BigResp = BigResp./max(TotResp,[],1); 
    RespTrials = RespTrials./permute(max(TotResp,[],1),[1,3,2]);
    Resp = Resp./max(TotResp,[],1); top = max(TotResp,[],1); TotResp = TotResp./max(TotResp,[],1); 
    BaseDR = cat(2,BaseDR,(squeeze(mean(mean(OptoMean(1:round(2*fr),:,:),1),3))-bot)./top);
    BaseOptoDR = cat(2,BaseOptoDR,(squeeze(mean(mean(OptoMean(round(4*fr):round(5*fr),:,[2,4,5]),1),3))-bot)./top);
    Reference = cat(2,Reference,Resp);
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
    R = cat(2,cat(1,data(i).Ret.Retinotopy.horzPref),cat(1,data(i).Ret.Retinotopy.vertPref));
    R(R(:,1) < 0 | R(:,1) > 8,1) = NaN; R(R(:,2) < 0 | R(:,2) > 8,2) = NaN;
    if data(i).monitor == 2; R(:,1) = R(:,1) + 2; end
    ASI = cat(1,ASI,data(i).Ori.gASI);
    DSI = cat(1,DSI,data(i).Ori.gDSI);
    ASIangle = cat(1,ASIangle,data(i).Ori.ASI_angle);
    DSIangle = cat(1,DSIangle,data(i).Ori.DSI_angle);
    SFpref = cat(1,SFpref,data(i).Ori.SFprefws');
    end
end

optoQIff = max(QI(:,1:2),[],2); optoQIbar = max(QI(:,3:4),[],2);
toss = ~((QIbar == 2 | QIff == 1) & (optoQIff >= QI_thresh | optoQIbar >= QI_thresh));
SBC(toss) = []; BigRespAll(:,toss) = []; Reference(:,toss) = []; QIbar(toss) = []; QIff(toss) = []; QI(toss,:) = []; 
FOV2 = FOV; FOV(toss) = []; 
optoQIff(toss) = []; optoQIbar(toss) = [];
ArousalCoeff(:,toss) = []; BigRespTrials(:,:,toss) = []; PupilAll(:,:,toss) = [];
BaseDR(toss) = []; BaseOptoDR(toss) = [];
RetPref(QIbar < 2,:) = NaN; RetPref(toss,:) = []; 
ASI(toss,:) = []; DSI(toss,:) = []; ASIangle(toss,:) = []; DSIangle(toss,:) = []; SFpref(toss) = []; 
SF = zeros(size(ASI)); % for j = 1:3; SF(SFpref == j,j) = 1; end
clearvars -except data ColourMap SBC ...
    BigRespAll Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
     ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

% define the 4 categories
QI_thresh = .15;
newFF = NaN(size(SBC)); 
newFF(QIff == 1 & QIbar < 2 & optoQIff >= QI_thresh & optoQIbar < QI_thresh) = 1; % FF
newFF(QIff == 1 & QIbar == 2 & optoQIff >= QI_thresh & optoQIbar >= QI_thresh) = 2; % both
newFF(QIff == 0 & QIbar == 2 & optoQIff < QI_thresh & optoQIbar >= QI_thresh) = 3; % bar
newFF(SBC == 1) = 4; % SBC

titles = {'FF','both','bar','SBC'};
cmap = cat(1,ColourMap.FF,ColourMap.Both,ColourMap.Bar,ColourMap.SBC);

% mouse vector
M{1} = [1,2,16]; % JR220
M{2} = [3,4,15];% JR227
M{3} = [5];% JR228
M{4} = [6,14];% JR230
M{5} = [7,13];% JR218
M{6} = [8,10];% JR219
M{7} = [9,11];% JR221
M{8} = [12,17];% JR212
M{9} = [18,19,20];% JR236
Mouse = M;

mouse = NaN(size(FOV));
for m = 1:length(M)
    mouse(sum(FOV == M{m},2)>0) = m;
end

%% fig 4c example trace (bar pref)

fr = 15.5;
Example = find(newFF == 3 & max(QI(:,3:4),[],2) > .25 & QI(:,11) >.25);

Stim = NaN(size(BigRespAll,1),1);
Lum = cat(1,zeros(ceil(5*fr),1),ones(round(2*fr),1)*-1,ones(round(2*fr),1),ones(round(2*fr),1)*-1,zeros(round(2*fr)+2,1));
Stim(1:204) = Lum; Stim(205:408) = Lum; Stim(409:1020) = 0; clear Lum
start = 1+1020+62*16; Stim(start:start+16) = 0; Stim(start+16:start+47) = -1; Stim(start+48:start+65) = 0;
start = 63+1020+62*16; Stim(start:start+16) = 0; Stim(start+16:start+47) = 1; Stim(start+48:start+65) = 0;
clear start

x = [204,62]; xpos = 0;
for j = 1:47
    if j <= 5
        xpos(j+1) = xpos(j)+x(1);
    else
        xpos(j+1) = xpos(j)+x(2);
    end
end
figure(1); clf
yyaxis left
patch([205+round(2*fr) 205+round(12*fr) 205+round(12*fr) 205+round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none'); hold on
patch([613+round(2*fr) 613+round(12*fr) 613+round(12*fr) 613+round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none');
patch([817+round(2*fr) 817+round(12*fr) 817+round(12*fr) 817+round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none');
for j = 1:length(xpos); 
    if ~(j == [1,2,3,4,5,22,23])
        patch([xpos(j)+(fr) xpos(j)+(3*fr) xpos(j)+(3*fr) xpos(j)+(fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
    end
end
plot(Stim,'-k','LineWidth',1);
box off; set(gca,'visible','off'); xlim([1 size(BigRespAll,1)]); ylim([-1.1 1.1])
set(gcf,'units','centimeters','position',[10 10 10 1]);

figure(2);
i = Example(74);
yyaxis left
patch([205+round(2*fr) 205+round(12*fr) 205+round(12*fr) 205+round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none'); hold on
patch([613+round(2*fr) 613+round(12*fr) 613+round(12*fr) 613+round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none');
patch([817+round(2*fr) 817+round(12*fr) 817+round(12*fr) 817+round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none');
for j = 1:length(xpos); 
    if j <=4
        patch([xpos(j)+(5*fr) xpos(j)+(11*fr) xpos(j)+(11*fr) xpos(j)+(5*fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
    else j >= 6
        patch([xpos(j)+(fr) xpos(j)+(3*fr) xpos(j)+(3*fr) xpos(j)+(fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
    end
end
plot(medfilt1(BigRespAll(:,i),7),'-','Color',cmap(3,:),'LineWidth',.8);
xlim([0 xpos(end)]); ylim([min(medfilt1(BigRespAll(:,i),7)) max(medfilt1(BigRespAll(:,i),7))]); box off
yticks([0 .5]); ylabel('\Delta F/F')
yyaxis right
B = medfilt1(BigRespAll(:,i),7); B = (B - min(B))./max(B - min(B));
plot(B,'-','Color',cmap(3,:),'LineWidth',.8);
xlim([0 xpos(end)]); ylim([min(B) max(B)]); box off
yticks([0 .5 1]); ylabel('dynamic range')
set(gca,'xColor','none');

set(gcf,'units','centimeters','position',[10 10 10 4]);

%% fig 4d and 5i example opto and ctrl for bar, both and FF with shading of areas chosen for response

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
     ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

lt = .5; ht = .5;
for b = 1:size(BigRespTrials,3)
    Low(1:204,b) = squeeze(mean(BigRespTrials(1:204,PupilAll(:,1,b)<lt,b),2));
    Low(205:408,b) = squeeze(mean(BigRespTrials(409:612,PupilAll(:,3,b)<lt,b),2));
    High(1:204,b) = squeeze(mean(BigRespTrials(1:204,PupilAll(:,1,b)>ht,b),2));
    High(205:408,b) = squeeze(mean(BigRespTrials(409:612,PupilAll(:,3,b)>ht,b),2));
end
DR = [.5 .7]; ctrl = [3:15];
C = Reference(ctrl,:); A = ArousalCoeff;
A(C < DR(1) | C > DR(2)) = NaN; A = nanmean(A,1);  
ctrl = [3:15]; opto = [17:29]; % all steps of the response
C = Reference(ctrl,:); O = Reference(opto,:)-Reference(ctrl,:);
O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); 
clear ctrl opto DR b

% set up stim trace
fr = 15.5; breaki = 10;
Stim = NaN(408+breaki,1);
Lum = cat(1,zeros(ceil(5*fr),1),ones(round(2*fr),1)*-1,ones(round(2*fr),1),ones(round(2*fr),1)*-1,zeros(round(2*fr)+2,1));
Stim(1:204) = Lum; Stim(205+breaki:end) = 0; clear Lum

xstart = [4,5.1,6.1,7.1,8.1,9.1,10.1,11.1];
xend = [4.9,5.9,6.9,7.9,8.9,9.9,10.9,11.9];

xpos = [1,204+breaki];
figure(1); clf
patch(205+breaki+[round(2*fr) round(12*fr) round(12*fr) round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none'); hold on
patch(1+[round(2*fr) round(12*fr) round(12*fr) round(2*fr)],[-1.5 -1.5 1.5 1.5],'red','FaceAlpha',.2,'EdgeColor','none');
for j = 1:length(xpos); 
    patch([xpos(j)+(5*fr) xpos(j)+(11*fr) xpos(j)+(11*fr) xpos(j)+(5*fr)],[-1.1 -1.1 1.1 1.1],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none');
end
plot(Stim,'-k','LineWidth',1);
box off; set(gca,'visible','off'); xlim([1 408+breaki]); ylim([-1.1 1.1])
set(gcf,'units','centimeters','position',[10 10 4 1]);
print(gcf,'-painters','\\anastasia\data\2p\jasmine\DRNoptoFigures_Revisions\Fig4dStim','-depsc')

% bar example
i = 3; Bar = 34; yl = [0 1.2]; %[7,8,12,18,23,25,28,34];
I = newFF == i & (~isnan(sum(sum(cat(3,Low,High),3),1)))';
B = BigRespTrials(:,:,I); P = PupilAll(:,:,I); A1 = A(I); O1 = O(I);
[~,idx] = sort(A(I));

figure(2); 
last = length(xstart);
for i = [1,205+breaki]
    for j = 1:length(xstart)
        patch([i+xstart(j)*fr i+xend(j)*fr i+xend(j)*fr i+xstart(j)*fr],[yl(1) yl(1) yl(2) yl(2)],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none'); hold on
    end
    last = last-1;
end
patch([xstart(1)*fr xend(end)*fr xend(end)*fr xstart(1)*fr],[.5 .5 .7 .7],[1 1 1],'EdgeColor','k','FaceAlpha',0); hold on
patch([xstart(1)*fr xend(end)*fr xend(end)*fr xstart(1)*fr]+204,[.5 .5 .7 .7],[1 1 1],'EdgeColor','k','FaceAlpha',0); hold on
M = nanmean(B(1:1020,:,idx(Bar)),2); M = medfilt1(M,7); M = M([1:204,409:612]);
S = nanstd(B([1:204,409:612],:,idx(Bar)),[],2)./sqrt(size(B,2));
shadedErrorBarJR(1:204,M(1:204),S(1:204),{'color',ColourMap.Ctrl,'LineWidth',1},ColourMap.Ctrl); hold on
shadedErrorBarJR(205+breaki:408+breaki,M([205:408]),S([205:408]),{'color',ColourMap.Ctrl,'LineWidth',1},ColourMap.Ctrl);
M = nanmean(B(1:1020,:,idx(Bar)),2); M = medfilt1(M,7); M = M([205:408,613:816]);
S = nanstd(B([205:408,613:816],:,idx(Bar)),[],2)./sqrt(size(B,2));
shadedErrorBarJR(1:204,M(1:204),S(1:204),{'color',ColourMap.Opto,'LineWidth',1},ColourMap.Opto); hold on
shadedErrorBarJR(205+breaki:408+breaki,M([205:408]),S([205:408]),{'color',ColourMap.Opto,'LineWidth',1},ColourMap.Opto);
xlim([0 408+breaki]); ylim(yl); box off
set(gca,'xcolor','none'); set(gca,'ycolor','none')
set(gcf,'units','centimeters','position',[10 10 4 3]);

% FF example
i=1; FF = 8; yl = [-.05 1.1];%[8,46,58,66,72,211];
lt = .5; ht = .5; t = 8;
I = newFF == i & (~isnan(sum(sum(cat(3,Low,High),3),1)))';
B = BigRespTrials(:,:,I); P = PupilAll(:,:,I); A1 = A(I); O1 = O(I);
[~,idx] = sort(A(I));

figure(4); clf
last = length(xstart);
for i = [1,205+breaki]
    for j = 1:length(xstart)
        patch([i+xstart(j)*fr i+xend(j)*fr i+xend(j)*fr i+xstart(j)*fr],[yl(1) yl(1) yl(2) yl(2)],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none'); hold on
    end
    last = last-1;
end
patch([xstart(1)*fr xend(end)*fr xend(end)*fr xstart(1)*fr],[.5 .5 .7 .7],[1 1 1],'EdgeColor','k','FaceAlpha',0); hold on
patch([xstart(1)*fr xend(end)*fr xend(end)*fr xstart(1)*fr]+204,[.5 .5 .7 .7],[1 1 1],'EdgeColor','k','FaceAlpha',0); hold on
M = nanmean(B(1:1020,:,idx(FF)),2); M = medfilt1(M,7); M = (M - min(M))./max(M - min(M)); M = M([1:204,409:612]);
S = nanstd(B([1:204,409:612],:,idx(FF)),[],2)./sqrt(size(B,2));
shadedErrorBarJR(1:204,M(1:204),S(1:204),{'color',ColourMap.Ctrl,'LineWidth',1},ColourMap.Ctrl); hold on
shadedErrorBarJR(205+breaki:408+breaki,M(205:408),S(205:408),{'color',ColourMap.Ctrl,'LineWidth',1},ColourMap.Ctrl);
M = nanmean(B(1:1020,:,idx(FF)),2); M = medfilt1(M,7); M = (M - min(M))./max(M - min(M)); M = M([205:408,613:816]);
S = nanstd(B([205:408,613:816],:,idx(FF)),[],2)./sqrt(size(B,2));
shadedErrorBarJR(1:204,M(1:204),S(1:204),{'color',ColourMap.Opto,'LineWidth',1},ColourMap.Opto); hold on
shadedErrorBarJR(205+breaki:408+breaki,M(205:408),S(205:408),{'color',ColourMap.Opto,'LineWidth',1},ColourMap.Opto);
xlim([0 408+breaki]); ylim(yl); box off
set(gca,'xcolor','none'); set(gca,'ycolor','none')
set(gcf,'units','centimeters','position',[10 10 4 3]);

%% fig 4e response vs suppression for FF, both, bar and SBC

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
     ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

lt = .5; ht = .5;
for b = 1:size(BigRespTrials,3)
    Low(1:204,b) = squeeze(mean(BigRespTrials(1:204,PupilAll(:,1,b)<lt,b),2));
    Low(205:408,b) = squeeze(mean(BigRespTrials(409:612,PupilAll(:,3,b)<lt,b),2));
    High(1:204,b) = squeeze(mean(BigRespTrials(1:204,PupilAll(:,1,b)>ht,b),2));
    High(205:408,b) = squeeze(mean(BigRespTrials(409:612,PupilAll(:,3,b)>ht,b),2));
end

ctrl = [2:15]; opto = [16:29]; % all steps of the response

Bstart = 0; Bstep = .1; Bend = 1;
bins = [Bstart:Bstep:Bend]; plotbins = [Bstart+(Bstep/2):Bstep:Bend];
figure(1); clf
patch([.5 .7 .7 .5],[-.2 -.2 -.02 -.02],[1 1 1],'EdgeColor','k'); hold on
for i = 1:max(newFF)
    VisIdx = newFF == i;
    C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
    for k = 1:length(bins)-1
        theseO = O(C >= bins(k) & C < bins(k+1));
        Supp(k) = nanmean(theseO);
        SuppStd(k) = nanstd(theseO)./sqrt(sum(~isnan(theseO)));
    end
    shadedErrorBarJR(plotbins,Supp,SuppStd,{'color',cmap(i,:),'LineWidth',1},cmap(i,:)); hold on
    plot(plotbins,Supp,'o','Color',cmap(i,:),'LineWidth',1,'MarkerSize',2); hold on
end
plot([0 1],[0 0],'--k','LineWidth',.5);
xlim([0 1]); ylim([-.35 .05]); axis square; box off
xticks([0 .5 1]); yticks([-.2 -.1 0]); yticklabels({'-.2','-.1','0'})
xlabel('R_{ctrl}'); ylabel('R_{opto} - R_{ctrl}')
set(gcf,'units','centimeters','position',[10 10 4 4]);

% inset
inset1 = axes('Position',[.37 .36 .2 .2]);
DR = [.5 .7]; ctrl = [3:15];
C = Reference(ctrl,:); A = ArousalCoeff;
A(C < DR(1) | C > DR(2)) = NaN; A = nanmean(A,1); 
FF = 8; Bar = 34;
ctrl = [2:15]; opto = [16:29];
for i = [1,3]
    I = newFF == i & (~isnan(sum(sum(cat(3,Low,High),3),1)))';
    R = Reference(:,I); Supp = NaN(size(plotbins));
    A1 = A(I); [~,idx] = sort(A(I));
    if i == 1; j = FF; elseif i == 3; j = Bar; end
    C = R(ctrl,idx(j)); O = R(opto,idx(j))-R(ctrl,idx(j));
    for k = 1:length(bins)-1
        Supp(k) = nanmean(O(C >= bins(k) & C < bins(k+1)));
    end
    plot(plotbins(~isnan(Supp)),Supp(~isnan(Supp)),'-o','Color',cmap(i,:),'LineWidth',1,'MarkerSize',2); hold on
end
plot([0 1],[0 0],'--k','LineWidth',.5);
xlim([0 1]); box off; ylim([-.35 .05]); axis square;
xticks([0 .5 1]); yticks([-.2 0]); yticklabels({'-.2','0'})

%% fig 4f cdf of suppression between .5 and .7 for FF, both, bar and SBC

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
     ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

Bstart = -.6; Bstep = .05; Bend = .2;
bins = [Bstart:Bstep:Bend]; plotbins = [Bstart+(Bstep/2):Bstep:Bend];
DR = [.5 .7]; ctrl = [2:15]; opto = [16:29]; % all steps of the response
clear p
LME = [];
for i = 1:max(newFF)
    VisIdx = newFF == i;
    thisM = mouse(VisIdx); thisF = FOV(VisIdx);
    C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
    O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); 
    thisM(isnan(O)) = []; thisF(isnan(O)) = []; O(isnan(O)) = []; Boutons(i) = length(O);
    figure(4); h = histogram(O,bins,'Normalization','cdf');
    [~,p(i)] = kstest(O);
    thisLME = NaN(Boutons(i),4); 
    thisLME(:,1) = O; thisLME(:,2) = i; thisLME(:,3) = thisM; thisLME(:,4) = thisF;
    LME = cat(1,LME,thisLME);
    figure(2); plot(plotbins,h.Values,'Color',cmap(i,:),'LineWidth',1); hold on
end
ylim([0 1]); box off; 
plot([0 0],[0 1],'--k','LineWidth',.5); xlim([min(plotbins) max(plotbins)])
xticks([-.4 -.2 0]);
xlabel({'R_{opto} - R_{ctrl}'}); ylabel('cumulative probability')
set(gcf,'units','centimeters','position',[15 10 2 4]);

%% LME
Suppression = LME(:,1);
MouseLME = nominal(LME(:,3));

% FF reference
Type = cellstr(num2str(LME(:,2)));
Type(LME(:,2) == 1) = {'1FF'}; Type(LME(:,2) == 2) = {'2FF+Bar'}; Type(LME(:,2) == 3) = {'3Bar'}; Type(LME(:,2) == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFF = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% FF+bAr reference
Type = cellstr(num2str(LME(:,2)));
Type(LME(:,2) == 1) = {'2FF'}; Type(LME(:,2) == 2) = {'1FF+Bar'}; Type(LME(:,2) == 3) = {'3Bar'}; Type(LME(:,2) == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFFBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% Bar reference
Type = cellstr(num2str(LME(:,2)));
Type(LME(:,2) == 1) = {'3FF'}; Type(LME(:,2) == 2) = {'2FF+Bar'}; Type(LME(:,2) == 3) = {'1Bar'}; Type(LME(:,2) == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% SBC reference
Type = cellstr(num2str(thisLME(:,2)));
Type(LME(:,2) == 1) = {'4FF'}; Type(LME(:,2) == 2) = {'2FF+Bar'}; Type(LME(:,2) == 3) = {'3Bar'}; Type(LME(:,2) == 4) = {'1SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeSBC = fitlme(Opto,'Suppression ~ Type + (1|MouseLME)')

%% fig 4G, opto for bar and FF

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
    dataRaw ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

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
        % bar opto
        VisIdx = newFF == 3 & F(:,i); 
        Bnum(t) = sum(VisIdx)/sum(F(:,i));
        C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
        O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); O(isnan(O)) = [];
        Bar(t,1) = mean(O);
        whichM(t) = j;
        if FFnum(t) < .05 | Bnum(t) < .05; Bar(t,:) = NaN; FF(t,:) = NaN; end
        t = t+1;
    end
end

figure(3); clf
for j = 1:length(M)
    B_m(j,:) = nanmean(Bar(whichM == j,:),1);
    F_m(j,:) = nanmean(FF(whichM == j,:),1);
end
bar(1,nanmean(F_m(:,1)),'FaceColor',cmap(1,:),'EdgeColor','none'); hold on
bar(2,nanmean(B_m(:,1)),'FaceColor',cmap(3,:),'EdgeColor','none');
for i = 1:length(M)
    plot([1,2],cat(2,F_m(i,1),B_m(i,1)),'-k');
end
[h,p] = ttest(F_m,B_m);
box off
plot([1 2],[.05 .05],'-k','LineWidth',.5); text(1.4,.06,'*');
xticks([1,2]); xticklabels({'FF','bar'})
ylim([-.25 .07]); yticks([-.2 0]); yticklabels({'-.2','0'})
xlim([.2 2.7])
set(gca,'YColor','k')
ylabel('R_{opto} - R_{ctrl}')
set(gcf,'units','centimeters','position',[15 10 2.5 4]); 

%% fig 4H SF pref for FF+bar
% and supp fig 4G

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
    dataRaw ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

ctrl = [2:15]; opto = [16:29]; DR = [.5 .7]; M = Mouse;
keeps = QI(:,11) > .15;
keeps = QI(:,11) > .15 & newFF == 2;% regress SF pref with FF+Bar boutons

Ms = NaN(size(FOV));
for m = 1:length(M)
    Ms(sum(FOV == M{m},2)>0) = m;
end
Ms = Ms(keeps); SFpreffy = SFpref(keeps);

C = Reference(ctrl,keeps); O = Reference(opto,keeps)-Reference(ctrl,keeps);
O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); 
Ms(isnan(O)) = []; SFpreffy(isnan(O)) = []; O(isnan(O)) = []; 

figure(1);
scatter(SFpreffy,O,2,cmap(2,:),'filled'); hold on
whichM = unique(Ms); b = NaN(length(whichM),2); C = NaN(length(whichM),1);
for m = 1:length(whichM)
    if sum(Ms == whichM(m)) > 10
    b(m,:) = regress(O(Ms == whichM(m))',cat(1,ones(size(O(Ms == whichM(m)))),SFpreffy(Ms == whichM(m))')');
    plot([1 3],[sum(b(m,:)) (b(m,2)*3) + b(m,1)],'Color',[.7 .7 .7],'LineWidth',.5)
    c = corrcoef(O(Ms == whichM(m) & ~isnan(SFpreffy))', SFpreffy(Ms == whichM(m) & ~isnan(SFpreffy)));
    C(m) = c(2,1);
    end
end
bb = regress(O',cat(1,ones(size(O)),SFpreffy')');
plot([1 3],[sum(bb) (bb(2)*3) + bb(1)],'-k','LineWidth',.5)
plot([.8 3.2],[0 0],'--k');
xlim([.8 3.2]); ylim([-.4 .1])
yticks([-.4 -.2 0])
axis square; box off
xlabel('SF pref'); ylabel('R_{opto}-R_{ctrl}')
set(gcf,'units','centimeters','position',[10 10 4 4]);

figure(2); % supp figure 4G
bar(1,nanmean(C),'EdgeColor',cmap(2,:),'FaceColor','none'); hold on
scatter(ones(1,length(C)),C,2,'k','filled')
axis square; box off; ylabel({'corr coeff', '(opto vs SFpref)'}); xticklabels({'FF+bar'})
ylim([-.3 1]); yticks([0 .5 1])
set(gcf,'units','centimeters','position',[15 10 2 4]);

%% ------- SUPPLEMENTARY PANELS
%% 4a big heatmap

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
     ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

ctrl = [3:15]; opto = [17:29]; % all steps of the response
DR = [.5 .7]; M = Mouse;
C = Reference(ctrl,:); O = Reference(opto,:)-Reference(ctrl,:);
O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); 

AllB = []; boutons = 0; Mid = [];
for i = 1:max(newFF)
    I = newFF == i & ~isnan(O)';
    for j = 1:length(M)
        I2 = I & sum(FOV == M{j},2)>0;
        thisB = cat(1,BigRespAll([1:204],I2),...
            BigRespAll([205:408],I2) - BigRespAll([1:204],I2),...
            BigRespAll([409:612],I2),...
            BigRespAll([613:816],I2) - BigRespAll([409:612],I2),...
            BigRespAll([817:end],I2));
        [~,idx] = sort(O(I2));
        AllB = cat(2,AllB,thisB(:,idx));
        Mid = cat(1,Mid,ones(length(idx),1).*j);
    end
    boutons(i+1) = sum(I) + boutons(i);
end

figure(1); clf
imagesc(smoothdata(AllB,1)'); caxis([-.3 .3]); colormap('bluewhitered'); colorbar('northoutside')
hold on
for i = 2:length(boutons)
    plot([1 size(AllB,1)],[boutons(i) boutons(i)],'-k');
end
xticks([])
ylabel('boutons')
set(gcf,'units','centimeters','position',[10 10 15 10.5]); 
set(gca,'units','centimeters','position',[1.5 .1 13.4 9]);

figure(2); clf
Type = ones(1,max(boutons));
for j = 2:4
    Type(boutons(j):boutons(j+1)) = j;
end
imagesc(Type'); colormap(cmap)
xticks([]); yticks([])
set(gcf,'units','centimeters','position',[20 10 .5 10]); 

figure(3); clf
imagesc(Mid); colormap('jet')
xticks([]); yticks([])
set(gcf,'units','centimeters','position',[25 10 .5 10]); 

%% supp figure S4B-C, retinotopy

t = 4;
for i = 1:length(data)
    if ~isempty(data(i).Opto)
        % which boutons have good QIbar
        keeps = sum([max(data(i).Ret.QI(:,1:8),[],2) >= QI_thresh,max(data(i).Ret.QI(:,9:16),[],2) >= QI_thresh],2) == 2;
        % what is their position in FOV (convert to um)
        XY = data(i).ROIxy*(100/267); XY(:,~keeps) = []; % XY(1,:) is A-P axis, XY(2,:) is M-L axis
        % what is their retinotopic preference
        R = cat(2,cat(1,data(i).Ret.Retinotopy.horzPref),cat(1,data(i).Ret.Retinotopy.vertPref));
        R(R(:,1) < 0 | R(:,1) > 8,1) = NaN; R(R(:,2) < 0 | R(:,2) > 8,2) = NaN;
        if data(i).monitor == 2; R(:,1) = R(:,1) + 2; end
        R(~keeps,:) = [];  % R(:,1) is elevation, which moves A-P, R(:,2) is azimuth which moves M-L
        % regress it?
        [r] = corrcoef(XY(1,~isnan(R(:,1))),R(~isnan(R(:,1)),1)'); ele(t) = r(2,1);
        [r] = corrcoef(XY(2,~isnan(R(:,2))),R(~isnan(R(:,2)),2)'); azi(t) = r(2,1);
        % make a pretty plot
        if t == 4% sum(t == [3,4,5,12])>0
            tosser = sum(isnan(R),2) == 0;
            bins = [0:.5:9]; R_ele = discretize(R(:,1),bins); R_azi = discretize(R(:,2),bins); cc = jet(length(bins-1));
            for j = 1:size(cc,1)
                figure(1); scatter(XY(1,tosser & R_ele == j),XY(2,tosser & R_ele == j),2,cc(j,:),'filled'); hold on; 
                set(gca,'YDir','reverse'); daspect([1 1 1]); axis equal; xlim([0 250]); ylim([0 150]); box on; xticks([]); yticks([]); title('elevation')
                set(gcf,'units','centimeters','position',[10 10 5 4]);
                set(gca,'units','centimeters','position',[0 0 5 4]);
                figure(2); scatter(XY(1,tosser & R_azi == j),XY(2,tosser & R_azi == j),2,cc(j,:),'filled'); hold on; 
                set(gca,'YDir','reverse'); daspect([1 1 1]); axis equal; xlim([0 250]); ylim([0 150]); box on; xticks([]); yticks([]); title('azimuth')
                set(gcf,'units','centimeters','position',[10 10 5 4]);
                set(gca,'units','centimeters','position',[0 0 5 4]);
            end
        end
        t = t+1;
    end
end

figure;
bar(1,nanmean(azi),'FaceColor','none'); hold on; scatter(ones(size(azi)),azi,5,'k','filled');
bar(2,nanmean(ele),'FaceColor','none'); hold on; scatter(ones(size(ele)).*2,ele,5,'k','filled');
ylabel({'correlation coefficient','bouton location vs. retinotopic preference'})
ylim([-.2 .8]); yticks([0 .5])
xticks([1 2]); xticklabels({'azimuth','elevation'}); xtickangle(45)
box off; axis square
set(gcf,'units','centimeters','position',[15 10 2.5 4]); 

%% supp fig 4D1 (new) cdf + bar plot dynamic range of suppression 

bins = [-.4:.05:.1]; plotbins = [-.375:.05:.1];
B = BaseOptoDR - BaseDR;
for i = 1:max(newFF)
    figure(1); h = histogram(B(newFF == i),bins,'Normalization','cdf');
    [~,p(i)] = kstest(B(newFF == i));
    figure(2); plot(plotbins,h.Values,'-','Color',cmap(i,:)); hold on
end
plot([0 0],[0 1],'--k');
box off; 
xlabel('Baseline_{opto}-Baseline_{ctrl}'); ylabel('cumulative probability')
xlim([-.4 .1]); xticks([-.2 0]); ylim([0 1]); yticks([0 .5 1])

B = BaseOptoDR - BaseDR;
clear FOVbase
f = unique(FOV);
for i = 1:max(newFF)
    for j = 1:length(f)
        FOVbase(i,j) = nanmean(B(newFF == i & FOV == f(j)));
    end
end
[~,~,stats] = kruskalwallis(FOVbase');
c = multcompare(stats);

figure(2); inset1 = axes('Position',[.25 .5 .2 .35]);
for j = 1:size(FOVbase,1)
    bar(j,nanmean(FOVbase(j,:),2),'FaceColor',cmap(j,:),'EdgeColor','none'); hold on
end
for j = 1:size(FOVbase,2)
    plot([1:size(FOVbase,1)],FOVbase(:,j),'-k','LineWidth',.5);
end
plot([0 size(FOVbase,1)+1],[0 0],'-k');
plot([1 3],[-.28 -.28],'-k'); text(1.6,-.30,'**');
box off; xticks([1:size(FOVbase,1)]); xticklabels(titles(1:size(FOVbase,1))); xtickangle(45)
set(gca, 'XAxisLocation', 'top'); set(gca, 'YAxisLocation', 'right');
xlim([0 size(FOVbase,1)+1]); ylabel('Baseline_{opto}-Baseline_{ctrl}')
ylim([-.3 0]); yticks([-.15 0])

set(gcf,'units','centimeters','position',[15 10 5 4]); 

%% LME
tossy = isnan(newFF);
Suppression = BaseOptoDR(~tossy)' - BaseDR(~tossy)';
MouseLME = nominal(mouse(~tossy));

% FF reference
Type = cellstr(num2str(newFF(~tossy)));
Type(newFF(~tossy) == 1) = {'1FF'}; Type(newFF(~tossy) == 2) = {'2FF+Bar'}; Type(newFF(~tossy) == 3) = {'3Bar'}; Type(newFF(~tossy) == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFF = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% FF+Bar reference
Type = cellstr(num2str(newFF(~tossy)));
Type(newFF(~tossy) == 1) = {'2FF'}; Type(newFF(~tossy) == 2) = {'1FF+Bar'}; Type(newFF(~tossy) == 3) = {'3Bar'}; Type(newFF(~tossy) == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFFBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% Bar reference
Type = cellstr(num2str(newFF(~tossy)));
Type(newFF(~tossy) == 1) = {'3FF'}; Type(newFF(~tossy) == 2) = {'2FF+Bar'}; Type(newFF(~tossy) == 3) = {'1Bar'}; Type(newFF(~tossy) == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% SBC reference
Type = cellstr(num2str(newFF(~tossy)));
Type(newFF(~tossy) == 1) = {'4FF'}; Type(newFF(~tossy) == 2) = {'2FF+Bar'}; Type(newFF(~tossy) == 3) = {'3Bar'}; Type(newFF(~tossy) == 4) = {'1SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeSBC = fitlme(Opto,'Suppression ~ Type + (1|MouseLME)')

%% supp fig 4D2 distribution and bar plot of baseline levels

bins = [0:.1:1]; plotbins = [.05:.1:1];
for i = 1:max(newFF)
    figure(1); h = histogram(BaseDR(newFF == i),bins,'Normalization','cdf');
    figure(2); plot(plotbins,h.Values,'-','Color',cmap(i,:)); hold on
end
axis square; box off; xticks([0 .5 1]); yticks([0 .5 1]);
ylim([0 1]); xlim([0 1]); xlabel('baseline level'); ylabel('cumulative probability')

% inset
clear FOVbase
f = unique(FOV);
for i = 1:max(newFF)
    for j = 1:length(f)
        FOVbase(i,j) = nanmean(BaseDR(newFF == i & FOV == f(j)));
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
plot([1 3],[1 1],'-k'); text(1.8,1.03,'***');
plot([2 4],[.7 .7],'-k'); text(2.8,.73,'***');
plot([3 4],[.8 .8],'-k'); text(2.2,.83,'***');
plot([1 2],[.9 .9],'-k'); text(1.2,.93,'*');
box off; xticks([1:size(FOVbase,1)]); xticklabels(titles(1:size(FOVbase,1))); xtickangle(45)
xlim([0 size(FOVbase,1)+1]); ylabel('baseline level')
ylim([0 1]); yticks([0 .5 1])
set(gca, 'YAxisLocation', 'right');
set(gcf,'units','centimeters','position',[10 10 5 4]); 

%% LME
tossy = isnan(newFF);
Suppression = BaseDR(~tossy)';
MouseLME = nominal(mouse(~tossy));

% FF reference
Type = cellstr(num2str(newFF(~tossy)));
Type(newFF(~tossy) == 1) = {'1FF'}; Type(newFF(~tossy) == 2) = {'2FF+Bar'}; Type(newFF(~tossy) == 3) = {'3Bar'}; Type(newFF(~tossy) == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFF = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% FF+Bar reference
Type = cellstr(num2str(newFF(~tossy)));
Type(newFF(~tossy) == 1) = {'2FF'}; Type(newFF(~tossy) == 2) = {'1FF+Bar'}; Type(newFF(~tossy) == 3) = {'3Bar'}; Type(newFF(~tossy) == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFFBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% Bar reference
Type = cellstr(num2str(newFF(~tossy)));
Type(newFF(~tossy) == 1) = {'3FF'}; Type(newFF(~tossy) == 2) = {'2FF+Bar'}; Type(newFF(~tossy) == 3) = {'1Bar'}; Type(newFF(~tossy) == 4) = {'4SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% SBC reference
Type = cellstr(num2str(newFF(~tossy)));
Type(newFF(~tossy) == 1) = {'4FF'}; Type(newFF(~tossy) == 2) = {'2FF+Bar'}; Type(newFF(~tossy) == 3) = {'3Bar'}; Type(newFF(~tossy) == 4) = {'1SBC'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeSBC = fitlme(Opto,'Suppression ~ Type + (1|MouseLME)')

%% fig supp4E-F mean suppression between .5 and .7 for each FOV and each mouse

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
    dataRaw ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

ctrl = [2:15]; opto = [16:29]; % all steps of the response
DR = [.5 .7]; M = Mouse;

figure(2);
patch([0 1 1 0],[0 0 .05 .05],[.7 .7 .7],'EdgeColor','none','FaceAlpha',.5); hold on
patch([0 .05 .05 0],[.05 .05 1 1],[.7 .7 .7],'EdgeColor','none','FaceAlpha',.5);

c = hsv(length(M)); t = 1; Mbad = zeros(1,length(M));
for j = 1:length(M)
    F = FOV == M{j};
    badcounter = 0;
    for i = 1:size(F,2)
        VisIdx = newFF == 1 & F(:,i); FFnum(t) = sum(VisIdx)/sum(F(:,i));
        C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
        O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); O(isnan(O)) = [];
        FF(t) = mean(O);
        VisIdx = newFF == 3 & F(:,i); Bnum(t) = sum(VisIdx)/sum(F(:,i));
        C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
        O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); O(isnan(O)) = [];
        Bar(t) = mean(O);
        if FFnum(t) < .05 | Bnum(t) < .05; Bar(t) = NaN; FF(t) = NaN; badcounter = badcounter+1; end % at least 5% of boutons
        figure(1); scatter(Bar(t),FF(t),2,c(j,:),'filled'); hold on
        figure(2); scatter(Bnum(t),FFnum(t),2,c(j,:),'filled'); hold on
        t = t+1;
    end
    if badcounter == size(F,2); Mbad(j) = 1; end
end

figure(1); % supp fig 4F
plot([-.3 0],[-.3 0],'--k'); 
ypos = -.02;
for j = 1:length(M)
    if Mbad(j) == 0
        scatter(-.3,ypos,2,c(j,:),'filled'); 
        text(-.275,ypos,strcat('M',num2str(j)));
        ypos = ypos - .02;
    end
end
xlabel('bar R_{opto}-R_{ctrl}'); ylabel('FF R_{opto}-R_{ctrl}');
axis square; box off
set(gca, 'XAxisLocation', 'top'); set(gca, 'YAxisLocation', 'right');
set(gcf,'units','centimeters','position',[10 10 4 4]);

figure(2); % supp fig 4E
ypos = .95;
for j = 1:length(M)
    if Mbad(j) == 0
        scatter(.8,ypos,2,c(j,:),'filled'); 
        text(.85,ypos,strcat('M',num2str(j)));
        ypos = ypos - .05;
    end
end
xlim([0 1]); ylim([0 1])
xlabel('fraction of bar boutons'); ylabel('fraction of FF boutons')
axis square; box off
set(gcf,'units','centimeters','position',[15 10 4 4]);

%% SUPP FIGURE 4H, DSI and ASI

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
     ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

% get DSI and ASI at pref SF
DSIprefSF = NaN(size(SFpref)); ASIprefSF = NaN(size(SFpref));
DSIprefSFangle = NaN(size(SFpref)); ASIprefSFangle = NaN(size(SFpref));
V = SFpref <= 1.5; 
DSIprefSF(V) = DSI(V,1); ASIprefSF(V) = ASI(V,1); 
DSIprefSFangle(V) = DSIangle(V,1); ASIprefSFangle(V) = ASIangle(V,1);
V = SFpref > 1.5 & SFpref < 2.5; 
DSIprefSF(V) = DSI(V,2); ASIprefSF(V) = ASI(V,2);
DSIprefSFangle(V) = DSIangle(V,2); ASIprefSFangle(V) = ASIangle(V,2);
V = SFpref >= 2.5; 
DSIprefSFangle(V) = DSIangle(V,3); ASIprefSFangle(V) = ASIangle(V,3);
DSIprefSF(V) = DSI(V,3); ASIprefSF(V) = ASI(V,3);

% get opto suppression at all DR levels for DS vs nonDS and AS vs nonAS in all categories?? together

ctrl = [2:15]; opto = [16:29]; % all steps of the response
bins = [0:.1:1]; plotbins = [.05:.1:1];
DS_thresh = .4; AS_thresh = .4;
Bary = sum(newFF == [2,3],2)>0;

figure(1); clf
VisIdx = DSIprefSF > DS_thresh & QI(:,11)>.15 & Bary;
C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
for k = 1:length(bins)-1
    theseO = O(C >= bins(k) & C < bins(k+1));
    Supp(k) = nanmean(theseO);
    SuppStd(k) = nanstd(theseO)./sqrt(sum(~isnan(theseO)));
end
shadedErrorBarJR(plotbins,Supp,SuppStd,{'color','k','LineWidth',1},'k'); hold on
plot(plotbins,Supp,'o','Color','k','LineWidth',1,'MarkerSize',2); hold on
    
VisIdx = DSIprefSF < DS_thresh & Bary & QI(:,11)>.15;
C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
clear Supp SuppStd
for k = 1:length(bins)-1
    theseO = O(C >= bins(k) & C < bins(k+1));
    Supp(k) = nanmean(theseO);
    SuppStd(k) = nanstd(theseO)./sqrt(sum(~isnan(theseO)));
end
shadedErrorBarJR(plotbins,Supp,SuppStd,{'--','color','k','LineWidth',1},'k'); hold on
plot(plotbins,Supp,'o','Color','k','LineWidth',1,'MarkerSize',2); hold on
axis square; box off
plot([0 1],[0 0],'--k');
ylim([-.15 .02]); yticks([-.1 0]); xticks([0 .5 1])
ylabel('R_{opto}-R_{ctrl}'); xlabel('R_{ctrl}')
set(gcf,'units','centimeters','position',[10 10 4 4]);

figure(2);
VisIdx = DSIprefSF < DS_thresh & ASIprefSF > AS_thresh & QI(:,11)>.15 & Bary;
C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
clear Supp SuppStd
for k = 1:length(bins)-1
    theseO = O(C >= bins(k) & C < bins(k+1));
    Supp(k) = nanmean(theseO);
    SuppStd(k) = nanstd(theseO)./sqrt(sum(~isnan(theseO)));
end
shadedErrorBarJR(plotbins,Supp,SuppStd,{'color','k','LineWidth',1},'k'); hold on
plot(plotbins,Supp,'o','Color','k','LineWidth',1,'MarkerSize',2); hold on
VisIdx = ASIprefSF < AS_thresh & Bary & QI(:,11)>.15;
C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
clear Supp SuppStd
for k = 1:length(bins)-1
    theseO = O(C >= bins(k) & C < bins(k+1));
    Supp(k) = nanmean(theseO);
    SuppStd(k) = nanstd(theseO)./sqrt(sum(~isnan(theseO)));
end
shadedErrorBarJR(plotbins,Supp,SuppStd,{'--','color','k','LineWidth',1},'k'); hold on
plot(plotbins,Supp,'o','Color','k','LineWidth',1,'MarkerSize',2); hold on
axis square; box off
plot([0 1],[0 0],'--k');
ylim([-.15 .02]); yticks([-.1 0]); xticks([0 .5 1])
ylabel('R_{opto}-R_{ctrl}'); xlabel('R_{ctrl}')
set(gcf,'units','centimeters','position',[15 10 4 4]);

