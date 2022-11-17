%% Figure 1

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
    'DefaultHistogramLineWidth',1.5,...
    'DefaultTextFontUnits','Points',...
    'DefaultTextFontSize',7,...
    'DefaultTextFontName','Arial');

%% load GRAB data

clear
% load Photometry.mat
PlotColours

%% fig 1b trace of GRAB response

i = 5;
figure(2); clf
set(gcf,'units','centimeters','position',[10 10 4 4]);
set(gca,'units','centimeters','position',[10 10 4 4]);
h = axes;
B = nanmean(GRAB(i).vis_resp(:,GRAB(i).Stim.Trials == GRAB(i).Stim.Blank),2);
C = nanstd(GRAB(i).vis_resp(:,GRAB(i).Stim.Trials == GRAB(i).Stim.Blank),[],2)./sqrt(sum(GRAB(i).Stim.Trials == GRAB(i).Stim.Blank));
shadedErrorBarJR([],B,C,{'color',ColourMap.Ctrl,'LineWidth',1.5},ColourMap.Ctrl); hold on
B = nanmean(GRAB(i).vis_resp(:,GRAB(i).Stim.Trials == GRAB(i).Stim.BlankOpto),2);
C = nanstd(GRAB(i).vis_resp(:,GRAB(i).Stim.Trials == GRAB(i).Stim.BlankOpto),[],2)./sqrt(sum(GRAB(i).Stim.Trials == GRAB(i).Stim.BlankOpto));
shadedErrorBarJR([],B,C,{'color',ColourMap.Opto,'LineWidth',1.5},ColourMap.Opto); hold on
patch([200 800 800 200],[.015 .015 .016 .016],'red','FaceAlpha',.2,'EdgeColor','none');
plot([0 200],[.005 .005],'-k','LineWidth',1); 
plot([2 2],[.005 .015],'-k','LineWidth',1);
text(80,.0065,'2 s','Fontsize',7,'FontName','Arial')
text(-80,.007,'1% \DeltaF/F','rotation',90,'Fontsize',7,'FontName','Arial')
plot([200 400],[.018 .018],'-','Color',ColourMap.Ctrl,'LineWidth',1.5);
text(450,.0185,'ctrl trials','Fontsize',7,'FontName','Arial')
plot([200 400],[.0205 .0205],'-','Color',ColourMap.Opto,'LineWidth',1.5);
text(450,.021,'opto trials','Fontsize',7,'FontName','Arial')
xlim([0 1501]); ylim([-.01 .022]); box off
set(h,'xcolor','none'); set(h,'ycolor','none')

%% fig 1C and 1D opto effect on GRAB, pupil and running

clearvars -except GRAB GRAB_data RGC_new RGC_glu ColourMap TTX TTX_glu
t = 1; u = 1; v= 1; nfr = 100;
pupCtrlAll = []; pupOptoAll = []; runCtrlAll = []; runOptoAll = [];
for i = 1:length(GRAB)
    if GRAB(i).frequency== 20 & strcmp(GRAB(i).genotype,'pet1')
        if isfield(GRAB(i).Stim,'Blank');
            V = GRAB(i).Stim.Trials == GRAB(i).Stim.Blank; VO = GRAB(i).Stim.Trials == GRAB(i).Stim.BlankOpto;
        elseif isfield(GRAB(i).Stim,'Vis');
            V = GRAB(i).Stim.Trials == GRAB(i).Stim.Vis; VO = GRAB(i).Stim.Trials == GRAB(i).Stim.VisOpto;
        end
        dffCtrl(u) = nanmean(mean(GRAB(i).vis_resp(4*nfr:8*nfr,V),1),2);
        dffOpto(u) = nanmean(mean(GRAB(i).vis_resp(4*nfr:8*nfr,VO),1),2);
        mouseG{u} = GRAB(i).mouse;
        if ~isempty(GRAB(i).pup_resp)
        pupCtrl(t) = nanmean(mean(GRAB(i).pup_resp(4*nfr:8*nfr,V),1)-mean(GRAB(i).pup_resp(1:2*nfr,V),1),2);
        pupCtrl(t) = pupCtrl(t)/max(GRAB(i).peak_pup);
        pupOpto(t) = nanmean(mean(GRAB(i).pup_resp(4*nfr:8*nfr,VO),1)-mean(GRAB(i).pup_resp(1:2*nfr,VO),1),2);
        pupOpto(t) = pupOpto(t)/max(GRAB(i).peak_pup);
        pupCtrlAll = cat(1,pupCtrlAll,cat(2,ones(sum(V),1).*t,(mean(GRAB(i).pup_resp(4*nfr:8*nfr,V),1)-mean(GRAB(i).pup_resp(1:2*nfr,V),1))'./max(GRAB(i).peak_pup)));
        pupOptoAll = cat(1,pupOptoAll,cat(2,ones(sum(VO),1).*t,(mean(GRAB(i).pup_resp(4*nfr:8*nfr,VO),1)-mean(GRAB(i).pup_resp(1:2*nfr,VO),1))'./max(GRAB(i).peak_pup)));
        mouseP{t} = GRAB(i).mouse;
        t = t+1;
        end
        if ~isempty(GRAB(i).run_resp)
        runCtrl(v) = nanmean(mean(GRAB(i).run_resp(4*nfr:8*nfr,V),1)-mean(GRAB(i).run_resp(1:2*nfr,V),1),2);
        runOpto(v) = nanmean(mean(GRAB(i).run_resp(4*nfr:8*nfr,VO),1)-mean(GRAB(i).run_resp(1:2*nfr,VO),1),2);
        runCtrlAll = cat(1,runCtrlAll,cat(2,ones(sum(V),1).*v,mean(GRAB(i).run_resp(4*nfr:8*nfr,V),1)'-mean(GRAB(i).run_resp(1:2*nfr,V),1)'));
        runOptoAll = cat(1,runOptoAll,cat(2,ones(sum(VO),1).*v,mean(GRAB(i).run_resp(4*nfr:8*nfr,VO),1)'-mean(GRAB(i).run_resp(1:2*nfr,VO),1)'));
        mouseR{v} = GRAB(i).mouse;
        v = v+1;
        end
        u = u+1;
    end
end

for i = 1:length(RGC_new)
    if RGC_new(i).frequency == 20 & strcmp(RGC_new(i).genotype,'pet1')
        if RGC_new(i).prepulse == 1; b = 0; else b = 3; end
        if ~isempty(RGC_new(i).pup_resp)
        pupCtrl(t) = nanmean(mean(RGC_new(i).pup_resp(4*nfr:(5+b)*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1)-mean(RGC_new(i).pup_resp(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1),2);
        pupCtrl(t) = pupCtrl(t)/max(RGC_new(i).peak_pup);
        pupOpto(t) = nanmean(mean(RGC_new(i).pup_resp(4*nfr:(5+b)*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1)-mean(RGC_new(i).pup_resp(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1),2);
        pupOpto(t) = pupOpto(t)/max(RGC_new(i).peak_pup);
        pupCtrlAll = cat(1,pupCtrlAll,cat(2,ones(sum(RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1).*t,(mean(RGC_new(i).pup_resp(4*nfr:(5+b)*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1)-mean(RGC_new(i).pup_resp(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1))'./max(RGC_new(i).peak_pup)));
        pupOptoAll = cat(1,pupOptoAll,cat(2,ones(sum(RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1).*t,(mean(RGC_new(i).pup_resp(4*nfr:(5+b)*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1)-mean(RGC_new(i).pup_resp(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1))'./max(RGC_new(i).peak_pup)));
        mouseP{t} = RGC_new(i).mouse;
        t = t+1;
        end
        if ~isempty(RGC_new(i).run_resp)
        runCtrl(v) = nanmean(mean(RGC_new(i).run_resp(4*nfr:(5+b)*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1)-mean(RGC_new(i).run_resp(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1),2);
        runOpto(v) = nanmean(mean(RGC_new(i).run_resp(4*nfr:(5+b)*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1)-mean(RGC_new(i).run_resp(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1),2);
        runCtrlAll = cat(1,runCtrlAll,cat(2,ones(sum(RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1).*v,mean(RGC_new(i).run_resp(4*nfr:(5+b)*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1)'-mean(RGC_new(i).run_resp(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1)'));
        runOptoAll = cat(1,runOptoAll,cat(2,ones(sum(RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1).*v,mean(RGC_new(i).run_resp(4*nfr:(5+b)*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1)'-mean(RGC_new(i).run_resp(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1)'));
        mouseR{v} = RGC_new(i).mouse;
        v = v+1;
        end
    end
end

for i = 1:length(RGC_glu)
    if RGC_glu(i).frequency == 20 & strcmp(RGC_glu(i).genotype,'pet1')
        b = 0;
        if ~isempty(RGC_glu(i).pup_resp)
        pupCtrl(t) = nanmean(mean(RGC_glu(i).pup_resp(4*nfr:(5+b)*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis),1)-mean(RGC_glu(i).pup_resp(1:2*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis),1),2);
        pupCtrl(t) = pupCtrl(t)/max(RGC_glu(i).peak_pup);
        pupOpto(t) = nanmean(mean(RGC_glu(i).pup_resp(4*nfr:(5+b)*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto),1)-mean(RGC_glu(i).pup_resp(1:2*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto),1),2);
        pupOpto(t) = pupOpto(t)/max(RGC_glu(i).peak_pup);
        pupCtrlAll = cat(1,pupCtrlAll,cat(2,ones(sum(RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis),1).*t,(mean(RGC_glu(i).pup_resp(4*nfr:(5+b)*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis),1)-mean(RGC_glu(i).pup_resp(1:2*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis),1))'./max(RGC_glu(i).peak_pup)));
        pupOptoAll = cat(1,pupOptoAll,cat(2,ones(sum(RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto),1).*t,(mean(RGC_glu(i).pup_resp(4*nfr:(5+b)*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto),1)-mean(RGC_glu(i).pup_resp(1:2*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto),1))'./max(RGC_glu(i).peak_pup)));
        mouseP{t} = RGC_glu(i).mouse;
        t = t+1;
        end
        if ~isempty(RGC_glu(i).run_resp)
        runCtrl(v) = nanmean(mean(RGC_glu(i).run_resp(4*nfr:(5+b)*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis),1)-mean(RGC_glu(i).run_resp(1:2*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis),1),2);
        runOpto(v) = nanmean(mean(RGC_glu(i).run_resp(4*nfr:(5+b)*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto),1)-mean(RGC_glu(i).run_resp(1:2*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto),1),2);
        runCtrlAll = cat(1,runCtrlAll,cat(2,ones(sum(RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis),1).*v,mean(RGC_glu(i).run_resp(4*nfr:(5+b)*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis),1)'-mean(RGC_glu(i).run_resp(1:2*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis),1)'));
        runOptoAll = cat(1,runOptoAll,cat(2,ones(sum(RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto),1).*v,mean(RGC_glu(i).run_resp(4*nfr:(5+b)*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto),1)'-mean(RGC_glu(i).run_resp(1:2*nfr,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto),1)'));
        mouseR{v} = RGC_glu(i).mouse;
        v = v+1;
        end
    end
end

% with each run for main
figure(1); clf % fig 1C
set(gcf,'units','centimeters','position',[10 10 2.3 4]);
bar(1,nanmean(dffCtrl),'FaceColor',ColourMap.Ctrl,'FaceAlpha',.8,'EdgeColor','none'); hold on
bar(2,nanmean(dffOpto),'FaceColor',ColourMap.Opto,'FaceAlpha',.8,'EdgeColor','none'); hold on
plot(cat(1,dffCtrl,dffOpto),'-','Color',[.7 .7 .7])
[~,p1] = ttest(dffCtrl,dffOpto);
M = unique(mouseG);
for i = 1:length(M)
    thisM = find(strcmp(mouseG,M{i}));
    plot(mean(cat(1,dffCtrl(thisM),dffOpto(thisM)),2),'-k','LineWidth',1)
end
box off; 
plot([1 2],[.017 .017],'-k'); text(1.3,.0175,'***')
yticks([0 .01]); yticklabels({'0','1'}); ylim([-.006 .018]); xtickangle(45)
xlim([0 3]); xticks([1 2]); xticklabels({'ctrl','opto'})
ylabel('% \DeltaF/F'); 

Diff = dffCtrl' - dffOpto';
Mouse = nominal(mouseG)';
Opto = table(Diff,Mouse);
lme = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')

figure(2); clf % fig 1D pupil
set(gcf,'units','centimeters','position',[15 10 2.3 4]);
bar(1,nanmean(pupCtrl),'FaceColor',ColourMap.Ctrl,'FaceAlpha',.8,'EdgeColor','none'); hold on
bar(2,nanmean(pupOpto),'FaceColor',ColourMap.Opto,'FaceAlpha',.8,'EdgeColor','none'); hold on
plot(cat(1,pupCtrl,pupOpto),'-','Color',[.7 .7 .7])
[~,p2] = ttest(pupCtrl,pupOpto);
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot(mean(cat(1,pupCtrl(thisM),pupOpto(thisM)),2),'-k','LineWidth',1)
end
box off; xlim([0 3]); xticks([1 2]); xticklabels({'ctrl','opto'});xtickangle(45)
ylim([-.1 .1]); yticks([-.1 0 .1])
ylabel('\Delta pupil (fraction of max)'); title('pupil')

Diff = pupCtrl' - pupOpto';
Mouse = nominal(mouseP)';
Opto = table(Diff,Mouse);
lmeP = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')

figure(3); clf % fig 1D running
set(gcf,'units','centimeters','position',[20 10 2.3 4]);
bar(1,nanmean(runCtrl),'FaceColor',ColourMap.Ctrl,'FaceAlpha',.8,'EdgeColor','none'); hold on
bar(2,nanmean(runOpto),'FaceColor',ColourMap.Opto,'FaceAlpha',.8,'EdgeColor','none'); hold on
plot(cat(1,runCtrl,runOpto),'-','Color',[.7 .7 .7])
[~,p3] = ttest(runCtrl,runOpto);
M = unique(mouseR);
for i = 1:length(M)
    thisM = find(strcmp(mouseR,M{i}));
    plot(mean(cat(1,runCtrl(thisM),runOpto(thisM)),2),'-k','LineWidth',1)
end
box off; xlim([0 3]); xticks([1 2]); xticklabels({'ctrl','opto'}); xtickangle(45)
ylim([-12 12])
ylabel('\Delta speed (cm/s)'); title('running')

Diff = runCtrl' - runOpto';
Mouse = nominal(mouseR)';
Opto = table(Diff,Mouse);
lmeR = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')

%% fig 1F-G trace of GCaMP with Chrimson and without Chrimson

i = 3;
figure(2); clf
set(gcf,'units','centimeters','position',[10 10 3 4]);
set(gca,'units','centimeters','position',[10 10 3 4]);
h = axes;
B = nanmean(RGC_new(i).vis_resp(:,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(1)),2);
C = nanstd(RGC_new(i).vis_resp(:,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(1)),[],2)./sqrt(sum(RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(1)));
B(1:100) = []; C(1:100) = []; k = round(max(B+C)*1000)./1000; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Ctrl,'LineWidth',1.5},ColourMap.Ctrl); hold on
B = nanmean(RGC_new(i).vis_resp(:,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(1)),2);
C = nanstd(RGC_new(i).vis_resp(:,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(1)),[],2)./sqrt(sum(RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(1)));
B(1:100) = []; C(1:100) = []; r = round(min(B-C)*1000)./1000; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Opto,'LineWidth',1.5},ColourMap.Opto); hold on
plot([0 200],[.005 .005],'-k','LineWidth',1); 
plot([2 2],[.005 .015],'-k','LineWidth',1);
text(80,.007,'2 s','Fontsize',7,'FontName','Arial')
text(-60,.007,'1% \DeltaF/F','rotation',90,'Fontsize',7,'FontName','Arial')

patch([100 700 700 100],[k+.001 k+.001 k+.002 k+.002],'red','FaceAlpha',.2,'EdgeColor','none');
patch([400 600 600 400],[k+.002 k+.002 k+.003 k+.003],'k','FaceAlpha',.2,'EdgeColor','none');
plot([200 400],[k+.005 k+.005],'-','Color',ColourMap.Ctrl,'LineWidth',1.5);
text(450,k+.0055,'ctrl trials','Fontsize',7,'FontName','Arial')
plot([200 400],[k+.008 k+.008],'-','Color',ColourMap.Opto,'LineWidth',1.5);
text(450,k+.0085,'opto trials','Fontsize',7,'FontName','Arial')
xlim([0 length(B)]); ylim([r k+.009]); box off
set(h,'xcolor','none'); set(h,'ycolor','none')

% ctrl trace
i = 19;
figure(3); clf
set(gcf,'units','centimeters','position',[20 10 3 4]);
set(gca,'units','centimeters','position',[20 10 3 4]);
h = axes;
B = nanmean(RGC_new(i).vis_resp(:,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(1)),2);
C = nanstd(RGC_new(i).vis_resp(:,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(1)),[],2)./sqrt(sum(RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(1)));
B(1:100) = []; C(1:100) = []; k = round(max(B+C)*1000)./1000; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Ctrl,'LineWidth',1.5},ColourMap.Ctrl); hold on
B = nanmean(RGC_new(i).vis_resp(:,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(1)),2);
C = nanstd(RGC_new(i).vis_resp(:,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(1)),[],2)./sqrt(sum(RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(1)));
B(1:100) = []; C(1:100) = []; r = round(min(B-C)*1000)./1000; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Opto,'LineWidth',1.5},ColourMap.Opto); hold on
plot([0 200],[.005 .005],'-k','LineWidth',1); 
plot([2 2],[.005 .015],'-k','LineWidth',1);
text(80,.007,'2 s','Fontsize',7,'FontName','Arial')
text(-60,.007,'1% \DeltaF/F','rotation',90,'Fontsize',7,'FontName','Arial')

patch([100 700 700 100],[k+.001 k+.001 k+.002 k+.002],'red','FaceAlpha',.2,'EdgeColor','none');
patch([400 600 600 400],[k+.002 k+.002 k+.003 k+.003],'k','FaceAlpha',.2,'EdgeColor','none');
xlim([0 length(B)]); ylim([r k+.009]); box off
set(h,'xcolor','none'); set(h,'ycolor','none')

%% fig 1h and 1I bar plot of GCaMP with Chrimson
% and supplementary figure s1H

clearvars -except GRAB GRAB_data RGC_new RGC_glu ColourMap TTX TTX_glu
t = 1; u = 1; nfr = 100;
for i = 1:length(RGC_new)
    if RGC_new(i).frequency == 20 
        if strcmp(RGC_new(i).genotype,'pet1')
            if RGC_new(i).prepulse == 1;
                ChR(t,1) = mean(mean(RGC_new(i).vis_resp(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),2),1);
                ChR(t,2) = mean(mean(RGC_new(i).vis_resp(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),2),1);
                S = contains(RGC_new(i).Stim.VisStim,'Sine_0_lSF_primetime');
                ChRv(t,1) = mean(mean(RGC_new(i).vis_resp(3*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S)),2),1);
                ChRv(t,2) = mean(mean(RGC_new(i).vis_resp(3*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S)),2),1);
                % correct for baseline drop
                ChRvS(t,1) = mean(mean(RGC_new(i).vis_resp(3*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S)),1)...
                    - mean(RGC_new(i).vis_resp(2.75*nfr:3*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S)),1),2);
                ChRvS(t,2) = mean(mean(RGC_new(i).vis_resp(3*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S)),1)...
                    - mean(RGC_new(i).vis_resp(2.75*nfr:3*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S)),1),2);
                mouseP{t} = RGC_new(i).mouse;
                t = t+1;
            else
                ChR(t,1) = mean(mean(RGC_new(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.Blank,RGC_new(i).Stim.Vis],2))),2),1);
                ChR(t,2) = mean(mean(RGC_new(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.BlankOpto,RGC_new(i).Stim.VisOpto],2))),2),1);
                S = contains(RGC_new(i).Stim.VisStim,'Sine_0_lSF_primetime');
                ChRv(t,1) = mean(mean(RGC_new(i).vis_resp(5*nfr:7*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S)),2),1);
                ChRv(t,2) = mean(mean(RGC_new(i).vis_resp(5*nfr:7*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S)),2),1);
                % correct baseline drop
                ChRvS(t,1) = mean(mean(RGC_new(i).vis_resp(5*nfr:7*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S)),1)...
                    - mean(RGC_new(i).vis_resp(4*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S)),1),2);
                ChRvS(t,2) = mean(mean(RGC_new(i).vis_resp(5*nfr:7*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S)),1)...
                    - mean(RGC_new(i).vis_resp(4*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S)),1),2);
                mouseP{t} = RGC_new(i).mouse;
                t = t+1;
            end
        elseif strcmp(RGC_new(i).genotype,'WT')
            if RGC_new(i).prepulse == 1;
                noChR(u,1) = mean(mean(RGC_new(i).vis_resp(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),2),1);
                noChR(u,2) = mean(mean(RGC_new(i).vis_resp(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),2),1);
                S = contains(RGC_new(i).Stim.VisStim,'Sine_0_lSF_primetime') | contains(RGC_new(i).Stim.VisStim,'NoiseBarHorz1of1');
                if sum(S) > 0 
                noChRv(u,1) = mean(mean(RGC_new(i).vis_resp(3*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S)),2),1);
                noChRv(u,2) = mean(mean(RGC_new(i).vis_resp(3*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S)),2),1);
                % correct baseline
                noChRvS(u,1) = mean(mean(RGC_new(i).vis_resp(3*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S)),1)...
                    - mean(RGC_new(i).vis_resp(2.75*nfr:3*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S)),1),2);
                noChRvS(u,2) = mean(mean(RGC_new(i).vis_resp(3*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S)),1)...
                    - mean(RGC_new(i).vis_resp(2.75*nfr:3*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S)),1),2);
                else; noChRv(u,:) = NaN; noChRvS(u,:) = NaN; end
                mouseC{u} = RGC_new(i).mouse;
                u = u+1;
            else
                noChR(u,1) = mean(mean(RGC_new(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.Blank,RGC_new(i).Stim.Vis],2))),2),1);
                noChR(u,2) = mean(mean(RGC_new(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.BlankOpto,RGC_new(i).Stim.VisOpto],2))),2),1);
                S = contains(RGC_new(i).Stim.VisStim,'Sine_0_lSF_primetime') | contains(RGC_new(i).Stim.VisStim,'NoiseBarHorz1of1');
                if sum(S) > 0 
                noChRv(u,1) = mean(mean(RGC_new(i).vis_resp(5*nfr:7*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S)),2),1);
                noChRv(u,2) = mean(mean(RGC_new(i).vis_resp(5*nfr:7*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S)),2),1);
                % baseline correct
                noChRvS(u,1) = mean(mean(RGC_new(i).vis_resp(5*nfr:7*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S)),1)...
                    - mean(RGC_new(i).vis_resp(4*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S)),1),2);
                noChRvS(u,2) = mean(mean(RGC_new(i).vis_resp(5*nfr:7*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S)),1)...
                    - mean(RGC_new(i).vis_resp(4*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S)),1),2);
                else; noChRv(u,:) = NaN; end
                mouseC{u} = RGC_new(i).mouse;
                u = u+1;
            end
        end
    end
end

% BASELINE with mouse averages for main fig 1H
figure(1); clf
% Pet1
bar(1,nanmean(ChR(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(2,nanmean(ChR(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot([1,2],ChR','-','Color',[.7 .7 .7])
[~,p1] = ttest2(ChR(:,1),ChR(:,2));
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot([1,2],nanmean(ChR(thisM,:),1)','-k','LineWidth',1)
end
Diff = ChR(:,1) - ChR(:,2);
Mouse = nominal(mouseP)';
Opto = table(Diff,Mouse);
lmePet1Baseline = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
% WT
bar(3,nanmean(noChR(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(4,nanmean(noChR(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot([3,4],noChR','-','Color',[.7 .7 .7])
[~,p2] = ttest2(noChR(:,1),noChR(:,2));
M = unique(mouseC);
for i = 1:length(M)
    thisM = find(strcmp(mouseC,M{i}));
    plot([3,4],nanmean(noChR(thisM,:),1)','-k','LineWidth',1)
end
Diff = noChR(:,1) - noChR(:,2);
Mouse = nominal(mouseC)';
Opto = table(Diff,Mouse);
lmeWTBaseline = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')

plot([1 2],[.005 .005],'-k'); text(1,.006,'Pet1')
plot([1 2],[.004 .004],'-k'); text(1.2,.0043,'***')
plot([3 4],[.005 .005],'-k'); text(3,.006,'WT')
box off; xlim([0 5]); xticks([1 2 3 4]); xticklabels({'ctrl','opto','ctrl','opto'}); xtickangle(45)
ylim([-.015 .006]); yticks([-.01 0]); yticklabels({'-1','0'}) % in percentage
ylabel('% \Delta F/F'); title('baseline')
set(gcf,'units','centimeters','position',[10 10 3.5 4]);

% visual response, fig 1I
figure(2); clf
% Pet1
set(gcf,'units','centimeters','position',[15 10 3.5 4]);
bar(1,nanmean(ChRv(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(2,nanmean(ChRv(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot([1,2],ChRv','-','Color',[.7 .7 .7])
[~,p3] = ttest2(ChRv(:,1),ChRv(:,2));
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot([1,2],nanmean(ChRv(thisM,:),1)','-k','LineWidth',1)
end
Diff = ChRv(:,1) - ChRv(:,2);
Mouse = nominal(mouseP)';
Opto = table(Diff,Mouse);
lmePet1Vis = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
% WT
bar(3,nanmean(noChRv(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(4,nanmean(noChRv(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot([3,4],noChRv','-','Color',[.7 .7 .7])
[~,p4] = ttest2(noChRv(:,1),noChRv(:,2));
M = unique(mouseC);
for i = 1:length(M)
    thisM = find(strcmp(mouseC,M{i}));
    plot([3,4],nanmean(noChRv(thisM,:),1)','-k','LineWidth',1)
end
Diff = noChRv(:,1) - noChRv(:,2);
Mouse = nominal(mouseC)';
Opto = table(Diff,Mouse);
lmeWTVis = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
plot([1 2],[.056 .056],'-k'); text(1,.059,'Pet1')
plot([1 2],[.053 .053],'-k'); text(1.2,.054,'**')
plot([3 4],[.056 .056],'-k'); text(3,.059,'WT')
box off; xlim([0 5]); xticks([1 2 3 4]); xticklabels({'ctrl','opto','ctrl','opto'}); xtickangle(45)
ylim([0 .06]); yticks([0 .02 .04]); yticklabels({'0','2','4'})
ylabel('% \Delta F/F'); title('vis resp')

% visual response, baseline adjusted, figure S1H
figure(3); clf
% Pet1
set(gcf,'units','centimeters','position',[20 10 3.5 4]);
bar(1,nanmean(ChRvS(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(2,nanmean(ChRvS(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot([1,2],ChRvS','-','Color',[.7 .7 .7])
[~,p5] = ttest2(ChRvS(:,1),ChRvS(:,2));
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot([1,2],nanmean(ChRvS(thisM,:),1)','-k','LineWidth',1)
end
Diff = ChRvS(:,1) - ChRvS(:,2);
Mouse = nominal(mouseP)';
Opto = table(Diff,Mouse);
lmePet1VisBase = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
% WT
bar(3,nanmean(noChRvS(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(4,nanmean(noChRvS(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot([3,4],noChRvS','-','Color',[.7 .7 .7])
[~,p6] = ttest2(noChRvS(:,1),noChRvS(:,2));
M = unique(mouseC);
for i = 1:length(M)
    thisM = find(strcmp(mouseC,M{i}));
    plot([3,4],nanmean(noChRvS(thisM,:),1)','-k','LineWidth',1)
end
Diff = noChRvS(:,1) - noChRvS(:,2);
Mouse = nominal(mouseC)';
Opto = table(Diff,Mouse);
lmeWTVisBase = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
plot([1 2],[.056 .056],'-k'); text(1,.059,'Pet1')
plot([3 4],[.056 .056],'-k'); text(3,.059,'WT')
box off; xlim([0 5]); xticks([1 2 3 4]); xticklabels({'ctrl','opto','ctrl','opto'}); xtickangle(45)
ylim([0 .06]); yticks([0 .02 .04]); yticklabels({'0','2','4'})
ylabel('% \Delta F/F'); title('vis resp')

%% fig 1J trace of GCaMP and opto with TTX and washout

clearvars -except GRAB GRAB_data RGC_new RGC_glu ColourMap TTX TTX_glu
i=2;

% washout
figure(3); clf
set(gcf,'units','centimeters','position',[20 10 3 4]);
set(gca,'units','centimeters','position',[20 10 3 4]);
h = axes;
B = nanmean(TTX.post(i).vis_resp(:,TTX.post(i).Stim.Trials == TTX.post(i).Stim.Vis),2);
C = nanstd(TTX.post(i).vis_resp(:,TTX.post(i).Stim.Trials == TTX.post(i).Stim.Vis),[],2)./sqrt(sum(TTX.post(i).Stim.Trials == TTX.post(i).Stim.Vis));
B(1:100) = []; C(1:100) = []; k = round(max(B+C)*1000)./1000; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Ctrl,'LineWidth',1.5},ColourMap.Ctrl); hold on
B = nanmean(TTX.post(i).vis_resp(:,TTX.post(i).Stim.Trials == TTX.post(i).Stim.VisOpto),2);
C = nanstd(TTX.post(i).vis_resp(:,TTX.post(i).Stim.Trials == TTX.post(i).Stim.VisOpto),[],2)./sqrt(sum(TTX.post(i).Stim.Trials == TTX.post(i).Stim.VisOpto));
B(1:100) = []; C(1:100) = []; r = round(min(B-C)*1000)./1000; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Opto,'LineWidth',1.5},ColourMap.Opto); hold on

patch([100 700 700 100],[k+.001 k+.001 k+.002 k+.002],'red','FaceAlpha',.2,'EdgeColor','none');
patch([400 600 600 400],[k+.002 k+.002 k+.003 k+.003],'k','FaceAlpha',.2,'EdgeColor','none');
xlim([0 length(B)]); ylim([r k+.009]); box off
set(h,'xcolor','none'); set(h,'ycolor','none')
title('washout')

% pre TTX
figure(1); clf
set(gcf,'units','centimeters','position',[10 10 3 4]);
set(gca,'units','centimeters','position',[10 10 3 4]);
h = axes;
B = nanmean(TTX.pre(i).vis_resp(:,TTX.pre(i).Stim.Trials == TTX.pre(i).Stim.Vis),2);
C = nanstd(TTX.pre(i).vis_resp(:,TTX.pre(i).Stim.Trials == TTX.pre(i).Stim.Vis),[],2)./sqrt(sum(TTX.pre(i).Stim.Trials == TTX.pre(i).Stim.Vis));
B(1:100) = []; C(1:100) = []; %k = round(max(B+C)*1000)./1000; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Ctrl,'LineWidth',1.5},ColourMap.Ctrl); hold on
B = nanmean(TTX.pre(i).vis_resp(:,TTX.pre(i).Stim.Trials == TTX.pre(i).Stim.VisOpto),2);
C = nanstd(TTX.pre(i).vis_resp(:,TTX.pre(i).Stim.Trials == TTX.pre(i).Stim.VisOpto),[],2)./sqrt(sum(TTX.pre(i).Stim.Trials == TTX.pre(i).Stim.VisOpto));
B(1:100) = []; C(1:100) = []; %r = round(min(B-C)*1000)./1000; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Opto,'LineWidth',1.5},ColourMap.Opto); hold on
plot([0 200],[.005 .005],'-k','LineWidth',1); 
plot([2 2],[.005 .015],'-k','LineWidth',1);
text(80,.007,'2 s','Fontsize',7,'FontName','Arial')
text(-60,.007,'1% \DeltaF/F','rotation',90,'Fontsize',7,'FontName','Arial')

patch([100 700 700 100],[k+.001 k+.001 k+.002 k+.002],'red','FaceAlpha',.2,'EdgeColor','none');
patch([400 600 600 400],[k+.002 k+.002 k+.003 k+.003],'k','FaceAlpha',.2,'EdgeColor','none');
xlim([0 length(B)]); ylim([r k+.009]); box off
set(h,'xcolor','none'); set(h,'ycolor','none')
title('pre-TTX')

% post TTX
figure(2); clf
set(gcf,'units','centimeters','position',[15 10 3 4]);
set(gca,'units','centimeters','position',[15 10 3 4]);
h = axes;
B = nanmean(TTX.drug(i).vis_resp(:,TTX.drug(i).Stim.Trials == TTX.drug(i).Stim.Vis),2);
C = nanstd(TTX.drug(i).vis_resp(:,TTX.drug(i).Stim.Trials == TTX.drug(i).Stim.Vis),[],2)./sqrt(sum(TTX.drug(i).Stim.Trials == TTX.drug(i).Stim.Vis));
B(1:100) = []; C(1:100) = []; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Ctrl,'LineWidth',1.5},ColourMap.Ctrl); hold on
B = nanmean(TTX.drug(i).vis_resp(:,TTX.drug(i).Stim.Trials == TTX.drug(i).Stim.VisOpto),2);
C = nanstd(TTX.drug(i).vis_resp(:,TTX.drug(i).Stim.Trials == TTX.drug(i).Stim.VisOpto),[],2)./sqrt(sum(TTX.drug(i).Stim.Trials == TTX.drug(i).Stim.VisOpto));
B(1:100) = []; C(1:100) = []; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Opto,'LineWidth',1.5},ColourMap.Opto); hold on

patch([100 700 700 100],[k+.001 k+.001 k+.002 k+.002],'red','FaceAlpha',.2,'EdgeColor','none');
patch([400 600 600 400],[k+.002 k+.002 k+.003 k+.003],'k','FaceAlpha',.2,'EdgeColor','none');
plot([200 400],[k+.005 k+.005],'-','Color',ColourMap.Ctrl,'LineWidth',1.5);
text(450,k+.0055,'ctrl trials','Fontsize',7,'FontName','Arial')
plot([200 400],[k+.008 k+.008],'-','Color',ColourMap.Opto,'LineWidth',1.5);
text(450,k+.0085,'opto trials','Fontsize',7,'FontName','Arial')
xlim([0 length(B)]); ylim([r k+.009]); box off
set(h,'xcolor','none'); set(h,'ycolor','none')
title('with TTX')

%% fig 1k bar plot of GCaMP with TTX

for i = 1:length(TTX.pre)
    Base(1,i) = mean(mean(TTX.pre(i).vis_resp(4*nfr:5*nfr,find(sum(TTX.pre(i).Stim.Trials == [TTX.pre(i).Stim.BlankOpto,TTX.pre(i).Stim.VisOpto],2))),2),1);
    Base(2,i) = mean(mean(TTX.drug(i).vis_resp(4*nfr:5*nfr,find(sum(TTX.drug(i).Stim.Trials == [TTX.drug(i).Stim.BlankOpto,TTX.drug(i).Stim.VisOpto],2))),2),1);
    Base(3,i) = mean(mean(TTX.post(i).vis_resp(4*nfr:5*nfr,find(sum(TTX.post(i).Stim.Trials == [TTX.post(i).Stim.BlankOpto,TTX.post(i).Stim.VisOpto],2))),2),1);
    Vis(1,i) = mean(mean(TTX.pre(i).vis_resp(5*nfr:7*nfr,TTX.pre(i).Stim.Trials == TTX.pre(i).Stim.Vis),2),1);
    Vis(2,i) = mean(mean(TTX.drug(i).vis_resp(5*nfr:7*nfr,TTX.drug(i).Stim.Trials == TTX.drug(i).Stim.Vis),2),1);
    Vis(3,i) = mean(mean(TTX.post(i).vis_resp(5*nfr:7*nfr,TTX.post(i).Stim.Trials == TTX.post(i).Stim.Vis),2),1);
end

figure(1); clf
set(gcf,'units','centimeters','position',[10 10 3 4]);
bar([1:3],mean(Base,2),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
for i = 1:size(Base,2); plot([1:3],Base(:,i),'-k'); end
xticklabels({'pre','TTX','washout'}); xtickangle(45)
title('baseline'); box off
ylabel('% \Delta F/F'); ylim([-.005 .002]); yticks([-.003 0]); yticklabels({'-.3','0'})

figure(2); clf
set(gcf,'units','centimeters','position',[15 10 3 4]);
bar([1:3],mean(Vis,2),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
for i = 1:size(Vis,2); plot([1:3],Vis(:,i),'-k'); end
xticklabels({'pre','TTX','washout'}); xtickangle(45)
title('visual response'); box off
ylabel('% \Delta F/F'); ylim([-.005 .02]); yticks([0 .01]); yticklabels({'0','1'})

%% fig 1m glusnfr example

i = 16;
figure(2); clf
set(gcf,'units','centimeters','position',[10 10 3 4]);
set(gca,'units','centimeters','position',[10 10 3 4]);
h = axes;
B = nanmean(RGC_glu(i).vis_resp(:,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis(1)),2);
C = nanstd(RGC_glu(i).vis_resp(:,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis(1)),[],2)./sqrt(sum(RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.Vis(1)));
B(1:100) = []; C(1:100) = []; k = round(max(B+C)*1000)./1000; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Ctrl,'LineWidth',1.5},ColourMap.Ctrl); hold on
B = nanmean(RGC_glu(i).vis_resp(:,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto(1)),2);
C = nanstd(RGC_glu(i).vis_resp(:,RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto(1)),[],2)./sqrt(sum(RGC_glu(i).Stim.Trials == RGC_glu(i).Stim.VisOpto(1)));
B(1:100) = []; C(1:100) = []; r = round(min(B-C)*1000)./1000; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Opto,'LineWidth',1.5},ColourMap.Opto); hold on
plot([0 200],[.005 .005],'-k','LineWidth',1); 
plot([2 2],[.005 .015],'-k','LineWidth',1);
text(80,.007,'2 s','Fontsize',7,'FontName','Arial')
text(-60,.007,'1% \DeltaF/F','rotation',90,'Fontsize',7,'FontName','Arial')

patch([100 700 700 100],[k+.001 k+.001 k+.002 k+.002],'red','FaceAlpha',.2,'EdgeColor','none');
patch([400 600 600 400],[k+.002 k+.002 k+.003 k+.003],'k','FaceAlpha',.2,'EdgeColor','none');
plot([200 400],[k+.005 k+.005],'-','Color',ColourMap.Ctrl,'LineWidth',1.5);
text(450,k+.0055,'ctrl trials','Fontsize',7,'FontName','Arial')
plot([200 400],[k+.008 k+.008],'-','Color',ColourMap.Opto,'LineWidth',1.5);
text(450,k+.0085,'opto trials','Fontsize',7,'FontName','Arial')
xlim([0 length(B)]); ylim([r k+.009]); box off
set(h,'xcolor','none'); set(h,'ycolor','none')

%% fig 1n glusnfr
% and supp fig 1P

clearvars -except GRAB GRAB_data RGC_new RGC_glu ColourMap TTX TTX_glu
t = 1; u = 1; nfr = 100;
for i = 1:length(RGC_glu)
    if strcmp(RGC_glu(i).genotype,'pet1')
        ChR(t,1) = mean(mean(RGC_glu(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_glu(i).Stim.Trials == [RGC_glu(i).Stim.Vis],2))),2),1);
        ChR(t,2) = mean(mean(RGC_glu(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_glu(i).Stim.Trials == [RGC_glu(i).Stim.VisOpto],2))),2),1);
        V_ChR(t,1) = mean(mean(RGC_glu(i).vis_resp(5*nfr:7*nfr,find(sum(RGC_glu(i).Stim.Trials == [RGC_glu(i).Stim.Vis],2))),2),1);
        V_ChR(t,2) = mean(mean(RGC_glu(i).vis_resp(5*nfr:7*nfr,find(sum(RGC_glu(i).Stim.Trials == [RGC_glu(i).Stim.VisOpto],2))),2),1);
        V_ChRS(t,1) = mean(mean(RGC_glu(i).vis_resp(5*nfr:7*nfr,find(sum(RGC_glu(i).Stim.Trials == [RGC_glu(i).Stim.Vis],2))),2),1)...
            - mean(mean(RGC_glu(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_glu(i).Stim.Trials == [RGC_glu(i).Stim.Vis],2))),2),1);
        V_ChRS(t,2) = mean(mean(RGC_glu(i).vis_resp(5*nfr:7*nfr,find(sum(RGC_glu(i).Stim.Trials == [RGC_glu(i).Stim.VisOpto],2))),2),1)...
            - mean(mean(RGC_glu(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_glu(i).Stim.Trials == [RGC_glu(i).Stim.VisOpto],2))),2),1);
        mouseP{t} = RGC_glu(i).mouse;
        t = t+1;
    end
end

% BASELINE 
figure(3); clf
bar(1,nanmean(ChR(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(2,nanmean(ChR(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot(ChR','-','Color',[.7 .7 .7]);
[~,p1] = ttest(ChR(:,1),ChR(:,2));
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot(mean(ChR(thisM,:),1)','-k','LineWidth',1)
end
Diff = ChR(:,1) - ChR(:,2);
Mouse = nominal(mouseP)';
Opto = table(Diff,Mouse);
lmePet1BaseGlu = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
plot([1,2],[.0022 .0022],'-k'); text(1.3,.0023,'***');
box off; xlim([.25 4.75]); xticks([1 2]); xticklabels({'ctrl','opto'}); xtickangle(45)
ylim([-.003 .0025]); yticks([-.002 0]); yticklabels({'-0.2','0'})
ylabel('% \Delta F/F'); title('baseline')
set(gcf,'units','centimeters','position',[20 10 2.3 4]);

% vis response
figure(4); clf
set(gcf,'units','centimeters','position',[25 10 2.3 4]);
bar(1,nanmean(V_ChR(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(2,nanmean(V_ChR(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot(V_ChR','-','Color',[.7 .7 .7])
[~,p2] = ttest(V_ChR(:,1),V_ChR(:,2));
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot(mean(V_ChR(thisM,:),1)','-k','LineWidth',1)
end
Diff = V_ChR(:,1) - V_ChR(:,2);
Mouse = nominal(mouseP)';
Opto = table(Diff,Mouse);
lmePet1VisGlu = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
plot([1,2],[.016 .016],'-k'); text(1.3,.0165,'***');
box off; xlim([.25 4.75]); xticks([1 2]); xticklabels({'ctrl','opto'}); xtickangle(45)
ylim([-.002 .018]); yticks([0 .01]); yticklabels({'0','1'})
ylabel('% \Delta F/F'); title('vis resp')

% vis response, baseline sub for supp fig 1P
figure(5); clf
set(gcf,'units','centimeters','position',[5 10 2.3 4]);
set(gca,'units','centimeters','position',[.8 .8 1 2.8]);
bar(1,nanmean(V_ChRS(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(2,nanmean(V_ChRS(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot(V_ChR','-','Color',[.7 .7 .7])
[~,p2] = ttest(V_ChRS(:,1),V_ChRS(:,2));
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot(mean(V_ChRS(thisM,:),1)','-k','LineWidth',1)
end
Diff = V_ChRS(:,1) - V_ChRS(:,2);
Mouse = nominal(mouseP)';
Opto = table(Diff,Mouse);
lmePet1VisGluSub = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
plot([1,2],[.016 .016],'-k'); text(1.3,.0165,'*');
box off; xticks([1 2]); xticklabels({'ctrl','opto'}); xtickangle(45)
ylim([-.002 .018]); yticks([0 .01]); yticklabels({'0','1'})
ylabel('% \Delta F/F'); title('vis resp')

%% fig 1o bar plot of GluSnFR with TTX

clearvars -except GRAB GRAB_data RGC_new RGC_glu ColourMap TTX TTX_glu
nfr = 100;
for i = 1:length(TTX_glu.pre)
    Base(1,i) = mean(mean(TTX_glu.pre(i).vis_resp(4*nfr:5*nfr,TTX_glu.pre(i).Stim.Trials == TTX_glu.pre(i).Stim.VisOpto),2),1);
    Base(2,i) = mean(mean(TTX_glu.drug(i).vis_resp(4*nfr:5*nfr,TTX_glu.drug(i).Stim.Trials == TTX_glu.drug(i).Stim.VisOpto),2),1);
    Base(3,i) = mean(mean(TTX_glu.post(i).vis_resp(4*nfr:5*nfr,TTX_glu.post(i).Stim.Trials == TTX_glu.post(i).Stim.VisOpto),2),1);
    Vis(1,i) = mean(mean(TTX_glu.pre(i).vis_resp(5*nfr:7*nfr,TTX_glu.pre(i).Stim.Trials == TTX_glu.pre(i).Stim.Vis),2),1);
    Vis(2,i) = mean(mean(TTX_glu.drug(i).vis_resp(5*nfr:7*nfr,TTX_glu.drug(i).Stim.Trials == TTX_glu.drug(i).Stim.Vis),2),1);
    Vis(3,i) = mean(mean(TTX_glu.post(i).vis_resp(5*nfr:7*nfr,TTX_glu.post(i).Stim.Trials == TTX_glu.post(i).Stim.Vis),2),1);
end

figure(1); clf
set(gcf,'units','centimeters','position',[10 10 2 4]);
bar([1:3],mean(Base,2),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
for i = 1:size(Base,2); plot([1:3],Base(:,i),'-k'); end
xticklabels({'pre','TTX','washout'}); xtickangle(45)
title('baseline'); box off
xlim([.25 5.25])
ylabel('% \Delta F/F'); ylim([-.003 .002]); yticks([-.002 0]); yticklabels({'-.3','0'})

figure(2); clf
set(gcf,'units','centimeters','position',[15 10 2 4]);
bar([1:3],mean(Vis,2),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
for i = 1:size(Vis,2); plot([1:3],Vis(:,i),'-k'); end
xticklabels({'pre','TTX','washout'}); xtickangle(45)
title('visual response'); box off
xlim([.25 5.25])
ylabel('% \Delta F/F'); ylim([-.005 .02]); yticks([0 .01]); yticklabels({'0','1'})

%% --------------- supplementary figure panels

%% fig S1C GRAB response in ctrl animals

clearvars -except GRAB GRAB_data RGC_new RGC_glu ColourMap TTX TTX_glu
 u = 1; nfr = 100;
for i = 1:length(GRAB)
    if GRAB(i).frequency== 20 & strcmp(GRAB(i).genotype,'WT')
        dffCtrl(u) = nanmean(mean(GRAB(i).vis_resp(4*nfr:8*nfr,GRAB(i).Stim.Trials == GRAB(i).Stim.Vis),1),2);
        dffOpto(u) = nanmean(mean(GRAB(i).vis_resp(4*nfr:8*nfr,GRAB(i).Stim.Trials == GRAB(i).Stim.VisOpto),1),2);
        mouseG{u} = GRAB(i).mouse;
        u = u+1;
    end
end

% with each run for main
figure(1); clf
set(gcf,'units','centimeters','position',[10 10 2.3 4]);
bar(1,nanmean(dffCtrl),'FaceColor',ColourMap.Ctrl,'FaceAlpha',.8,'EdgeColor','none'); hold on
bar(2,nanmean(dffOpto),'FaceColor',ColourMap.Opto,'FaceAlpha',.8,'EdgeColor','none'); hold on
plot(cat(1,dffCtrl,dffOpto),'-','Color',[.7 .7 .7]); box off; 
[~,p1] = ttest(dffCtrl,dffOpto);
M = unique(mouseG);
for i = 1:length(M)
    thisM = find(strcmp(mouseG,M{i}));
    plot(mean(cat(1,dffCtrl(thisM),dffOpto(thisM)),2),'-k','LineWidth',1)
end
Diff = dffCtrl' - dffOpto';
Mouse = nominal(mouseG)';
Opto = table(Diff,Mouse);
lmeWTGRAB = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
yticks([-.01 0 .01]); yticklabels({'-1','0','1'}); ylim([-.01 .01]); xtickangle(45)
xlim([0 3]); xticks([1 2]); xticklabels({'ctrl','opto'})
ylabel('% \DeltaF/F'); 

%% fluoxetine fig S1D

clear pre salineR fluoR
runs = [6,12,13];
saline = [19.66,27.75,10.75];
fluo = [32,40.16,22];
Fs = 100;
fluoR = NaN(30*60*Fs,3);
t=1;
for i = runs
    GRAB = GRAB_data(i).GRAB;
    pre(:,t) = GRAB(1:10*60*Fs); % 10 minutes pre saline
    F0 = mean(pre(:,t));
    pre(:,t) = (pre(:,t)-F0)/F0;
    salineR(:,t) = (GRAB((saline(t)+1)*60*Fs:(saline(t)+11)*60*Fs)-F0)/F0; % 10 minutes with saline
    R = GRAB((fluo(t))*60*Fs:end);
    fluoR(1:length(R),t) = (R-F0)/F0;
    Fluo(t) = prctile(fluoR(:,t),90);
    t = t+1;
end

Trace = cat(1,pre,NaN(1000,3),salineR,NaN(1000,3),fluoR);

figure(1); clf; h = axes;
plot(Trace,'-','Color',[.5 .5 .5],'LineWidth',.5); hold on
xlim([15 length(Trace)]); ylim([-.1 .35]);
plot([2*60*Fs 8*60*Fs],[.05 .05],'-k','LineWidth',1); 
plot([2+(2*60*Fs) 2+(2*60*Fs)],[.05 .15],'-k','LineWidth',1);
text(3*60*Fs,.08,'5 mins','Fontsize',7,'FontName','Arial')
text(20,.06,'10% \DeltaF/F','rotation',90,'Fontsize',7,'FontName','Arial')
set(gcf,'units','centimeters','position',[10 10 5 4]);  set(gca,'units','centimeters','position',[.3 .3 4.5 3.5]);
set(h,'xcolor','none'); set(h,'ycolor','none')

%% fig S1E for all runs, dff on mean of trace, std of dff distributions

framerate = 100; clear GRAB_std GRAB_std_M 
for i = 1:length(GRAB_data)
    if isfield(GRAB_data(i).injections,'saline')
        saline = GRAB_data(i).injections.saline*60*framerate;
    else
        saline = length(GRAB_data(i).GRAB);
    end
    F0 = nanmean(GRAB_data(i).GRAB(1:saline));
    GRAB_data(i).GRAB_dff = (GRAB_data(i).GRAB - F0)./F0;
    GRAB_std(i,1) = std(GRAB_data(i).GRAB_dff);
    mins = length(GRAB_data(i).GRAB)/(60*framerate);
    for j = 1:(floor(mins-20)/5)+1 % subsample to make the comparisons more accurate across mice and runs
        x(1) = 1+((j-1)*5)*60*framerate; x(2) = x(1)+(20*60*framerate);
        GRAB_std_subS{i}(j) = std(GRAB_data(i).GRAB_dff(x(1):x(2)));
    end
    GRAB_std(i,2) = mean(GRAB_std_subS{i});
    mouse_std{i} = GRAB_data(i).mouse;
end
M = unique(mouse_std);
for i = 1:length(M); GRAB_std_M(i,:) = mean(GRAB_std(strcmp(mouse_std,M{i}),:),1); end

figure(2);
bar(1,mean(GRAB_std(:,2),1),'FaceColor','none'); hold on
scatter(ones(1,length(GRAB_std)).*.8,GRAB_std(:,2),2,[.5 .5 .5]);
scatter(ones(1,length(M)).*1.2,GRAB_std_M(:,2),2,'k');
yticks([0 .02 .04]); yticklabels({'0','2','4'}); xticks([])
box off; ylabel('%\DeltaF/F')
set(gcf,'units','centimeters','position',[10 10 2.5 4]);  set(gca,'units','centimeters','position',[1 .8 1 2.5]);
GRAB_std(:,1) = []; GRAB_std_M(:,1) = []; % only keep the subsampled one!

% LME
Diff = GRAB_std;
Mouse = nominal(mouse_std)';
Opto = table(Diff,Mouse);
lme = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')

%% fig s1g trace of GCaMP with Chrimson baseline suppression

i = 3;
figure(2); clf
set(gcf,'units','centimeters','position',[10 10 3 4]);
set(gca,'units','centimeters','position',[10 10 3 4]);
h = axes;
B = nanmean(RGC_new(i).vis_resp(:,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),2);
C = nanstd(RGC_new(i).vis_resp(:,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank(1)),[],2)./sqrt(sum(RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank(1)));
B(1:100) = []; C(1:100) = []; k = round(max(B+C)*1000)./1000; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Ctrl,'LineWidth',1.5},ColourMap.Ctrl); hold on
B = nanmean(RGC_new(i).vis_resp(:,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto(1)),2);
C = nanstd(RGC_new(i).vis_resp(:,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto(1)),[],2)./sqrt(sum(RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto(1)));
B(1:100) = []; C(1:100) = []; r = round(min(B-C)*1000)./1000; 
shadedErrorBarJR([],B,C,{'color',ColourMap.Opto,'LineWidth',1.5},ColourMap.Opto); hold on
plot([0 200],[.005 .005],'-k','LineWidth',1); 
plot([2 2],[.005 .015],'-k','LineWidth',1);
text(80,.007,'2 s','Fontsize',7,'FontName','Arial')
text(-60,.007,'1% \DeltaF/F','rotation',90,'Fontsize',7,'FontName','Arial')

patch([100 700 700 100],[k+.001 k+.001 k+.002 k+.002],'red','FaceAlpha',.2,'EdgeColor','none');
patch([400 600 600 400],[k+.002 k+.002 k+.003 k+.003],'k','FaceAlpha',.2,'EdgeColor','none');
plot([200 400],[k+.005 k+.005],'-','Color',ColourMap.Ctrl,'LineWidth',1.5);
text(450,k+.0055,'ctrl trials','Fontsize',7,'FontName','Arial')
plot([200 400],[k+.008 k+.008],'-','Color',ColourMap.Opto,'LineWidth',1.5);
text(450,k+.0085,'opto trials','Fontsize',7,'FontName','Arial')
xlim([0 length(B)]); ylim([r k+.009]); box off
set(h,'xcolor','none'); set(h,'ycolor','none')

%% fig s1I, correlation of opto effect at baseline with pupil size

clearvars -except GRAB GRAB_data RGC_new RGC_glu ColourMap TTX TTX_glu

t = 1; u = 1; nfr = 100;
for i = 1:length(RGC_new)
    if RGC_new(i).frequency == 20 & ~isempty(RGC_new(i).pup_resp)
        if strcmp(RGC_new(i).genotype,'pet1')
            pup = RGC_new(i).pup_resp./max(RGC_new(i).peak_pup);
            closed = sum(isnan(pup(2*nfr:(2+RGC_new(i).prepulse+3)*nfr,:)),1) > nfr;
            pup = nanmean(pup(1:2*nfr,:),1); pup(closed) = NaN;
            if RGC_new(i).prepulse == 1;
                BaseOpto = mean(RGC_new(i).vis_resp(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1)...
                    - mean(RGC_new(i).vis_resp(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1);
                Base = mean(RGC_new(i).vis_resp(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1)...
                    - mean(RGC_new(i).vis_resp(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1);
                Popto = pup(RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto); Pctrl = pup(RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank);
            else
                BaseOpto = mean(RGC_new(i).vis_resp(1:2*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.BlankOpto,RGC_new(i).Stim.VisOpto],2))),1);...
                    - mean(RGC_new(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.BlankOpto,RGC_new(i).Stim.VisOpto],2))),1);
                Base = mean(RGC_new(i).vis_resp(1:2*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.Blank,RGC_new(i).Stim.Vis],2))),1);...
                    - mean(RGC_new(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.Blank,RGC_new(i).Stim.Vis],2))),1);
                Popto = pup(find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.BlankOpto,RGC_new(i).Stim.VisOpto],2)));
                Pctrl = pup(find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.Blank,RGC_new(i).Stim.Vis],2)));
            end
            [p,r] = corrcoef(BaseOpto(~isnan(Popto)),Popto(~isnan(Popto)));
            r_all_opto(t) = r(2,1); p_all_opto(t) = p(2,1);
            [p,r] = corrcoef(Base(~isnan(Pctrl)),Pctrl(~isnan(Pctrl)));
            r_all_ctrl(t) = r(2,1); p_all_ctrl(t) = p(2,1);
            mouseP{t} = RGC_new(i).mouse;
            t = t+1;
        end
    end
end

% with mouse averages for main
figure(1); clf
bar(1,nanmean(p_all_ctrl),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(2,nanmean(p_all_opto),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot([1,2],cat(1,p_all_ctrl,p_all_opto),'-','Color',[.7 .7 .7])
[~,p1] = ttest(p_all_ctrl,p_all_opto);
p_all = cat(1,p_all_ctrl,p_all_opto);
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot([1,2],nanmean(p_all(:,thisM),2)','-k','LineWidth',1)
end
box off; xlim([0 5]); xticks([1 2]); xticklabels({'ctrl','opto'}); xtickangle(45)
ylim([-.5 .5]); yticks([-.4 -.2 0 .2 .4]);
ylabel('Correlation coefficient'); 
set(gcf,'units','centimeters','position',[10 10 3.5 4]);

figure(2);
bar(1,nanmean(p_all_opto),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
scatter(ones(size(p_all_opto))*.8,p_all_opto,2,[.7 .7 .7]);
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    scatter(1.2,nanmean(p_all_opto(thisM)),2,'k');
end
box off; xlim([0 5]); xticks([1]); xticklabels({'Opto'}); xtickangle(45)
ylim([-.3 .3]); yticks([-.2 0 .2]);
ylabel('Correlation coefficient'); 
set(gcf,'units','centimeters','position',[10 10 3.5 4]);

%% fig s1J, GCaMP matching for pupil size; quantify baseline and vis response opto suppression
clearvars -except GRAB GRAB_data RGC_new RGC_glu ColourMap TTX TTX_glu

t = 1; u = 1; nfr = 100;
for i = 1:length(RGC_new)
    if RGC_new(i).frequency == 20 & ~isempty(RGC_new(i).pup_resp)
        if strcmp(RGC_new(i).genotype,'pet1')
            PupStd = NaN(1,length(RGC_new(i).Stim.onsets)); StdP = NaN(1,length(RGC_new(i).Stim.onsets));
            pup = RGC_new(i).pup_resp./max(RGC_new(i).peak_pup);
            pup = pup - mean(pup(1:2*nfr,:),1);
            for v = 1:length(RGC_new(i).Stim.onsets)
                PupStd(v) = mean(pup(3*nfr:5*nfr,v)./std(pup(1:2*nfr,v))); StdP(v) = std(pup(1:2*nfr,v));
            end
            PupStd = abs(PupStd)<=1; TrialFract(t) = sum(PupStd)/length(PupStd); StdPerct(t) = nanmean(StdP);
            if RGC_new(i).prepulse == 1;
                ChR(t,1) = mean(mean(RGC_new(i).vis_resp(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank & PupStd'),2),1);
                ChR(t,2) = mean(mean(RGC_new(i).vis_resp(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto & PupStd'),2),1);
                S = contains(RGC_new(i).Stim.VisStim,'Sine_0_lSF_primetime');
                ChRv(t,1) = mean(mean(RGC_new(i).vis_resp(3*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S) & PupStd'),2),1);
                ChRv(t,2) = mean(mean(RGC_new(i).vis_resp(3*nfr:5*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S) & PupStd'),2),1);
                mouseP{t} = RGC_new(i).mouse;
                t = t+1;
            else
                ChR(t,1) = mean(mean(RGC_new(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.Blank,RGC_new(i).Stim.Vis],2) & PupStd')),2),1);
                ChR(t,2) = mean(mean(RGC_new(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.BlankOpto,RGC_new(i).Stim.VisOpto],2) & PupStd')),2),1);
                S = contains(RGC_new(i).Stim.VisStim,'Sine_0_lSF_primetime');
                ChRv(t,1) = mean(mean(RGC_new(i).vis_resp(5*nfr:7*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S) & PupStd'),2),1);
                ChRv(t,2) = mean(mean(RGC_new(i).vis_resp(5*nfr:7*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S) & PupStd'),2),1);
                mouseP{t} = RGC_new(i).mouse;
                t = t+1;
            end
        end
    end
end

% baseline with mouse averages for main
figure(1); clf
bar(1,nanmean(ChR(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(2,nanmean(ChR(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot([1,2],ChR','-','Color',[.7 .7 .7])
[~,p1] = ttest(ChR(:,1),ChR(:,2));
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot([1,2],nanmean(ChR(thisM,:),1)','-k','LineWidth',1)
end
Diff = ChR(:,1) - ChR(:,2);
Mouse = nominal(mouseP)';
Opto = table(Diff,Mouse);
lmePet1VisPup = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
plot([1 2],[.005 .005],'-k'); text(1,.006,'Pet1')
plot([1 2],[.004 .004],'-k'); text(1.2,.0043,'***')
box off; xlim([0 5]); xticks([1 2]); xticklabels({'ctrl','opto'}); xtickangle(45)
ylim([-.018 .006]); yticks([-.01 0]); yticklabels({'-1','0'}) % in percentage
ylabel('% \Delta F/F'); title('baseline')
set(gcf,'units','centimeters','position',[10 10 3.5 4]);

% visual response
figure(2); clf
set(gcf,'units','centimeters','position',[15 10 3.5 4]);
bar(1,nanmean(ChRv(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(2,nanmean(ChRv(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot([1,2],ChRv','-','Color',[.7 .7 .7])
[~,p3] = ttest(ChRv(:,1),ChRv(:,2));
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot([1,2],nanmean(ChRv(thisM,:),1)','-k','LineWidth',1)
end
Diff = ChRv(:,1) - ChRv(:,2);
Mouse = nominal(mouseP)';
Opto = table(Diff,Mouse);
lmePet1VisPup = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
plot([1 2],[.056 .056],'-k'); text(1,.059,'Pet1')
plot([1 2],[.053 .053],'-k'); text(1.2,.054,'***')
box off; xlim([0 5]); xticks([1 2]); xticklabels({'ctrl','opto'}); xtickangle(45)
ylim([0 .06]); yticks([0 .02 .04]); yticklabels({'0','2','4'})
ylabel('% \Delta F/F'); title('vis resp')

%% fig S1K, effect of opto on pupil displacement

clearvars -except GRAB GRAB_data RGC_new RGC_glu ColourMap TTX TTX_glu
t = 1; u = 1; nfr = 100;
for i = 1:length(RGC_new)
    clear X Y x y 
    if RGC_new(i).frequency == 20 & ~isempty(RGC_new(i).pup_resp)
        if strcmp(RGC_new(i).genotype,'pet1')
            for v = 1:length(RGC_new(i).Stim.onsets)
                X(:,v) = RGC_new(i).pupX_100Hz(RGC_new(i).Stim.onsets(v) - 2*nfr:RGC_new(i).Stim.onsets(v) + (3+RGC_new(i).prepulse)*nfr);
                Y(:,v) = RGC_new(i).pupY_100Hz(RGC_new(i).Stim.onsets(v) - 2*nfr:RGC_new(i).Stim.onsets(v) + (3+RGC_new(i).prepulse)*nfr);
            end
            if RGC_new(i).prepulse == 1;
                x = nanmean(X(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1) - nanmean(X(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1);
                y = nanmean(Y(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1) - nanmean(Y(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1);
                deltaDisp(t,1) = mean(sqrt((x.^2)+(y.^2)));
                x = nanmean(X(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1) - nanmean(X(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1);
                y = nanmean(Y(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1) - nanmean(Y(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1);
                deltaDisp(t,2) = mean(sqrt((x.^2)+(y.^2)));
                mouseP{t} = RGC_new(i).mouse;
                t = t+1;
            else
                T = find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.Blank,RGC_new(i).Stim.Vis],2));
                x = nanmean(X(4*nfr:5*nfr,T),1) - nanmean(X(1:2*nfr,T),1);
                y = nanmean(Y(4*nfr:5*nfr,T),1) - nanmean(Y(1:2*nfr,T),1);
                deltaDisp(t,1) = mean(sqrt((x.^2)+(y.^2)));
                T = find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.BlankOpto,RGC_new(i).Stim.VisOpto],2));
                x = nanmean(X(4*nfr:5*nfr,T),1) - nanmean(X(1:2*nfr,T),1);
                y = nanmean(Y(4*nfr:5*nfr,T),1) - nanmean(Y(1:2*nfr,T),1);
                deltaDisp(t,2) = mean(sqrt((x.^2)+(y.^2)));
                mouseP{t} = RGC_new(i).mouse;
                t = t+1;
            end
        end
    end
end

% displacement
figure(3); clf
set(gcf,'units','centimeters','position',[20 10 3.5 4]);
bar(1,nanmean(deltaDisp(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(2,nanmean(deltaDisp(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot([1,2],deltaDisp','-','Color',[.7 .7 .7])
[~,p4] = ttest(deltaDisp(:,1),deltaDisp(:,2));
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot([1,2],nanmean(deltaDisp(thisM,:),1)','-k','LineWidth',1)
end
Diff = deltaDisp(:,1) - deltaDisp(:,2);
Mouse = nominal(mouseP)';
Opto = table(Diff,Mouse);
lmePupDisp = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
plot([1 2],[12 12],'-k'); text(1,12.8,'Pet1')
box off; xticks([1 2]); xticklabels({'ctrl','opto'}); xtickangle(45)
ylim([0 13]); yticks([0 5 10]);
ylabel('\Delta pupil position'); 
set(gca,'units','centimeters','position',[.8 .8 1.5 3]);

%% fig S1L, correlation of opto effect at baseline with pupil displacement

clearvars -except GRAB GRAB_data RGC_new RGC_glu ColourMap TTX TTX_glu

t = 1; u = 1; nfr = 100;
for i = 1:length(RGC_new)
    if RGC_new(i).frequency == 20 & ~isempty(RGC_new(i).pup_resp)
        if strcmp(RGC_new(i).genotype,'pet1')
            clear X Y 
            for v = 1:length(RGC_new(i).Stim.onsets)
                X(:,v) = RGC_new(i).pupX_100Hz(RGC_new(i).Stim.onsets(v) - 2*nfr:RGC_new(i).Stim.onsets(v) + (3+RGC_new(i).prepulse)*nfr);
                Y(:,v) = RGC_new(i).pupY_100Hz(RGC_new(i).Stim.onsets(v) - 2*nfr:RGC_new(i).Stim.onsets(v) + (3+RGC_new(i).prepulse)*nfr);
            end
            pupS = RGC_new(i).pup_resp./max(RGC_new(i).peak_pup);
            closed = sum(isnan(pupS(2*nfr:(2+RGC_new(i).prepulse+3)*nfr,:)),1) > nfr;
            if RGC_new(i).prepulse == 1
                BaseOpto = mean(RGC_new(i).vis_resp(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1)...
                    - mean(RGC_new(i).vis_resp(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1);
                Base = mean(RGC_new(i).vis_resp(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1)...
                    - mean(RGC_new(i).vis_resp(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1);
                x = nanmean(X(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1) - nanmean(X(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1);
                y = nanmean(Y(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1) - nanmean(Y(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank),1);
                Pctrl = sqrt((x.^2)+(y.^2));
                x = nanmean(X(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1) - nanmean(X(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1);
                y = nanmean(Y(4*nfr:6*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1) - nanmean(Y(1:2*nfr,RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto),1);
                Popto = sqrt((x.^2)+(y.^2));
            else
                BaseOpto = mean(RGC_new(i).vis_resp(1:2*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.BlankOpto,RGC_new(i).Stim.VisOpto],2))),1);...
                    - mean(RGC_new(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.BlankOpto,RGC_new(i).Stim.VisOpto],2))),1);
                Base = mean(RGC_new(i).vis_resp(1:2*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.Blank,RGC_new(i).Stim.Vis],2))),1);...
                    - mean(RGC_new(i).vis_resp(4*nfr:5*nfr,find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.Blank,RGC_new(i).Stim.Vis],2))),1);
                T = find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.Blank,RGC_new(i).Stim.Vis],2));
                x = nanmean(X(4*nfr:5*nfr,T),1) - nanmean(X(1:2*nfr,T),1);
                y = nanmean(Y(4*nfr:5*nfr,T),1) - nanmean(Y(1:2*nfr,T),1);
                Pctrl = sqrt((x.^2)+(y.^2));
                T = find(sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.BlankOpto,RGC_new(i).Stim.VisOpto],2));
                x = nanmean(X(4*nfr:5*nfr,T),1) - nanmean(X(1:2*nfr,T),1);
                y = nanmean(Y(4*nfr:5*nfr,T),1) - nanmean(Y(1:2*nfr,T),1);
                Popto = sqrt((x.^2)+(y.^2));
            end
            [p,r] = corrcoef(BaseOpto(~isnan(Popto)),Popto(~isnan(Popto)));
            r_all_opto(t) = r(2,1); p_all_opto(t) = p(2,1);
            [p,r] = corrcoef(Base(~isnan(Pctrl)),Pctrl(~isnan(Pctrl)));
            r_all_ctrl(t) = r(2,1); p_all_ctrl(t) = p(2,1);
            mouseP{t} = RGC_new(i).mouse;
            t = t+1;
        end
    end
end

% with mouse averages for main

figure(2);
bar(1,nanmean(p_all_opto),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
scatter(ones(size(p_all_opto))*.8,p_all_opto,2,[.7 .7 .7]);
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    scatter(1.2,nanmean(p_all_opto(thisM)),2,'k');
end
box off; xticks([1]); xticklabels({'Opto'}); xtickangle(45)
ylim([-.5 .5]); yticks([-.4 -.2 0 .2 .4]);
ylabel('Correlation coefficient'); 
set(gcf,'units','centimeters','position',[10 10 3.5 4]);
set(gca,'units','centimeters','position',[1.8 .8 1 3]);

%% supp fig 1M, for pupil movement trial matching

clearvars -except GRAB GRAB_data RGC_new RGC_glu ColourMap TTX TTX_glu
t = 1; u = 1; nfr = 100;
for i = 1:length(RGC_new)
    if RGC_new(i).frequency == 20 & ~isempty(RGC_new(i).pup_resp)
        if strcmp(RGC_new(i).genotype,'pet1')
            clear X Y 
            for v = 1:length(RGC_new(i).Stim.onsets)
                X(:,v) = RGC_new(i).pupX_100Hz(RGC_new(i).Stim.onsets(v) - 2*nfr:RGC_new(i).Stim.onsets(v) + (3+RGC_new(i).prepulse)*nfr);
                Y(:,v) = RGC_new(i).pupY_100Hz(RGC_new(i).Stim.onsets(v) - 2*nfr:RGC_new(i).Stim.onsets(v) + (3+RGC_new(i).prepulse)*nfr);
            end
            pupS = RGC_new(i).pup_resp./max(RGC_new(i).peak_pup);
            closed = sum(isnan(pupS(2*nfr:(2+RGC_new(i).prepulse+3)*nfr,:)),1) > nfr;
            if RGC_new(i).prepulse == 1; x2 = 6; else; x2 = 5; end
            x = nanmean(X(4*nfr:x2*nfr,:),1) - nanmean(X(1:2*nfr,:),1); y = nanmean(Y(4*nfr:x2*nfr,:),1) - nanmean(Y(1:2*nfr,:),1);
            Pctrl = sqrt((x.^2)+(y.^2)); 
            % baseline ctrl
            if RGC_new(i).prepulse == 1; 
                T = RGC_new(i).Stim.Trials == RGC_new(i).Stim.Blank;
            else
                T = sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.Blank,RGC_new(i).Stim.Vis],2)>0;
            end
            tot = sum(T); T = T' & ~closed & Pctrl < nanstd(Pctrl); Tnum(t,1) = sum(T)/tot;
            ChR(t,1) = mean(mean(RGC_new(i).vis_resp(4*nfr:x2*nfr,T),2),1);
            % baseline opto
            if RGC_new(i).prepulse == 1; 
                T = RGC_new(i).Stim.Trials == RGC_new(i).Stim.BlankOpto;
            else
                T = sum(RGC_new(i).Stim.Trials == [RGC_new(i).Stim.BlankOpto,RGC_new(i).Stim.VisOpto],2)>0;
            end
            tot = sum(T); T = T' & ~closed & Pctrl < nanstd(Pctrl); Tnum(t,2) = sum(T)/tot;
            ChR(t,2) = mean(mean(RGC_new(i).vis_resp(4*nfr:x2*nfr,T),2),1);
            % vis ctrl
            if RGC_new(i).prepulse == 1; x1 = 3; else; x1 = 5; end
            S = contains(RGC_new(i).Stim.VisStim,'Sine_0_lSF_primetime');
            T = RGC_new(i).Stim.Trials == RGC_new(i).Stim.Vis(S);
            tot = sum(T); T = T' & ~closed & Pctrl < nanstd(Pctrl); Tnum(t,3) = sum(T)/tot;
            ChRv(t,1) = mean(mean(RGC_new(i).vis_resp(x1*nfr:(x1+2)*nfr,T),2),1);
            % vis opto
            T = RGC_new(i).Stim.Trials == RGC_new(i).Stim.VisOpto(S);
            tot = sum(T); T = T' & ~closed & Pctrl < nanstd(Pctrl); Tnum(t,4) = sum(T)/tot;
            ChRv(t,2) = mean(mean(RGC_new(i).vis_resp(x1*nfr:(x1+2)*nfr,T),2),1);
            mouseP{t} = RGC_new(i).mouse;
            t = t+1;
        end
    end
end

% baseline with mouse averages
figure(1); clf
bar(1,nanmean(ChR(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(2,nanmean(ChR(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot([1,2],ChR','-','Color',[.7 .7 .7])
[~,p1] = ttest(ChR(:,1),ChR(:,2));
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot([1,2],nanmean(ChR(thisM,:),1)','-k','LineWidth',1)
end
Diff = ChR(:,1) - ChR(:,2);
Mouse = nominal(mouseP)';
Opto = table(Diff,Mouse);
lmePupDispVis = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
plot([1 2],[.007 .007],'-k'); text(1,.009,'Pet1')
plot([1 2],[.005 .005],'-k'); text(1.2,.006,'***')
box off; xticks([1 2]); xticklabels({'ctrl','opto'}); xtickangle(45)
ylim([-.02 .01]); yticks([-.02 -.01 0]); yticklabels({'-2','-1','0'}) % in percentage
ylabel('% \Delta F/F'); title('baseline')
set(gcf,'units','centimeters','position',[10 10 3.5 4]);
set(gca,'units','centimeters','position',[1 1 1.5 3]);

% visual response
figure(2); clf
set(gcf,'units','centimeters','position',[15 10 3.5 4]);
bar(1,nanmean(ChRv(:,1)),'FaceColor',ColourMap.Ctrl,'EdgeColor','none','FaceAlpha',.8); hold on
bar(2,nanmean(ChRv(:,2)),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
plot([1,2],ChRv','-','Color',[.7 .7 .7])
[~,p3] = ttest(ChRv(:,1),ChRv(:,2));
M = unique(mouseP);
for i = 1:length(M)
    thisM = find(strcmp(mouseP,M{i}));
    plot([1,2],nanmean(ChRv(thisM,:),1)','-k','LineWidth',1)
end
Diff = ChRv(:,1) - ChRv(:,2);
Mouse = nominal(mouseP)';
Opto = table(Diff,Mouse);
lmePupDispVis = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
plot([1 2],[.046 .046],'-k'); text(1,.049,'Pet1')
plot([1 2],[.043 .043],'-k'); text(1.2,.044,'***')
box off; xticks([1 2]); xticklabels({'ctrl','opto'}); xtickangle(45)
ylim([0 .05]); yticks([0 .02 .04]); yticklabels({'0','2','4'})
ylabel('% \Delta F/F'); title('vis resp')
set(gca,'units','centimeters','position',[1 1 1.5 3]);

