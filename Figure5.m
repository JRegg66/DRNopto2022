%% Figure 5

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

%% load photometry data

clear
% load Photometry.mat
clearvars -except GRAB_data GRAB
PlotColours

%% fig 5B, example trace of slow changes 

i = 4; nfr = 100;
close all
figure(1); h = axes;
plot(medfilt1(GRAB_data(i).GRAB - mean(GRAB_data(i).GRAB),nfr),'k','LineWidth',.5); hold on
xlim([1 length(GRAB_data(i).GRAB)])
plot([2*60*nfr 8*60*nfr],[.1 .1],'-k','LineWidth',1); 
text(3*60*nfr,.12,'5 mins','Fontsize',7,'FontName','Arial')
set(gcf,'units','centimeters','position',[10 15 20 5]);  set(gca,'units','centimeters','position',[.3 .3 19.5 4.5]);
set(h,'xcolor','none'); set(h,'ycolor','none')

figure(2); h = axes;
plot(medfilt1(GRAB_data(i).Pupil ./ max(GRAB_data(i).Pupil),nfr),'k','LineWidth',.5); hold on
xlim([1 length(GRAB_data(i).GRAB)])
set(gcf,'units','centimeters','position',[10 10 20 5]);  set(gca,'units','centimeters','position',[.3 .3 19.5 4.5]);
set(h,'xcolor','none'); set(h,'ycolor','none')

%% fluoxetine 

clear pre salineR fluoR
runs = [6,12,13];
saline = [19.66,27.75,10.75];
fluo = [32,40.16,22];
Fs = 100;
fluoR = NaN(30*60*Fs,3);
t=1;
for i = runs
    GRAB = GRAB_data(i).GRAB;
    R = GRAB((fluo(t))*60*Fs:end);
    fluoR(1:length(R),t) = (R-F0)/F0;
    Fluo(t) = prctile(fluoR(:,t),90);
    t = t+1;
end

%% for all runs, dff on mean of trace, std of dff distributions

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

%% grab opto data
% and redo dff

nfr = 100;
for i = 1:length(GRAB)
    F0 = mean(GRAB(i).trace);
    GRABy = (GRAB(i).trace - F0)./F0;
    for j = 1:length(GRAB(i).Stim.onsets)
        GRAB(i).vis_resp(:,j) = GRABy(GRAB(i).Stim.onsets(j)-(2*nfr):GRAB(i).Stim.onsets(j)+(13*nfr)) - mean(GRABy(GRAB(i).Stim.onsets(j)-(2*nfr):GRAB(i).Stim.onsets(j)));
    end
end

%%  dff opto response on mean of trace, report amplitude of responses
clearvars -except GRAB GRAB_data ColourMap Fluo GRAB_std GRAB_std_M GRAB_Diff

u = 1; nfr = 100;
for i = 1:length(GRAB)
    if GRAB(i).frequency== 20 & strcmp(GRAB(i).genotype,'pet1')
        if isfield(GRAB(i).Stim,'Blank');
            VO = GRAB(i).Stim.Trials == GRAB(i).Stim.BlankOpto;
        elseif isfield(GRAB(i).Stim,'Vis');
            VO = GRAB(i).Stim.Trials == GRAB(i).Stim.VisOpto;
        end
        Opto(u) = nanmean(mean(GRAB(i).vis_resp(4*nfr:8*nfr,VO),1),2);
        mouseG{u} = GRAB(i).mouse;
        u = u+1;
    end
end

% with each run for main
M = unique(mouseG);
for i = 1:length(M)
    thisM = find(strcmp(mouseG,M{i})); Opto_M(i) = mean(Opto(thisM));
end

%% ensure and shock response amplitude

% ensure
t = 1; clear G_ensure; Fs = 2500; nfr = 100;
for i = 1:length(GRAB_data)
    if ~isempty(GRAB_data(i).ensure)
        s = GRAB_data(i).ensure.ensure; clear GRAB_ensure
        for j = 1:length(s)
            G = GRAB_data(i).GRAB_dff;
            F = G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
            F0 = mean(G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
            GRAB_ensure(j,1:length(F)) = (F-F0); %./F0;
        end
        G_ensure(t,:) = nanmean(nanmean(GRAB_ensure(:,5*nfr:8*nfr),1),2); 
        mouse{t} = GRAB_data(i).mouse;
        t = t+1; clear GRAB_ensure F F0 run_ensure
    end
end

M = unique(mouse);
for i = 1:length(M); thisM = find(strcmp(mouse,M{i})); G_ensure_M(i) = nanmean(G_ensure(thisM)); end

% shock
t = 1; clear G_shock; Fs = 2500; nfr = 100;
for i = 1:length(GRAB_data)
    if ~isempty(GRAB_data(i).shock)
        s = GRAB_data(i).shock.shock; clear GRAB_shock
        for j = 1:length(s)
            G = GRAB_data(i).GRAB_dff;
            F = G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
            F0 = mean(G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
            GRAB_shock(j,1:length(F)) = F - F0; %(F-F0)./F0;
        end
        G_shock(t,:) = nanmean(nanmean(GRAB_shock(:,5*nfr:8*nfr),1),2); 
        mouse{t} = GRAB_data(i).mouse;
        t = t+1; clear GRAB_shock
    end
end

M = unique(mouse);
for i = 1:length(M); thisM = find(strcmp(mouse,M{i})); G_shock_M(i) = nanmean(G_shock(thisM)); end

%% fig 5C, correlation between pupil and GRAB
% and supp fig S5D
clearvars -except GRAB_data GRAB Fluo GRAB_std GRAB_std_M Opto Opto_M G_shock G_shock_M G_ensure G_ensure_M
close all
fr = 100;
Ysmooth = fr; t = 1;
coeffP = NaN(length(GRAB_data),2); coeffR = NaN(length(GRAB_data),2); GRAB_Diff = NaN(length(GRAB_data),2);
coeffR_Rsq = NaN(length(GRAB_data),2); coeffP_Rsq = NaN(length(GRAB_data),2);
for tr = 1:length(GRAB_data)
    if sum(isnan(GRAB_data(tr).Pupil)) < length(GRAB_data(tr).Pupil)

        G = smooth(GRAB_data(tr).GRAB_dff,Ysmooth);
        ch405 = smooth(GRAB_data(tr).ch405_F_corr,Ysmooth);
        Y = cat(2,G,ch405);
        % cut off fluoxetine administration
        Run = GRAB_data(tr).Run;
        if isfield(GRAB_data(tr).injections,'saline')
            saline = GRAB_data(tr).injections.fluoxetine*60*fr;
        else
            saline = [];
        end
        Run(saline+1:end) = []; Y(saline+1:end,:) = [];

        % running vs no running
        figure(tr); subplot(3,2,[1:2]) 
        plot([1:length(Run)]/(fr*60),Run); hold on; 
        plot([1:length(Run)]/(fr*60),abs(Run)>1);
        u = abs(Run)>1;
        drops = find(diff(abs(Run)>1)<0);
        scatter(drops/(100*60),u(drops));
        ups = find(diff(abs(Run)>1)>0);

        % does it start up or down?
        % if it starts down keep it if it's at least 1 minute
        minTime = 1*60*fr; % 1 minute
        if ups(1) < drops(1)
            if ups(1) > minTime;
                drops = cat(2,1,drops);
            else
                ups(1) = [];
            end
        end

        % does it end up or down?
        if drops(end) > ups(end)
            if (length(Run) - drops(end)) > minTime;
                ups = cat(2,ups,length(Run));
            else
                drops(end) = [];
            end
        end

        % only keep drops/ups pairs with at least 1 minute in between
        DownTime = ups - drops;
        drops(DownTime < minTime) = [];
        ups(DownTime < minTime) = [];
        DownTime(DownTime < minTime) = [];

        % should be at least 5 mins per session of no running
        minTotTime = 5*60*fr;
        if sum(DownTime) < minTotTime; drops = []; ups = []; end
        noRun = zeros(size(Run));
        for j = 1:length(drops)
            noRun(drops(j):ups(j)) = 1;
        end
        figure(tr); subplot(3,2,[1:2])
        plot([1:length(Run)]/(fr*60),noRun);

        % there should be a range of pupil values in that no running period
        % smooth pupil
        P = GRAB_data(tr).Pupil'./max(GRAB_data(tr).Pupil); 
        P(saline+1:end) = [];
        v = ~isnan(P);
        G = griddedInterpolant(find(v), P(v), 'previous');
        idx = find(~v); bp = G(idx);
        G.Method = 'next'; bn = G(idx);
        P(idx) = (bp + bn) / 2; 
        idx2 = find(isnan(P)); if ~isempty(idx2); P(idx2) = P(idx2(end)+1); end
        P = smooth(P,Ysmooth); 

        figure(tr); subplot(3,2,3:4)
        plot([1:length(Run)]/(fr*60),P); hold on;
        plot([1:length(Run)]/(fr*60),noRun);
        figure(tr); subplot(3,2,5:6)
        plot([1:length(Run)]/(fr*60),mat2gray(Y(:,1))); hold on;
        plot([1:length(Run)]/(fr*60),noRun);

        % only stationary
        if prctile(P(noRun == 1),90) > .5 & prctile(P(noRun == 1),10) < .5
            for j = 1:size(Y,2)
                [b,r] = corrcoef(Y(noRun == 1,j),P(noRun == 1));
                [~,~,~,~,stats] = regress(Y(noRun == 1,j),cat(2,ones(size(P(noRun == 1))),P(noRun == 1)));
                coeffP(tr,j) = b(2,1); coeffP_pval(tr,j) = r(2,1);
                coeffP_Rsq(tr,j) = stats(1);
                if j == 1
                GRAB_Diff(tr,1) = nanmean(Y(noRun == 1 & P' > .5,j));
                GRAB_Diff(tr,2) = nanmean(Y(noRun == 1 & P' < .5,j));
                end
            end
        end
        
        % all data
        if prctile(P(noRun == 0),90) > .5 & prctile(P(noRun == 0),10) < .5
            for j = 1:size(Y,2)
                [b,r] = corrcoef(Y(noRun == 0,j),P(noRun == 0));
                [~,~,~,~,stats] = regress(Y(noRun == 0,j),cat(2,ones(size(P(noRun == 0))),P(noRun == 0)));
                coeffR(tr,j) = b(2,1); coeffR_pval(tr,j) = r(2,1);
                coeffR_Rsq(tr,j) = stats(1);
            end
        end
        
        if sum(tr == [1,6,7,16,17,19]) > 0
            figure(30+t);
            set(gcf,'units','centimeters','position',[5 5 10 3]);
            h = axes;
            for j = 1:length(drops);
                patch([drops(j) ups(j) ups(j) drops(j)],[0 0 1 1],[.7 .7 .7],'EdgeColor','none'); hold on
            end
            plot(P,'k','LineWidth',.5);
            x(1) = fr*60*.5; x(2) = (fr*60*5)+x(1);
            plot([x(1) x(2)],[.1 .1],'-k','LineWidth',.5);
            xlim([1 length(P)]); ylim([0 1]); set(h,'xcolor','none'); set(h,'ycolor','none')

            figure(30+t+1); t = t+2;
            set(gcf,'units','centimeters','position',[5 15 10 3]);
            h = axes;
            y(1) = min(GRAB_data(tr).GRAB_dff); y(2) = max(GRAB_data(tr).GRAB_dff);
            for j = 1:length(drops);
                patch([drops(j) ups(j) ups(j) drops(j)],[y(1) y(1) y(2) y(2)],[.7 .7 .7],'EdgeColor','none'); hold on
            end
            plot(GRAB_data(tr).GRAB_dff,'k','LineWidth',.5);
            xlim([1 length(P)]); ylim([y(1) y(2)]); set(h,'xcolor','none'); box off % set(h,'ycolor','none')
        end
        clearvars -except tr GRAB_data GRAB coeffP fr Ysmooth coeffR mouse GRAB_Diff...
            Fluo GRAB_std GRAB_std_M Opto Opto_M G_shock G_shock_M...
            G_ensure G_ensure_M t coeffR_Rsq coeffP_Rsq coeffP_pval coeffR_pval
    end
    mouse{tr} = GRAB_data(tr).mouse;
end

close all
M = unique(mouse);
for j = 1:length(M); 
    coeffP_M(j,:) = nanmean(coeffP(strcmp(mouse,M{j}),:),1);
    coeffR_M(j,:) = nanmean(coeffR(strcmp(mouse,M{j}),:),1); 
    coeffP_RsqM(j,:) = nanmean(coeffP_Rsq(strcmp(mouse,M{j}),:),1);
    coeffR_RsqM(j,:) = nanmean(coeffR_Rsq(strcmp(mouse,M{j}),:),1); 
    GRAB_Diff_M(j,:) = nanmean(GRAB_Diff(strcmp(mouse,M{j}),:),1); 
end

figure(22); % look at coefficients, stationary, fig S5D1
bar(nanmean(coeffP),'FaceColor','none'); hold on
for i = 1:size(coeffP,2)
    scatter(ones(size(coeffP,1),1).*(i-.2),coeffP(:,i)',2,[.5 .5 .5]);
    scatter(ones(size(coeffP_M,1),1).*(i+.2),coeffP_M(:,i)',2,'k');
end
box off; ylim([-.6 .4])
xticklabels({'GRAB','405'}); xtickangle(45)
set(gcf,'units','centimeters','position',[15 10 3 4]);
set(gca,'units','centimeters','position',[1 .8 2.5 3]);

figure(23); % look at coefficients, stationary, fig 5C
bar(nanmean(coeffP(:,1)),'FaceColor','none'); hold on
for i = 1
    scatter(ones(size(coeffP,1),1).*(i-.2),coeffP(:,i)',2,[.5 .5 .5]);
    scatter(ones(size(coeffP_M,1),1).*(i+.2),coeffP_M(:,i)',2,'k');
end
box off; ylim([-.6 .4])
xticklabels({'GRAB'}); xtickangle(45)
set(gcf,'units','centimeters','position',[15 10 3 4]);
set(gca,'units','centimeters','position',[1 .8 1 3]);

figure(24); % look at coefficients, locomotion, fig S5D2
bar(nanmean(coeffR),'FaceColor','none'); hold on
for i = 1:size(coeffR,2)
    scatter(ones(size(coeffR,1),1).*(i-.2),coeffR(:,i)',2,[.5 .5 .5]);
    scatter(ones(size(coeffR_M,1),1).*(i+.2),coeffR_M(:,i)',2,'k');
end
box off; ylim([-.6 .4])
xticklabels({'GRAB','405'}); xtickangle(45)
set(gcf,'units','centimeters','position',[15 10 3 4]);
set(gca,'units','centimeters','position',[1 .8 2.5 3]);

figure(24); % look at Rsquare, stationary, fig S5D1
bar(nanmean(coeffP_Rsq),'FaceColor','none'); hold on
for i = 1:size(coeffP,2)
    scatter(ones(size(coeffP_Rsq,1),1).*(i-.2),coeffP_Rsq(:,i)',2,[.5 .5 .5]);
    scatter(ones(size(coeffP_RsqM,1),1).*(i+.2),coeffP_RsqM(:,i)',2,'k');
end
box off; ylim([0 .4])
xticklabels({'GRAB','405'}); xtickangle(45)
set(gcf,'units','centimeters','position',[15 10 3 4]);
set(gca,'units','centimeters','position',[1 .8 2.5 3]);

figure(25); % look at Rsquare, stationary, fig 5C
bar(nanmean(coeffP_Rsq(:,1)),'FaceColor','none'); hold on
for i = 1
    scatter(ones(size(coeffP_Rsq,1),1).*(i-.2),coeffP_Rsq(:,i)',2,[.5 .5 .5]);
    scatter(ones(size(coeffP_RsqM,1),1).*(i+.2),coeffP_RsqM(:,i)',2,'k');
end
box off; ylim([0 .4])
xticklabels({'GRAB'}); xtickangle(45)
set(gcf,'units','centimeters','position',[15 10 3 4]);
set(gca,'units','centimeters','position',[1 .8 1 3]);

figure(26); % look at Rsquare, locomotion, fig S5D2
bar(nanmean(coeffR_Rsq),'FaceColor','none'); hold on
for i = 1:size(coeffR_Rsq,2)
    scatter(ones(size(coeffR_Rsq,1),1).*(i-.2),coeffR_Rsq(:,i)',2,[.5 .5 .5]);
    scatter(ones(size(coeffR_RsqM,1),1).*(i+.2),coeffR_RsqM(:,i)',2,'k');
end
box off; ylim([0 .4])
xticklabels({'GRAB','405'}); xtickangle(45)
set(gcf,'units','centimeters','position',[15 10 3 4]);
set(gca,'units','centimeters','position',[1 .8 2.5 3]);

%% fig 5D, compare amplitudes!!

figure(1); siz = 2;
% fluoxetine
bar(1,mean(Fluo),'FaceColor','none'); hold on
scatter(ones(size(Fluo)).*1.2,Fluo,siz,'k');
% spontaneous
bar(2,mean(GRAB_std),'FaceColor','none'); hold on
scatter(ones(size(GRAB_std)).*1.8,GRAB_std,siz,[.5 .5 .5]);
scatter(ones(size(GRAB_std_M)).*2.2,GRAB_std_M,siz,'k');
% opto evoked
bar(3,mean(Opto),'FaceColor','none','EdgeColor',ColourMap.Opto);
scatter(ones(size(Opto)).*2.8,Opto,siz,[.5 .5 .5]);
scatter(ones(size(Opto_M)).*3.2,Opto_M,siz,'k');
% arousal
G = GRAB_Diff(:,1) - GRAB_Diff(:,2); G_M = GRAB_Diff_M(:,1) - GRAB_Diff_M(:,2);
bar(4,nanmean(G),'FaceColor','none');
scatter(ones(size(G)).*3.8,G,siz,[.5 .5 .5]);
scatter(ones(size(G_M)).*4.2,G_M,siz,'k');
% shock
bar(5,mean(G_shock),'FaceColor','none','EdgeColor','m');
scatter(ones(size(G_shock)).*4.8,G_shock,siz,[.5 .5 .5]);
scatter(ones(size(G_shock_M)).*5.2,G_shock_M,siz,'k');
% ensure
bar(6,mean(G_ensure),'FaceColor','none','EdgeColor',[1 .65 0]);
scatter(ones(size(G_ensure)).*5.8,G_ensure,siz,[.5 .5 .5]);
scatter(ones(size(G_ensure_M)).*6.2,G_ensure_M,siz,'k');

xticks([1:6]); xticklabels({'fluoxetine','std. spontaneous','Opto evoked','arousal','shock','ensure'}); xtickangle(45)
ylabel('%\DeltaF/F')
box off; ylim([-.04 .35]); yticks([0 .1 .2 .3]); yticklabels({'0','10','20','30'})

% redo it for the inset
inset1 = axes('Position',[.6 .6 .35 .35]);
% spontaneous
bar(1,mean(GRAB_std),'FaceColor','none'); hold on
scatter(ones(size(GRAB_std)).*.8,GRAB_std,siz,[.5 .5 .5]); scatter(ones(size(GRAB_std_M)).*1.2,GRAB_std_M,siz,'k');
% opto evoked
bar(2,mean(Opto),'FaceColor','none','EdgeColor',ColourMap.Opto);
scatter(ones(size(Opto)).*1.8,Opto,siz,[.5 .5 .5]); scatter(ones(size(Opto_M)).*2.2,Opto_M,siz,'k');
% arousal
bar(3,nanmean(G),'FaceColor','none');
scatter(ones(size(G)).*2.8,G,siz,[.5 .5 .5]); scatter(ones(size(G_M)).*3.2,G_M,siz,'k');
% shock
bar(4,mean(G_shock),'FaceColor','none','EdgeColor','m');
scatter(ones(size(G_shock)).*3.8,G_shock,siz,[.5 .5 .5]); scatter(ones(size(G_shock_M)).*4.2,G_shock_M,siz,'k');
% ensure
bar(5,mean(G_ensure),'FaceColor','none','EdgeColor',[1 .65 0]);
scatter(ones(size(G_ensure)).*4.8,G_ensure,siz,[.5 .5 .5]); scatter(ones(size(G_ensure_M)).*5.2,G_ensure_M,siz,'k');
xticks([1:5]); xticklabels({'spont','Opto','arousal','shock','ensure'}); xtickangle(45)
ylabel('%\DeltaF/F')
box off; ylim([-.04 .05]); yticks([-.02 0 .02 .04]); yticklabels({'-2','0','2','4'})
set(gcf,'units','centimeters','position',[10 10 7 5]);

%% ----- SUPPLEMENTARY PANELS OF GRAB PHOTOMETRY

%% fig S5A, motion correction  

i = 20; close all
clear G_R GRAB_shock F F0 Raw_shock Raw405_shock
Fs = 2500;
if ~isempty(GRAB_data(i).shock)
s = GRAB_data(i).shock.shock;
for j = 1:length(s)
    G = GRAB_data(i).ch465_F_smooth; G0 = nanmean(G); G = (G - G0)./G0;
    F = G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
    F0 = mean(G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
    Raw_shock(j,1:length(F)) = (F-F0);%./F0;
    G = GRAB_data(i).ch405_F_corr; G0 = nanmean(G); G = (G - G0)./G0;
    F = G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
    F0 = mean(G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
    Raw405_shock(j,1:length(F)) = (F-F0);%./F0;
    G = GRAB_data(i).GRAB_dff;
    F = G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
    F0 = mean(G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
    GRAB_shock(j,1:length(F)) = (F-F0);%./F0;
end

yl = [-.025 .01];
figure(1); h = axes;
patch([390 410 410 390],[yl(1) yl(1) yl(2) yl(2)],'m','FaceAlpha',.2,'EdgeColor','none'); hold on; 
shadedErrorBarJR([],mean(Raw_shock,1),std(Raw_shock,[],1)./sqrt(size(Raw_shock,1)),{'color','k','LineWidth',.5},'k');
shadedErrorBarJR([],mean(Raw405_shock,1),std(Raw405_shock,[],1)./sqrt(size(Raw405_shock,1)),{'--','color','k','LineWidth',.5},'k');
shadedErrorBarJR([],mean(GRAB_shock,1),std(GRAB_shock,[],1)./sqrt(size(GRAB_shock,1)),{'color',[0 .65 0],'LineWidth',.5},[0 .5 0]);
box off; ylim([yl]); 
plot([100 300],[.005 .005],'-k','LineWidth',1); 
plot([102 102],[.005 .01],'-k','LineWidth',1);
text(180,.0065,'2 s','Fontsize',7,'FontName','Arial')
text(20,.005,'.5% \DeltaF/F','rotation',90,'Fontsize',7,'FontName','Arial')
set(gcf,'units','centimeters','position',[10 10 4 4]);  set(gca,'units','centimeters','position',[.3 .3 3 3]);
set(h,'xcolor','none'); set(h,'ycolor','none')
end
%
clear G_R GRAB_ensure F F0 Raw_ensure Raw405_ensure; Fs = 2500;
if ~isempty(GRAB_data(i).ensure)
s = GRAB_data(i).ensure.ensure; yl = [-.02 .01];
for j = 1:length(s)
    G = GRAB_data(i).ch465_F_smooth; G0 = nanmean(G); G = (G - G0)./G0;
    F = G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
    F0 = mean(G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
    Raw_ensure(j,1:length(F)) = (F-F0);%./F0;
    G = GRAB_data(i).ch405_F_corr; G0 = nanmean(G); G = (G - G0)./G0;
    F = G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
    F0 = mean(G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
    Raw405_ensure(j,1:length(F)) = (F-F0);%./F0;
    G = GRAB_data(i).GRAB_dff;
    F = G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
    F0 = mean(G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
    GRAB_ensure(j,1:length(F)) = (F-F0);%./F0;
end

figure(2); h = axes;
patch([490 510 510 490],[yl(1) yl(1) yl(2) yl(2)],[1 .65 0],'FaceAlpha',.2,'EdgeColor','none'); hold on
shadedErrorBarJR([],mean(Raw_ensure,1),std(Raw_ensure,[],1)./sqrt(size(Raw_ensure,1)),{'color','k','LineWidth',.5},'k');
shadedErrorBarJR([],mean(Raw405_ensure,1),std(Raw405_ensure,[],1)./sqrt(size(Raw405_ensure,1)),{'--','color','k','LineWidth',.5},'k');
shadedErrorBarJR([],mean(GRAB_ensure,1),std(GRAB_ensure,[],1)./sqrt(size(GRAB_ensure,1)),{'color',[0 .65 0],'LineWidth',.5},[0 .5 0]);
box off; ylim([yl]); 
plot([100 300],[.005 .005],'-k','LineWidth',1); 
plot([102 102],[.005 .01],'-k','LineWidth',1);
text(180,.0065,'2 s','Fontsize',7,'FontName','Arial')
text(20,.005,'.5% \DeltaF/F','rotation',90,'Fontsize',7,'FontName','Arial')
set(gcf,'units','centimeters','position',[5 10 4 4]);  set(gca,'units','centimeters','position',[.3 .3 3 3]);
set(h,'xcolor','none'); set(h,'ycolor','none')
end

%% supp fig S5B, response to ensure 
% show change in running, pupil and photometry
% ensure response corrected

t = 1; clear G_R Run_R Lick_R Pup_R; Fs = 2500; nfr = 100;
for i = 1:length(GRAB_data)
    if ~isempty(GRAB_data(i).ensure)
        s = GRAB_data(i).ensure.ensure; clear GRAB_ensure Run_ensure Pup_ensure Lick_ensure
        for j = 1:length(s)
            G = GRAB_data(i).GRAB_dff;
            F = G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
            F0 = mean(G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
            GRAB_ensure(j,1:length(F)) = (F-F0);%./F0;
            Run = GRAB_data(i).Run;
            F = Run(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
            F0 = mean(Run(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
            Run_ensure(j,1:length(F)) = F-F0;
            Pup = GRAB_data(i).Pupil./max(GRAB_data(i).Pupil);
            F = Pup(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
            F0 = mean(Pup(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
            Pup_ensure(j,1:length(F)) = F-F0;
            Lick = GRAB_data(i).Lick;
            F = Lick(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
            F0 = mean(Lick(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
            Lick_ensure(j,1:length(F)) = F-F0;
        end
        G_R(t,:) = nanmean(GRAB_ensure,1); 
        Run_R(t,:) = nanmean(Run_ensure,1);
        Pup_R(t,:) = nanmean(Pup_ensure,1);
        Lick_R(t,:) = nanmean(Lick_ensure,1);
        mouse{t} = GRAB_data(i).mouse;
        t = t+1; clear GRAB_ensure F F0 run_ensure
    end
end

close all
figure(1);
h = axes;
plot(G_R','LineWidth',.5,'Color',[.3 .3 .3]); box off; hold on
xlim([1 size(G_R,2)]); ylim([-.01 .01])
patch([490 510 510 490],[-.01 -.01 .01 .01],'m','FaceAlpha',.2,'EdgeColor','none');
plot([100 300],[.005 .005],'-k','LineWidth',1); 
plot([102 102],[.005 .01],'-k','LineWidth',1);
text(180,.0065,'2 s','Fontsize',7,'FontName','Arial')
text(20,.005,'.5% \DeltaF/F','rotation',90,'Fontsize',7,'FontName','Arial')
set(gcf,'units','centimeters','position',[5 10 4 4]); set(gca,'units','centimeters','position',[.5 .5 3 3]);
set(h,'xcolor','none'); set(h,'ycolor','none')

figure(2);
bar(1,nanmean(nanmean(G_R(:,5*nfr:8*nfr),1),2),'FaceColor','none','EdgeColor',[1 0 1]); hold on;
scatter(ones(size(G_R,1),1)*.8,nanmean(G_R(:,5*nfr:8*nfr),2),2,[.5 .5 .5]);
M = unique(mouse);
for i = 1:length(M)
    thisM = find(strcmp(mouse,M{i}));
    scatter(1.2,nanmean(nanmean(G_R(thisM,5*nfr:8*nfr),1),2),2,'k');
end
ylim([-.01 .01]); yticks([-.01 0 .01]); yticklabels({'-1','0','1'}); ylabel({'% \DeltaF/F'}); 
xticks([1]); xticklabels({'ensure'}); xtickangle(45)
hold on; box off
set(gcf,'units','centimeters','position',[10 10 2 4]);
set(gca,'units','centimeters','position',[1 .8 1 3]);
% LME
Diff = nanmean(G_R(:,5*nfr:8*nfr),2);
Mouse = nominal(mouse)';
Opto = table(Diff,Mouse);
lme = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')

figure(3);
h = axes;
plot(Run_R','LineWidth',.5,'Color',[.3 .3 .3]); box off; hold on
xlim([1 size(Run_R,2)]); ylim([-15 4])
patch([490 510 510 490],[-15 -15 4 4],'m','FaceAlpha',.2,'EdgeColor','none');
plot([100 300],[1 1],'-k','LineWidth',1); 
plot([102 102],[1 3],'-k','LineWidth',1);
text(180,1.5,'2 s','Fontsize',7,'FontName','Arial')
text(20,1,'2 cm/s','rotation',90,'Fontsize',7,'FontName','Arial')
set(gcf,'units','centimeters','position',[15 10 4 4]); set(gca,'units','centimeters','position',[.5 .5 3 3]);
set(h,'xcolor','none'); set(h,'ycolor','none')

figure(4);
bar(1,nanmean(nanmean(Run_R(:,5*nfr:8*nfr),1),2),'FaceColor','none','EdgeColor',[1 0 1]); hold on;
scatter(ones(size(Run_R,1),1)*.8,nanmean(Run_R(:,5*nfr:8*nfr),2),2,[.5 .5 .5]);
M = unique(mouse);
for i = 1:length(M)
    thisM = find(strcmp(mouse,M{i}));
    scatter(1.2,nanmean(nanmean(Run_R(thisM,5*nfr:8*nfr),1),2),2,'k');
end
ylim([-10 2]); yticks([-10 -5 0]); ylabel({'\Delta speed (cm/s)'}); 
xticks([1]); xticklabels({'ensure'}); xtickangle(45)
hold on; box off
set(gcf,'units','centimeters','position',[20 10 2 4]);
set(gca,'units','centimeters','position',[1 .8 1 3]);
% LME
Diff = nanmean(Run_R(:,5*nfr:8*nfr),2);
Mouse = nominal(mouse)';
Opto = table(Diff,Mouse);
lme = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')

figure(5);
h = axes;
plot(Pup_R','LineWidth',.5,'Color',[.3 .3 .3]); box off; hold on
xlim([1 size(Pup_R,2)]); ylim([-.1 .1])
patch([490 510 510 490],[-.1 -.1 .1 .1],'m','FaceAlpha',.2,'EdgeColor','none');
plot([100 300],[.05 .05],'-k','LineWidth',1); 
plot([102 102],[.05 .1],'-k','LineWidth',1);
text(180,.065,'2 s','Fontsize',7,'FontName','Arial')
text(20,.065,'5%','rotation',90,'Fontsize',7,'FontName','Arial')
set(gcf,'units','centimeters','position',[25 10 4 4]); set(gca,'units','centimeters','position',[.5 .5 3 3]);
set(h,'xcolor','none'); set(h,'ycolor','none')

figure(6);
bar(1,nanmean(nanmean(Pup_R(:,5*nfr:8*nfr),1),2),'FaceColor','none','EdgeColor',[1 0 1]); hold on;
scatter(ones(size(Pup_R,1),1)*.8,nanmean(Pup_R(:,5*nfr:8*nfr),2),2,[.5 .5 .5]);
M = unique(mouse);
for i = 1:length(M)
    thisM = find(strcmp(mouse,M{i}));
    scatter(1.2,nanmean(nanmean(Pup_R(thisM,5*nfr:8*nfr),1),2),2,'k');
end
ylim([-.05 .05]); yticks([-.05 0 .05]); yticklabels({'-5','0','5'}); ylabel({'% \Delta normalized pupil size'}); 
xticklabels({'ensure'}); xtickangle(45)
hold on; box off
set(gcf,'units','centimeters','position',[30 10 2 4]);
set(gca,'units','centimeters','position',[1 .8 1 3]);
% LME
Diff = nanmean(Pup_R(:,5*nfr:8*nfr),2);
Mouse = nominal(mouse)';
Opto = table(Diff,Mouse);
lme = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')

%% supp fig s5C, response to shock
% show change in running, pupil and photometry
% ensure response corrected

clear G_R Run_R Lick_R Pup_R GRAB_ensure Run_ensure Pup_ensure Lick_ensure mouse
t = 1; Fs = 2500; nfr = 100;
for i = 1:length(GRAB_data)
    if ~isempty(GRAB_data(i).shock)
        if length(GRAB_data(i).shock.shock) > 1
        s = GRAB_data(i).shock.shock; clear GRAB_ensure Run_ensure Pup_ensure Lick_ensure
        for j = 1:length(s)
            G = GRAB_data(i).GRAB_dff;
            F = G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
            F0 = mean(G(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
            GRAB_ensure(j,1:length(F)) = (F-F0);%./F0;
            Run = GRAB_data(i).Run;
            F = Run(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
            F0 = mean(Run(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
            Run_ensure(j,1:length(F)) = F-F0;
            Pup = GRAB_data(i).Pupil./max(GRAB_data(i).Pupil);
            F = Pup(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)+(10*Fs));
            F0 = mean(Pup(GRAB_data(i).TimeStamp-5 > s(j)-(Fs*5) & GRAB_data(i).TimeStamp-5 < s(j)));
            Pup_ensure(j,1:length(F)) = F-F0;
        end
        G_R(t,:) = nanmean(GRAB_ensure,1); 
        Run_R(t,:) = nanmean(Run_ensure,1);
        Pup_R(t,:) = nanmean(Pup_ensure,1);
        mouse{t} = GRAB_data(i).mouse;
        t = t+1;
        end
    end
end
clear GRAB_ensure Run_ensure Pup_ensure Lick_ensure
close all
figure(1);
h = axes;
plot(G_R','LineWidth',.5,'Color',[.3 .3 .3]); box off; hold on
xlim([1 size(G_R,2)]); ylim([-.01 .01])
patch([490 510 510 490],[-.01 -.01 .01 .01],[1 .65 0],'FaceAlpha',.2,'EdgeColor','none');
plot([100 300],[.005 .005],'-k','LineWidth',1); 
plot([102 102],[.005 .01],'-k','LineWidth',1);
text(180,.0065,'2 s','Fontsize',7,'FontName','Arial')
text(20,.005,'.5% \DeltaF/F','rotation',90,'Fontsize',7,'FontName','Arial')
set(gcf,'units','centimeters','position',[5 10 4 4]); set(gca,'units','centimeters','position',[.5 .5 3 3]);
set(h,'xcolor','none'); set(h,'ycolor','none')

figure(2);
bar(1,nanmean(nanmean(G_R(:,5*nfr:8*nfr),1),2),'FaceColor','none','EdgeColor',[1 .65 0]); hold on;
scatter(ones(size(G_R,1),1)*.8,nanmean(G_R(:,5*nfr:8*nfr),2),2,[.5 .5 .5]);
M = unique(mouse);
for i = 1:length(M)
    thisM = find(strcmp(mouse,M{i}));
    scatter(1.2,nanmean(nanmean(G_R(thisM,5*nfr:8*nfr),1),2),2,'k');
end
ylim([-.01 .01]); yticks([-.01 0 .01]); yticklabels({'-1','0','1'}); ylabel({'% \DeltaF/F'}); 
xticks([1]); xticklabels({'shock'}); xtickangle(45)
hold on; box off
set(gcf,'units','centimeters','position',[10 10 2 4]);
set(gca,'units','centimeters','position',[1 .8 1 3]);
% LME
Diff = nanmean(G_R(:,5*nfr:8*nfr),2);
Mouse = nominal(mouse)';
Opto = table(Diff,Mouse);
lmeGRAB = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')

figure(3);
h = axes;
plot(Run_R','LineWidth',.5,'Color',[.3 .3 .3]); box off; hold on
xlim([1 size(Run_R,2)]); ylim([-2 12])
patch([490 510 510 490],[-2 -2 12 12],[1 .65 0],'FaceAlpha',.2,'EdgeColor','none');
plot([100 300],[1 1],'-k','LineWidth',1); 
plot([102 102],[1 3],'-k','LineWidth',1);
text(180,1.5,'2 s','Fontsize',7,'FontName','Arial')
text(20,1,'2 cm/s','rotation',90,'Fontsize',7,'FontName','Arial')
set(gcf,'units','centimeters','position',[15 10 4 4]); set(gca,'units','centimeters','position',[.5 .5 3 3]);
set(h,'xcolor','none'); set(h,'ycolor','none')

figure(4);
bar(1,nanmean(nanmean(Run_R(:,5*nfr:8*nfr),1),2),'FaceColor','none','EdgeColor',[1 .65 0]); hold on;
scatter(ones(size(Run_R,1),1)*.8,nanmean(Run_R(:,5*nfr:8*nfr),2),2,[.5 .5 .5]);
M = unique(mouse);
for i = 1:length(M)
    thisM = find(strcmp(mouse,M{i}));
    scatter(1.2,nanmean(nanmean(Run_R(thisM,5*nfr:8*nfr),1),2),2,'k');
end
ylim([-2 10]); yticks([0 5 10]); ylabel({'\Delta speed (cm/s)'}); 
xticks([1]); xticklabels({'shock'}); xtickangle(45)
hold on; box off
set(gcf,'units','centimeters','position',[20 10 2 4]);
set(gca,'units','centimeters','position',[1 .8 1 3]);
% LME
Diff = nanmean(Run_R(:,5*nfr:8*nfr),2);
Mouse = nominal(mouse)';
Opto = table(Diff,Mouse);
lmeRun = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')

figure(5);
h = axes;
plot(Pup_R','LineWidth',.5,'Color',[.3 .3 .3]); box off; hold on
xlim([1 size(Pup_R,2)]); ylim([-.1 .1])
patch([490 510 510 490],[-.1 -.1 .1 .1],[1 .65 0],'FaceAlpha',.2,'EdgeColor','none');
plot([100 300],[.05 .05],'-k','LineWidth',1); 
plot([102 102],[.05 .1],'-k','LineWidth',1);
text(180,.065,'2 s','Fontsize',7,'FontName','Arial')
text(20,.065,'5%','rotation',90,'Fontsize',7,'FontName','Arial')
set(gcf,'units','centimeters','position',[25 10 4 4]); set(gca,'units','centimeters','position',[.5 .5 3 3]);
set(h,'xcolor','none'); set(h,'ycolor','none')

figure(6);
bar(1,nanmean(nanmean(Pup_R(:,5*nfr:8*nfr),1),2),'FaceColor','none','EdgeColor',[1 .65 0]); hold on;
scatter(ones(size(Pup_R,1),1)*.8,nanmean(Pup_R(:,5*nfr:8*nfr),2),2,[.5 .5 .5]);
M = unique(mouse);
for i = 1:length(M)
    thisM = find(strcmp(mouse,M{i}));
    scatter(1.2,nanmean(nanmean(Pup_R(thisM,5*nfr:8*nfr),1),2),2,'k');
end
ylim([-.05 .05]); yticks([-.05 0 .05]); yticklabels({'-5','0','5'}); ylabel({'% \Delta normalized pupil size'}); 
xticklabels({'shock'}); xtickangle(45)
hold on; box off
set(gcf,'units','centimeters','position',[30 10 2 4]);
set(gca,'units','centimeters','position',[1 .8 1 3]);
% LME
Diff = nanmean(Pup_R(:,5*nfr:8*nfr),2);
Mouse = nominal(mouse)';
Opto = table(Diff,Mouse);
lmePup = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')

%% supp fig S5E, GRAB response amplitude on high and low arousal trials

% load Photometry.mat

clearvars -except GRAB 
u = 1; nfr = 100;
for i = 1:length(GRAB)
    if GRAB(i).frequency== 20 & ~isempty(GRAB(i).pup_resp)
        Arousal = nanmean(GRAB(i).pup_resp(1:2*nfr,:),1)./max(GRAB(i).peak_pup); 
        H_thresh = prctile(GRAB(i).pup_100Hz./max(GRAB(i).peak_pup),75);
        L_thresh = prctile(GRAB(i).pup_100Hz./max(GRAB(i).peak_pup),25);
        if isfield(GRAB(i).Stim,'Blank')
            dffOpto(u,1) = nanmean(mean(GRAB(i).vis_resp(4*nfr:8*nfr,GRAB(i).Stim.Trials == GRAB(i).Stim.BlankOpto & Arousal' < L_thresh),1),2);
            dffOpto(u,2) = nanmean(mean(GRAB(i).vis_resp(4*nfr:8*nfr,GRAB(i).Stim.Trials == GRAB(i).Stim.BlankOpto & Arousal' > H_thresh),1),2);
            numT(u,1) = sum(GRAB(i).Stim.Trials == GRAB(i).Stim.BlankOpto & Arousal' < L_thresh); numT(u,2) = sum(GRAB(i).Stim.Trials == GRAB(i).Stim.BlankOpto & Arousal' > H_thresh); 
            mouseG{u} = GRAB(i).mouse;
            u = u+1;
        elseif isfield(GRAB(i).Stim,'Vis')
            dffOpto(u,1) = nanmean(mean(GRAB(i).vis_resp(4*nfr:8*nfr,GRAB(i).Stim.Trials == GRAB(i).Stim.VisOpto & Arousal' < L_thresh),1),2);
            dffOpto(u,2) = nanmean(mean(GRAB(i).vis_resp(4*nfr:8*nfr,GRAB(i).Stim.Trials == GRAB(i).Stim.VisOpto & Arousal' > H_thresh),1),2);
            numT(u,1) = sum(GRAB(i).Stim.Trials == GRAB(i).Stim.VisOpto & Arousal' < L_thresh); numT(u,2) = sum(GRAB(i).Stim.Trials == GRAB(i).Stim.VisOpto & Arousal' > H_thresh);
            mouseG{u} = GRAB(i).mouse;
            u = u+1;
        end
    end
end

% with each run for main
figure(1); clf
set(gcf,'units','centimeters','position',[10 10 2.3 4]);
bar(1,nanmean(dffOpto(:,1)),'FaceColor',ColourMap.LowArousal,'FaceAlpha',.8,'EdgeColor','none'); hold on
bar(2,nanmean(dffOpto(:,2)),'FaceColor',ColourMap.HighArousal,'FaceAlpha',.8,'EdgeColor','none'); hold on
plot(dffOpto','-','Color',[.7 .7 .7])
[~,p1] = ttest(dffOpto(:,1),dffOpto(:,2));
M = unique(mouseG);
for i = 1:length(M)
    thisM = find(strcmp(mouseG,M{i}));
    plot(mean(dffOpto(thisM,:),1),'-k','LineWidth',1)
    [~,p_ind(i)] = ttest(dffOpto(thisM,1),dffOpto(thisM,2));
end
box off; 
Diff = dffOpto(:,1) - dffOpto(:,2);
Mouse = nominal(mouseG)';
Opto = table(Diff,Mouse);
lmeGRABarousal = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')
yticks([0 .01 .02]); yticklabels({'0','1','2'}); ylim([0 .027]); xtickangle(45)
xlim([0 3]); xticks([1 2]); xticklabels({'ctrl','opto'})
ylabel('% \DeltaF/F'); 

%% supp fig S5F correlation of opto response amplitude in GRAB with pupil size

clearvars -except GRAB ColourMap

t = 1; u = 1; nfr = 100;
for i = 1:length(GRAB)
    if ~isempty(GRAB(i).pup_resp)
        pup = GRAB(i).pup_resp./max(GRAB(i).peak_pup);
        closed = sum(isnan(pup(2*nfr:8*nfr,:)),1) > nfr;
        pup = nanmean(pup(1:2*nfr,:),1); pup(closed) = NaN;
        if isfield(GRAB(i).Stim,'Blank')
            BaseOpto = mean(GRAB(i).vis_resp(4*nfr:8*nfr,GRAB(i).Stim.Trials == GRAB(i).Stim.BlankOpto),1);
            Base = mean(GRAB(i).vis_resp(4*nfr:8*nfr,GRAB(i).Stim.Trials == GRAB(i).Stim.Blank),1);
            Popto = pup(GRAB(i).Stim.Trials == GRAB(i).Stim.BlankOpto);
            Pctrl = pup(GRAB(i).Stim.Trials == GRAB(i).Stim.Blank);
        elseif isfield(GRAB(i).Stim,'Vis')
            BaseOpto = mean(GRAB(i).vis_resp(4*nfr:8*nfr,GRAB(i).Stim.Trials == GRAB(i).Stim.VisOpto),1);
            Base = mean(GRAB(i).vis_resp(4*nfr:8*nfr,GRAB(i).Stim.Trials == GRAB(i).Stim.Vis),1);
            Popto = pup(GRAB(i).Stim.Trials == GRAB(i).Stim.VisOpto);
            Pctrl = pup(GRAB(i).Stim.Trials == GRAB(i).Stim.Vis);
        end
        [p,r] = corrcoef(BaseOpto(~isnan(Popto)),Popto(~isnan(Popto)));
        r_all_opto(t) = r(2,1); p_all_opto(t) = p(2,1);
        [p,r] = corrcoef(Base(~isnan(Pctrl)),Pctrl(~isnan(Pctrl)));
        r_all_ctrl(t) = r(2,1); p_all_ctrl(t) = p(2,1);
        mouseP{t} = GRAB(i).mouse;
        t = t+1;
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
% ylim([-.3 .3]); yticks([-.2 0 .2]);
ylabel('Correlation coefficient'); 
set(gcf,'units','centimeters','position',[10 10 3.5 4]); set(gca,'units','centimeters','position',[1 .8 1 3]);

%% ---------- ZOOMED IN 2P DATA
%% load data
clear
% load Zoom_2p_Data_Nov2022.mat
load('\\anastasia\data\2p\jasmine\RawDataStruct_Dec2021.mat');
PlotColours

%% pull out data

clearvars -except data dataRaw ColourMap
fr = 15.5; QI_thresh = .15; 
BigRespAll = []; QIbar = []; QIff = []; SBC = []; Reference = []; QI = [];  FOV = [];
PupilTrials = []; BaseTrials = []; ArousalCoeff = []; BigRespTrials = []; PupilAll = []; 
OOidx = []; SFpref = []; 
HighLow = [];

for i = 10:31
    if ~isempty(data(i).Opto)
    boutons = size(data(i).OptoMean,2);
    FOV = cat(1,FOV,ones(boutons,1)*i);
    clear OptoTraces OptoMean
    OptoQI = nan(boutons,length(data(i).OptoTraces));
    OptoMean = nan(204,boutons,5); OptoStd = nan(204,boutons,5);
    OptoTraces = nan(204,40,boutons,5); OptoTracesRaw = nan(204,40,boutons,5); Pupil = nan(40,5);
    for k = 1:length(data(i).OptoTraces)
        open = find(data(i).Arousal.Opto.FractClosedEye{k} < .15);
        OptoTraces(:,1:length(open),:,k) = data(i).OptoTraces{k}(:,open,:);
        OptoTracesRaw(:,1:length(open),:,k) = dataRaw(i).OptoTraces{k}(:,open,:);
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
    PupilTrials = cat(2,PupilTrials,Pt); BaseTrials = cat(2,BaseTrials,Rt);
    t = 1; AC = NaN(14,boutons); Pt = NaN(14,40,boutons); Rt = NaN(14,40,boutons); ACsig = NaN(14,boutons);
    HminL = NaN(14,boutons); HminL_prct = NaN(14,boutons);
    for k = [1,3]
        for j = 2:length(xstart)
            clear R
            R = squeeze(mean(OptoTraces(round(xstart(j)*fr):round(xend(j)*fr),:,:,k),1)); 
            R = R - minR; R = R./maxR;
            for b = 1:boutons
%                 c = regress(R(:,b),cat(2,ones(size(Pupil(:,k))),Pupil(:,k)));
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
    HighLow = cat(2,HighLow,HminL); 
    % SF pref new
    MC = data(i).Ori.MeanCurve - minR'; MC = MC./maxR';
    newSF = cat(1,newSF,squeeze(mean(MC(:,:,1),2))-squeeze(mean(MC(:,:,3),2)));
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
    Reference = cat(2,Reference,Resp);clear TotResp
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
    OOidx = cat(1,OOidx,data(i).Ret.OnOffPref);
    SFpref = cat(1,SFpref,data(i).Ori.SFprefws');
    end
end

optoQIff = max(QI(:,1:2),[],2); optoQIbar = max(QI(:,3:4),[],2);
toss = ~((QIbar == 2 | QIff == 1) & (optoQIff >= QI_thresh | optoQIbar >= QI_thresh));
SBC(toss) = []; BigRespAll(:,toss) = [];  Reference(:,toss) = []; QIbar(toss) = []; QIff(toss) = []; QI(toss,:) = []; 
FOV2 = FOV; FOV(toss) = []; 
optoQIff(toss) = []; optoQIbar(toss) = []; BaseTrials(:,toss) = []; PupilTrials(:,toss) = []; 
ArousalCoeff(:,toss) = []; BigRespTrials(:,:,toss) = []; PupilAll(:,:,toss) = []; 
 HighLow(:,toss) = [];
OOidx(QIff == 0) = NaN; OOidx(toss) = []; RetPref(QIbar < 2,:) = NaN; SFpref(toss) = []; 
SF = zeros(size(ASI)); % for j = 1:3; SF(SFpref == j,j) = 1; end
clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
    dataRaw ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    SFpref SF  HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
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

%% fig 5E, example traces

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

% bar
i = 3; Bar = 34; yl = [0 1.2]; %[7,8,12,18,23,25,28,34];
I = newFF == i & (~isnan(sum(sum(cat(3,Low,High),3),1)))';
B = BigRespTrials(:,:,I); P = PupilAll(:,:,I); A1 = A(I); O1 = O(I);
[~,idx] = sort(A(I));

figure(3);
last = length(xstart);
for i = [1,205+breaki]
    for j = 2:length(xstart)
        patch([i+xstart(j)*fr i+xend(j)*fr i+xend(j)*fr i+xstart(j)*fr],[yl(1) yl(1) yl(2) yl(2)],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none'); hold on
    end
    last = last-1;
end
patch([xstart(2)*fr xend(end)*fr xend(end)*fr xstart(2)*fr],[.5 .5 .7 .7],[1 1 1],'EdgeColor','k','FaceAlpha',0); hold on
patch([xstart(2)*fr xend(end)*fr xend(end)*fr xstart(2)*fr]+204,[.5 .5 .7 .7],[1 1 1],'EdgeColor','k','FaceAlpha',0); hold on
M = cat(1,mean(B([1:204],P(:,1,idx(Bar))<lt,idx(Bar)),2),mean(B([409:612],P(:,3,idx(Bar))<lt,idx(Bar)),2));
S = cat(1,std(B([1:204],P(:,1,idx(Bar))<lt,idx(Bar)),[],2)./sqrt(sum(P(:,1,idx(Bar))<lt)),std(B([409:612],P(:,3,idx(Bar))<lt,idx(Bar)),[],2)./sqrt(sum(P(:,3,idx(Bar))<lt)));
shadedErrorBarJR(1:204,M(1:204),S(1:204),{'color',ColourMap.LowArousal,'LineWidth',1},ColourMap.LowArousal); hold on
shadedErrorBarJR(205+breaki:408+breaki,M(205:408),S(205:408),{'color',ColourMap.LowArousal,'LineWidth',1},ColourMap.LowArousal);
M = cat(1,mean(B([1:204],P(:,1,idx(Bar))>ht,idx(Bar)),2),mean(B([409:612],P(:,3,idx(Bar))>ht,idx(Bar)),2)); 
S = cat(1,std(B([1:204],P(:,1,idx(Bar))>ht,idx(Bar)),[],2)./sqrt(sum(P(:,1,idx(Bar))>ht)),std(B([409:612],P(:,3,idx(Bar))>ht,idx(Bar)),[],2)./sqrt(sum(P(:,3,idx(Bar))>ht)));
shadedErrorBarJR(1:204,M(1:204),S(1:204),{'color',ColourMap.HighArousal,'LineWidth',1},ColourMap.HighArousal); hold on
shadedErrorBarJR(205+breaki:408+breaki,M(205:408),S(205:408),{'color',ColourMap.HighArousal,'LineWidth',1},ColourMap.HighArousal);
xlim([0 408+breaki]); ylim(yl); box off
set(gca,'xcolor','none'); set(gca,'ycolor','none')
set(gcf,'units','centimeters','position',[10 10 4 3]);

% FF example
i=1; FF = 8; yl = [-.05 1.1];%[8,46,58,66,72,211];
lt = .5; ht = .5; t = 8;
I = newFF == i & (~isnan(sum(sum(cat(3,Low,High),3),1)))';
B = BigRespTrials(:,:,I); P = PupilAll(:,:,I); A1 = A(I); O1 = O(I);
[~,idx] = sort(A(I));

figure(5);
last = length(xstart);
for i = [1,205+breaki]
    for j = 2:length(xstart)
        patch([i+xstart(j)*fr i+xend(j)*fr i+xend(j)*fr i+xstart(j)*fr],[yl(1) yl(1) yl(2) yl(2)],[.5 .5 .5],'FaceAlpha',.2,'EdgeColor','none'); hold on
    end
    last = last-1;
end
patch([xstart(2)*fr xend(end)*fr xend(end)*fr xstart(2)*fr],[.5 .5 .7 .7],[1 1 1],'EdgeColor','k','FaceAlpha',0); hold on
patch([xstart(2)*fr xend(end)*fr xend(end)*fr xstart(2)*fr]+204,[.5 .5 .7 .7],[1 1 1],'EdgeColor','k','FaceAlpha',0); hold on
M = nanmean(B(1:1020,:,idx(FF)),2); M = medfilt1(M,7); minM = min(M); maxM = max(M - min(M));
M = cat(1,mean(B([1:204],P(:,1,idx(FF))<lt,idx(FF)),2),mean(B([409:612],P(:,3,idx(FF))<lt,idx(FF)),2)); M = (M-minM)./maxM;
S = cat(1,std(B([1:204],P(:,1,idx(FF))<lt,idx(FF)),[],2)./sqrt(sum(P(:,1,idx(FF))<lt)),std(B([409:612],P(:,3,idx(FF))<lt,idx(FF)),[],2)./sqrt(sum(P(:,3,idx(FF))<lt)));
shadedErrorBarJR(1:204,M(1:204),S(1:204),{'color',ColourMap.LowArousal,'LineWidth',1},ColourMap.LowArousal); hold on
shadedErrorBarJR(205+breaki:408+breaki,M(205:408),S(205:408),{'color',ColourMap.LowArousal,'LineWidth',1},ColourMap.LowArousal);
M = cat(1,mean(B([1:204],P(:,1,idx(FF))>ht,idx(FF)),2),mean(B([409:612],P(:,3,idx(FF))>ht,idx(FF)),2)); M = (M-minM)./maxM;
S = cat(1,std(B([1:204],P(:,1,idx(FF))>ht,idx(FF)),[],2)./sqrt(sum(P(:,1,idx(FF))>ht)),std(B([409:612],P(:,3,idx(FF))>ht,idx(FF)),[],2)./sqrt(sum(P(:,3,idx(FF))>ht)));
shadedErrorBarJR(1:204,M(1:204),S(1:204),{'color',ColourMap.HighArousal,'LineWidth',1},ColourMap.HighArousal); hold on
shadedErrorBarJR(205+breaki:408+breaki,M(205:408),S(205:408),{'color',ColourMap.HighArousal,'LineWidth',1},ColourMap.HighArousal);
xlim([0 408+breaki]); ylim(yl); box off
set(gca,'xcolor','none'); set(gca,'ycolor','none')
set(gcf,'units','centimeters','position',[10 10 4 3]);

%% fig 5F arousal supp across dynamic range

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
    dataRaw ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

QI_thresh = .15;
HighLow_All = HighLow;
HighLow_All(8:13,optoQIbar < QI_thresh) = NaN;
HighLow_All(1:7,optoQIff < QI_thresh) = NaN;

CC_All = ArousalCoeff;
CC_All(8:13,optoQIbar < QI_thresh) = NaN;
CC_All(1:7,optoQIff < QI_thresh) = NaN;

lt = .5; ht = .5;
for b = 1:size(BigRespTrials,3)
    if sum(PupilAll(:,1,b)>ht) > 5 & sum(PupilAll(:,1,b)<lt) > 5
        Low(1:204,b) = squeeze(mean(BigRespTrials(1:204,PupilAll(:,1,b)<lt,b),2));
        Low(205:408,b) = squeeze(mean(BigRespTrials(409:612,PupilAll(:,3,b)<lt,b),2));
        High(1:204,b) = squeeze(mean(BigRespTrials(1:204,PupilAll(:,1,b)>ht,b),2));
        High(205:408,b) = squeeze(mean(BigRespTrials(409:612,PupilAll(:,3,b)>ht,b),2));
    else
        Low(1:408,b) = NaN;
        High(1:408,b) = NaN;
    end
end

ctrl = [3:15]; ctrlA = 1:13; % all steps of the response
Bstart = -.5; Bstep = .05; Bend = .4;
bins = [Bstart:Bstep:Bend]; plotbins = [Bstart+(Bstep/2):Bstep:Bend];
DR = [.5 .7]; LME = []; clear p Oeffect p1
for i = 1:max(newFF)
    VisIdx = newFF == i & (~isnan(sum(sum(cat(3,Low,High),3),1)))';
    thisM = mouse(VisIdx); thisF = FOV(VisIdx);
    C = Reference(ctrl,VisIdx); O = HighLow_All(ctrlA,VisIdx);
    O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); 
    thisM(isnan(O)) = []; thisF(isnan(O)) = []; O(isnan(O)) = []; Boutons(i) = length(O);
    figure(4); h = histogram(O,bins,'Normalization','cdf');
    [~,p1(i)] = kstest(O); Oeffect{i} = O;
    thisLME = NaN(Boutons(i),4); 
    thisLME(:,1) = O; thisLME(:,2) = i; thisLME(:,3) = thisM; thisLME(:,4) = thisF;
    LME = cat(1,LME,thisLME);
    figure(2); plot(plotbins,h.Values,'Color',cmap(i,:)); hold on
end
ylim([0 1]); box off; xlim([min(plotbins) max(plotbins)])
xticks([-.4 -.2 0 .2])
plot([0 0],[0 1],'--k');
xlabel({'R_{high} - R_{low}', '(Dynamic Range .5 to .7)'}); ylabel('cumulative probability')
set(gcf,'units','centimeters','position',[10 10 2 4]);

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

%% fig 5G scatterhistogram of arousal vs opto suppression through dynamic range

DR = [.5 .7];
ctrl = [3:15];
C = Reference(ctrl,:); A = HighLow_All;
A(C < DR(1) | C > DR(2)) = NaN; A = nanmean(A,1); 
ctrl = [3:15]; opto = [17:29]; % all steps of the response
C = Reference(ctrl,:); O = Reference(opto,:)-Reference(ctrl,:);
O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); 
clear ctrl opto DR lt ht b bar

figure;
x = A(~isnan(sum(sum(cat(3,Low,High),3),1)));
y = O(~isnan(sum(sum(cat(3,Low,High),3),1)));
[c,d] = corrcoef(x,y);
g = newFF(~isnan(sum(sum(cat(3,Low,High),3),1)));
scatterhist(x,y,'Group',g,'Kernel','on','Location','NorthEast',...
    'Direction','out','Marker','.','MarkerSize',.2,'Color',cmap(1:max(newFF),:),...
    'LineStyle',{'-'},'LineWidth',.5);
axis square; box off

xlabel('R_{high} - R_{low}'); ylabel('R_{opto}-R_{ctrl}')
xticks([-.5 0 .5]); yticks([-.5 0]); 
hold on; 
plot([-1 1],[0 0],'--k','LineWidth',.5);
plot([0 0],[-1 .2],'--k','LineWidth',.5);
legend(titles(1:max(newFF)))
xlim([-.5 .5]); ylim([-.8 .2])
set(gcf,'units','centimeters','position',[10 10 5 4.5]);

%% fig 5H, arousal for bar and FF

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
    dataRaw ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

HighLow_All = HighLow;
HighLow_All(8:13,optoQIbar < QI_thresh) = NaN;
HighLow_All(1:7,optoQIff < QI_thresh) = NaN;

lt = .5; ht = .5;
for b = 1:size(BigRespTrials,3)
    if sum(PupilAll(:,1,b)>ht) > 5 & sum(PupilAll(:,1,b)<lt) > 5
        Low(1:204,b) = squeeze(mean(BigRespTrials(1:204,PupilAll(:,1,b)<lt,b),2));
        Low(205:408,b) = squeeze(mean(BigRespTrials(409:612,PupilAll(:,3,b)<lt,b),2));
        High(1:204,b) = squeeze(mean(BigRespTrials(1:204,PupilAll(:,1,b)>ht,b),2));
        High(205:408,b) = squeeze(mean(BigRespTrials(409:612,PupilAll(:,3,b)>ht,b),2));
    else
        Low(1:408,b) = NaN;
        High(1:408,b) = NaN;
    end
end

ctrlA = [3:15];
DR = [.5 .7];

c = hsv(length(M)); t = 1;
for j = 1:length(M)
    F = FOV == M{j};
    for i = 1:size(F,2)
        % FF arousal
        VisIdx = newFF == 1 & F(:,i) & (~isnan(sum(sum(cat(3,Low,High),3),1)))'; 
        FFnum(t) = sum(VisIdx)/sum(F(:,i) & (~isnan(sum(sum(cat(3,Low,High),3),1)))');
        C = Reference(ctrlA,VisIdx); O = HighLow_All(:,VisIdx);
        O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); O(isnan(O)) = [];
        FF(t,2) = mean(O);
        % bar arousal
        VisIdx = newFF == 3 & F(:,i) & (~isnan(sum(sum(cat(3,Low,High),3),1)))'; 
        Bnum(t) = sum(VisIdx)/sum(F(:,i) & (~isnan(sum(sum(cat(3,Low,High),3),1)))');
        C = Reference(ctrlA,VisIdx); O = HighLow_All(:,VisIdx);
        O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); O(isnan(O)) = [];
        Bar(t,2) = mean(O);
        whichM(t) = j;
        if FFnum(t) < .05 | Bnum(t) < .05; Bar(t,:) = NaN; FF(t,:) = NaN; end
        t = t+1;
    end
end

figure(4); clf
for j = 1:length(M)
    B_m(j,:) = nanmean(Bar(whichM == j,:),1);
    F_m(j,:) = nanmean(FF(whichM == j,:),1);
end
bar(1,nanmean(F_m(:,2)),'FaceColor',cmap(1,:),'EdgeColor','none'); hold on
bar(2,nanmean(B_m(:,2)),'FaceColor',cmap(3,:),'EdgeColor','none');
[~,p] = ttest(F_m(:,2),B_m(:,2));
for i = 1:length(M)
    plot([1,2],cat(2,F_m(i,2),B_m(i,2)),'-k');
end
box off
plot([1 2],[.08 .08],'-k','LineWidth',.5); text(1.5,.085,'*');
xticks([1,2]); xticklabels({'FF','bar'})
ylim([-.15 .1]); yticks([-.1 0]); yticklabels({'-.1','0'}); xlim([.2 2.7])
set(gca,'YColor','k')
ylabel('R_{high} - R_{low}')
set(gcf,'units','centimeters','position',[20 10 2.5 4]); 

%% ------ SUPPLEMENTARY PANELS 2P IMAGING

%% supp fig 5h raw baseline

clear AC ACsig
for b = 1:size(BaseTrials,2)
    B = BaseTrials(~isnan(BaseTrials(:,b)),b); P = PupilTrials(~isnan(PupilTrials(:,b)),b);
    [c,d] = corrcoef(B(1:round(length(B)/2)),P(1:round(length(B)/2)));
    AC(1,b) = c(2,1); ACsig(1,b) = d(2,1);
    [c,d] = corrcoef(B(round(length(B)/2)+1:end),P(round(length(B)/2)+1:end));
    AC(2,b) = c(2); ACsig(2,b) = d(2,1);
end

for j = 1:max(newFF)
    for i = 1:2
        figure(4); h = histogram(AC(i,newFF == j),[-1:.1:1]);
        figure(1); plot([-.95:.1:1],h.Values,'-','Color',cmap(j,:)); hold on
    end
end
plot([0 0],[1 1500],'--k');
xlabel({'corr coef', 'BaseRawF v. Pupil'}); ylabel('number of boutons')
box off
set(gcf,'units','centimeters','position',[10 10 2 4]); 

for j = 1:max(newFF)
    for i = 1:2
        figure(4); h = histogram(AC(i,newFF == j & ACsig(i,:)'<.05),[-1:.1:1]);
        figure(2); plot([-.95:.1:1],h.Values,'-','Color',cmap(j,:)); hold on
    end
end
plot([0 0],[1 1500],'--k');
xlabel('corr coef BaseRawF v. Pupil'); ylabel('number of boutons')
box off

figure(3);
for i = 1:max(newFF)
    subplot(1,2,1);
    bar(i,mean(mean(AC(:,newFF == i))),'EdgeColor',cmap(i,:),'FaceColor','none','LineWidth',2); hold on
end

%% supp fig 5I suppression by run and by mouse

ctrl = [3:15]; ctrlA = 1:13;
DR = [.5 .7]; 

t = 1;
figure(2); clf; c = hsv(length(M)); Mbad = zeros(1,length(M));
for j = 1:length(M)
    F = FOV == M{j}; badcounter = 0;
    for i = 1:size(F,2)
        VisIdx = newFF == 1 & F(:,i) & (~isnan(sum(sum(cat(3,Low,High),3),1)))'; FFnum(t) = sum(VisIdx)/sum(F(:,i) & (~isnan(sum(sum(cat(3,Low,High),3),1)))');
        C = Reference(ctrl,VisIdx); O = HighLow_All(ctrlA,VisIdx);
        O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); O(isnan(O)) = [];
        FF(t) = mean(O); FFnum2(t) = length(O)/sum(F(:,i) & (~isnan(sum(sum(cat(3,Low,High),3),1)))');
        VisIdx = newFF == 3 & F(:,i) & (~isnan(sum(sum(cat(3,Low,High),3),1)))'; Bnum(t) = sum(VisIdx)/sum(F(:,i) & (~isnan(sum(sum(cat(3,Low,High),3),1)))');
        C = Reference(ctrl,VisIdx); O = HighLow_All(ctrlA,VisIdx);
        O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); O(isnan(O)) = [];
        Bar(t) = mean(O); Bnum2(t) = length(O)/sum(F(:,i) & (~isnan(sum(sum(cat(3,Low,High),3),1)))');
        if FFnum(t) < .05 | Bnum(t) < .05; Bar(t) = NaN; FF(t) = NaN; end % at least 5% of boutons
        if isnan(Bar(t)) | isnan(FF(t)); badcounter = badcounter+1; end
        figure(2); scatter(Bar(t),FF(t),2,c(j,:),'filled'); hold on
        t = t + 1;
    end
    if badcounter == size(F,2); Mbad(j) = 1; end
end

figure(2);
ypos = -.07;
for j = 1:length(M)
    if Mbad(j) == 0
        scatter(-.045,ypos,2,c(j,:),'filled'); 
        text(-.035,ypos,strcat('M',num2str(j)));
        ypos = ypos - .015;
    end
end
plot([-.2 .05],[-.2 .05],'--k','LineWidth',.5); 
plot([-.2 .1],[0 0],'--k'); plot([0 0],[-.2 .05],'--k','LineWidth',.5); 
xlim([-.15 .05]); ylim([-.15 .05]); 
xticks([-.1 0]); yticks([-.1 0])
xlabel('bar R_{high} - R_{Low}'); ylabel('FF R_{high} - R_{low}');
axis square; box off
set(gca, 'XAxisLocation', 'top'); set(gca, 'YAxisLocation', 'right');
set(gcf,'units','centimeters','position',[15 10 4 4]); 

%% fig S5J, SF preference

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
    dataRaw ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

QI_thresh = .15;
HighLow_All = HighLow;
HighLow_All(8:13,optoQIbar < QI_thresh) = NaN;
HighLow_All(1:7,optoQIff < QI_thresh) = NaN;

CC_All = ArousalCoeff;
CC_All(8:13,optoQIbar < QI_thresh) = NaN;
CC_All(1:7,optoQIff < QI_thresh) = NaN;

lt = .5; ht = .5;
for b = 1:size(BigRespTrials,3)
    if sum(PupilAll(:,1,b)>ht) > 5 & sum(PupilAll(:,1,b)<lt) > 5
        Low(1:204,b) = squeeze(mean(BigRespTrials(1:204,PupilAll(:,1,b)<lt,b),2));
        Low(205:408,b) = squeeze(mean(BigRespTrials(409:612,PupilAll(:,3,b)<lt,b),2));
        High(1:204,b) = squeeze(mean(BigRespTrials(1:204,PupilAll(:,1,b)>ht,b),2));
        High(205:408,b) = squeeze(mean(BigRespTrials(409:612,PupilAll(:,3,b)>ht,b),2));
    else
        Low(1:408,b) = NaN;
        High(1:408,b) = NaN;
    end
end

ctrl = [3:15]; ctrlA = 1:13;

Bstart = -.5; Bstep = .05; Bend = .4;
bins = [Bstart:Bstep:Bend]; plotbins = [Bstart+(Bstep/2):Bstep:Bend];
DR = [.5 .7]; LME = []; clear p Oeffect p1
SS = discretize(SFpref,[1,1.5,2.5,3]);
n = (newFF == 2);
cc = hsv(max(SS));
for i = 1:max(SS)
    VisIdx = SS == i & n & (~isnan(sum(sum(cat(3,Low,High),3),1)))';
    thisM = mouse(VisIdx); thisF = FOV(VisIdx);
    C = Reference(ctrl,VisIdx); O = HighLow_All(ctrlA,VisIdx);
    O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); 
    thisM(isnan(O)) = []; thisF(isnan(O)) = []; O(isnan(O)) = []; Boutons(i) = length(O);
    figure(4); h = histogram(O,bins,'Normalization','cdf');
    [~,p1(i)] = kstest(O); Oeffect{i} = O;
    thisLME = NaN(Boutons(i),4); 
    thisLME(:,1) = O; thisLME(:,2) = i; thisLME(:,3) = thisM; thisLME(:,4) = thisF;
    LME = cat(1,LME,thisLME);
    figure(2); plot(plotbins,h.Values,'Color',cc(i,:)); hold on
end
ylim([0 1]); box off; xlim([min(plotbins) max(plotbins)])
xticks([-.4 -.2 0 .2])
plot([0 0],[0 1],'--k');
xlabel({'R_{high} - R_{low}', '(Dynamic Range .5 to .7)'}); ylabel('cumulative probability')
legend({'SF1','SF2','SF3'},'location','best')
set(gcf,'units','centimeters','position',[10 10 2 4]);

%% LME for on vs off arousal suppression
Suppression = LME(:,1);
MouseLME = nominal(LME(:,3));

% SF1 reference
Type = cellstr(num2str(LME(:,2)));
Type(LME(:,2) == 1) = {'1SF1'}; Type(LME(:,2) == 2) = {'2SF2'}; Type(LME(:,2) == 3) = {'3SF3'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFF = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% SF2 reference
Type = cellstr(num2str(LME(:,2)));
Type(LME(:,2) == 1) = {'2SF1'}; Type(LME(:,2) == 2) = {'1SF2'}; Type(LME(:,2) == 3) = {'3SF3'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFFBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% SF3 reference
Type = cellstr(num2str(LME(:,2)));
Type(LME(:,2) == 1) = {'2SF1'}; Type(LME(:,2) == 2) = {'3SF2'}; Type(LME(:,2) == 3) = {'1SF3'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFFBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

%% fig S5K, On Off arousal suppression!

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
    dataRaw ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

QI_thresh = .15;
HighLow_All = HighLow;
HighLow_All(8:13,optoQIbar < QI_thresh) = NaN;
HighLow_All(1:7,optoQIff < QI_thresh) = NaN;

CC_All = ArousalCoeff;
CC_All(8:13,optoQIbar < QI_thresh) = NaN;
CC_All(1:7,optoQIff < QI_thresh) = NaN;

lt = .5; ht = .5;
for b = 1:size(BigRespTrials,3)
    if sum(PupilAll(:,1,b)>ht) > 5 & sum(PupilAll(:,1,b)<lt) > 5
        Low(1:204,b) = squeeze(mean(BigRespTrials(1:204,PupilAll(:,1,b)<lt,b),2));
        Low(205:408,b) = squeeze(mean(BigRespTrials(409:612,PupilAll(:,3,b)<lt,b),2));
        High(1:204,b) = squeeze(mean(BigRespTrials(1:204,PupilAll(:,1,b)>ht,b),2));
        High(205:408,b) = squeeze(mean(BigRespTrials(409:612,PupilAll(:,3,b)>ht,b),2));
    else
        Low(1:408,b) = NaN;
        High(1:408,b) = NaN;
    end
end

ctrl = [3:15]; ctrlA = 1:13;

Bstart = -.5; Bstep = .05; Bend = .4;
bins = [Bstart:Bstep:Bend]; plotbins = [Bstart+(Bstep/2):Bstep:Bend];
DR = [.5 .7]; LME = []; clear p Oeffect p1
OO = discretize(OOidx,[-1,-.9,0,.9,1]);
n = (newFF == 2 | newFF == 1);
cc = [ColourMap.On;ColourMap.OnOff;ColourMap.OnOff;ColourMap.Off]; 
for i = [1,4]%1:max(OO)
    VisIdx = OO == i & n & (~isnan(sum(sum(cat(3,Low,High),3),1)))';
    thisM = mouse(VisIdx); thisF = FOV(VisIdx);
    C = Reference(ctrl,VisIdx); O = HighLow_All(ctrlA,VisIdx);
    O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); 
    thisM(isnan(O)) = []; thisF(isnan(O)) = []; O(isnan(O)) = []; Boutons(i) = length(O);
    figure(4); h = histogram(O,bins,'Normalization','cdf');
    [~,p1(i)] = kstest(O); Oeffect{i} = O;
    thisLME = NaN(Boutons(i),4); 
    thisLME(:,1) = O; thisLME(:,2) = i; thisLME(:,3) = thisM; thisLME(:,4) = thisF;
    LME = cat(1,LME,thisLME);
    figure(2); plot(plotbins,h.Values,'Color',cc(i,:)); hold on
end
ylim([0 1]); box off; xlim([min(plotbins) max(plotbins)])
xticks([-.4 -.2 0 .2])
plot([0 0],[0 1],'--k');
xlabel({'R_{high} - R_{low}', '(Dynamic Range .5 to .7)'}); ylabel('cumulative probability')
legend({'Off','On',' '},'location','best'); 
set(gcf,'units','centimeters','position',[10 10 2 4]);

%% LME for on vs off arousal suppression
Suppression = LME(:,1);
MouseLME = nominal(LME(:,3));

% ON reference
Type = cellstr(num2str(LME(:,2)));
Type(LME(:,2) == 1) = {'1OFF'}; Type(LME(:,2) == 4) = {'2ON'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFF = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% OFF reference
Type = cellstr(num2str(LME(:,2)));
Type(LME(:,2) == 1) = {'2OFF'}; Type(LME(:,2) == 4) = {'1ON'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFFBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

%% fig S5L, on off opto suppression

clearvars -except data OOidx ColourMap SBC ...
    BigRespAll BigStd Reference QIbar QIff QI optoQIff optoQIbar PupilTrials BaseTrials ...
    dataRaw ArousalCoeff BigRespTrials PupilAll ArousalCoefR FOV...
    ASI DSI ASIangle DSIangle SFpref SF RetPref HighLow HighLowPrct newSF FOV2 toss BigRespAllnorm BaseDR BaseOptoDR...
    newFF mouse M Mouse cmap titles

Bstart = -.6; Bstep = .05; Bend = .2;
bins = [Bstart:Bstep:Bend]; plotbins = [Bstart+(Bstep/2):Bstep:Bend];
DR = [.5 .7]; ctrl = [2:15]; opto = [16:29]; % all steps of the response
clear p
OO = discretize(OOidx,[-1,-.9,0,.9,1]);
n = (newFF == 2 | newFF == 1);
cc = [ColourMap.On;ColourMap.OnOff;ColourMap.OnOff;ColourMap.Off]; 
LME = [];
for i = [1,4]%1:max(OO)
    VisIdx = OO == i & n;
    thisM = mouse(VisIdx); thisF = FOV(VisIdx);
    C = Reference(ctrl,VisIdx); O = Reference(opto,VisIdx)-Reference(ctrl,VisIdx);
    O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); 
    thisM(isnan(O)) = []; thisF(isnan(O)) = []; O(isnan(O)) = []; Boutons(i) = length(O);
    figure(4); h = histogram(O,bins,'Normalization','cdf');
    [~,p(i)] = kstest(O);
    thisLME = NaN(Boutons(i),4); 
    thisLME(:,1) = O; thisLME(:,2) = i; thisLME(:,3) = thisM; thisLME(:,4) = thisF;
    LME = cat(1,LME,thisLME);
    figure(2); plot(plotbins,h.Values,'Color',cc(i,:),'LineWidth',1); hold on
end
ylim([0 1]); box off; 
plot([0 0],[0 1],'--k','LineWidth',.5); xlim([min(plotbins) max(plotbins)])
xticks([-.4 -.2 0]);
xlabel({'R_{opto} - R_{ctrl}'}); ylabel('cumulative probability')
legend({'Off','On',' '},'location','best'); 
set(gcf,'units','centimeters','position',[10 10 2 4]);

%% LME for on vs off opto suppression
Suppression = LME(:,1);
MouseLME = nominal(LME(:,3));

% ON reference
Type = cellstr(num2str(LME(:,2)));
Type(LME(:,2) == 4) = {'1ON'}; Type(LME(:,2) == 1) = {'2OFF'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFF = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

% OFF reference
Type = cellstr(num2str(LME(:,2)));
Type(LME(:,2) == 4) = {'2ON'}; Type(LME(:,2) == 1) = {'1OFF'}; 
Type = nominal(Type);
Opto = table(Suppression,Type,MouseLME);
lmeFFBar = fitlme(Opto,'Suppression ~ 1 + Type + (1|MouseLME)')

%% pull out data by arousal

clearvars -except data ColourMap
fr = 15.5; QI_thresh = .15; 
BigRespAll = []; QIbar = []; QIff = []; SBC = []; Reference = []; QI = []; FOV = [];
TrialNum = []; TrialNumFOV = []; RefHigh = []; RefLow = []; RefMid = []; TrialNumLow = []; TrialNumHigh = []; TrialNumMid = []; QIarousal = [];

for i = 10:31
    if ~isempty(data(i).Opto)
    boutons = size(data(i).OptoMean,2);
    FOV = cat(1,FOV,ones(boutons,1)*i);
    clear OptoTraces OptoMean
    OptoQI = nan(boutons,length(data(i).OptoTraces));
    OptoMean = nan(204,boutons,5); 
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
    % response matrix
    BigResp = cat(1,reshape(permute(OptoMean,[1,3,2]),204*5,boutons),...
            reshape(permute(data(i).RetMean(1:62,:,:),[1,3,2]),62*18,boutons),...
            reshape(permute(data(i).OriMean(1:62,:,:,:),[1,3,4,2]),62*24,boutons));
    % high arousal
    clear OptoTraces OptoMean OptoQIhigh
    OptoMean = nan(204,boutons,5); Trialhigh = ones(5,boutons);
    for k = 1:length(data(i).OptoTraces)
        OptoTraces{k} = data(i).OptoTraces{k}(:,data(i).Arousal.Opto.PupTrial{k}>=.5 & data(i).Arousal.Opto.FractClosedEye{k} < .15,:);
        for b = 1:boutons
            OptoTraces{k}(:,:,b) = smoothdata(OptoTraces{k}(:,:,b),1,'movmean',7);
        end
        OptoMean(:,:,k) = squeeze(nanmean(OptoTraces{k},2));
        Trialhigh(k,:) = Trialhigh(k,:)*size(OptoTraces{k},2);
        meanresp = squeeze(nanmean(OptoTraces{k}(round(3*fr):end,:,:),2));
        varmeanresp = squeeze(var(meanresp,[],1));
        varresp = squeeze(var(OptoTraces{k}(round(3*fr):end,:,:),[],1,'omitnan'));
        meanvarresp = squeeze(nanmean(varresp,1));
        OptoQIhigh(:,k) = varmeanresp./meanvarresp;
        clear meanresp varmeanresp varresp meanvarresp
    end
    TrialNumHigh = cat(2,TrialNumHigh,Trialhigh); 
    RespHigh = (mean(OptoMean(1:round(2*fr),:,1),1));
    xstart = [4,5.1,6,7.1,8,9.1,10,11.1]; xend = [5,6,7,8,9,10,11,12];
    for k = 1:2
        for j = 1:length(xstart) % ctrl lum
            RespHigh = cat(1,RespHigh,(mean(OptoMean(round(xstart(j)*fr):round(xend(j)*fr),:,0+k),1)));
        end
        for j = 2:length(xstart)-1 % ctrl bar
            RespHigh = cat(1,RespHigh,(mean(OptoMean(round(xstart(j)*fr):round(xend(j)*fr),:,2+k),1)));
        end
    end
    % low arousal
    Triallow = ones(5,boutons);
    clear OptoTraces OptoMean OptoQIlow
    OptoMean = nan(204,boutons,5);
    for k = 1:length(data(i).OptoTraces)
        OptoTraces{k} = data(i).OptoTraces{k}(:,data(i).Arousal.Opto.PupTrial{k}<.5 & data(i).Arousal.Opto.FractClosedEye{k} < .15,:);
        for b = 1:boutons
            OptoTraces{k}(:,:,b) = smoothdata(OptoTraces{k}(:,:,b),1,'movmean',7);
        end
        OptoMean(:,:,k) = squeeze(nanmean(OptoTraces{k},2));
        Triallow(k,:) = Triallow(k,:)*size(OptoTraces{k},2);
        meanresp = squeeze(nanmean(OptoTraces{k}(round(3*fr):end,:,:),2));
        varmeanresp = squeeze(var(meanresp,[],1));
        varresp = squeeze(var(OptoTraces{k}(round(3*fr):end,:,:),[],1,'omitnan'));
        meanvarresp = squeeze(nanmean(varresp,1));
        OptoQIlow(:,k) = varmeanresp./meanvarresp;
        clear meanresp varmeanresp varresp meanvarresp
    end
    TrialNumLow = cat(2,TrialNumLow,Triallow); 
    RespLow = (mean(OptoMean(1:round(2*fr),:,1),1));
    xstart = [4,5.1,6,7.1,8,9.1,10,11.1]; xend = [5,6,7,8,9,10,11,12];
    for k = 1:2
        for j = 1:length(xstart) % ctrl lum
            RespLow = cat(1,RespLow,(mean(OptoMean(round(xstart(j)*fr):round(xend(j)*fr),:,0+k),1)));
        end
        for j = 2:length(xstart)-1 % ctrl bar
            RespLow = cat(1,RespLow,(mean(OptoMean(round(xstart(j)*fr):round(xend(j)*fr),:,2+k),1)));
        end
    end
    % mid arousal
    Trialmid = ones(5,boutons);
    clear OptoTraces OptoMean OptoQImid
    OptoMean = nan(204,boutons,5);
    for k = 1:length(data(i).OptoTraces)
        OptoTraces{k} = data(i).OptoTraces{k}(:,data(i).Arousal.Opto.PupTrial{k}>=.25 & data(i).Arousal.Opto.PupTrial{k}<.75 & data(i).Arousal.Opto.FractClosedEye{k} < .15,:);
        for b = 1:boutons
            OptoTraces{k}(:,:,b) = smoothdata(OptoTraces{k}(:,:,b),1,'movmean',7);
        end
        OptoMean(:,:,k) = squeeze(nanmean(OptoTraces{k},2));
        Trialmid(k,:) = Trialmid(k,:)*size(OptoTraces{k},2);
        meanresp = squeeze(nanmean(OptoTraces{k}(round(3*fr):end,:,:),2));
        varmeanresp = squeeze(var(meanresp,[],1));
        varresp = squeeze(var(OptoTraces{k}(round(3*fr):end,:,:),[],1,'omitnan'));
        meanvarresp = squeeze(nanmean(varresp,1));
        OptoQImid(:,k) = varmeanresp./meanvarresp;
        clear meanresp varmeanresp varresp meanvarresp
    end
    TrialNumMid = cat(2,TrialNumMid,Trialmid); 
    RespMid = (mean(OptoMean(1:round(2*fr),:,1),1));
    xstart = [4,5.1,6,7.1,8,9.1,10,11.1]; xend = [5,6,7,8,9,10,11,12];
    for k = 1:2
        for j = 1:length(xstart) % ctrl lum
            RespMid = cat(1,RespMid,(mean(OptoMean(round(xstart(j)*fr):round(xend(j)*fr),:,0+k),1)));
        end
        for j = 2:length(xstart)-1 % ctrl bar
            RespMid = cat(1,RespMid,(mean(OptoMean(round(xstart(j)*fr):round(xend(j)*fr),:,2+k),1)));
        end
    end
    %normalize
    % BigResp = BigResp - min(TotResp,[],1); 
    Resp = Resp - min(TotResp,[],1); 
    RespHigh = RespHigh - min(TotResp,[],1); RespLow = RespLow - min(TotResp,[],1); RespMid = RespMid - min(TotResp,[],1);
    TotResp = TotResp - min(TotResp,[],1);
    % BigResp = BigResp./max(TotResp,[],1); 
    Resp = Resp./max(TotResp,[],1); 
    RespHigh = RespHigh./max(TotResp,[],1); RespLow = RespLow./max(TotResp,[],1); RespMid = RespMid./max(TotResp,[],1);
    TotResp = TotResp./max(TotResp,[],1);
    Reference = cat(2,Reference,Resp); RefHigh = cat(2,RefHigh,RespHigh); RefLow = cat(2,RefLow,RespLow); RefMid = cat(2,RefMid,RespMid);
    BigRespAll = cat(2,BigRespAll,BigResp); clear TotResp
    QIbar = cat(1,QIbar,sum([max(data(i).Ret.QI(:,1:8),[],2) >= QI_thresh,max(data(i).Ret.QI(:,9:16),[],2) >= QI_thresh],2));
    QIff = cat(1,QIff,max(data(i).Ret.QI(:,17:18),[],2) >= QI_thresh);
    Ret = cat(2,squeeze(mean(data(i).RetMean(16:47,:,1:16),1)),reshape(squeeze(mean(data(i).OriMean(16:47,:,:,:),1)),boutons,24));
    bad = cat(2,data(i).Ret.QI(:,1:16),reshape(data(i).Ori.QI,boutons,24)) < QI_thresh;
    Ret(bad) = NaN;
    SBC = cat(1,SBC,(sum(~isnan(Ret),2) == sum(Ret<0,2)) & sum(~isnan(Ret),2)>0);
    QI = cat(1,QI,cat(2,OptoQI(:,1:4),...
        max(data(i).Ret.QI(:,1:8),[],2),max(data(i).Ret.QI(:,9:16),[],2),data(i).Ret.QI(:,17:18),...
        squeeze(max(data(i).Ori.QI,[],2))));
    QIarousal =cat(1,QIarousal,cat(2,OptoQIhigh,OptoQIlow,OptoQImid));
    end
end

optoQIff = max(QI(:,1:2),[],2); optoQIbar = max(QI(:,3:4),[],2);
toss = ~((QIbar == 2 | QIff == 1) & (optoQIff >= QI_thresh | optoQIbar >= QI_thresh));
SBC(toss) = []; FOV(toss) = []; BigRespAll(:,toss) = []; Reference(:,toss) = []; QIbar(toss) = []; QIff(toss) = []; QI(toss,:) = []; QIarousal(toss,:) = [];
optoQIff(toss) = []; optoQIbar(toss) = []; RefHigh(:,toss) = []; RefMid(:,toss) = []; RefLow(:,toss) = []; TrialNumMid(:,toss) = []; TrialNumHigh(:,toss) = []; TrialNumLow(:,toss) = [];
clearvars -except data ColourMap SBC BigRespAll Reference QIbar QIff QI optoQIff optoQIbar QIarousal...
    RefHigh RefLow RefMid TrialNumMid TrialNumHigh TrialNumLow FOV

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

%% fig 5g response vs suppression and cdf for responses below .5 pupil
% with pupil inset, line plot of resp vs suppression for supp

QI_thresh = .25; Tmin = 8;
M = (max(QIarousal(:,6:7),[],2) >= QI_thresh | max(QIarousal(:,8:9),[],2) >= QI_thresh)' & min(TrialNumLow(1:4,:),[],1) > Tmin; % criteria for keeping mid boutons
R = RefMid(:,M); F = newFF(M); mouseT = mouse(M);

ctrl = [2:15]; opto = [16:29]; % all steps of the response

% fig 4e cdf of suppression between .5 and .7 for FF, both, bar and SBC
DR = [.5 .7]; 
Bstart = -.8; Bstep = .05; Bend = .4; LME = [];
bins = [Bstart:Bstep:Bend]; plotbins = [Bstart+(Bstep/2):Bstep:Bend];
for i = 1:max(newFF)
    VisIdx = F == i;
    thisM = mouseT(VisIdx);
    C = R(ctrl,VisIdx); O = R(opto,VisIdx)-R(ctrl,VisIdx);
    O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); 
    thisM(isnan(O)) = []; O(isnan(O)) = []; Boutons(i) = length(O);
    [~,p2(i)] = kstest(O);
    thisLME = NaN(Boutons(i),3); 
    thisLME(:,1) = O; thisLME(:,2) = i; thisLME(:,3) = thisM;
    LME = cat(1,LME,thisLME);
    figure(4); h = histogram(O,bins,'Normalization','cdf');
    figure(2); plot(plotbins,h.Values,'Color',cmap(i,:),'LineWidth',.5); hold on
end
ylim([0 1]); xlim([min(plotbins) max(plotbins)]); box off; 
plot([0 0],[0 1],'--k');
xlabel('R_{opto} - R_{ctrl}'); ylabel('cumulative probability')
title('pupil < 50% of max')
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

%% fig 5G response vs suppression and cdf for responses above .5 pupil
% with pupil inset, line plot of resp vs suppression for supp

QI_thresh = .25; Tmin = 8;
M = (max(QIarousal(:,1:2),[],2) >= QI_thresh | max(QIarousal(:,3:4),[],2) >= QI_thresh)' & min(TrialNumHigh(1:4,:),[],1) > Tmin; % criteria for keeping mid boutons
R = RefMid(:,M); F = newFF(M); mouseT = mouse(M);

% ctrl = [2:2:10,11:15]; opto = [16:2:22,24:29]; % only steady state
ctrl = [2:15]; opto = [16:29]; % all steps of the response
% fig 4e cdf of suppression between .5 and .7 for FF, both, bar and SBC
DR = [.5 .7]; 
Bstart = -.8; Bstep = .05; Bend = .4;
bins = [Bstart:Bstep:Bend]; plotbins = [Bstart+(Bstep/2):Bstep:Bend]; LME = [];
for i = 1:max(newFF)
    VisIdx = F == i;
    thisM = mouseT(VisIdx);
    C = R(ctrl,VisIdx); O = R(opto,VisIdx)-R(ctrl,VisIdx);
    O(C < DR(1) | C > DR(2)) = NaN; O = nanmean(O,1); 
    thisM(isnan(O)) = []; O(isnan(O)) = []; Boutons(i) = length(O);
    [~,p3(i)] = kstest(O);
    thisLME = NaN(Boutons(i),3); 
    thisLME(:,1) = O; thisLME(:,2) = i; thisLME(:,3) = thisM;
    LME = cat(1,LME,thisLME);
    figure(4); h = histogram(O,bins,'Normalization','cdf');
    figure(2); plot(plotbins,h.Values,'Color',cmap(i,:),'LineWidth',.5); hold on
end
ylim([0 1]); xlim([min(plotbins) max(plotbins)]); box off;
plot([0 0],[0 1],'--k');
xlabel('R_{opto} - R_{ctrl}'); ylabel('cumulative probability')
title('pupil > 50% of max')
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
