%% Figure 2
% plotting script
% set default plot settings
set(groot,...
    'DefaultFigureRenderer','painters',...
    'DefaultFigureColor','w',...
    'DefaultAxesLineWidth',.5,...
    'DefaultAxesXColor','k',...
    'DefaultAxesYColor','k',...
    'DefaultAxesFontUnits','points',...
    'DefaultAxesFontSize',8,...
    'DefaultAxesFontName','Arial',...
    'DefaultAxesBox','off',...
    'DefaultAxesTickDir','in',...
    'DefaultAxesTickDirMode','manual',...
    'DefaultLineLineWidth',1,...
    'DefaultHistogramLineWidth',1.5,...
    'DefaultTextFontUnits','Points',...
    'DefaultTextFontSize',8,...
    'DefaultTextFontName','Arial');

%% load data
clear
PlotColours
% load WholeLGN_Data_Nov2022.mat

%% pull out the data

clearvars -except data ColourMap
fr = 15.5; QI_thresh = .15;
QIbar = []; QIff = []; QI = []; Base = []; XY = []; FOV = []; resp = []; allresp = []; respStd = []; SBC = [];

for i = 1:6
    boutons = size(data(i).OptoMean,2);
    OptoQI = nan(boutons,length(data(i).OptoTraces));
    OptoMean = nan(204,boutons,3);
    for k = 1:length(data(i).OptoTraces)
        OptoTraces{k} = data(i).OptoTraces{k}(:,data(i).Arousal.Opto.FractClosedEye{k} < .15,:);
        for b = 1:boutons
            OptoTraces{k}(:,:,b) = OptoTraces{k}(:,:,b);
        end
        OptoMean(:,:,k) = squeeze(nanmean(OptoTraces{k},2));
        meanresp = squeeze(nanmean(OptoTraces{k}(round(3*fr):end,:,:),2));
        varmeanresp = squeeze(var(meanresp,[],1));
        varresp = squeeze(var(OptoTraces{k}(round(3*fr):end,:,:),[],1,'omitnan'));
        meanvarresp = squeeze(nanmean(varresp,1));
        OptoQI(:,k) = varmeanresp./meanvarresp;
        clear meanresp varmeanresp varresp meanvarresp
    end
    AR = NaN(204,61,boutons); AR(:,1:size(OptoTraces{3},2),:) = OptoTraces{3};
    allresp = cat(3,allresp,AR);
    resp = cat(2,resp,OptoMean(:,:,3));
    S = squeeze(nanstd(OptoTraces{3},[],2))./sqrt(size(OptoTraces{3},2));
    respStd = cat(2,respStd,S);
    Base = cat(2,Base,squeeze(mean(mean(OptoMean(round(4*fr):round(5*fr),:,[2:3]),1),3)));
    QIbar = cat(1,QIbar,sum([max(data(i).Ret.QI(:,1:8),[],2) > QI_thresh,max(data(i).Ret.QI(:,9:16),[],2) > QI_thresh],2));
    QIff = cat(1,QIff,max(data(i).Ret.QI(:,17:18),[],2) > QI_thresh);
    QI = cat(1,QI,cat(2,OptoQI(:,1:2),...
        max(data(i).Ret.QI(:,1:8),[],2),max(data(i).Ret.QI(:,9:16),[],2),data(i).Ret.QI(:,17:18)));
    Ret = squeeze(mean(data(i).RetMean(16:47,:,1:16),1));
    bad = data(i).Ret.QI(:,1:16) < QI_thresh;
    Ret(bad) = NaN;
    S = cat(2,sum(Ret(:,1:8) < 0 & ~isnan(Ret(:,1:8)),2)>0 & sum(Ret > 0 & ~isnan(Ret),2) == 0,...
        sum(Ret(:,9:16) < 0 & ~isnan(Ret(:,9:16)),2)>0 & sum(Ret > 0 & ~isnan(Ret),2) == 0);
    SBC = cat(1,SBC,sum(S,2)>=2);
    XY = cat(1,XY,data(i).ROIxyReg);
    FOV = cat(1,FOV, ones(boutons,1).*i);
end

toss = ~((QIbar == 2 | QIff == 1)); 
QIbar(toss) = []; QIff(toss) = []; QI(toss,:) = []; resp(:,toss) = [];
Base(toss) = []; XY(toss,:) = []; allresp(:,:,toss) = []; respStd(:,toss) = []; SBC(toss) = [];
clearvars -except data QIbar QIff QI Base XY ColourMap FOV resp fr toss allresp respStd toss SBC

%% fig 2d example boutons of JR219 during baseline opto

pos = 1; F = FOV(~toss);
R = resp(:,F == pos); S = respStd(:,F == pos);
B = mean(R(round(4*fr):round(11*fr),:));
[~,idx] = sort(B);
for i = 1:size(R,2)
    R2(:,i) = smooth(R(:,i),7.5);
    S2(:,i) = smooth(S(:,i),7.5);
end

figure(3); clf
patch([31 186 186 31],[-.9 -.9 .1 .1],'red','FaceAlpha',.15,'EdgeColor','none'); hold on
for i = 1:length(chosen)
    shadedErrorBarJR([],R2(:,idx(chosen(i)))-(.2*(i-1)),S2(:,idx(chosen(i))),{'color','k','LineWidth',1},'k'); hold on
    scatter(-10,0-(.2*(i-1)),10,H(i,:),'filled');
end
xlim([-22 204]); box off; ylim([-.9 .1])
set(gcf,'units','centimeters','position',[10 10 3 4]); 
set(gcf,'units','centimeters','position',[10 10 3 4]); 

% scalebar
plot([-20 11],[-.9 -.9],'-k','LineWidth',1); 
plot([-19 -19],[-.9 -.8],'-k','LineWidth',1);
text(-5,-.86,'2 s','Fontsize',7,'FontName','Arial')
text(-30,-.9,'10% \DeltaF/F','rotation',90,'Fontsize',7,'FontName','Arial')

set(gca,'visible','off');

%% fig 2e cdf of base suppression compare 6 mice
% indicate trace of JR219

bins = [-.2:.01:.05]; plotbins = [-.19:.01:.05]; F = FOV(~toss);

for i = 1:6
    figure(2); h = histogram(Base(F == i),bins,'Normalization','cdf');
    [~,p(i)] = kstest(Base(F == i));
    if i == 1
        figure(1); plot(plotbins,h.Values,'-r','LineWidth',.5); hold on
    else
        figure(1); plot(plotbins,h.Values,'-k','LineWidth',.5);
    end
end
plot([0 0],[0 1],'--k','LineWidth',.5); box off; xlim([min(plotbins) max(plotbins)])
xlabel('\Delta F/F'); ylabel('cumulative probability')
set(gcf,'units','centimeters','position',[10 10 3 4]); 

Mouse = nominal(FOV); Base = Base';
Opto = table(Base,Mouse);
lmeBaseOpto = fitlme(Opto,'Base ~ 1 + (1|Mouse)')

%% fig 2g scatter of Base supp of all boutons from 6 mice
% indicate bouton of the trace!

clear cmap2
edges = [-.1:.025:.1]; ref = edges(1:end-1)+.0125;
cmap = ones(length(edges)-1,3);
b = [0:1/(size(cmap,1)/2):1];
cmap(1:size(cmap,1)/2,3) = fliplr(b(2:end));
cmap((size(cmap,1)/2)+1:end,1) = b(2:end);
cmap = [0 0 .5; 0 0 1; 1 1 1; 1 0 0; .5 0 0]; % bluewhitered standard!
for i = 1:3
    cmap2(:,i) = interp1([1:5],cmap(:,i),[1:5/length(edges):5]);
end

A = Base;
A(A >= .1) = .1; A(A <= -.1) = -.1; 
A = discretize(A,edges);
figure(1);
for j = 1:max(A)
    scatter(XY(A == j,1),XY(A == j,2),.15,cmap2(j,:),'filled'); hold on
end
    
figure(1); set(gca,'YDir','reverse'); daspect([1 1 1]); set(gca,'visible','off');
set(gcf,'units','centimeters','position',[10 10 4.5 3]); 
set(gcf,'units','centimeters','position',[10 10 4.5 3]); 

figure(2); % this is an imperfect colormap!!
imagesc(ref'); colormap(bluewhitered); ylabel('baseline \DeltaF/F')
set(gca,'xtick',[])
yticks([1 length(ref)])
yticklabels({'-.1', '.1'})
set(gcf,'units','centimeters','position',[20 10 .5 2]); 
set(gcf,'units','centimeters','position',[20 10 .1 2]); 

%% ------------ SUPPLEMENTARY PANELS

%% supp fig 2a and 2b - align ipsi patches

scalebarLength = 1.25*120; % 120 pixels is 200um
[optimizer, metric] = imregconfig('monomodal');
t = 1; clear Anatomy
for i = [1,2,4,5]
    if i > 1
        tform = imregtform(data(i).anatomy(:,:,3), data(1).anatomy(:,:,3), 'rigid', optimizer, metric);
        for j = 1:4
            Anatomy(:,:,j,t) = imwarp(data(i).anatomy(:,:,j),tform,'OutputView',imref2d(size(data(i).anatomy(:,:,3))));
        end
    else
        Anatomy(:,:,:,t) = data(i).anatomy;
    end
    t = t+1;
end
Anatomy2 = Anatomy;
%%
Anatomy = Anatomy2;
for i = 1:4 
    A = Anatomy(:,:,1,i); A(A>20000) = 20000; A(A<2500) = NaN; Anatomy(:,:,1,i) = A; 
    A = Anatomy(:,:,2,i); A(A>3200) = 3200;  A(A<2750) = NaN; Anatomy(:,:,2,i) = A; 
end

figure(3);
imshowpair(mat2gray(nanmean(Anatomy(:,:,1,:),4)),mat2gray(nanmean(Anatomy(:,:,2,:),4))); hold on
plot([30 30+scalebarLength],[512-30 512-30],'-w','Linewidth',1); text(30,512-70,'250\mum','Color','w')
set(gcf,'units','centimeters','position',[20 10 4 3]);
set(gca,'units','centimeters','position',[0 0 4 2.9]);

% AP averaging
clear AP_*
AP_gcamp = squeeze(nanmean(Anatomy(:,:,1,:),1)); AP_ctb = squeeze(nanmean(Anatomy(:,:,2,:),1)); 
bins = [0:60:size(AP_gcamp,1)]; if max(bins) < size(AP_gcamp,1); bins(end+1) = size(AP_gcamp,1); end
plotbins = bins(1:end-1)+(diff(bins)/2);
AP = [1:size(AP_gcamp,1)]; AP = discretize(AP,bins);
for i = 1:size(AP_gcamp,2)
    for j = 1:max(AP)
        AP_contra(j,i) = nanmean(AP_gcamp(AP == j,i)); 
        AP_ipsi(j,i) = nanmean(AP_ctb(AP == j,i)); 
    end
    AP_contra(:,i) = mat2gray(AP_contra(:,i));
    AP_ipsi(:,i) = mat2gray(AP_ipsi(:,i));
end

% ML averaging
clear ML_*
ML_gcamp = squeeze(nanmean(Anatomy(:,:,1,:),2)); ML_ctb = squeeze(nanmean(Anatomy(:,:,2,:),2)); 
binsML = [0:60:size(ML_gcamp,1)]; if max(bins) < size(ML_gcamp,1); bins(end+1) = size(ML_gcamp,1); end
plotbinsML = binsML(1:end-1)+(diff(binsML)/2);
ML = [1:size(ML_gcamp,1)]; ML = discretize(ML,binsML);
for i = 1:size(ML_gcamp,2)
    for j = 1:max(ML)
        ML_contra(j,i) = nanmean(ML_gcamp(ML == j,i)); 
        ML_ipsi(j,i) = nanmean(ML_ctb(ML == j,i)); 
    end
    ML_contra(:,i) = mat2gray(ML_contra(:,i));
    ML_ipsi(:,i) = mat2gray(ML_ipsi(:,i));
end

figure(1); clf % AP average
for i = 1:size(AP_ipsi,2)
    plot(plotbins*(5/3),AP_contra(:,i),'Color','g'); hold on
    plot(plotbins*(5/3),AP_ipsi(:,i),'Color','m'); hold on
end
xlim([min(plotbins) max(plotbins)]*(5/3))
yticks([0 1])
box off; axis square
xlabel('dist. anterior edge (\mum)')
ylabel('norm. fluorescence (a.u.)')
set(gcf,'units','centimeters','position',[10 10 4 4]); 

figure(2); clf % AP average
for i = 1:size(ML_ipsi,2)
    plot(plotbinsML*(5/3),ML_contra(:,i),'Color','g'); hold on
    plot(plotbinsML*(5/3),ML_ipsi(:,i),'Color','m'); hold on
end
xlim([min(plotbinsML) max(plotbinsML)]*(5/3))
yticks([0 1])
box off; axis square
xlabel('dist. medial edge (\mum)')
ylabel('norm. fluorescence (a.u.)')
set(gcf,'units','centimeters','position',[15 10 4 4]); 

%% supp fig 2c, 2d and 2e - Allen data

AP_2p = (size(Anatomy,2)/60)*100; % 60pix is 100um
ML_2p = (size(Anatomy,1)/60)*100; % 60pix is 100um

slices = [82:120];
x = [270,325]; y = [280, 325]; z = [6,9];

% vglut2 ipsi and contra transverse AVERAGE
clear Kmod Fmod F_stack K_stack
% link path to these files
% F = 'vglut2_rottransverse.tif';
% K = 'vglut2_rottransverse_flp.tif';
% M = 'dLGNmask\dLGNmaskrot_transverseslice';
t = 1;
for i = slices
    F_stack(:,:,t) = rgb2gray(imread(F,i));
    K_stack(:,:,t) = rgb2gray(imread(K,i));
    Mpath = strcat(M,num2str(i),'.tiff');
    M_stack(:,:,t) = sum(imread(Mpath),3)>409;
    t = t+1;
end

Kmod = mean(mat2gray(K_stack(x(1):x(2),y(1):y(2),z(1):z(2))),3)'; AP_allen = size(Kmod,2)*25; ML_allen = size(Kmod,1)*25; % 25um per pix
Fmod = mean(mat2gray(F_stack(x(1):x(2),y(1):y(2),z(1):z(2))),3)'; 
Gmod = zeros(size(Kmod));
RGBmod = cat(3,Kmod,Fmod,Gmod);
figure(4); imshow(RGBmod); hold on
scalebarLength = 10;
plot([2 2+scalebarLength],[size(RGBmod,1)-2 size(RGBmod,1)-2],'-w','Linewidth',1); text(5,size(RGBmod,1)-5,'250\mum','Color','w')
set(gcf,'units','centimeters','position',[5 10 4 4]); 
set(gca,'units','centimeters','position',[0 0 4 4]); 

% extended view
Y = 280+(4*4); X = 270+(6*4);
C_extended = mean(mat2gray(K_stack(X-(1500/25):X+(1500/25),Y-(1500/25):Y+(1500/25),z(1):z(2))),3)';
I_extended = mean(mat2gray(F_stack(X-(1500/25):X+(1500/25),Y-(1500/25):Y+(1500/25),z(1):z(2))),3)';
G_extended = zeros(size(C_extended));
M_extended = sum(mat2gray(M_stack(X-(1500/25):X+(1500/25),Y-(1500/25):Y+(1500/25),z(1):z(2))),3)'>0;
RGB_entended = cat(3,C_extended,I_extended,G_extended);
B = bwboundaries(imdilate(M_extended,strel('disk',2)));
figure(3); imshow(RGB_entended); hold on
plot(B{1}(:,2),B{1}(:,1),'-w')
scalebarLength = 40;
plot([2 2+scalebarLength],[size(RGB_entended,1)-2 size(RGB_entended,1)-2],'-w','Linewidth',1); 
text(5,size(RGB_entended,1)-5,'1 mm','Color','w')
set(gcf,'units','centimeters','position',[10 10 4 4]); 
set(gca,'units','centimeters','position',[0 0 4 4]); 

% AP averaging
clear AP_*
AP_gcamp = squeeze(nanmean(Fmod,1)); AP_ctb = squeeze(nanmean(Kmod,1)); 
bins = [0:4:size(AP_gcamp,2)]; if max(bins) < size(AP_gcamp,2); bins(end+1) = size(AP_gcamp,2); end
plotbins = bins(1:end-1)+(diff(bins)/2);
AP = [1:size(AP_gcamp,2)]; AP = discretize(AP,bins);
for j = 1:max(AP)
    AP_contra(j) = nanmean(AP_gcamp(AP == j)); 
    AP_ipsi(j) = nanmean(AP_ctb(AP == j)); 
end
AP_contra = mat2gray(AP_contra);
AP_ipsi = mat2gray(AP_ipsi);

% ML averaging
clear ML_*
ML_gcamp = squeeze(nanmean(Fmod,2)); ML_ctb = squeeze(nanmean(Kmod,2)); 
binsML = [0:4:size(ML_gcamp,1)]; if max(bins) < size(ML_gcamp,1); bins(end+1) = size(ML_gcamp,1); end
plotbinsML = binsML(1:end-1)+(diff(binsML)/2);
ML = [1:size(ML_gcamp,1)]; ML = discretize(ML,binsML);
for j = 1:max(ML)
    ML_contra(j) = nanmean(ML_gcamp(ML == j)); 
    ML_ipsi(j) = nanmean(ML_ctb(ML == j)); 
end
ML_contra = mat2gray(ML_contra);
ML_ipsi = mat2gray(ML_ipsi);

figure(1); clf; % AP average
plot(plotbins*25,AP_contra,'Color','g'); hold on
plot(plotbins*25,AP_ipsi,'Color','m'); hold on
xlim([min(plotbins) max(plotbins)]*25)
yticks([0 1])
box off; axis square
xlabel('dist. anterior edge (\mum)')
ylabel('norm. fluorescence (a.u.)')
set(gcf,'units','centimeters','position',[15 10 4 4]); 

figure(2); clf; % AP average
plot(plotbinsML*25,ML_contra,'Color','g'); hold on
plot(plotbinsML*25,ML_ipsi,'Color','m'); hold on
xlim([min(plotbinsML) max(plotbinsML)]*25)
yticks([0 1])
box off; axis square
xlabel('dist. medial edge (\mum)')
ylabel('norm. fluorescence (a.u.)')
set(gcf,'units','centimeters','position',[20 10 4 4]); 

%% supp fig 2f pupil effect across 2p runs

fr = 15.5;
for i = 1:6
    sub = data(i).Arousal.Opto.PupTraces{3}./nanstd(data(i).Arousal.normArea);
    sub = sub - nanmean(sub(1:30,:),1);
    P(i) = nanmean(nanmean(sub(round(3*fr):round(12*fr))));
end

figure;
bar(nanmean(P),'FaceColor','none','EdgeColor',ColourMap.Opto,'LineWidth',1);
[~,p] = ttest(P);
hold on
scatter(ones(size(P)),P,2,'k');
box off
xticks([]); %ylim([-.25 .25])
ylabel('\Delta pupil')
set(gcf,'units','centimeters','position',[10 10 3 4]); 
set(gca,'units','centimeters','position',[1 .8 1 3]);

% LME
Diff = P';
Mouse = nominal([1:6])';
Opto = table(Diff,Mouse);
lme = fitlme(Opto,'Diff ~ 1 + (1|Mouse)')

%% fig supp 2f2 correlation of opto effect at baseline with pupil size

clearvars -except data ColourMap

fr = 15.5; QI_thresh = .15;
for i = 1:6
    sub{1} = mean(data(i).Arousal.Opto.PupTraces{2}(1:2*fr,:)./nanstd(data(i).Arousal.normArea),1); sub{1}(data(i).Arousal.Opto.FractClosedEye{2} >= .15) = NaN;
    sub{2} = mean(data(i).Arousal.Opto.PupTraces{3}(1:2*fr,:)./nanstd(data(i).Arousal.normArea),1); sub{2}(data(i).Arousal.Opto.FractClosedEye{3} >= .15) = NaN;
    Popto = cat(2,sub{1},sub{2});
    b{1} = squeeze(mean(data(i).OptoTraces{2}(1:2*fr,:,:),1) - mean(data(i).OptoTraces{2}(4*fr:5*fr,:,:),1));
    b{2} = squeeze(mean(data(i).OptoTraces{3}(1:2*fr,:,:),1) - mean(data(i).OptoTraces{3}(4*fr:5*fr,:,:),1));
    BaseOpto = cat(1,b{1},b{2});
    QIbar = sum([max(data(i).Ret.QI(:,1:8),[],2) > QI_thresh,max(data(i).Ret.QI(:,9:16),[],2) > QI_thresh],2);
    QIff = max(data(i).Ret.QI(:,17:18),[],2) > QI_thresh;
    bads = ~((QIbar == 2 | QIff == 1)); BaseOpto(:,bads) = []; clear bads QIbar QIff
    for j = 1:size(BaseOpto,2)
        [p,r] = corrcoef(BaseOpto(~isnan(Popto),j),Popto(~isnan(Popto)));
        P(j) = p(2,1);
    end
    p_all_opto(i) = nanmean(P);
    clear P sub b BaseOpto Popto
end

% with mouse averages for main
figure(2);
bar(1,nanmean(p_all_opto),'FaceColor',ColourMap.Opto,'EdgeColor','none','FaceAlpha',.8); hold on
scatter(ones(size(p_all_opto)),p_all_opto,2,'k');
box off; xticks([1]); xticklabels({'Opto'}); xtickangle(45)
ylim([-.15 .15]); yticks([-.1 0 .1]);
ylabel('Correlation coefficient'); 
set(gcf,'units','centimeters','position',[10 10 3 4]);
set(gca,'units','centimeters','position',[1 .8 1 3]);

%% supp fig 2g example bouton individual traces

pos = 1; F = FOV(~toss);
R = resp(:,F == pos); S = respStd(:,F == pos); % A = allresp(:,:,F == pos);
B = mean(R(round(4*fr):round(11*fr),:));
[~,idx] = sort(B); chosen = [1,10,35,500,1000];

for b = 1:length(chosen)
    i = idx(chosen(b));
    figure(b); clf
    B = allresp(:,:,i); B(:,isnan(mean(B,1))) = [];
    for j = 1:size(B,2)
        B(:,j) = smooth(B(:,j),7.5);
    end
    imagesc(B'); caxis([-.5 .5]); colormap('bluewhitered')
    ylabel('trials'); xlabel('time (s)')
    xticks([31 109 187]); xticklabels({'0','5','10'})
    set(gcf,'units','centimeters','position',[10 10 3 4]); 
end

% make a colorbar
i = idx(chosen(1));
figure(6); clf
B = vis(:,:,i); B(:,isnan(mean(B,1))) = [];
for j = 1:size(B,2)
    B(:,j) = smooth(B(:,j),7.5);
end
imagesc(B'); caxis([-.5 .5]); colormap('bluewhitered'); colorbar
ylabel('trials'); xlabel('time (s)')
xticks([31 109 187]); xticklabels({'0','5','10'})
set(gcf,'units','centimeters','position',[10 10 3 4]); 

%% supp fig 2H scatter of spatial map for each mouse

clear cmap2
edges = [-.1:.025:.1]; ref = edges(1:end-1)+.0125;
cmap = ones(length(edges)-1,3);
b = [0:1/(size(cmap,1)/2):1];
cmap(1:size(cmap,1)/2,3) = fliplr(b(2:end));
cmap((size(cmap,1)/2)+1:end,1) = b(2:end);
cmap = [0 0 .5; 0 0 1; 1 1 1; 1 0 0; .5 0 0]; % bluewhitered standard!
for i = 1:3
    cmap2(:,i) = interp1([1:5],cmap(:,i),[1:5/length(edges):5]);
end

F = FOV'; F(toss) = [];
A = Base;
A(A >= .1) = .1; A(A <= -.1) = -.1; 
A = discretize(A,edges);
for m = 1:6
    figure(m);
    for j = 1:max(A)
        scatter(XY(A == j & F == m,1),XY(A == j & F == m,2),.2,cmap2(j,:),'filled'); hold on
    end 
    figure(m); set(gca,'YDir','reverse'); daspect([1 1 1]); set(gca,'visible','off'); ylim([-50 400]); xlim([0 600])
    set(gcf,'units','centimeters','position',[10 10 4.5 3]); 
end
figure(7); % this is a colormap!!
imagesc(ref'); colormap(bluewhitered); ylabel('baseline \DeltaF/F')
set(gca,'xtick',[])
yticks([1 length(ref)])
yticklabels({'-.1', '.1'})
set(gcf,'units','centimeters','position',[20 10 .5 2]); 

%% supp 2I anterior/posterior axis difference in suppression

F = FOV; F(toss) = []; cc = hsv(6);
clear B_AP
for i = 1:max(F)
    AP = XY(F == i,1)'; B = Base(F == i);
    AP = discretize(AP,[0:60:550]);
    for j = 1:max(AP)
        B_AP(j,i) = mean(B(AP == j));
    end
end

figure(1); clf
ypos = -.03;
for i = 1:6
    plot([30:60:520]*(5/3),B_AP(:,i),'Color',cc(i,:)); hold on
    plot([700 750],[ypos ypos],'-','Color',cc(i,:));
    text(770,ypos,strcat('M',num2str(i)));
    ypos = ypos-.007;
end
plot([30 520]*(5/3),[0 0],'--k');
xlim([30 520]*(5/3)); ylim([-.08 .01])
yticks([-.05 0])
box off
xlabel('\mum')
ylabel('\DeltaF/F')
set(gcf,'units','centimeters','position',[10 10 4 4]); 

F = FOV; F(toss) = []; cc = hsv(6);
clear B_AP
for i = 1:max(F)
    ML = XY(F == i,2)'; B = Base(F == i);
    ML = discretize(ML,[0:60:400]);
    for j = 1:max(ML)
        B_ML(j,i) = mean(B(ML == j));
    end
end

figure(2); clf
ypos = -.03;
for i = 1:6
    plot([30:60:370]*(5/3),B_ML(:,i),'Color',cc(i,:)); hold on
end
plot([30 370]*(5/3),[0 0],'--k');
xlim([30 370]*(5/3)); ylim([-.08 .01])
yticks([-.05 0])
box off
xlabel('\mum')
ylabel('\DeltaF/F')
set(gcf,'units','centimeters','position',[10 10 4 4]); 

