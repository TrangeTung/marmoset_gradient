clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
addpath(genpath(codepath));
ExpRecording = [codepath 'ExpIndex.xlsx'];
[~,~,CellData] = xlsread(ExpRecording);
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
WholeMaskNII = fullfile(codepath,'marmoset_atlas','Trange_Template_mask_V3.nii');

AtlasNII = fullfile(codepath,'marmoset_atlas','Trange_atlas_RikenBMA_cortex.nii');
AtlasMask = spm_read_vols_4D(spm_vol(AtlasNII));

AtlasExcel = fullfile(codepath,'marmoset_atlas','Trange_modifid_RikenBMA_atlas_labels.xlsx');
[~,~,CellData] = xlsread(AtlasExcel);
AtlasTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
atlas_label = AtlasTable.Labels;
Abbre = AtlasTable.Abbre;

xhdr=spm_vol(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.nii.gz'));
Labels = spm_read_vols(xhdr);
AtlasTXT = loadtxt(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.txt'));

%% ION Group
%
a_=0;

a=65;b=86;c=64;d=512; % 4D-EPI size
lmask = AtlasMask>20;
Xmask = imresize3(lmask,round([a,b,c]/2),'nearest');
NUM=5522;

filepath = fullfile(WholePath,'DiffusionGradient_Structure','.');
load(fullfile(filepath,'G.mat'));
Gx = G.gradients{1};
G_Ref = nan(116,5); G_Ref([1:50,52:116],:)=Gx ;
G_Ref = fillmissing(G_Ref,'linear');

Vf_ION = zeros(116,345,10);

ArousalExcel = fullfile('E:\rsfMRI_marmoset\Figure1','ArousalIndex.xlsx');
[~,~,ArousalTable] = xlsread(ArousalExcel);
Vec = {};
for idx = 1:318
    idx
    path = ArousalTable{idx,2};
    
    kk = ArousalTable{idx,3};
    
    GEPath = fullfile(path,'Functions','Seed-based',num2str(kk));
    V=xlsread(fullfile(GEPath,'TimeSeries_Riken_new.xlsx')); 
    V(isnan(V))=0;
    V = (V-mean(V,2))./std(V,0,2);
    %V(isnan(V))=0;
    
    Q = V(21:end,:);
    Qm = mean(Q,2);
    zeroidx = find(isnan(Qm)| std(Q,0,2)<0.1);
    if ~isempty(zeroidx)
        Q(zeroidx,:)=nan;
        Q = fillmissing(Q,'nearest');
    end
    
    if idx==118
        1;
    end
    
    Vec{idx} = corr(Q');;%
end

filepath = fullfile(WholePath,'DiffusionGradient_NS');
load(fullfile(filepath,'G.mat'));
NS_Gd = G.aligned{1};

G = GradientMaps('kernel','cs','approach','dm','n_components',5,'align','pa');
G = G.fit(Vec,'Reference',G_Ref);
CC_SC = zeros(a,5);
CC_NS = zeros(a,5);
for al = 1:318
    Gdx = G.aligned{al};
    CX = corr(G_Ref,Gdx);
    CC_SC(al,:) = diag(CX);
    CX = corr(Gdx,NS_Gd);
    CC_NS(al,:) = diag(CX);
end

CC = CC_NS;
EYE = cell2mat(ArousalTable(:,1));
%EYE = atanh(EYE*2-1);
bins=10; clear M
[no,xo]  = hist(EYE,bins);
Group =  zeros(size(EYE,1),1);
for i=1:bins
    loc = find((EYE>xo(i)-1/bins/2)& (EYE<xo(i)+1/bins/2));
    Group(loc) = i;
    M(i,:) = nanmedian(CC(loc,:),1);   
    SEM(i,:) = nanstd(CC(loc,:))/sqrt(size(CC(loc,:),1));
end
figure;plot(M)


F = figure('Position', [680 58 1007 920]);
EYE = cell2mat(ArousalTable(:,1));
clear RR
for loop=1:4

    subplot(2,2,loop);
    plot(EYE*100,CC(1:318,loop),'o','MarkerSize',10,...
        'MarkerEdgeColor','none','MarkerFaceColor','k');ylim([0 1])
    xlabel('EyeOpenRatio (%)');ylabel('Similarity (C.C.)');
    NANidx = ~isnan(EYE) & ~isnan(CC(:,loop));
    [p,s] = polyfit(EYE(NANidx)*100,CC(NANidx,loop),2);
    y_ = polyval(p,EYE(NANidx));
    [b,bint,r,rint,stats] = regress(y_,[ones(size(CC(NANidx,loop),1),1),CC(NANidx,loop)]);
    RR(loop,:) = [stats(1),stats(3),s.normr];
stats(3)
    x = 0.01:0.01:1;
    y = polyval(p,x*100);
    [Y,DELTA]=polyconf(p,x,s,'predopt','curve');
    hold on;
    patch([x,fliplr(x)]*100,[y+DELTA,fliplr(y-DELTA)],'red',...
         'EdgeColor','none','FaceAlpha',.5);
     plot(x*100,y,'linewidth',2);
    set(gca,'fontsize',15,'linewidth',1.5,'box','off',...
        'tickdir','out','ticklength',[0.05 0.025]);
end


F = figure('Position', [680 58 1007 920]);
EYE = cell2mat(ArousalTable(:,1));
EYE = atanh(EYE*2-1);
clear RR
for loop=1:4

    subplot(2,2,loop);
    plot(EYE,CC(1:318,loop),'o','MarkerSize',10,...
        'MarkerEdgeColor','none','MarkerFaceColor','k');ylim([0 1])
    xlabel('EyeOpenRatio (z)');ylabel('Similarity (C.C.)');
    NANidx = ~isnan(EYE) & ~isnan(CC(:,loop)) & ~(EYE==Inf) & ~(EYE==-Inf);
    [p,s] = polyfit(EYE(NANidx),CC(NANidx,loop),2);
    y_ = polyval(p,EYE(NANidx));
    [b,bint,r,rint,stats] = regress(y_,[ones(size(CC(NANidx,loop),1),1),CC(NANidx,loop)]);
    RR(loop,:) = [stats(1),stats(3)];
    x = -5:0.01:5;
    y = polyval(p,x);
    [Y,DELTA]=polyconf(p,x,s,'predopt','curve');
    hold on;
    patch([x,fliplr(x)],[y+DELTA,fliplr(y-DELTA)],'red',...
         'EdgeColor','none','FaceAlpha',.5);
     plot(x,y,'linewidth',2);
    set(gca,'fontsize',15,'linewidth',1.5,'box','off',...
        'tickdir','out','ticklength',[0.05 0.025]);
end


F = figure('Position', [680 787 1027 191]);
for loop=1:4
subplot(1,4,loop);
raincloud_plot(CC(:,loop),'band_width',0.05)
    set(gca,'fontsize',15,'linewidth',1.5,'box','off',...
        'tickdir','out','ticklength',[0.05 0.025]);
    xlim([-0.1 1])
end

    
F = figure('Position', [680 787 1027 191]);
subplot(1,4,1);
EYE = cell2mat(ArousalTable(:,1));
raincloud_plot(EYE(:),'band_width',0.01)
    set(gca,'fontsize',15,'linewidth',1.5,'box','off',...
        'tickdir','out','ticklength',[0.05 0.025]);
subplot(1,4,2);
EYE = cell2mat(ArousalTable(:,1));
EYE = atanh(EYE*2-1);
raincloud_plot(EYE(NANidx),'band_width',0.1)
    set(gca,'fontsize',15,'linewidth',1.5,'box','off',...
        'tickdir','out','ticklength',[0.05 0.025]);
subplot(1,4,3);
raincloud_plot(CC(:),'band_width',0.01)
    set(gca,'fontsize',15,'linewidth',1.5,'box','off',...
        'tickdir','out','ticklength',[0.05 0.025]);
    

EYE = cell2mat(ArousalTable(:,1));
% EYE = atanh(EYE*2-1);
clear RR
for loop=1:4

%     NANidx = ~isnan(EYE) & ~isnan(CC(:,loop)) & ~(EYE==Inf) & ~(EYE==-Inf);
    
    [p,s] = polyfit(EYE(NANidx)*100,CC(NANidx,loop),1);
    y_ = polyval(p,EYE(NANidx));
    [b,bint,r,rint,stats] = regress(y_,[ones(size(CC(NANidx,loop),1),1),CC(NANidx,loop)]);
    RR(loop,1) = stats(1);%s.normr;
%     RR(loop,1) = sqrt(sum((y_ - CC(NANidx,loop) ).^2));
    
    [p,s] = polyfit(EYE(NANidx)*100,CC(NANidx,loop),2);
    y_ = polyval(p,EYE(NANidx));
    [b,bint,r,rint,stats] = regress(y_,[ones(size(CC(NANidx,loop),1),1),CC(NANidx,loop)]);
    RR(loop,2) = stats(1);%s.normr;
%     RR(loop,3) = sqrt(sum((y_(NANidx) - CC(NANidx,loop) ).^2));
    
    [p,s] = polyfit(EYE(NANidx)*100,CC(NANidx,loop),3);
    y_ = polyval(p,EYE(NANidx));
    [b,bint,r,rint,stats] = regress(y_,[ones(size(CC(NANidx,loop),1),1),CC(NANidx,loop)]);
    RR(loop,3) = stats(1);%s.normr;
%     RR(loop,3) = sqrt(sum((y_(NANidx) - CC(NANidx,loop) ).^2));
    
end
