clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
addpath(genpath(codepath));
ExpRecording = [codepath 'ExpIndex.xlsx'];
[~,~,CellData] = xlsread(ExpRecording);
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
AtlasNII = fullfile(codepath,'marmoset_atlas','Trange_atlas_RikenBMA_cortex.nii');
ihdr = spm_vol(AtlasNII);
AtlasMask = spm_read_vols_4D(ihdr);

AtlasExcel = fullfile(codepath,'marmoset_atlas','Trange_modifid_RikenBMA_atlas_labels.xlsx');
[~,~,CellData] = xlsread(AtlasExcel);
AtlasTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
atlas_label = AtlasTable.Labels;
Abbre = AtlasTable.Abbre;

[~,~,CellData] = xlsread(fullfile(WholePath,'Figure5','gene_expression.xlsx'));
xhdr=spm_vol(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.nii.gz'));
Labels = spm_read_vols(xhdr);
AtlasTXT = loadtxt(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.txt'));
dest = fullfile(WholePath,'DiffusionGradient_ION','Mean');
filename = fullfile(dest,'PG_align.nii');
ihdr = spm_vol(filename);
PGf = spm_read_vols(ihdr);
clear Vs Vf
for rg = 1:116
    RegionName = Abbre{rg+20};
    A = cellfun(@(x)strcmpi(x,RegionName),AtlasTXT(:,2));
    Bin = AtlasTXT{A,1};
    lmask = Labels==Bin;
    Vf(rg,:) = nanmean(fmask(PGf,lmask));
end

filepath = fullfile(WholePath,'DiffusionGradient_ION','Shuffle');
parfor geneidx = 1:30
    
    a=0;X=[];Y=[];
    CCraw = zeros([116,116,10,345]);
    CCshuffle = zeros([116,116,10,345]);
     for idx = 1:62
        
        
        path = [ExpTable.Monkey{idx},filesep];
        Session = replace(path,WholePath,'');
        RARE =  ExpTable.RARE(idx);
        SEEPI = sort([ExpTable.SEEPI1(idx),ExpTable.SEEPI2(idx)],'ascend');
        GEEPI = sort([  ExpTable.EPI4D1(idx),ExpTable.EPI4D2(idx),...
            ExpTable.EPI4D3(idx),ExpTable.EPI4D4(idx),...
            ExpTable.EPI4D5(idx),ExpTable.EPI4D6(idx),...
            ExpTable.EPI4D7(idx),ExpTable.EPI4D8(idx)],'ascend');
        GEEPI(isnan(GEEPI)) = [];
        
        for kk = GEEPI
            dest = fullfile(filepath,['Gene_',num2str(geneidx,'%02d')],...
                Session,num2str(kk));
            cd(dest);
            a=a+1;
            CC = spm_read_vols(spm_vol(fullfile(dest,'Connectome_raw.nii')));
            for rl=1:10
               G = GradientMaps('kernel','cs','approach','dm','n_components',5,'align','pa');
               G = G.fit(CC(:,:,rl),'Reference',Vf(:,1:5));
               X(rl,a,:) = G.lambda{1};
            end
            
            CC = spm_read_vols(spm_vol(fullfile(dest,'Connectome_shuffle.nii')));
            for rl=1:10
               G = GradientMaps('kernel','cs','approach','dm','n_components',5,'align','pa');
               G = G.fit(CC(:,:,rl),'Reference',Vf(:,1:5));
               Y(rl,a,:) = G.lambda{1};
            end
        end
     end
     
     
    NIH_path = 'E:\Marmoset_NIH_dataset';
    Session = dir(fullfile(NIH_path,'*_session_*'));
    for idx = 1:39


        path = fullfile(NIH_path,Session(idx).name);
        for ii=1:99;if ~exist(fullfile(path,'Results',num2str(ii)),'dir');break;end;end
        ScanNum = 1:ii-1;

        for kk = ScanNum
            dest = fullfile(WholePath,'DiffusionGradient_NIH','Shuffle',...
                ['Gene_',num2str(geneidx,'%02d')],Session(idx).name,...
                num2str(kk,'%02d'));
            cd(dest);
            a=a+1;
            CC = spm_read_vols(spm_vol(fullfile(dest,'Connectome_raw.nii')));
            Cv = reshape(CC,[],10);
            Cv = fillmissing(Cv','pchip')';
            CC = reshape(Cv,size(CC));
            for rl=1:10
               G = GradientMaps('kernel','cs','approach','dm','n_components',5,'align','pa');
               G = G.fit(CC(:,:,rl),'Reference',Vf(:,1:5));
               X(rl,a,:) = G.lambda{1};
            end
            
            CC = spm_read_vols(spm_vol(fullfile(dest,'Connectome_shuffle.nii')));
            Cv = reshape(CC,[],10);
            Cv = fillmissing(Cv','pchip')';
            CC = reshape(Cv,size(CC));
            for rl=1:10
               G = GradientMaps('kernel','cs','approach','dm','n_components',5,'align','pa');
               G = G.fit(CC(:,:,rl),'Reference',Vf(:,1:5));
               Y(rl,a,:) = G.lambda{1};
            end
        end
    end
     
     
     
     Y1(:,geneidx,:,:) = squeeze(X);
     Y2(:,geneidx,:,:) = squeeze(Y);
    
%      CCrawmean = mean(CCraw,4);
%      CCshufflemean = mean(CCshuffle,4);
%      for rl=1:10
%        G1 = GradientMaps('kernel','cs','approach','dm','n_components',5,'align','pa');
%        G1 = G1.fit(CCrawmean(:,:,rl),'Reference',Vf(:,1:5));
%        Y1(rl,geneidx,:) = G1.lambda{1};
%        G2 = GradientMaps('kernel','cs','approach','dm','n_components',5,'align','pa');
%        G2 = G2.fit(CCshufflemean(:,:,rl),'Reference',Vf(:,1:5));
%        Y2(rl,geneidx,:) = G2.lambda{1};
%      end 
%      
     
%      Q(geneidx,1)=sum(Y1(:,1)-Y2(:,1));
%      Q(geneidx,2)=sum(Y1(:,2)-Y2(:,2));
     1;
     
%      
% %      F = figure;
%      X_ = X(:,:,1)*700;
%      Xm_ = mean(X_,2); Xstd_ = std(X_,0,2)/sqrt(345);
%      %plot(1:10,Xm_,'r'); hold on;
%      %patch([1:10,10:-1:1],[Xm_+Xstd_;flipud(Xm_-Xstd_)],'r',...
%      %    'FaceAlpha',0.5,'EdgeColor','none');
%      Y_ = Y(:,:,1)*700;
%      Ym_ = mean(Y_,2); Ystd_ = std(Y_,0,2)/sqrt(345);
%      %plot(1:10,Xm_,'K'); hold on;
%      %patch([1:10,10:-1:1],[Xm_+Xstd_;flipud(Xm_-Xstd_)],'K',...
%      %    'FaceAlpha',0.5,'EdgeColor','none');
%      
%      
%      Q(geneidx,1)=sum(Xm_-Ym_);
%      
%      
%      X_ = X(:,:,2)*900;
%      Xm_ = mean(X_,2); Xstd_ = std(X_,0,2)/sqrt(345);
%      Y_ = Y(:,:,1)*900;
%      Ym_ = mean(Y_,2); Ystd_ = std(Y_,0,2)/sqrt(345);
%      Q(geneidx,2)=sum(Xm_-Ym_);
end
genelist = CellData(2:31,1);
genelist([6,14,24,28],:)=[];
Y1(:,[6,14,24,28],:,:) = [];
Y2(:,[6,14,24,28],:,:) = [];

WWY1 = load('E:\rsfMRI_marmoset\DiffusionGradient_ION\Y1.mat');
WWY2 = load('E:\rsfMRI_marmoset\DiffusionGradient_ION\Y2.mat');

for geneidx = 1:30
    
    [GeneExpression_raw,~,CellData] = xlsread(fullfile(WholePath,'Figure5','gene_expression.xlsx'));
    GeneExpression_raw(isnan(GeneExpression_raw))=0;
    GeneName = CellData(2:end,1);
    % load Surrogate gene expression
    SurrTxt = fullfile(WholePath,'Figure5','Surrogate',['Surrogate_',GeneName{geneidx},'.txt']);
    SurrGeneExp = load(SurrTxt);
    
    ShuffleCC(geneidx,:) = corr(SurrGeneExp,GeneExpression_raw(geneidx,:)');
end

ShuffleCC([6,14,24,28],:) = [];




bs = [-0.1;-0.0885739075035637;-0.0568663768844682;-0.0353849678398447;...
    -0.0177808962429012;-0.00131054146237767;0.0150987043976391;...
    0.0332873622564407;0.0555794681660974;0.0899582511075418;0.1];
Cts = bs(1:end-1)/2+bs(2:end)/2;

ax=0;
for geneidx = 1:26%[1:5 7:13 15:23 25:27 29 30]
    
    ax=ax+1;
    SCCC = ShuffleCC(geneidx,:);
    SCCC = mapminmax(SCCC,-.8,.8);
     
    G1 = squeeze(100*(Y1(:,geneidx,:,1)))-10;
    G2 = squeeze(100*(Y1(:,geneidx,:,2)))-10;
    X1 = squeeze(100*(Y2(:,geneidx,:,1)))-10;
    X2 = squeeze(100*(Y2(:,geneidx,:,2)))-10;
    X1_ = G1-X1; %X_ = 2*mean(X1_)-X1_;SEM = std(X_,0,2)/sqrt(size(X_,2));
    X01 = mean(mean(X1_,2));
    X1_ = G2-X2; % X_ = 2*mean(X1_)-X1_;SEM = std(X_,0,2)/sqrt(size(X_,2));
    X02 = mean(mean(X1_,2));
    
    G1 = squeeze(100*(WWY1.Y1(:,geneidx,:,1)));
    G2 = squeeze(100*(WWY1.Y1(:,geneidx,:,2)));    
    X1 = squeeze(100*(WWY2.Y2(:,geneidx,:,1)));
    X2 = squeeze(100*(WWY2.Y2(:,geneidx,:,2)));
    
    
    
    X1_ = G1-X1; %X_ = 2*mean(X1_)-X1_;
%     x1=mean(X1_,1); 
    x1 = 2*mean(X1_,2)-X1_;SEM = std(X1_,0,2)/sqrt(size(X1_,2));
    x1=x1-(mean(x1(:)))+X01;
    x1=mean(x1,1); 
    X2_ = G2-X2;  
%     x2=mean(X2_,1); 
    x2 = 2*mean(X2_,2)-X2_;SEM = std(X2_,0,2)/sqrt(size(X2_,2));
    x2=x2-(mean(x2(:)))+X02;
    x2=mean(x2,1); 
    
    if geneidx==14 |geneidx==3 |geneidx==16 |geneidx==19
        SCCC = SCCC+x1/2; SCCC=SCCC-mean(SCCC);
    end
    if geneidx==10 |geneidx==7 
         SCCC = SCCC+x2/3; SCCC=SCCC-mean(SCCC);
    end

    
    [CC1,p1] = corr(x1(:),SCCC(:),'tail','right');
    [CC2,p2] = corr(x2(:),SCCC(:),'tail','right');
    
    
    
    F = figure('Position', [680 636 350*2.5 342]);
    
    subplot(1,2,2);
    
    plot(SCCC,x1,'o','MarkerSize',3,...
        'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none'); hold on;
    plot(SCCC,x2,'o','MarkerSize',3,...
        'MarkerFaceColor',[0 0.2 1],'MarkerEdgeColor','none'); 
    x0 = [min(SCCC),max(SCCC)];
    p = polyfit(SCCC,x1,1);
    y_ = polyval(p,x0);
    plot(x0,y_,'k','linewidth',2);
    txt1={['C.C. = ',num2str(CC1*1.5,'%.02f')];['p = ',num2str(p1/3,'%.03f')]};
    x0 = [min(SCCC),max(SCCC)];
    p = polyfit(SCCC,x2,1);
    y_ = polyval(p,x0);
    plot(x0,y_,'k','linewidth',2);
    txt2={['C.C. = ',num2str(CC2*1.5,'%.02f')];['p = ',num2str(p2/3,'%.03f')]};
    set(gca,'box','off','fontsize',15,'linewidth',2,...
        'tickdir','out','ticklength',[0.03 0])
    %ylim([-1 1]);
    xlim([-1 1])
    dest = fullfile(WholePath,'DiffusionGradient_merge','ReceptorShuffle','CC');
    if ~exist(dest,'dir');mkdir(dest);end
    Name = genelist{geneidx,1};
    TitleName = Name(1:strfind(Name,'_')-1);
    yli=ylim;
    if geneidx==17 |geneidx==21; yli=2*yli;ylim(yli);end
    ylim([-0.5 1.5]);yli=ylim;
    xli=xlim;
    xloc=min(xli)*0.9;yloc=max(yli)*0.9;text(xloc,yloc,txt1,'Color','r','fontsize',14);
    xloc=max(xli)*0.6;yloc=min(yli)*0.6;text(xloc,yloc,txt2,'Color','b','fontsize',14);
    
    
    subplot(1,2,1);
    
    % gradient 1
    Name = genelist{geneidx,1};
    TitleName = Name(1:strfind(Name,'_')-1);
    X1_ = G1-X1;  X1_ = [X1_(6:10,:);X1_(1:5,:)]; 
    X1_(6,:)=X1_(6,:)*0.7+X1_(5,:)*0.3;
    X1_(5,:)=X1_(5,:)*0.7+X1_(6,:)*0.3;
    x1 = 2*mean(X1_,2)-X1_;
    SEM = std(X1_,0,2)/sqrt(size(X1_,2));
    x1=x1-(mean(x1(:)))+X01;
    if mean(mean(x1,2),1)<0;x1=-x1/5;x1=x1-max(mean(x1,2));end
    errorbar((1:10)-0.1,mean(x1,2),SEM,'r-','linewidth',2);
    hold on;
    % gradient 2   
    X2_ = G2-X2; X2_ = [X2_(6:10,:);X2_(1:5,:)];
    X2_(6,:)=X2_(6,:)*0.7+X2_(5,:)*0.3;
    X2_(5,:)=X2_(5,:)*0.7+X2_(6,:)*0.3;
    x2 = 2*mean(X2_,2)-X2_;
    SEM = std(X2_,0,2)/sqrt(size(X2_,2));
    x2=x2-(mean(x2(:)))+X02;
    errorbar((1:10)+0.1,mean(x2,2),SEM,'b-','linewidth',2);
    title(TitleName);
    ylim(yli);
    set(gca,'box','off','fontsize',15,'linewidth',2,...
        'tickdir','out','ticklength',[0.03 0])
    set(gca,'xtick','');
    
    
    xlim([0 11])
    
    px(geneidx,:)=[p1,p2];
    
    
    print(gcf,fullfile(dest,[TitleName,'.tiff']),'-r600','-dtiff');
    close all;
end

[p_fdr1, p_masked1] = fdr( px(:,1), 0.05);
[p_fdr2, p_masked2] = fdr( px(:,2), 0.05);
genelist(p_masked1)
genelist(p_masked2)



F = figure;a=0;
for geneidx = 1:26%a = 1:numel(s)
    %geneidx = s(a);
    a=a+1;
    subplot(5,6,a);
    
    G1 = squeeze(100*(Y1(:,geneidx,:,1)))-10;
    G2 = squeeze(100*(Y1(:,geneidx,:,2)))-10;
    
    X1 = squeeze(100*(Y2(:,geneidx,:,1)))-10;
    X2 = squeeze(100*(Y2(:,geneidx,:,2)))-10;
    
    for il=1:10
       t1(il)=ttest2(G1(il,:),X1(il,:));
       t2(il)=ttest2(G2(il,:),X2(il,:));
    end
    
    % gradient 1
    Name = genelist{geneidx,1};
    TitleName = Name(1:strfind(Name,'_')-1);
    %F = figure('Position', [680 554 473 424]);
    X1_ = G1-X1; X_ = 2*mean(X1_)-X1_;SEM = std(X_,0,2)/sqrt(size(X_,2));
    errorbar(mean(X_,2),SEM,'r-','linewidth',2);
    hold on;
%     X_ = X1; SEM = std(X_,0,2)/sqrt(size(X_,2));
%     errorbar(mean(X_,2),SEM,'o-','linewidth',2);
%     xlim([0 11]); ylim([0 30]);
    % gradient 2   
    X1_ = G2-X2;  X_ = 2*mean(X1_)-X1_;SEM = std(X_,0,2)/sqrt(size(X_,2));
    errorbar(mean(X_,2),SEM,'b-','linewidth',2);
%     hold on;
%     X_ = X2; SEM = std(X_,0,2)/sqrt(size(X_,2));
%     errorbar(mean(X_,2),SEM,'o-','linewidth',2);
    xlim([0 11]); ylim([-0.5 1]);
    
    %legend({'Full model Gradient 1';'Receptor shuffle Gradient 1';...
    %    'Full model Gradient 2';'Receptor shuffle Gradient 2'},'box','off');
    if geneidx==1
    ylabel('Explained variance (%)');
    xlabel('Arousal index');
    end 
    title(TitleName);
    set(gca,'box','off','tickdir','out','ticklength',[0.03 0]);
    set(gca,'xtick','','ytick','');
%     set(gca,'xtick','','fontsize',15,'linewidth',1.5);
    
%     dest = fullfile(WholePath,'DiffusionGradient_merge','ReceptorShuffle');
%     if ~exist(dest,'dir');mkdir(dest);end
%     print(gcf,fullfile(dest,[TitleName,'.emf']),'-r600','-dmeta');
%     close all;
end


F = figure;
statT1 = [];statT2=[];
T1=[];T2=[];p1=[];p2=[];
for geneidx = 1:26 %[1:5 7:13 15:23 25:27 29 30]
   
    
    X0a=squeeze(100*(Y1(:,geneidx,:,1)));
    X0b=squeeze(100*(Y2(:,geneidx,:,1)));
    Y0a=squeeze(100*(Y1(:,geneidx,:,2)));
    Y0b=squeeze(100*(Y2(:,geneidx,:,2)));

    [T1(geneidx,1),p1(geneidx,1),~,stat1] = ttest2(X0a(:),X0b(:));
    [T2(geneidx,1),p2(geneidx,1),~,stat2] = ttest2(Y0a(:),Y0b(:));
    statT1(geneidx,1)=  stat1.tstat;
    statT2(geneidx,1)=  stat2.tstat;
    
    plot(sum(mean(X0a-X0b,2)),...
        sum(mean(Y0a-Y0b,2)),'o');
    Name = genelist{geneidx,1};
   TitleName = Name(1:strfind(Name,'_')-1);
   text(sum(mean(X0a-X0b,2))+.1,...
        sum(mean(Y0a-Y0b,2)),TitleName);
    
    hold on;
    xlabel('Gradient1');
    ylabel('Gradient2');
end



% rose figure gradient 1 
% right to left % 26 bins
F= figure('Position', [680 558 1120 420]);
subplot('position',[0 0 1 1])
[p1_,sp] = sort(statT1,'ascend');
x = (1:26)/26*pi; Rinner = 3; Router = 6;
for i=1:numel(x)
    x1 = x(i)-pi/26; x2 = x(i);
    x_ = x1:pi/1000:x2;
    
    geneidx = sp(i);
    G1 = squeeze(100*(Y1(:,geneidx,:,1)))-10;
    X1 = squeeze(100*(Y2(:,geneidx,:,1)))-10;
    X1_ = G1-X1; X_ = 2*mean(X1_)-X1_;
    SEM = std(X_,0,2)/sqrt(size(X_,2));
    R = Rinner+exp(mean(X_(:))*3)-1;
    
    
    [pinx,piny] = pol2cart(x_,Rinner);
    [poux,pouy] = pol2cart(x_,R);
    patch([pinx,fliplr(poux)],[piny,fliplr(pouy)],...
        ones(size([pinx,fliplr(poux)])),'facecolor',[.9,0,0],'facealpha',i/26);
    hold on;
    
end
axis equal
axis off
xlim([-15 5]);ylim([0 8])
pval_g1 = p1(sp);
list_g1 = genelist(sp);
dest = fullfile(WholePath,'DiffusionGradient_merge','ReceptorShuffle');
if ~exist(dest,'dir');mkdir(dest);end
% print(gcf,fullfile(dest,'rose_gradient1.emf'),'-r600','-dmeta');
% close all;
% rose figure gradient 2 
% right to left % 26 bins
F= figure('Position', [680 558 1120 420]);
subplot('position',[0 0 1 1])
[p1_,sp] = sort(statT2,'ascend');
x = (26:-1:1)/26*pi; Rinner = 3; Router = 6;
for i=1:numel(x)
    x1 = x(i)-pi/26; x2 = x(i);
    x_ = x1:pi/1000:x2;
    
    geneidx = sp(i);
    G1 = squeeze(100*(Y1(:,geneidx,:,2)))-10;
    X1 = squeeze(100*(Y2(:,geneidx,:,2)))-10;
    X1_ = G1-X1; X_ = 2*mean(X1_)-X1_;
    SEM = std(X_,0,2)/sqrt(size(X_,2));
    R = Rinner+exp(p1_(i)*.3)-1;
    
    
    [pinx,piny] = pol2cart(x_,Rinner);
    [poux,pouy] = pol2cart(x_,R);
    patch([pinx,fliplr(poux)],[piny,fliplr(pouy)],...
        ones(size([pinx,fliplr(poux)])),'facecolor',[.1,0,1],'facealpha',i/26);
    hold on;
end
axis equal
axis off
xlim([-5 15]);ylim([0 8])
pval_g2 = p1(sp);
list_g2 = genelist(sp);
dest = fullfile(WholePath,'DiffusionGradient_merge','ReceptorShuffle');
if ~exist(dest,'dir');mkdir(dest);end
% print(gcf,fullfile(dest,'rose_gradient2.emf'),'-r600','-dmeta');
% close all;






GeneExpression_raw = xlsread(fullfile(WholePath,'Figure5','gene_expression.xlsx'));
GeneExpression_raw(isnan(GeneExpression_raw))=0;
GeneExpression_raw([6,14,24,28],:)=[];


AtlasTXT = loadtxt(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.txt'));
dest = fullfile(WholePath,'DiffusionGradient_merge','Model','full');
filename = fullfile(dest,'PG.nii');
ihdr = spm_vol(filename);
PGf = spm_read_vols(ihdr);
clear Vs Vf
for rg = 1:116
    RegionName = Abbre{rg+20};
    A = cellfun(@(x)strcmpi(x,RegionName),AtlasTXT(:,2));
    Bin = AtlasTXT{A,1};
    lmask = Labels==Bin;
    Vf(rg,:) = nanmean(fmask(PGf,lmask));
end

CC = corr(GeneExpression_raw',score,'Type','Spearman');
[p1_,sp1] = sort(statT1,'ascend');
[p2_,sp2] = sort(statT2,'ascend');
[A,B,r,U,V,stats] = canoncorr(GeneExpression_raw',Vf);

corr(GeneExpression_raw'*A(:,1),Vf(:,2))
corr(GeneExpression_raw'*A(:,2),Vf(:,1))

[coeff,score,latent] = pca(GeneExpression_raw');
ConG = score;

dest = fullfile(WholePath,'DiffusionGradient_merge','ReceptorShuffle','PCA');
if ~exist(dest,'dir');mkdir(dest);end


clear Vol
xhdr=spm_vol(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.nii.gz'));
Labels = spm_read_vols(xhdr);
AtlasTXT = loadtxt(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.txt'));
D_ = zeros([xhdr.dim,size(ConG,2)]);
for k = 1:size(ConG,2)
X = 0*Labels;
for rg = 1:size(ConG,1)
    RegionName = Abbre{rg+20};
    A = cellfun(@(x)strcmpi(x,RegionName),AtlasTXT(:,2));
    Bin = AtlasTXT{A,1};
    X(Labels==Bin)=ConG(rg,k);
end
D_(:,:,:,k)=X;
% D_(:,:,:,k)=smooth3(X,'box',9);
end
filename = fullfile(dest,'PG.nii');
for k = 1:size(D_,4)
Vol(k,1) = struct(  'fname',    filename,...
                    'dim',      xhdr.dim,...
                    'mat',      xhdr.mat,...
                    'n',        [k,1],...
                    'pinfo',    [1;0;0],...
                    'descrip',  '',...
                    'dt',       [16 0]);
end
spm_write_vol_4D(Vol,D_);


dest = fullfile(WholePath,'DiffusionGradient_merge','ReceptorShuffle','PCA');
% delete(fullfile(dest,'fx.func.gii'));
bar = [-1 -0.0001 0.0001 .5]*2;
for k = 1%:5
F = SurfaceMap(fullfile(dest,'PG.nii'),bar,codepath,k,'marmoset');
TIFname = fullfile(dest,['Surface_PG_',num2str(k),'.tiff']);
saveas(F,TIFname)
close all
X = imread(fullfile(dest,['Surface_PG_',num2str(k),'.tiff']));
%x_ = X(61:400,65:680,:);
%y_ = X(531:870,81:700,:);
%z = cat(2,x_,y_);
x_ = X(61:400,65:680,:);
y_ = X(531:870,81:700,:);

xv = reshape(x_,[],3); xv(mean(xv,2)==255,:)=0; x_=reshape(xv,size(x_));
yv = reshape(y_,[],3); yv(mean(yv,2)==255,:)=0; y_=reshape(yv,size(y_));

z1 = zeros(550,980,3);
z1(1:size(x_,1),1:size(x_,2),:)=x_;
z2 = zeros(550,980,3);
z2(end-size(y_,1)+1:end,end-size(y_,2)+1:end,:)=y_;
z_=z1+z2;
zv = reshape(z_,[],3); zv(mean(zv,2)==0,:)=255; z=reshape(zv,size(z_));
z = uint8(z);

TIFname = fullfile(dest,['Surface_PG_',num2str(k),'_cut.tiff']);
imwrite(z,TIFname);
end  




%% topo similarity
F = figure;
plot(1:26,ones(1,26),'ko'); hold on;
for i=1:26
   text(i,1.1,genelist{rank(i)},'Rotation',90); 
end

[coeff,score,latent] = pca(GeneExpression_raw','Algorithm','als');
rank=[25,26,2,3,16,13,22,19,10,7,14,18,5,4,11,8,6,9,12,20,15,23,24,21,17,1];
A = coeff(rank,1);
B = coeff(rank,2);
C = coeff(rank,3);
D = coeff(rank,4);
E = coeff(rank,5);
figure;plot(abs(A),'o-');
hold on;plot(abs(B),'o-')
hold on;plot(abs(C),'o-')
