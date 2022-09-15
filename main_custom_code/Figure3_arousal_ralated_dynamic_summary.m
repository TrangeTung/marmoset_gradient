clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));

ExpRecording = [codepath 'ExpIndex.xlsx'];
[~,~,CellData] = xlsread(ExpRecording);
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
WholeMaskNII = fullfile(codepath,'marmoset_atlas','Trange_Template_mask_V3.nii');
lmask = spm_read_vols(spm_vol(WholeMaskNII));


bins = 10; W1=[];  W2=[]; Wr=[];
for bl=1:bins
    dest_ = fullfile(WholePath,'DiffusionGradient_ION','Dynamic_VIweighted',num2str(bl,'%02d'));
    load(fullfile(dest_,'G.mat'));
    lambda = G.lambda{1};
    W1(:,bl) = lambda;
end


for bl=1:bins
    dest_ = fullfile(WholePath,'DiffusionGradient_ION','Dynamic_VIremoval',num2str(bl,'%02d'));
    load(fullfile(dest_,'G.mat'));
    lambda = G.lambda{1};
    Wr(:,bl) = lambda;
end

F = figure( 'Position', [680 464 441 514]);
X = W1(1:4,:);
Xr = Wr(1:4,:);
plot(X','linewidth',2);
hold on;plot([0,0],[0.04 0.06],'k');
plot(Xr','color',[.5 .5 .5],'linewidth',2);
xlim([-1 11])





%% topographic flow
ICApath = fullfile(WholePath,'SubMit','ICA_final','processed');
filename = fullfile(ICApath,'normICAv1_all_vtsplit_premotor.nii.gz');
ihdr = spm_vol(filename);
img = spm_read_vols(ihdr);
ICAmask = img>0;
ICAname = {'ventralsoma';'dorsalsoma';'frontalpole';'parahip';'OPFC';...
    'AUD';'frontalpaietal';'DMN';'visualMT';'visual0';'visual1';...
    'visual2';'visual3';'mACC';'premotor'};

% gradient
bins = 10; V_ION=[];
for bl=1:bins
    dest = fullfile(WholePath,'DiffusionGradient_ION','Dynamic_VIweighted',num2str(bl,'%02d'));
    filename = fullfile(dest,'PG_align.nii');
    ihdr = spm_vol(filename);
    PG = spm_read_vols(ihdr);
    [a,b,c,d] = size(PG);
    
    for i=1:numel(ICAname)
        X = ICAmask(:,:,:,i);
        Xmask = imresize3(X,[a,b,c],'nearest');
        V_ION(i,bl,:) = nanmean(fmask(PG,Xmask),1); 
    end
end

Map = [0:255;zeros(1,256);255:-1:0]';
MyMap = 1/255*interp1(1:256,Map,10:25:256);
F = figure( 'Position', [680 607 409 371]);
for iy = 1:9
    for ix=1%:numel(ICAname)
       X = V_ION(ix,iy:iy+1,2);
       Y = V_ION(ix,iy:iy+1,1);
       quiver(X(1),Y(1),X(2)-X(1),Y(2)-Y(1),'color',MyMap(iy,:),...
           'linewidth',1,'MaxHeadSize',1);
       hold on;
       if iy==5
          text(X(1),Y(1),ICAname{ix},'HorizontalAlignment','center') ;
       end
       
    end 
end
xlabel('Gradient 2');ylabel('Gradient 1');
set(gca,'linewidth',2,'fontsize',15,'box','off');


%% gradient similarity
Wmask = imresize3(lmask,[147,200,135],'nearest');
% structure gradient
dest = fullfile(WholePath,'DiffusionGradient_Structure');
filename = fullfile(dest,'StructurePG.nii');
ihdr = spm_vol(filename);
PG = spm_read_vols(ihdr);
[a,b,c,d] = size(PG);
V = fmask(PG,Wmask);

bins = 10; V_ION=[];
for bl=1:bins
    dest = fullfile(WholePath,'DiffusionGradient_ION','Dynamic_VIweighted',num2str(bl,'%02d'));
    filename = fullfile(dest,'PG_align.nii');
    ihdr = spm_vol(filename);
    PG = spm_read_vols(ihdr);
    [a,b,c,d] = size(PG);
    
    V_ION(:,:,bl) = fmask(PG,Wmask);
end
bins = 10; V_NIH=[];
for bl=1:bins
    dest = fullfile(WholePath,'DiffusionGradient_NIH','Dynamic_VIweighted',num2str(bl,'%02d'));
    filename = fullfile(dest,'PG_align.nii');
    ihdr = spm_vol(filename);
    PG = spm_read_vols(ihdr);
    [a,b,c,d] = size(PG);
    
    V_NIH(:,:,bl) = fmask(PG,Wmask);
end
bins = 10;
for bl=1:bins
    X = V(:,1:4);
    Y_ION = squeeze(V_ION(:,1:4,bl));
    Cx_ION = corr(X,Y_ION);
    Sim_ION(:,bl) = diag(Cx_ION);
    Y_NIH = squeeze(V_NIH(:,1:4,bl));
    Cx_NIH = corr(X,Y_NIH);
    Sim_NIH(:,bl) = diag(Cx_NIH);    
end
F = figure( 'Position',  [680 678 752 300]);
subplot(1,2,1);
plot(Sim_ION(1,:),'k','linewidth',2); hold on
plot(Sim_NIH(1,:),'color',[.5 .5 .5],'linewidth',2); 
ylabel('Spatial correlation (C.C.)');
set(gca,'linewidth',2,'fontsize',15,'xtick','','box','off');
subplot(1,2,2);
plot(Sim_ION(2,:),'k','linewidth',2); hold on
plot(-Sim_NIH(2,:),'color',[.5 .5 .5],'linewidth',2); 
ylabel('Spatial correlation (C.C.)');
set(gca,'linewidth',2,'fontsize',15,'xtick','','box','off');



%% Human
ICAname = {'Visual';'Som/Motor';'Sal/VentAttn';'DorsAttn';'Limbic';...
    'Control';'DMN'};
ICA_labels = {[1:81,501:581];[82:172,582:672];[173:233,673:733];...
    [234:288,734:788];[289:317,789:817];[318:374,848:874];[375:500,875:1000]};

V_HCP = [];
for bl=1:10
    
dest = fullfile(filepath,'Dynamic_VIweighted',num2str(bl,'%02d'));
load(fullfile(dest,'G.mat'));
Gs = G.aligned{1};
for i=1:numel(ICAname)
    X = mean(Gs(ICA_labels{i},:),1);
    V_HCP(i,bl,:) = X;
end
end



Map = [0:255;zeros(1,256);255:-1:0]';
MyMap = 1/255*interp1(1:256,Map,10:25:256);
F = figure( 'Position', [680 607 409 371]);
for iy = 1:9
    for ix=1:numel(ICAname)
       X = V_HCP(ix,iy:iy+1,2)*50;
       Y = V_HCP(ix,iy:iy+1,1)*50;
       quiver(X(1),Y(1),X(2)-X(1),Y(2)-Y(1),'color',MyMap(iy+1,:),...
           'linewidth',2,'MaxHeadSize',0.5);
       hold on;
       if iy==5
          text(X(1),Y(1),ICAname{ix},'HorizontalAlignment','center') ;
       end
       
    end 
end
xlabel('Gradient 2');ylabel('Gradient 1');
set(gca,'linewidth',2,'fontsize',15,'box','off');
