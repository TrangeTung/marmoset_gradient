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
lmask = AtlasMask>20;
a=65;b=86;c=64;d=512; % 4D-EPI size
lmask = AtlasMask>20;
Xmask = imresize3(lmask,round([a,b,c]/2),'nearest');
NUM=5522;


%{
for idx = 5:5:62
    idx
    path = [ExpTable.Monkey{idx},filesep];
    RARE =  ExpTable.RARE(idx);
    SEEPI = sort([ExpTable.SEEPI1(idx),ExpTable.SEEPI2(idx)],'ascend');
    GEEPI = sort([  ExpTable.EPI4D1(idx),ExpTable.EPI4D2(idx),...
        ExpTable.EPI4D3(idx),ExpTable.EPI4D4(idx),...
        ExpTable.EPI4D5(idx),ExpTable.EPI4D6(idx),...
        ExpTable.EPI4D7(idx),ExpTable.EPI4D8(idx)],'ascend');
    GEEPI(isnan(GEEPI)) = [];
   
    
    fclose('all');
    spm('defaults', 'FMRI');
    set(spm('CreateIntWin','off'),'Visible','on');

    
    all_func = MY_find_images_in_all_scans(path,'Results',{GEEPI(:)},'^snSRGrsmUw2dseq','.nii',[1 Inf],'separate_cells');
    realign_mlb = MY_get_default_realign_batch_struct(all_func);
    realign_mlb{1}.spm.spatial.realign.estwrite.roptions = [];
    F = spm_figure('GetWin');
    disp('Start to process realignment !')
    spm_jobman('run',realign_mlb);
       
end

%}



%% Fig.A Head motion Distribution
a=0; clear *FD*
for idx = 1:62
    idx
    path = [ExpTable.Monkey{idx},filesep];
    RARE =  ExpTable.RARE(idx);
    SEEPI = sort([ExpTable.SEEPI1(idx),ExpTable.SEEPI2(idx)],'ascend');
    GEEPI = sort([  ExpTable.EPI4D1(idx),ExpTable.EPI4D2(idx),...
        ExpTable.EPI4D3(idx),ExpTable.EPI4D4(idx),...
        ExpTable.EPI4D5(idx),ExpTable.EPI4D6(idx),...
        ExpTable.EPI4D7(idx),ExpTable.EPI4D8(idx)],'ascend');
    GEEPI(isnan(GEEPI)) = [];
    
    for kk=GEEPI
        GEPath = fullfile(path,'Results',num2str(kk));
        
        rp = load(fullfile(GEPath,'rp_smUw2dseq.txt'));
        rp(:,1:3) = rp(:,1:3)/6;
        rp(:,4:6) = rp(:,4:6)*15;
        drp = rp - [rp(1,:);rp(1:end-1,:)]; % mm
        FD = sqrt(drp(:,1).^2+drp(:,2).^2+drp(:,3).^2+drp(:,4).^2+drp(:,5).^2+drp(:,6).^2)*1000; % um
        
        a=a+1;
        maxFD_no(:,a) = max(FD);
        meanFD_no(:,a) = mean(FD);
        
        rp = load(fullfile(GEPath,'rp_snSRGrsmUw2dseq.txt'));
        rp(:,1:3) = rp(:,1:3)/1;
        rp(:,4:6) = rp(:,4:6)*15;
        drp = rp - [rp(1,:);rp(1:end-1,:)]; % mm
        FD = sqrt(drp(:,1).^2+drp(:,2).^2+drp(:,3).^2+drp(:,4).^2+drp(:,5).^2+drp(:,6).^2)*1000; % um
        
        
        maxFD_22(:,a) = max(FD);
        meanFD_22(:,a) = mean(FD);
    end
end

F = figure;
subplot(2,1,1);
histogram(maxFD_no,50); xlim([0 1000])
set(gca,'box','off','tickdir','out');
subplot(2,1,2);
histogram(maxFD_22,50); xlim([0 100])
set(gca,'box','off','tickdir','out');

%% Fig.B Between Group similarity

a1=0; a2=0; clear *IDX *Gd
NUM=5522;
for idx = 1:62
    idx
    path = [ExpTable.Monkey{idx},filesep];
    RARE =  ExpTable.RARE(idx);
    SEEPI = sort([ExpTable.SEEPI1(idx),ExpTable.SEEPI2(idx)],'ascend');
    GEEPI = sort([  ExpTable.EPI4D1(idx),ExpTable.EPI4D2(idx),...
        ExpTable.EPI4D3(idx),ExpTable.EPI4D4(idx),...
        ExpTable.EPI4D5(idx),ExpTable.EPI4D6(idx),...
        ExpTable.EPI4D7(idx),ExpTable.EPI4D8(idx)],'ascend');
    GEEPI(isnan(GEEPI)) = [];
    
    for kk=GEEPI
        GEPath = fullfile(path,'Results',num2str(kk));
        
        rp = load(fullfile(GEPath,'rp_smUw2dseq.txt'));
        rp(:,1:3) = rp(:,1:3)/6;
        rp(:,4:6) = rp(:,4:6)*15;
        drp = rp - [rp(1,:);rp(1:end-1,:)]; % mm
        FD = sqrt(drp(:,1).^2+drp(:,2).^2+drp(:,3).^2+drp(:,4).^2+drp(:,5).^2+drp(:,6).^2)*1000; % um
        
        a=a+1;
        
        dest = fullfile(WholePath,'DiffusionGradient_ION');
        dest_ = fullfile(dest,replace(ExpTable.Monkey{idx},WholePath,''),num2str(kk));
        if ~exist(dest_,'dir');mkdir(dest_);end
        load(fullfile(dest_,'rawdata.mat'));
        nanidx = rawdata.nanidx;
        V = rawdata.data;
        V_ = nan(NUM,size(V,2));
        V_(~nanidx,:)=V;
        %V_ = fillmissing(V_,'linear');
        FC = corr(V_');
        
        if max(FD)>250 % N=65
            a1=a1+1;
            if a1==1;LargeGd=zeros([size(FC),65]);end
            LargeGd(:,:,a1)=FC;
        end
        
        if max(FD)<50 % N=90
            a2=a2+1;
            if a2==1;LittleGd=zeros([size(FC),50]);end
            LittleGd(:,:,a2)=FC;
        end
    end
end

dest_ = fullfile(WholePath,'DiffusionGradient_ION','Mean');
load(fullfile(dest_,'Gradients.mat'));
Gdf = Gd;

% Large Motion
FCz = mean(LargeGd,3); % clear LargeGd
NANidx = isnan(nanmean(FCz,1));
G = GradientMaps('kernel','cs','approach','dm','n_components',20,'align','ja');
G = G.fit(FCz(~NANidx,~NANidx),'Reference',Gdf(~NANidx,:));

Gd = G.aligned{1}*50;
Gx = zeros(numel(NANidx),size(Gd,2));
Gx(~NANidx,:) = Gd;
Gd=Gx;

M = funmask(Gd,Xmask);
xhdr=spm_vol(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.nii.gz'));
M_ = zeros([xhdr.dim,size(M,4)]);
for m=1:size(M,4)
    X = M(:,:,:,m);
    X = imresize3(X,xhdr.dim,'linear');
    X = smooth3(X,'box',7);
    M_(:,:,:,m) = X;
end

dest_ = fullfile(WholePath,'DiffusionGradient_ION','Mean_LargeMotion');
if ~exist(dest_,'dir');mkdir(dest_);end
save(fullfile(dest_,'G.mat'),'G');
save(fullfile(dest_,'Gradients.mat'),'Gd');

filename = fullfile(dest_,'PG_align.nii');
clear Vol
for k = 1:size(M,4)
    Vol(k,1) = struct(  'fname',    filename,...
                        'dim',      xhdr.dim,...
                        'mat',      xhdr.mat,...
                        'n',        [k,1],...
                        'pinfo',    [1;0;0],...
                        'descrip',  '',...
                        'dt',       [16 0]);
end
spm_write_vol_4D(Vol,M_); 

delete(fullfile(dest_,'fx.func.gii'));
bar = [-1 -0.0001 0.0001 1]*5;
for k = 1:5
    F = SurfaceMap(fullfile(dest_,'PG_align.nii'),bar/k,codepath,k,'marmoset');
    TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'.tiff']);
    saveas(F,TIFname)
    close all
    X = imread(fullfile(dest_,['Surface_PG_align_',num2str(k),'.tiff']));
    x_ = X(61:400,65:680,:);
    y_ = X(531:870,81:700,:);
    z = cat(2,x_,y_);
    TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'_cut.tiff']);
    imwrite(z,TIFname);
    %pause(3);
    
    
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
    TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'_cut1.tiff']);
    imwrite(uint8(z),TIFname);
end




% Little Motion
FCz = mean(LittleGd,3);  %clear LittleGd
NANidx = isnan(nanmean(FCz,1));
G = GradientMaps('kernel','cs','approach','dm','n_components',20,'align','ja');
G = G.fit(FCz(~NANidx,~NANidx),'Reference',Gdf(~NANidx,:));

Gd = G.aligned{1}*50;
Gx = zeros(numel(NANidx),size(Gd,2));
Gx(~NANidx,:) = Gd;
Gd = Gx;

M = funmask(Gd,Xmask);
xhdr=spm_vol(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.nii.gz'));
M_ = zeros([xhdr.dim,size(M,4)]);
for m=1:size(M,4)
    X = M(:,:,:,m);
    X = imresize3(X,xhdr.dim,'linear');
    X = smooth3(X,'box',7);
    M_(:,:,:,m) = X;
end

dest_ = fullfile(WholePath,'DiffusionGradient_ION','Mean_LittleMotion');
if ~exist(dest_,'dir');mkdir(dest_);end
save(fullfile(dest_,'G.mat'),'G');
save(fullfile(dest_,'Gradients.mat'),'Gd');

filename = fullfile(dest_,'PG_align.nii');
clear Vol
for k = 1:size(M,4)
    Vol(k,1) = struct(  'fname',    filename,...
                        'dim',      xhdr.dim,...
                        'mat',      xhdr.mat,...
                        'n',        [k,1],...
                        'pinfo',    [1;0;0],...
                        'descrip',  '',...
                        'dt',       [16 0]);
end
spm_write_vol_4D(Vol,M_); 

delete(fullfile(dest_,'fx.func.gii'));
bar = [-1 -0.0001 0.0001 1]*5;
for k = 1:5
    bar = [-1 -0.0001 0.0001 1]*5;
    F = SurfaceMap(fullfile(dest_,'PG_align.nii'),bar/k,codepath,k,'marmoset');
    TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'.tiff']);
    saveas(F,TIFname)
    close all
    X = imread(fullfile(dest_,['Surface_PG_align_',num2str(k),'.tiff']));
    x_ = X(61:400,65:680,:);
    y_ = X(531:870,81:700,:);
    z = cat(2,x_,y_);
    TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'_cut.tiff']);
    imwrite(z,TIFname);
    %pause(3);
    
    
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
    TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'_cut1.tiff']);
    imwrite(uint8(z),TIFname);
end



% cross group similarity
dest_ = fullfile(WholePath,'DiffusionGradient_ION','Mean_LargeMotion');
load(fullfile(dest_,'Gradients.mat'));
Gd_Large = Gd;
dest_ = fullfile(WholePath,'DiffusionGradient_ION','Mean_LittleMotion');
load(fullfile(dest_,'Gradients.mat'));
Gd_Little = Gd;

F = figure('Position', [110 659 1689 319]);
for lp = 1:4
   subplot(1,4,lp);
   plot(Gd_Large(:,lp),Gd_Little(:,lp),'o',...
       'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
   [b,bint,r,rint,stats] = regress(Gd_Large(:,lp),[ones(size(Gd,1),1),Gd_Little(:,lp)]);
   [stats(1),stats(3)]
   xlim([-10 10]);ylim([-10 10]);
   set(gca,'box','off','fontsize',15,'tickdir','out');
end
