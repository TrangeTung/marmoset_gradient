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

%% ION Group
%{
a_=0;

a=65;b=86;c=64;d=512; % 4D-EPI size
lmask = AtlasMask>20;
Xmask = imresize3(lmask,round([a,b,c]/2),'nearest');
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

    for kk = GEEPI
        a_=a_+1;
        %{
        a_=a_+1;
        fprintf(['session: ',num2str(idx),';folder: ',num2str(kk),'\n'])
        
        GEPath = fullfile(path,'Results',num2str(kk));
        NIIfile = fullfile(GEPath,'snSRGrsmUw2dseq.nii');
        
        hdr = spm_vol(NIIfile);
        img = spm_read_vols_4D(hdr); 
        [a,b,c,d] = size(img);
        
        lmask = AtlasMask>20;
        for ii=1:d
           X = img(:,:,:,ii);
           X = imresize3(X,round([a,b,c]/2),'linear');
           if ii==1;I=zeros([round([a,b,c]/2),d]);end
           I(:,:,:,ii)=X;
        end
        I(isnan(I))=0;
        I = (I+flip(I,1))/2;
        Xmask = imresize3(lmask,round([a,b,c]/2),'nearest');
        S = fmask(I,Xmask);
        V = (S-mean(S,2))./std(S,0,2);
        nanidx = isnan(mean(V,2));
        V(nanidx,:)=[];
        
        rawdata.nanidx = nanidx;
        rawdata.data = V;
        
        dest = fullfile(WholePath,'DiffusionGradient');
        dest_ = fullfile(dest,replace(ExpTable.Monkey{idx},WholePath,''),num2str(kk));
        if ~exist(dest_,'dir');mkdir(dest_);end
        save(fullfile(dest_,'rawdata.mat'),'rawdata');
        clear I S V rawdata
        %}
        
        %
        
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
        
        if a_==1;Vec=zeros([size(FC),345]);end
        Vec(:,:,a_)=FC;
        %
    end
end


%
FCz = mean(Vec,3);  clear Vec
NANidx = isnan(nanmean(FCz,1));
G = GradientMaps('kernel','cs','approach','dm','n_components',20);
G = G.fit(FCz(~NANidx,~NANidx));

Gd = G.gradients{1}*50;
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

dest_ = fullfile(WholePath,'DiffusionGradient_ION','Mean');
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
    pause(3);
end
%


%}


%% NIH Group

%{


a_=0;

a=65;b=86;c=64;d=512; % 4D-EPI size
lmask = AtlasMask>20;
Xmask = imresize3(lmask,round([a,b,c]/2),'nearest');
NUM=5522;
NANps = [1:18,91:127,222:271];

NIH_path = 'E:\Marmoset_NIH_dataset';
Session = dir(fullfile(NIH_path,'*_session_*'));

for idx = 1:39
    
    
    path = fullfile(NIH_path,Session(idx).name);
    for ii=1:99;if ~exist(fullfile(path,'Results',num2str(ii)),'dir');break;end;end
    ScanNum = 1:ii-1;
    
    for nl = ScanNum
        
        fprintf(['marmoset: ',num2str(idx),';scan: ',num2str(nl),'\n'])
               
        %{
        GEPath = fullfile(path,'Results',num2str(nl));
        cd(GEPath)
        hdr = spm_vol(fullfile(GEPath,'sncRrsmUw2dseq.nii'));
        img = spm_read_vols_4D(hdr);
        
        a=65;b=86;c=64;d=512;
        lmask = AtlasMask>20;
        for ii=1:d
           X = img(:,:,:,ii);
           X = imresize3(X,round([a,b,c]/2),'linear');
           if ii==1;I=zeros([round([a,b,c]/2),d]);end
           I(:,:,:,ii)=X;
        end
        I(isnan(I))=0;
        I = (I+flip(I,1))/2;
        Xmask = imresize3(lmask,round([a,b,c]/2),'nearest');
        S = fmask(I,Xmask);
        V = (S-mean(S,2))./std(S,0,2);
        nanidx = isnan(mean(V,2));
        V(nanidx,:)=[];
        
        
        rawdata.nanidx = nanidx;
        rawdata.data = V;
        
        dest = fullfile(WholePath,'DiffusionGradient_NIH');
        dest_ = fullfile(dest,Session(idx).name,num2str(nl,'%02d'));
        if ~exist(dest_,'dir');mkdir(dest_);end
        save(fullfile(dest_,'rawdata.mat'),'rawdata');
        clear I S V rawdata        
        %}
        
        a_=a_+1;
        dest = fullfile(WholePath,'DiffusionGradient_NIH');
        dest_ = fullfile(dest,Session(idx).name,num2str(nl,'%02d'));
        if ~exist(dest_,'dir');mkdir(dest_);end
        load(fullfile(dest_,'rawdata.mat'));
        nanidx = rawdata.nanidx;
        V = rawdata.data;
        V_ = nan(NUM,size(V,2));
        V_(~nanidx,:)=V;
        %V_ = fillmissing(V_,'linear');
        FC = corr(V_');
        
        if a_==1;Vec=zeros([size(FC),278]);end
        Vec(:,:,a_)=FC;

    end
end

filepath = fullfile(WholePath,'DiffusionGradient_ION','Mean');
load(fullfile(filepath,'Gradients.mat'));
G_Ref = Gd;


FCz = mean(Vec,3); clear Vec
NANidx_ = find(isnan(nanmean(FCz,1)));
G = GradientMaps('kernel','cs','approach','dm','n_components',20,'align','pa');
NANnum = union(NANidx_,NANps);
NANidx = zeros(NUM,1);NANidx(NANnum)=1;
data = FCz(~NANidx,~NANidx);
G = G.fit(data,'reference',G_Ref(~NANidx,:));

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

dest_ = fullfile(WholePath,'DiffusionGradient_NIH','Mean');
if ~exist(dest_,'dir');mkdir(dest_);end
save(fullfile(dest_,'G.mat'),'G');


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

% delete(fullfile(dest_,'fx.func.gii'));
bar = [-1 -0.0001 0.0001 1]*6;
for k = 1:5
    F = SurfaceMap(fullfile(dest_,'PG_align.nii'),bar/k,codepath,k,'marmoset');
    TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'.tiff']);
    saveas(F,TIFname)
    close all
    X = imread(fullfile(dest_,['Surface_PG_align_',num2str(k),'.tiff']));

    x_ = X(61:400,65:680,:);
    y_ = X(531:870,81:700,:);
    z = cat(2,x_,y_);
    TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'_cut1.tiff']);
    imwrite(z,TIFname);

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
    TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'_cut.tiff']);
    imwrite(uint8(z),TIFname);
end



%}


%% Radar plot
%
% ica label
ICApath = fullfile(WholePath,'SubMit','ICA_final','processed');
filename = fullfile(ICApath,'normICAv1_all_vtsplit_premotor.nii.gz');
ihdr = spm_vol(filename);
img = spm_read_vols(ihdr);
ICAmask = img>0;
ICAname = {'ventralsoma';'dorsalsoma';'frontalpole';'parahip';'OPFC';...
    'AUD';'frontalpaietal';'DMN';'visualMT';'visual0';'visual1';...
    'visual2';'visual3';'mACC';'premotor'};
ICAcolor = {[061,139,198];[080,185,231];[235,051,066];[212,179,130];[207,062,079];...
    [062,181,174];[236,091,083];[245,173,051];[051,134,093];[225,230,051];...
    [171,194,119];[078,189,097];[057,166,097];[233,051,153];[067,084,160]};


% gradient
dest = fullfile(WholePath,'DiffusionGradient_ION','Mean');
filename = fullfile(dest,'PG_align.nii');
ihdr = spm_vol(filename);
PG_ION = spm_read_vols(ihdr);
[a,b,c,d] = size(PG_ION);
% clear V_ION
% for i=1:numel(ICAname)
%     X = ICAmask(:,:,:,i);
%     Xmask = imresize3(X,[a,b,c],'nearest');
%     V_ION{i,:} = (fmask(PG,Xmask)); 
% end
dest = fullfile(WholePath,'DiffusionGradient_NIH','Mean');
filename = fullfile(dest,'PG_align.nii');
ihdr = spm_vol(filename);
PG_NIH = spm_read_vols(ihdr);
[a,b,c,d] = size(PG_NIH);
clear V_
for i=1:numel(ICAname)
    X = ICAmask(:,:,:,i);
    Xmask = imresize3(X,[a,b,c],'nearest');
    V_{i,:} = (fmask(PG_NIH,Xmask)+fmask(PG_ION,Xmask))/2; 
end


dest = fullfile(WholePath,'DiffusionGradient_Structure','.');
filename = fullfile(dest,'StructurePG.nii');
ihdr = spm_vol(filename);
PG = spm_read_vols(ihdr);
[a,b,c,d] = size(PG);
clear V_
for i=1:numel(ICAname)
    X = ICAmask(:,:,:,i);
    Xmask = imresize3(X,[a,b,c],'nearest');
    V_{i,:} = fmask(PG,Xmask); 
end


Q = cell2mat(V_);

F = figure('Position', [438 628 1336 350]);
for gl=1:4

X=[];V=[];R=[];
for i=1:numel(ICAname);R(i) = median(V_{i}(:,gl));end
[~,Reorder] = sort(R,'descend');

Y = sort(Q(:,gl),'ascend'); YN = numel(Y);
LOW5 = Y(round(YN*20/100));
TOP5 = Y(round(YN*80/100));
for i=1:numel(Reorder)
    X = V_{Reorder(i)}(:,gl);
    if mean(X)<LOW5 || mean(X)>TOP5;h(gl,i)=1;else;h(gl,i)=0;end
    V=cat(1,V,[X(:),repmat(i,numel(X),1)]);
end
MyMap = cell2mat(ICAcolor)/255;
subplot(1,4,gl);
boxplot(V(:,1),V(:,2),'Symbol','','Colors' ,MyMap(Reorder,:),...
    'LabelOrientation','horizontal','Orientation','horizontal',...
'Labels',ICAname(Reorder),'PlotStyle','compact');
set(gca,'box','off','linewidth',1,'fontsize',12,'tickdir','out')


end

Reorder = [2 1 6 3 4 5 7:15];
F = draw_radar([V_ION(Reorder,1),V_NIH(Reorder,1)/3*2],[-4 6],[-0 0],ICAname(Reorder));
saveas(F,fullfile(dest,'Gradient1.tiff'))
F = draw_radar([V_ION(Reorder,2),V_NIH(Reorder,2)/3*2],[-4 6]/2,[-0 0],ICAname(Reorder));
saveas(F,fullfile(dest,'Gradient2.tiff'))
F = draw_radar([V_ION(Reorder,3),V_NIH(Reorder,3)/3*2],[-2 4]/2,[-0 0],ICAname(Reorder));
saveas(F,fullfile(dest,'Gradient3.tiff'))
F = draw_radar([V_ION(Reorder,4),V_NIH(Reorder,4)/3*2],[-1 2],[-0 0],ICAname(Reorder));
saveas(F,fullfile(dest,'Gradient4.tiff'))
F = draw_radar([V_ION(Reorder,5),V_NIH(Reorder,5)/3*2],[-1 1],[-0 0],ICAname(Reorder));
saveas(F,fullfile(dest,'Gradient5.tiff'))
close all


filepath = fullfile(WholePath,'DiffusionGradient_ION','Mean');
load(fullfile(filepath,'Gradients.mat'));

xx = Gd(:,1);yy=Gd(:,2);

F=figure;
for iii=1:numel(xx)
    r=0;g=0;b=0;
    if yy(iii)>0;r=min([1*yy(iii)/4,1]);end
    if  xx(iii)>1;g=min([1,1*(xx(iii))/6]);end
    if yy(iii)<0;b=min([1,-1*(yy(iii))/6]);end
    plot(yy(iii),xx(iii),'o','markersize',3,'markerfacecolor',[r,g,b],'markeredgecolor','none');
    hold on;
end
xlabel('Gradient 2');
ylabel('Gradient 1');
set(gca,'box','off','fontsize',15,'linewidth',1.5);
%}    