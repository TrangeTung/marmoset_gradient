% clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
addpath(genpath(codepath));
AtlasNII = fullfile(codepath,'marmoset_atlas','Trange_atlas_RikenBMA_cortex.nii');
AtlasMask = spm_read_vols_4D(spm_vol(AtlasNII));
%% Gradient Fingerprints

RegionName = {'pSomMot';'AUD';'pVis';'MTcomplex';'sVis';'RSC';
    'PFC';'lParietal';'lTemporal';'mParietal';'Insular'};
Marmoset_AtlasName = fullfile(codepath,'MBM_v3.0.1','atlas_MBM_cortex_vM.nii.gz');
Marmoset_Labels = {[16:18];[20:23];[46:48];[49:50];[51:53];[41:45];...
    [1:7];[36:40];[27:30];[8:11];[24:26]};
Human_AtlasName = fullfile('E:\HumanAtlas\functional',...
    'Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_1mm.nii.gz');
Human_Labels = {[67:108];[134:147];[1:66];[146:181];[359:377];[471:478];
    [303:356];[485:500];[289:302];[389:410];[227:235,260:267]};

% Human LandMarkers2
%{
ihdr = spm_vol(Human_AtlasName);
img = spm_read_vols(ihdr);
Lx = zeros([size(img),numel(Human_Labels)]);
for hl=1:numel(Human_Labels)
    ls = Human_Labels{hl};
    X=zeros(size(img));
    for sl=1:numel(ls);X(img==ls(sl))=1;end
    Lx(:,:,:,hl)=X;
end
path_out=fullfile(codepath,'Human_LandMarkers2.nii');
for k = 1:size(Lx,4)
    Vol(k,1) = struct(  'fname',    path_out,...
                        'dim',      ihdr.dim,...
                        'mat',      ihdr.mat,...
                        'n',        [k,1],...
                        'pinfo',    [1;0;0],...
                        'descrip',  '',...
                        'dt',       [16 0]);
end
spm_write_vol_4D(Vol,Lx);

NIIfile = fullfile(codepath,'Human_LandMarkers2.nii');
GIIname = fullfile(codepath,['L.Human_LandMarkers2.func.gii']);
system([fullfile(codepath,'workbench','bin_windows64','wb_command.exe'),...
        ' -volume-to-surface-mapping ',...
        NIIfile,' ',...
        fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','lh.inflated.gii'),' ',...
        GIIname,' ',...
        '-ribbon-constrained ',...
        fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','lh.white.gii'),' ',...
        fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','lh.pial.gii')]);
%}
% Marmoset LandMarkers2
%{
ihdr = spm_vol(Marmoset_AtlasName);
img = spm_read_vols(ihdr);
Lx = zeros([size(img),numel(Marmoset_Labels)]);
for hl=1:numel(Marmoset_Labels)
    ls = Marmoset_Labels{hl};
    X=zeros(size(img));
    for sl=1:numel(ls);X(img==ls(sl))=1;end
    Lx(:,:,:,hl)=X;
end
path_out=fullfile(codepath,'Marmoset_LandMarkers2.nii');
for k = 1:size(Lx,4)
    Vol(k,1) = struct(  'fname',    path_out,...
                        'dim',      ihdr.dim,...
                        'mat',      ihdr.mat,...
                        'n',        [k,1],...
                        'pinfo',    [1;0;0],...
                        'descrip',  '',...
                        'dt',       [16 0]);
end
spm_write_vol_4D(Vol,Lx);
NIIfile = fullfile(codepath,'Marmoset_LandMarkers2.nii');
GIIname = fullfile(codepath,['Marmoset_LandMarkers2.func.gii']);
system([fullfile(codepath,'workbench','bin_windows64','wb_command.exe'),...
            ' -volume-to-surface-mapping ',...
            NIIfile,' ',...
            fullfile(codepath,'MBM_v3.0.1','surfFS.lh.graymid.surf.gii'),' ',...
            GIIname,' ',...
            '-ribbon-constrained ',...
            fullfile(codepath,'MBM_v3.0.1','surfFS.lh.white.surf.gii'),' ',...
            fullfile(codepath,'MBM_v3.0.1','surfFS.lh.pial.surf.gii')]);
%}


GIIL = gifti(fullfile(codepath,['L.Human_LandMarkers2.func.gii']));
LMsk = [GIIL.cdata];

load(fullfile(WholePath,'NC_Resubmit','Cross-Species','Child_Adolescent_Gradient.mat'));
ChildG = zeros(20484,2);
ChildG(~child_NANidx,:) = child_G.gradients{1}(:,1:2);
ChildV = ChildG(1:10242,:)'*LMsk;

AdolescentG = zeros(20484,2);
AdolescentG(~adolescent_NANidx,:) = adolescent_G.gradients{1}(:,1:2);
AdolescentV = AdolescentG(1:10242,:)'*LMsk;


dest = fullfile(WholePath,'NC_Resubmit','Cross-Species','adult');
NIIDir = dir(fullfile(dest,'*.nii.gz'));
for nl=1:2%numel(NIIDir)
    
    GIInameL = fullfile(dest,['L.',replace(NIIDir(nl).name,'.nii.gz',''),'.func.gii']);
    GIInameR = fullfile(dest,['R.',replace(NIIDir(nl).name,'.nii.gz',''),'.func.gii']);
    GIIL = gifti(GIInameL);
    GIIR = gifti(GIInameR);
    HumanGII = [GIIL.cdata;GIIR.cdata];
    AdultG(:,nl) = HumanGII;
    
end
AdultV = AdultG(1:10242,:)'*LMsk;
    

GIIname = fullfile(WholePath,'DiffusionGradient_ION','Mean','fx.func.gii');
MarmosetGII = gifti(GIIname);
GIIL = gifti(fullfile(codepath,['Marmoset_LandMarkers2.func.gii']));
MarmosetV = MarmosetGII.cdata(:,1:2)'* GIIL.cdata;



X = [MarmosetV(1,:);MarmosetV(2,:);
    ChildV(1,:);ChildV(2,:);
    AdultV(1,:);AdultV(2,:)];
Names={'Marmoset-G1';'Marmoset-G2';
    'Children-G1';'Children-G2';
    'Adults-G1';'Adults-G2';};
[M,Q]=community_louvain(corr(X'),.8);
[Ms,s] = sort(M);

F = figure;
h = MY_Circos_plot(abs(corr(X(s,:)')));


F = figure('Position', [680 585 1419 393]);
subplot(1,3,1); imagesc(corr(X'));
colormap('jet');caxis([0 1])
set(gca,'yticklabel',Names);
set(gca,'xticklabel',Names,'xticklabelRotation',90);
subplot(1,3,2); imagesc(corr(X(s,:)'));
colormap('jet');caxis([0 1])
set(gca,'yticklabel',Names(s));
set(gca,'xticklabel',Names(s),'xticklabelRotation',90);
subplot(1,3,3);
colormap('jet');caxis([0 1])
colorbar

F = figure('Position', [680 585 141 393]);
X = (X-mean(X))./std(X,0,2);
imagesc(X);
set(gca,'xticklabel',RegionName);
set(gca,'yticklabel',Names,'xticklabelRotation',90);
colormap('autumn');caxis([-1 1])
colorbar('EastOutside')


figure;imagesc(corr(X'));caxis([0 1])
figure;imagesc(abs(corr(X')));caxis([0 1])

figure;plot(-MarmosetV(1,:),ChildV(1,:),'.')
figure;plot(-MarmosetV(1,:),AdultV(2,:),'.')






% marmoset
dest = fullfile(WholePath,'DiffusionGradient_ION','Mean');




% human child
dest = fullfile(WholePath,'NC_Resubmit','Cross-Species');
load(fullfile(dest,'child_SparseFCz.mat'));
X = std(SparseFCz,0,2);
NANidx = X==0;
Input = SparseFCz(~NANidx,~NANidx)/2+SparseFCz(~NANidx,~NANidx)'/2;
G = GradientMaps('kernel','cs','approach','dm','n_components',10);
G = G.fit(Input);
child_G = G;
child_NANidx = NANidx;

dest = fullfile(WholePath,'NC_Resubmit','Cross-Species');
load(fullfile(dest,'adolescent_SparseFCz.mat'));
X = std(SparseFCz,0,2);
NANidx = X==0;
Input = SparseFCz(~NANidx,~NANidx)/2+SparseFCz(~NANidx,~NANidx)'/2;
G = GradientMaps('kernel','cs','approach','dm','n_components',10);
G = G.fit(Input);
adolescent_G = G;
adolescent_NANidx = NANidx;

C = Cm(1:size(Cm,1)/2,1:size(Cm,1)/2);
X = std(C); clear STDd Cm
NANidx = X<0.03;
Input = C(~NANidx,~NANidx)/2+C(~NANidx,~NANidx)'/2;
G = GradientMaps('kernel','cs','approach','dm','n_components',10);
G = G.fit(Input);

% human adult
dest = fullfile(WholePath,'NC_Resubmit','Cross-Species','adult');
NIIDir = dir(fullfile(dest,'*.nii.gz'));
for nl=1:numel(NIIDir)
    
    NIIfile = fullfile(dest,NIIDir(nl).name);
    GIInameL = fullfile(dest,['L.',replace(NIIDir(nl).name,'.nii.gz',''),'.func.gii']);
    GIInameR = fullfile(dest,['R.',replace(NIIDir(nl).name,'.nii.gz',''),'.func.gii']);

    if ~exist(GIInameL,'file')
        system([fullfile(codepath,'workbench','bin_windows64','wb_command.exe'),...
                ' -volume-to-surface-mapping ',...
                NIIfile,' ',...
                fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','lh.inflated.gii'),' ',...
                GIInameL,' ',...
                '-ribbon-constrained ',...
                fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','lh.white.gii'),' ',...
                fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','lh.pial.gii')]);

        system([fullfile(codepath,'workbench','bin_windows64','wb_command.exe'),...
                ' -volume-to-surface-mapping ',...
                NIIfile,' ',...
                fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','rh.inflated.gii'),' ',...
                GIInameR,' ',...
                '-ribbon-constrained ',...
                fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','rh.white.gii'),' ',...
                fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','rh.pial.gii')]);
    end
end

dest = fullfile(WholePath,'NC_Resubmit','Cross-Species','adult');
NIIfile = fullfile(codepath,'PlanB_code','Human_LandMarkers.nii');
GIInameL = fullfile(dest,['L.Human_LandMarkers.func.gii']);
GIInameR = fullfile(dest,['R.Human_LandMarkers.func.gii']);

if ~exist(GIInameL,'file')
    system([fullfile(codepath,'workbench','bin_windows64','wb_command.exe'),...
            ' -volume-to-surface-mapping ',...
            NIIfile,' ',...
            fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','lh.inflated.gii'),' ',...
            GIInameL,' ',...
            '-ribbon-constrained ',...
            fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','lh.white.gii'),' ',...
            fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','lh.pial.gii')]);

    system([fullfile(codepath,'workbench','bin_windows64','wb_command.exe'),...
            ' -volume-to-surface-mapping ',...
            NIIfile,' ',...
            fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','rh.inflated.gii'),' ',...
            GIInameR,' ',...
            '-ribbon-constrained ',...
            fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','rh.white.gii'),' ',...
            fullfile(codepath,'HCP_MMP_v1.0','fsaverage5','rh.pial.gii')]);
end

GIIL = gifti(fullfile(dest,['L.Human_LandMarkers.func.gii']));
LMsk = [GIIL.cdata];

ChildG = zeros(20484,5);
ChildG(~child_NANidx,:) = child_G.gradients{1}(:,1:5);
ChildV = ChildG(1:10242,:)'*LMsk;

AdolescentG = zeros(20484,5);
AdolescentG(~adolescent_NANidx,:) = adolescent_G.gradients{1}(:,1:5);
AdolescentV = AdolescentG(1:10242,:)'*LMsk;

dest = fullfile(WholePath,'NC_Resubmit','Cross-Species','adult');
NIIDir = dir(fullfile(dest,'*.nii.gz'));
for nl=1:numel(NIIDir)
    
    GIInameL = fullfile(dest,['L.',replace(NIIDir(nl).name,'.nii.gz',''),'.func.gii']);
    GIInameR = fullfile(dest,['R.',replace(NIIDir(nl).name,'.nii.gz',''),'.func.gii']);
    GIIL = gifti(GIInameL);
    GIIR = gifti(GIInameR);
    HumanGII = [GIIL.cdata;GIIR.cdata];
    AdultG(:,nl) = HumanGII;
    
end
AdultV = AdultG(1:10242,:)'*LMsk;
    
    
    
    
dest = fullfile(codepath,'PlanB_code');
NIIfile = fullfile(dest,'Marmoset_LandMarkers.nii');
GIInameL = fullfile(dest,['L.Marmoset_LandMarkers.func.gii']);
GIInameR = fullfile(dest,['R.Marmoset_LandMarkers.func.gii']);
system([fullfile(codepath,'workbench','bin_windows64','wb_command.exe'),...
        ' -volume-to-surface-mapping ',...
        NIIfile,' ',...
        fullfile(codepath,'MBM_v3.0.1','surfFS.lh.graymid.surf.gii'),' ',...
        GIInameL,' ',...
        '-ribbon-constrained ',...
        fullfile(codepath,'MBM_v3.0.1','surfFS.lh.white.surf.gii'),' ',...
        fullfile(codepath,'MBM_v3.0.1','surfFS.lh.pial.surf.gii')]);

system([fullfile(codepath,'workbench','bin_windows64','wb_command.exe'),...
        ' -volume-to-surface-mapping ',...
        NIIfile,' ',...
        fullfile(codepath,'MBM_v3.0.1','surfFS.lh.graymid.surf.gii'),' ',...
        GIInameR,' ',...
        '-ribbon-constrained ',...
        fullfile(codepath,'MBM_v3.0.1','surfFS.lh.white.surf.gii'),' ',...
        fullfile(codepath,'MBM_v3.0.1','surfFS.lh.pial.surf.gii')]);
                
GIIname = fullfile(WholePath,'DiffusionGradient_ION','Mean','fx.func.gii');
MarmosetGII = gifti(GIIname);
GIIL = gifti(GIInameL);
GIIR = gifti(GIInameR);
MarmosetV = MarmosetGII.cdata(:,1:5)'* GIIL.cdata;


X1 = [MarmosetV(1,:);ChildV(1,:);AdultV(2,:);];
Y1 = (X1-mean(X1,2)) ./std(X1,0,2);
X2 = [MarmosetV(2,:);ChildV(2,:);AdultV(1,:);];
Y2 = (X2-mean(X2,2)) ./std(X2,0,2);
Z_ = corr([Y1;Y2]');
figure;imagesc(Z_)

for loop=1:20
    loop
Cx = nchoosek(1:31,loop+10);
err = zeros(size(Cx,1),1);
parfor cl=1:size(Cx,1)
    Z1 = Y1(:,Cx(cl,:));
    Z2 = Y2(:,Cx(cl,:));
    Z_ = corr([Z1;Z2]');
    Ref = [[ones(3,3),zeros(3,3)];[zeros(3,3),ones(3,3)]];
    Z0 = abs((Z_)-Ref);

    err(cl) = sum(Z0(:));
end
ERRpin(loop)=find(err==min(err));
ERR(loop)=min(err);
end


loop=4;
cl=ERRpin(loop);
Cx = nchoosek(1:31,loop+10);
Z1 = Y1(:,Cx(cl,:));
Z2 = Y2(:,Cx(cl,:));
Z_ = corr([Z1;Z2]');
figure;imagesc(Z_)









dest = fullfile(codepath,'PlanB_code');
NIIfile = fullfile(dest,'Marmoset_LandMarkers.nii');
fMRI_4D = spm_read_vols(spm_vol(NIIfile));
a=65;b=86;c=64;d=512; % 4D-EPI size
lmask = AtlasMask>20;
Xmask = imresize3(lmask,round([a,b,c]/2),'nearest');
for bs=1:size(fMRI_4D,4)
Makr(:,:,:,bs) = imresize3(fMRI_4D(:,:,:,bs),round([a,b,c]/2),'linear');
end
Vd = fmask(Makr,Xmask);
MarmosetV = Gd'*Vd;


X1 = [MarmosetV(1:2,:);ChildV(1:2,:);AdultV(1:2,:);];
Y1 = (X1-mean(X1)) ./std(X1,0,2);

M=link_communities(corr(Y1));
[M,Q]=community_louvain(XX,.8);
M = clique_communities(abs(corr(Y1')));

 x=0;
 for gm=0.1:0.05:1
[Ci,Q]=modularity_und(abs(corr(Y1')),0.80);
x=x+1;
dif1(x) = sum(abs(Ci-[1;2;1;2;2;1]));
dif2(x) = sum(abs(Ci-[2;1;2;1;1;2]));

 end

X2 = [MarmosetV(2,:);ChildV(2,:);AdultV(1,:);];
Y2 = (X2-mean(X2)) ./std(X2,0,2);
for loop=1:20
    loop
Cx = nchoosek(1:31,loop+10);
err = zeros(size(Cx,1),1);
parfor cl=1:size(Cx,1)
    Z1 = Y1(:,Cx(cl,:));
    Z2 = Y2(:,Cx(cl,:));
    Z_ = corr([Z1;Z2]');
    Ref = [[ones(3,3),zeros(3,3)];[zeros(3,3),ones(3,3)]];
    Z0 = abs((Z_)-Ref);

    err(cl) = sum(Z0(:));
end
ERRpin(loop)=find(err==min(err));
ERR(loop)=min(err);
end

