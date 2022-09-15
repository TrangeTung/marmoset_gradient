clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
addpath(genpath(codepath));

AtlasExcel = fullfile(codepath,'marmoset_atlas','Trange_modifid_RikenBMA_atlas_labels.xlsx');
[~,~,CellData] = xlsread(AtlasExcel);
AtlasTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
atlas_label = AtlasTable.Labels;
Abbre = AtlasTable.Abbre;

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



filename = fullfile(codepath,'MarmosetBrainConnectivity.xlsx');
FLNe_log = xlsread(filename,'log(FLNe)');
FLNe_log(isnan(FLNe_log))=-6;
FLNe_log(51,:)=[];


G = GradientMaps('kernel','cs','approach','dm','n_components',5,'align','pa');
G = G.fit(FLNe_log,'Reference',Vf([1:50 52:116],1:5));
Gd = G.aligned{1}*100;



Gx = zeros(116,size(Gd,2));
Gx([1:50,52:end],:) = Gd;
ConG = fillmissing(Gx,'linear');

dest = fullfile(WholePath,'DiffusionGradient_Structure');
if ~exist(dest,'dir');mkdir(dest);end
%save(fullfile(dest,'G.mat'),'G');

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
    D_(:,:,:,k)=smooth3(X,'box',9);
end
filename = fullfile(dest,'StructurePG.nii');
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

delete(fullfile(dest,'fx.func.gii'));
bar = [-2 -0.0001 0.0001 2]/2;
for k = 1:5
    F = SurfaceMap(fullfile(dest,'StructurePG.nii'),bar,codepath,k,'marmoset');
    TIFname = fullfile(dest,['Surface_PG_',num2str(k),'.tiff']);
    saveas(F,TIFname)
    close all
    X = imread(fullfile(dest,['Surface_PG_',num2str(k),'.tiff']));
    x_ = X(61:400,65:680,:);
    y_ = X(531:870,81:700,:);
    z = cat(2,x_,y_);
    TIFname = fullfile(dest,['Surface_PG_',num2str(k),'_cut.tiff']);
    imwrite(z,TIFname);
end


%% Radar plot
% ica label
ICApath = fullfile(WholePath,'SubMit','ICA_final','processed');
filename = fullfile(ICApath,'normICAv1_all_vtsplit_premotor.nii.gz');
ihdr = spm_vol(filename);
img = spm_read_vols(ihdr);
ICAmask = img>0;
ICAname = {'ventralsoma';'dorsalsoma';'frontalpole';'parahip';'OPFC';...
    'AUD';'frontalpaietal';'DMN';'visualMT';'visual0';'visual1';...
    'visual2';'visual3';'mACC';'premotor'};

% gradient
dest = fullfile(WholePath,'DiffusionGradient_Structure');
filename = fullfile(dest,'StructurePG.nii');
ihdr = spm_vol(filename);
PG = spm_read_vols(ihdr);
[a,b,c,d] = size(PG);
clear V
for i=1:numel(ICAname)
    X = ICAmask(:,:,:,i);
    Xmask = imresize3(X,[a,b,c],'nearest');
    V(i,:) = nanmean(fmask(PG,Xmask),1); 
end

Reorder = [2 1 6 3 4 5 7:15];
F = draw_radar(V(Reorder,1),[-4 6]/2,[-0 0],ICAname(Reorder));
saveas(F,fullfile(dest,'Gradient1.tiff'))
F = draw_radar(V(Reorder,2),[-4 6]/4,[-0 0],ICAname(Reorder));
saveas(F,fullfile(dest,'Gradient2.tiff'))
F = draw_radar(V(Reorder,3),[-2 4]/4,[-0 0],ICAname(Reorder));
saveas(F,fullfile(dest,'Gradient3.tiff'))
F = draw_radar(V(Reorder,4),[-2 2]/2,[-0 0],ICAname(Reorder));
saveas(F,fullfile(dest,'Gradient4.tiff'))
F = draw_radar(V(Reorder,5),[-1 1],[-0 0],ICAname(Reorder));
saveas(F,fullfile(dest,'Gradient5.tiff'))
close all



%% structure function similarity
xhdr=spm_vol(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.nii.gz'));
Labels = spm_read_vols(xhdr);
AtlasTXT = loadtxt(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.txt'));

dest = fullfile(WholePath,'DiffusionGradient_Structure');
filename = fullfile(dest,'StructurePG.nii');
ihdr = spm_vol(filename);
PGs = spm_read_vols(ihdr);
dest = fullfile(WholePath,'DiffusionGradient_ION','Mean');
filename = fullfile(dest,'PG_align.nii');
ihdr = spm_vol(filename);
PGf_ION = spm_read_vols(ihdr);
dest = fullfile(WholePath,'DiffusionGradient_NIH','Mean');
filename = fullfile(dest,'PG_align.nii');
ihdr = spm_vol(filename);
PGf_NIH = spm_read_vols(ihdr);


clear Vs Vf*
for rg = 1:116
    RegionName = Abbre{rg+20};
    A = cellfun(@(x)strcmpi(x,RegionName),AtlasTXT(:,2));
    Bin = AtlasTXT{A,1};
    lmask = Labels==Bin;
    Vs(rg,:) = nanmean(fmask(PGs,lmask));
    Vf_NIH(rg,:) = nanmean(fmask(PGf_NIH,lmask));
    Vf_ION(rg,:) = nanmean(fmask(PGf_ION,lmask));
end

CC_ION = corr(Vs,Vf_ION);
CC_NIH = corr(Vs,Vf_NIH);

