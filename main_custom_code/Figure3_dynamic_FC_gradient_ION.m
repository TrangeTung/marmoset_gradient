clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
cd(codepath); 
addpath(genpath(codepath));
ExpRecording = [codepath 'ExpIndex.xlsx'];
[~,~,CellData] = xlsread(ExpRecording);
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));

AtlasExcel = fullfile(codepath,'marmoset_atlas','Trange_modifid_RikenBMA_atlas_labels.xlsx');
[~,~,CellData] = xlsread(AtlasExcel);
AtlasTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
atlas_label = AtlasTable.Labels;
Abbre = AtlasTable.Abbre;


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
PGf = spm_read_vols(ihdr);

clear Vs Vf
for rg = 1:116
    RegionName = Abbre{rg+20};
    A = cellfun(@(x)strcmpi(x,RegionName),AtlasTXT(:,2));
    Bin = AtlasTXT{A,1};
    lmask = Labels==Bin;
    Vf(rg,:) = nanmean(fmask(PGs,lmask));
end

a=0;
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
        
        a=a+1;
        GEPath = fullfile(path,'Functions','Seed-based',num2str(kk));
        V=xlsread(fullfile(GEPath,'TimeSeries_Riken_new.xlsx')); V(isnan(V))=0;
        V = (V-mean(V,2))./std(V,0,2); V(isnan(V))=0;
        FC = corr(V(21:end,:)');
        
        
        FCx(:,:,a)=FC;
    end
    

end

FCz = nanmedian(FCx,3);
NANidx = (std(FCz,1)==0);
G = GradientMaps('kernel','cs','approach','dm');
G = G.fit(FCz(~NANidx,~NANidx));

Gd = G.gradients{1}*50;
Gx = zeros(116,size(Gd,2));
Gx(~NANidx,:) = Gd;
Gd = fillmissing(Gx,'pchip');
ConG = Gd;

dest = fullfile(WholePath,'DiffusionGradient_ION','Dynamic');
if ~exist(dest,'dir');mkdir(dest);end
% save(fullfile(dest,'G.mat'),'G');

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

delete(fullfile(dest,'fx.func.gii'));
bar = [-2 -0.0001 0.0001 2];
for k = 1:5
    F = SurfaceMap(fullfile(dest,'PG.nii'),bar,codepath,k,'marmoset');
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