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

AtlasExcel = fullfile(codepath,'marmoset_atlas','Trange_modifid_RikenBMA_atlas_labels.xlsx');
[~,~,CellData] = xlsread(AtlasExcel);
AtlasTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));

predictedFC_full_static       = zeros([116 116 345]);
predictedFC_SC_static         = zeros([116 116 345]);
predictedFC_GeM_static        = zeros([116 116 345]);
predictedFC_SCshuffle_static  = zeros([116 116 345]);
predictedFC_GeMshuffle_static = zeros([116 116 345]);
a_=0;
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
        dest = fullfile(WholePath,'DiffusionGradient_ION');
        dest_ = fullfile(dest,replace(ExpTable.Monkey{idx},WholePath,''),num2str(kk),'Model');
        predictedFC_full_static       (:,:,a_)= csvread(fullfile(dest_,'predictedFC_full_static.csv'));
        predictedFC_SC_static         (:,:,a_)= csvread(fullfile(dest_,'predictedFC_SC_static.csv'));
        predictedFC_GeM_static        (:,:,a_)= csvread(fullfile(dest_,'predictedFC_GeM_static.csv'));
        predictedFC_SCshuffle_static  (:,:,a_)= csvread(fullfile(dest_,'predictedFC_SCshuffle_static.csv'));
        predictedFC_GeMshuffle_static (:,:,a_)= csvread(fullfile(dest_,'predictedFC_GeMshuffle_static.csv'));
    end
end

Type = {'full';'SC';'GeM';'SCshuffle';'GeMshuffle'};
for tl=[3 5]%1:numel(Type)
    
    eval(['Vec=predictedFC_',Type{tl},'_static;']);
    FCz = mean(Vec,3);
    if tl>=4;FCz = mean(predictedFC_full_static,3)-mean(Vec,3);end
    clear Vec
    
    
    G = GradientMaps('kernel','cs','approach','dm','n_components',5,'align','pa');
    G = G.fit(FCz,'Reference',Vf(:,1:5));
    ConG = G.aligned{1};
    ConG = fillmissing(ConG,'linear');

    dest = fullfile(WholePath,'DiffusionGradient_ION','Model',Type{tl});
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
    
    
    dest = fullfile(WholePath,'DiffusionGradient_ION','Model',Type{tl});
    %delete(fullfile(dest,'fx.func.gii'));
    bar = [-1 -0.0001 0.0001 1]*3;
    for k = 1:2
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
end