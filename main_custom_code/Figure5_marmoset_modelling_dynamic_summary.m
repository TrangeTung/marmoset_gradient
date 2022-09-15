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

%% dynamic FC summary
%{
Type = {'full';'SC';'GeM';'SCshuffle';'GeMshuffle'};
for tl=1:numel(Type)

bins = 10;  Cut = 0;

for bl=1:bins

    bs = [-3;-0.0885739075035637;-0.0568663768844682;-0.0353849678398447;...
      -0.0177808962429012;-0.00131054146237767;0.0150987043976391;...
      0.0332873622564407;0.0555794681660974;0.0899582511075418;3];
    bands = [bs(bl),bs(bl+1)];
    
    
    a_=0;  Vec = [];
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
            dest_ = fullfile(dest,replace(ExpTable.Monkey{idx},WholePath,''),num2str(kk));
            load(fullfile(dest_,'VigilanceIndex.mat'));
            %load(fullfile(dest_,'rawdata.mat'));

            
            locs = find(VI>=bands(1)&VI<=bands(2));
            locs(locs<Cut | locs>numel(VI)-Cut+1)=[];
            
            dest_ = fullfile(dest,replace(ExpTable.Monkey{idx},WholePath,''),num2str(kk),'Model');
            filename = fullfile(dest_,['predictedFC_',Type{tl},'_dynamic.nii']);
            FC = spm_read_vols(spm_vol(filename));
            FC_ = FC(:,:,locs);
            
            Vec = cat(3,Vec,FC_);
            
            
        end
    end
    FCz = mean(Vec,3);
    dest = fullfile(WholePath,'DiffusionGradient_ION','DynamicModel',Type{tl},num2str(bl,'%02d'));
    if ~exist(dest,'dir');mkdir(dest);end
    csvwrite(fullfile(dest,'Connectome.csv'),FCz);
    dlmwrite(fullfile(dest,'number.txt'),size(Vec,3));
    clear Vec
end

end
%}

%% dynamic gradients
%
Type = {'full';'SC';'GeM';'SCshuffle';'GeMshuffle'};
for tl=1:numel(Type)
    
    bins = 10;  Cut = 16;
    
    for bl=1:bins
        
        dest = fullfile(WholePath,'DiffusionGradient_ION','DynamicModel',Type{tl},num2str(bl,'%02d'));
        data = csvread(fullfile(dest,'Connectome.csv'));
        FCz = data;
        if tl>=4
            dest_ = fullfile(WholePath,'DiffusionGradient_ION','DynamicModel',Type{1},num2str(bl,'%02d'));
            data_ = csvread(fullfile(dest_,'Connectome.csv'));
            FCz = data_-data;
            
        end
        
        G = GradientMaps('kernel','cs','approach','dm','n_components',5,'align','pa');
        G = G.fit(FCz,'Reference',Vf(:,1:5));
        ConG = G.aligned{1};
        
        dest = fullfile(WholePath,'DiffusionGradient_ION','DynamicModel',Type{tl},num2str(bl,'%02d'));
        if ~exist(dest,'dir');mkdir(dest);end
        save(fullfile(dest,'G.mat'),'G');
        
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
        
        
        
        dest = fullfile(WholePath,'DiffusionGradient_ION','DynamicModel',Type{tl},num2str(bl,'%02d'));
        delete(fullfile(dest,'fx.func.gii'));
        bar = [-1 -0.0001 0.0001 1]*3;
        for k = 1:5
            F = SurfaceMap(fullfile(dest,'PG.nii'),bar,codepath,k,'marmoset');
            TIFname = fullfile(dest,['Surface_PG_',num2str(k),'.tiff']);
            saveas(F,TIFname)
            close all
            X = imread(fullfile(dest,['Surface_PG_',num2str(k),'.tiff']));
            x_ = X(61:400,65:680,:);
            y_ = X(531:870,81:700,:);
            z = cat(2,x_,y_);
            TIFname = fullfile(dest,['Surface_PG_',num2str(k),'_cut0.tiff']);
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
            z = uint8(z);
            
            TIFname = fullfile(dest,['Surface_PG_',num2str(k),'_cut.tiff']);
            imwrite(z,TIFname);
        end
    end
end
%}