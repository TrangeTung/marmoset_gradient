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

genelist = fullfile(codepath,'Gene_list.xlsx');
receptor_ = {'histamine';'dopamine';'serotonin';'cholinergic';'adrenoceptor'};
gene_base = 'H:\Marmoset_Gene\MarmosetGeneHistology\';


AtlasExcel = fullfile(codepath,'marmoset_atlas','Trange_modifid_RikenBMA_atlas_labels.xlsx');
[~,~,CellData] = xlsread(AtlasExcel);
AtlasTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));



filename = fullfile(codepath,'MarmosetBrainConnectivity.xlsx');
FLNe_log = xlsread(filename,'log(FLNe)');
FLNe_log(isnan(FLNe_log))=-6;
SCcorr = corr(FLNe_log');
SCcorr=fillmissing(SCcorr,'linear');
SCcorr=fillmissing(SCcorr','linear');
SC = SCcorr;



GeM_raw = xlsread(fullfile(WholePath,'Figure5','gene_expression_correlation_matrix_edge_complex.xlsx'));
GeM = GeM_raw; 
    

GeM_Glutamate = xlsread(fullfile(WholePath,'Figure3','Glutamate_Gene_Expression.xlsx'));
GeM_Glutamate(isnan(GeM_Glutamate))=0;
GeM = corr(GeM_Glutamate(:,21:end));
GeM(isnan(GeM)) = 0;
dest = fullfile(WholePath,'DiffusionGradient_ION_Glutamate');


for idx = 1:62
%     idx
    path = [ExpTable.Monkey{idx},filesep];
    RARE =  ExpTable.RARE(idx);
    SEEPI = sort([ExpTable.SEEPI1(idx),ExpTable.SEEPI2(idx)],'ascend');
    GEEPI = sort([  ExpTable.EPI4D1(idx),ExpTable.EPI4D2(idx),...
        ExpTable.EPI4D3(idx),ExpTable.EPI4D4(idx),...
        ExpTable.EPI4D5(idx),ExpTable.EPI4D6(idx),...
        ExpTable.EPI4D7(idx),ExpTable.EPI4D8(idx)],'ascend');
    GEEPI(isnan(GEEPI)) = [];
    
    for kk = GEEPI
        
          
        fprintf(['Session ',num2str(idx),'; Folder ',num2str(kk),'\n']);

        GEPath = fullfile(path,'Functions','Seed-based',num2str(kk));
        V=xlsread(fullfile(GEPath,'TimeSeries_Riken_new.xlsx')); V(isnan(V))=0;
        V = (V-mean(V,2))./std(V,0,2); V(isnan(V))=0;
        
        Q = V(21:end,:);
        FC = corr(Q');
        FC(isnan(FC))=0;

        %% full model static
        DesignMatrix = [SC(:),GeM(:),ones([numel(GeM(:)),1])*1];
        lamda = 10:-1:-20;
        V = FC(:);
        DM = DesignMatrix;
        fa = 10^lamda(13);
        beta = pinv(DM'*DM+fa*eye(size(DM,2)))*(DM'*V);
        VecP = (beta'*DM')';
        predictedFC_full_static = reshape(VecP,size(FC));
        
        %% single variable static
        sloop = 1; VecP = (beta(sloop,:)'*DM(:,sloop)')';
        predictedFC_SC_static = reshape(VecP,size(FC));
        sloop = 2; VecP = (beta(sloop,:)'*DM(:,sloop)')';
        predictedFC_GeM_static = reshape(VecP,size(FC));
        
        
        %% shuffle static
        sloop = 1; [~,~,~,~,VecP] = MY_gpuArray_GLM_shuffle_static(V(:)',DesignMatrix,fa,sloop);
        predictedFC_SCshuffle_static = reshape(VecP,size(FC));
        sloop = 2; [~,~,~,~,VecP] = MY_gpuArray_GLM_shuffle_static(V(:)',DesignMatrix,fa,sloop);
        predictedFC_GeMshuffle_static = reshape(VecP,size(FC));

        
        %% full model Dynamic
        load(fullfile(GEPath,'DyCC_full.mat'));
        D  = DyCC_full;%(lottery,lottery,:);
        V = reshape(D,[],size(D,3));
        
        lamda = 10:-1:-20;
        DM = DesignMatrix;
        fa = 10^lamda(13);
        beta = pinv(DM'*DM+fa*eye(size(DM,2)))*(DM'*V);
        VecP = (beta'*DM')';
        predictedFC_full_dynamic = reshape(VecP,size(D));
        
        %% single variable Dynamic
        sloop = 1; VecP = (beta(sloop,:)'*DM(:,sloop)')';
        predictedFC_SC_dynamic = reshape(VecP,size(D));
        sloop = 2; VecP = (beta(sloop,:)'*DM(:,sloop)')';
        predictedFC_GeM_dynamic = reshape(VecP,size(D));
        
        
        %% shuffle Dynamic
        sloop = 1; [~,~,~,~,VecP] = MY_gpuArray_GLM_shuffle(V',DesignMatrix,fa,sloop);
        predictedFC_SCshuffle_dynamic = reshape(VecP,size(D));
        sloop = 2; [~,~,~,~,VecP] = MY_gpuArray_GLM_shuffle(V',DesignMatrix,fa,sloop);
        predictedFC_GeMshuffle_dynamic = reshape(VecP,size(D));
        
        %% static save
        dest_ = fullfile(dest,replace(ExpTable.Monkey{idx},WholePath,''),num2str(kk),'Model');
        if ~exist(dest_,'dir');mkdir(dest_);end
        csvwrite(fullfile(dest_,'predictedFC_full_static.csv'),predictedFC_full_static);
        csvwrite(fullfile(dest_,'predictedFC_SC_static.csv'),predictedFC_SC_static);
        csvwrite(fullfile(dest_,'predictedFC_GeM_static.csv'),predictedFC_GeM_static);
        csvwrite(fullfile(dest_,'predictedFC_SCshuffle_static.csv'),predictedFC_SCshuffle_static);
        csvwrite(fullfile(dest_,'predictedFC_GeMshuffle_static.csv'),predictedFC_GeMshuffle_static);
        %% dynamic save
        hdr = ihdr; hdr.dim=size(predictedFC_full_dynamic);
        hdr.fname = fullfile(dest_,'predictedFC_full_dynamic.nii'); 
        spm_write_vol(hdr,predictedFC_full_dynamic);
        hdr.fname = fullfile(dest_,'predictedFC_SC_dynamic.nii'); 
        spm_write_vol(hdr,predictedFC_SC_dynamic);
        hdr.fname = fullfile(dest_,'predictedFC_GeM_dynamic.nii'); 
        spm_write_vol(hdr,predictedFC_GeM_dynamic);
        hdr.fname = fullfile(dest_,'predictedFC_SCshuffle_dynamic.nii'); 
        spm_write_vol(hdr,predictedFC_SCshuffle_dynamic);
        hdr.fname = fullfile(dest_,'predictedFC_GeMshuffle_dynamic.nii'); 
        spm_write_vol(hdr,predictedFC_GeMshuffle_dynamic);
        
        clear predictedFC_* D V VecP
        
    end
end
    


