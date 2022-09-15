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
atlas_label = AtlasTable.Labels;
Abbre = AtlasTable.Abbre;

%{
filename = fullfile(codepath,'marmoset_connectivity_atlas',...
    'marmoset_brain_connectivity_1_0_fln_matrix.txt');
X = importdata(filename);
target = unique(X.textdata(:,1));
Abbre = AtlasTable.Abbre(21:end);
atlas_label = AtlasTable.Labels;
for  rloop = 1:numel(target)
    A = cellfun(@(x)strcmp(lower(x),lower(target{rloop})),Abbre);
    lottery(rloop) = find(A==true);
end

dest = fullfile(WholePath,'Figure5_1');
SC_ = xlsread(fullfile(dest,'FLNe_matrix_edge_complex.xlsx'));
SCl = log10(SC_); SCl(isnan(SCl))=1;SCl(SCl==-Inf)=-7;
SC = corr(SCl);
%SC = (SCl+7)/7;
%}


filename = fullfile(codepath,'MarmosetBrainConnectivity.xlsx');
FLNe_log = xlsread(filename,'log(FLNe)');
FLNe_log(isnan(FLNe_log))=-6;
SCcorr = corr(FLNe_log');
SCcorr=fillmissing(SCcorr,'linear');
SCcorr=fillmissing(SCcorr','linear');
SC = SCcorr;


GeM_raw = xlsread(fullfile(WholePath,'Figure5','gene_expression_correlation_matrix_edge_complex.xlsx'));
GeM = GeM_raw;%(lottery,lottery);

%% Dynamic CC region-wise
a=0;
%{
NUM=5522;
NIH_path = 'E:\Marmoset_NIH_dataset';
Session = dir(fullfile(NIH_path,'*_session_*'));
for idx = 2:2:39
    
    
    path = fullfile(NIH_path,Session(idx).name);
    for ii=1:99;if ~exist(fullfile(path,'Results',num2str(ii)),'dir');break;end;end
    ScanNum = 1:ii-1;
    
    for nl = ScanNum
        a=a+1;
        fprintf(['marmoset: ',num2str(idx),';scan: ',num2str(nl),'\n'])
        
        GEPath = fullfile(path,'Results',num2str(nl));
        hdr = spm_vol(fullfile(GEPath,'sncRrsmUw2dseq.nii'));
        fMRI_4D = spm_read_vols_4D(hdr); 
        fMRI_4D(isnan(fMRI_4D))=0;
        fMRI_4D = (fMRI_4D+flip(fMRI_4D,1))/2;
        
        V = zeros(116,512);
        for ix=1:116
           ilmask = AtlasMask==(AtlasTable.Labels(ix+20));
           X = mean(fmask(fMRI_4D,ilmask),1);
           V(ix,:) = (X-mean(X))/std(X);
        end
        clear fMRI_4D
        Excel = fullfile(GEPath,'TimeSeries_Riken_new.xlsx');
        xlswrite(Excel,V,'Sheet1','A1');
        
        GEPath = fullfile(path,'Results',num2str(nl));
        V = xlsread(fullfile(GEPath,'TimeSeries_Riken_new.xlsx'));
        
        dest = fullfile(WholePath,'DiffusionGradient_NIH');
        dest_ = fullfile(dest,Session(idx).name,num2str(nl,'%02d'));
        load(fullfile(dest_,'VigilanceIndex.mat'));
        AI = VI;
        data_ = V(1:end,:);
        NANidx = isnan(mean(data_,2));
        data = data_(~NANidx,:);
        
        data_reg = data*0;
        Reg = [ones(size(AI,1),1),AI];
        for i=1:size(data,1)
            [beta,~,residual] = regress(data(i,:)',Reg);
            data_reg(i,:) = beta(1)+residual;
        end
        
        %% sliding window tapered 15
        [ Ct ] = tapered_sliding_window(data',15, 1);
        Ct_15 = zeros([size(data_,1),size(data_,1),size(data_,2)]);
        Ct_15(~NANidx,~NANidx,:) = Ct;
        %% sliding window tapered 30
        [ Ct ] = tapered_sliding_window(data',30, 1);
        Ct_30 = zeros([size(data_,1),size(data_,1),size(data_,2)]);
        Ct_30(~NANidx,~NANidx,:) = Ct;
        %% sliding window tapered 15 + Arousal Regress
        [ Ct ] = tapered_sliding_window(data_reg',15, 1);
        Ct_reg_15 = zeros([size(data_,1),size(data_,1),size(data_,2)]);
        Ct_reg_15(~NANidx,~NANidx,:) = Ct;
        %% sliding window tapered 30 + Arousal Regress
        [ Ct ] = tapered_sliding_window(data_reg',30, 1);
        Ct_reg_30 = zeros([size(data_,1),size(data_,1),size(data_,2)]);
        Ct_reg_30(~NANidx,~NANidx,:) = Ct;
        %% DCC
        DyCC_full = DCCsimple(data');
        %% DCC + Arousal Regress
        DyCC = DCCsimple(data_reg');
        DyCC_reg_full = zeros([size(data_reg,1),size(data_reg,1),size(data_reg,2)]);
        DyCC_reg_full(~NANidx,~NANidx,:) = DyCC;
        
        cd(dest_);
        save('DyCC.mat','DyCC_full','DyCC_reg_full','Ct_15','Ct_30','Ct_reg_15','Ct_reg_30');
        clear *reg* Ct_*  Dy*      
        
        
        
    end
end
%}

%% model
%
a=0;
NUM=5522;
NIH_path = 'E:\Marmoset_NIH_dataset';
Session = dir(fullfile(NIH_path,'*_session_*'));
for idx = 1:39
    
    
    path = fullfile(NIH_path,Session(idx).name);
    for ii=1:99;if ~exist(fullfile(path,'Results',num2str(ii)),'dir');break;end;end
    ScanNum = 1:ii-1;
    
    for nl = ScanNum
        a=a+1;
        fprintf(['marmoset: ',num2str(idx),';scan: ',num2str(nl),'\n'])
        
        GEPath = fullfile(path,'Results',num2str(nl));
        V=xlsread(fullfile(GEPath,'TimeSeries_Riken_new.xlsx')); V(isnan(V))=0;
        V = (V-mean(V,2))./std(V,0,2); V(isnan(V))=0;
        
        Q = V(1:end,:);
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
        dest = fullfile(WholePath,'DiffusionGradient_NIH');
        dest_ = fullfile(dest,Session(idx).name,num2str(nl,'%02d'));
        load(fullfile(dest_,'DyCC.mat'));
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
        dest = fullfile(WholePath,'DiffusionGradient_NIH');
        dest_ = fullfile(dest,Session(idx).name,num2str(nl,'%02d'),'Model');
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


%}
