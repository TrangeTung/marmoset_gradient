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
GeM0 = GeM_raw;%(lottery,lottery);




for geneidx = 1:30
    
    [GeneExpression_raw,~,CellData] = xlsread(fullfile(WholePath,'Figure5','gene_expression.xlsx'));
    GeneExpression_raw(isnan(GeneExpression_raw))=0;
    GeneName = CellData(2:end,1);
    % load Surrogate gene expression
    SurrTxt = fullfile(WholePath,'Figure5','Surrogate',['Surrogate_',GeneName{geneidx},'.txt']);
    SurrGeneExp = load(SurrTxt);
    
    for idx = 1:62
        
        fprintf(['Gene:',num2str(geneidx),';Session:',num2str(idx),'\n']);
        
        path = [ExpTable.Monkey{idx},filesep];
        RARE =  ExpTable.RARE(idx);
        SEEPI = sort([ExpTable.SEEPI1(idx),ExpTable.SEEPI2(idx)],'ascend');
        GEEPI = sort([  ExpTable.EPI4D1(idx),ExpTable.EPI4D2(idx),...
            ExpTable.EPI4D3(idx),ExpTable.EPI4D4(idx),...
            ExpTable.EPI4D5(idx),ExpTable.EPI4D6(idx),...
            ExpTable.EPI4D7(idx),ExpTable.EPI4D8(idx)],'ascend');
        GEEPI(isnan(GEEPI)) = [];
        
        for kk = GEEPI
            
            parfor shf = 1:1000
                %SHFlist = randperm(116);
                %GeneExpression = GeneExpression_raw;
                %GeneExpression(geneidx,:) = GeneExpression(geneidx,SHFlist);
                GeneExpression = GeneExpression_raw;
                GeneExpression(geneidx,:) = SurrGeneExp(:,shf);
                
                GEPath = fullfile(path,'Functions','Seed-based',num2str(kk));
                lamda = 10:-1:-20;
                X = load(fullfile(GEPath,'DyCC_full.mat'));
                D  = X.DyCC_full;%(lottery,lottery,:);
                V = reshape(D,[],size(D,3));
                
                
                
                
                %% full model Dynamic raw
                GeM = GeM0;
                DesignMatrix = [SC(:),GeM(:),ones([numel(GeM(:)),1])*1];
                DM = DesignMatrix;
                fa = 10^lamda(13);
                beta = pinv(DM'*DM+fa*eye(size(DM,2)))*(DM'*V);
                VecP = (beta'*DM')';
                predictedFC_full_dynamic_raw = reshape(VecP,size(D));
                
                
                %% full model Dynamic shuffle
                GeM = corr(GeneExpression);
                DesignMatrix = [SC(:),GeM(:),ones([numel(GeM(:)),1])*1];
                DM = DesignMatrix;
                fa = 10^lamda(13);
                beta = pinv(DM'*DM+fa*eye(size(DM,2)))*(DM'*V);
                VecP = (beta'*DM')';
                predictedFC_full_dynamic_shuffle = reshape(VecP,size(D));
                
                dest = fullfile(WholePath,'DiffusionGradient_ION','Shuffle_ReSubmit',...
                    ['Gene_',num2str(geneidx,'%02d')],['Shuffle_',num2str(shf,'%04d')]);
                dest_ = fullfile(dest,replace(ExpTable.Monkey{idx},WholePath,''),num2str(kk));
                if ~exist(dest_,'dir');mkdir(dest_);end
                
                predictedFC_full_dynamic_diff = predictedFC_full_dynamic_raw-predictedFC_full_dynamic_shuffle;
                
                
                
                
                
                destx = fullfile(WholePath,'DiffusionGradient_ION');
                dest_x = fullfile(destx,replace(ExpTable.Monkey{idx},WholePath,''),num2str(kk));
                Y=load(fullfile(dest_x,'VigilanceIndex.mat'));
                VI=Y.VI;
                
                bins = 10;  Cut = 0;
                Vec1=zeros(116,116,10);
                for bl=1:bins
                    bs = [-3;-0.0885739075035637;-0.0568663768844682;-0.0353849678398447;...
                    -0.0177808962429012;-0.00131054146237767;0.0150987043976391;...
                    0.0332873622564407;0.0555794681660974;0.0899582511075418;3];
                    bands = [bs(bl),bs(bl+1)];
                    
                    locs = find(VI>=bands(1)&VI<=bands(2));
                    locs(locs<Cut | locs>numel(VI)-Cut+1)=[];
                    
                    FC1=predictedFC_full_dynamic_diff;
                    FC1_ = FC1(:,:,locs);
                    Vec1(:,:,bl) = mean(FC1_,3);
                    
                end
                
                hdr = ihdr; hdr.dim=size(Vec1);
                hdr.fname = fullfile(dest_,'Connectome_difference.nii');
                spm_write_vol(hdr,Vec1); Vec1=[];
                
                
                predictedFC_full_dynamic_diff=[];
                predictedFC_full_dynamic_raw=[];
                predictedFC_full_dynamic_shuffle=[];
                FC1=[];FC1_=[];
                
                Excel = fullfile(dest_,'GeneExpression_shuffle.csv');
                csvwrite(Excel,GeneExpression);
            end
        end
        
        
        
    end
end