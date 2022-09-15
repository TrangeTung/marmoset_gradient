%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fistly
% Download code packages from the links ("download_URL.txt" in each folder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% codepath = '....\..\marmoset_gradient\';
% filepath = '....\..\examplefile\';
addpath(genpath(codepath));

%% structural gradients
% Variable "Gradients" --> structual gradients based on retrograde SC
% Region names (1~116 regions) were list in the first column of 
% "MarmosetBrainConnectivity.xlsx"
tic;
filename = fullfile(filepath,'MarmosetBrainConnectivity.xlsx');
FLNe_log = xlsread(filename,'log(FLNe)');
FLNe_log(isnan(FLNe_log))=-6;
FLNe_log(51,:)=[]; % exclude this region due to no 
G = GradientMaps('kernel','cs','approach','dm','n_components',5);
G = G.fit(FLNe_log);
Gd = G.gradients{1};
Gradients = zeros(116,5);
Gradients([1:50 52:116],:) = Gd;
toc;
%% GLM model and reduced model
% Results stored in the folder "/filepath../Model/"
% predicted static results saved as .csv file
% predicted dynamic results saved as .nii file
% Region names (1~116 regions) were list in the first column of 
% "MarmosetBrainConnectivity.xlsx"
% Results of reduced model was the difference between results of 
% full_GLM and shuffled GLM 
tic;
AtlasNII = fullfile(filepath,'Trange_atlas_RikenBMA_cortex.nii');
ihdr = spm_vol(AtlasNII);
filename = fullfile(filepath,'MarmosetBrainConnectivity.xlsx');
FLNe_log = xlsread(filename,'log(FLNe)');
FLNe_log(isnan(FLNe_log))=-6;
SCcorr = corr(FLNe_log');
SCcorr=fillmissing(SCcorr,'linear');
SCcorr=fillmissing(SCcorr','linear');
SC = SCcorr;
GeM_raw = xlsread(fullfile(filepath,'gene_expression_correlation_matrix_edge_complex.xlsx'));
GeM = GeM_raw; 
V=xlsread(fullfile(filepath,'TimeSeries_Riken_new.xlsx')); V(isnan(V))=0;
V = (V-mean(V,2))./std(V,0,2); V(isnan(V))=0;
Q = V(21:end,:);
FC = corr(Q');
FC(isnan(FC))=0;
% full model static
DesignMatrix = [SC(:),GeM(:),ones([numel(GeM(:)),1])*1];
lamda = 10:-1:-20;
V = FC(:);
DM = DesignMatrix;
fa = 10^lamda(13);
beta = pinv(DM'*DM+fa*eye(size(DM,2)))*(DM'*V);
VecP = (beta'*DM')';
predictedFC_full_static = reshape(VecP,size(FC));
% single variable static
sloop = 1; VecP = (beta(sloop,:)'*DM(:,sloop)')';
predictedFC_SC_static = reshape(VecP,size(FC));
sloop = 2; VecP = (beta(sloop,:)'*DM(:,sloop)')';
predictedFC_GeM_static = reshape(VecP,size(FC));
% shuffle static
sloop = 1; [~,~,~,~,VecP] = MY_gpuArray_GLM_shuffle_static(V(:)',DesignMatrix,fa,sloop);
predictedFC_SCshuffle_static = reshape(VecP,size(FC));
sloop = 2; [~,~,~,~,VecP] = MY_gpuArray_GLM_shuffle_static(V(:)',DesignMatrix,fa,sloop);
predictedFC_GeMshuffle_static = reshape(VecP,size(FC));
% full model Dynamic
load(fullfile(filepath,'DyCC_full.mat'));
D  = DyCC_full;
V = reshape(D,[],size(D,3));
lamda = 10:-1:-20;
DM = DesignMatrix;
fa = 10^lamda(13);
beta = pinv(DM'*DM+fa*eye(size(DM,2)))*(DM'*V);
VecP = (beta'*DM')';
predictedFC_full_dynamic = reshape(VecP,size(D));
% single variable Dynamic
sloop = 1; VecP = (beta(sloop,:)'*DM(:,sloop)')';
predictedFC_SC_dynamic = reshape(VecP,size(D));
sloop = 2; VecP = (beta(sloop,:)'*DM(:,sloop)')';
predictedFC_GeM_dynamic = reshape(VecP,size(D));
% shuffle Dynamic
sloop = 1; [~,~,~,~,VecP] = MY_gpuArray_GLM_shuffle(V',DesignMatrix,fa,sloop);
predictedFC_SCshuffle_dynamic = reshape(VecP,size(D));
sloop = 2; [~,~,~,~,VecP] = MY_gpuArray_GLM_shuffle(V',DesignMatrix,fa,sloop);
predictedFC_GeMshuffle_dynamic = reshape(VecP,size(D));
% static save
dest_ = fullfile(filepath,'Model');
if ~exist(dest_,'dir');mkdir(dest_);end
csvwrite(fullfile(dest_,'predictedFC_full_static.csv'),predictedFC_full_static);
csvwrite(fullfile(dest_,'predictedFC_SC_static.csv'),predictedFC_SC_static);
csvwrite(fullfile(dest_,'predictedFC_GeM_static.csv'),predictedFC_GeM_static);
csvwrite(fullfile(dest_,'predictedFC_SCshuffle_static.csv'),predictedFC_SCshuffle_static);
csvwrite(fullfile(dest_,'predictedFC_GeMshuffle_static.csv'),predictedFC_GeMshuffle_static);
% dynamic save
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
toc;