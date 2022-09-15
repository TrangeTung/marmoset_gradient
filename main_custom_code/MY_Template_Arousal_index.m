%% preparation
% According to your experimental parameters and certain file direction, Change these bellow.
clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));
ExpRecording = [codepath 'ExpIndex.xlsx'];
[~,~,CellData] = xlsread(ExpRecording);
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
WholeMaskNII = fullfile(codepath,'marmoset_atlas','Trange_Template_mask_V3.nii');
lmask = spm_read_vols(spm_vol(WholeMaskNII));

ihdr = spm_vol(fullfile(WholePath,'Figure1','Arousal','convCC_Ovl.nii'));
ArousalTemplate = spm_read_vols(ihdr);
Vref = fmask(ArousalTemplate,lmask);

%
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
    
    for kk=GEEPI
        
        GEPath = fullfile(path,'Results',num2str(kk));
        filename = fullfile(GEPath,'snSRGrsmUw2dseq.nii');
        hdr = spm_vol(filename);
        fMRI_4D = spm_read_vols_4D(hdr);
        Vec = fmask(fMRI_4D,lmask);
        clear fMRI_4D
        
        V = (Vec-mean(Vec,2))./std(Vec,0,2);
        V(isnan(V))=0;
        
        AI = corr(V,Vref);
        a=a+1;
        if a==1;ArousalIndex = zeros(numel(AI),345);end
        ArousalIndex(:,a)=AI;
    end
    
end
Excel = fullfile(WholePath,'Figure1','MRI_based_Arousal_index.xlsx');
xlswrite(Excel,ArousalIndex,'Sheet1','A1');
