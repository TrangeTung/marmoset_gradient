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

X=[];

Q = numel(y)/10;
Seg = y(round(Q*(1:9)));

NUM=5522;HM1=cell(10,1);
for idx = 1:62
    
    path = [ExpTable.Monkey{idx},filesep];
    Session = replace(path,WholePath,'');
    RARE =  ExpTable.RARE(idx);
    SEEPI = sort([ExpTable.SEEPI1(idx),ExpTable.SEEPI2(idx)],'ascend');
    GEEPI = sort([  ExpTable.EPI4D1(idx),ExpTable.EPI4D2(idx),...
        ExpTable.EPI4D3(idx),ExpTable.EPI4D4(idx),...
        ExpTable.EPI4D5(idx),ExpTable.EPI4D6(idx),...
        ExpTable.EPI4D7(idx),ExpTable.EPI4D8(idx)],'ascend');
    GEEPI(isnan(GEEPI)) = [];
    
    
    for kk = GEEPI
        
        fprintf(['session: ',num2str(idx),';folder: ',num2str(kk),'\n']);
        
        GEPath = fullfile(path,'Results',num2str(kk));
        rp = load(fullfile(GEPath,'rp_smUw2dseq.txt'));
        rp(:,1:3) = rp(:,1:3)/6; rp(:,4:6) = rp(:,4:6)*16;
        rp_ = rp - [rp(1,:);rp(1:end-1,:)];
        FD = sqrt(rp_(:,1).^2+rp_(:,2).^2+rp_(:,3).^2+...
            rp_(:,4).^2+rp_(:,5).^2+rp_(:,6).^2)*1000; % um
        
        GEPath = fullfile(WholePath,'DiffusionGradient_ION',Session,num2str(kk));
        load(fullfile(GEPath,'VigilanceIndex.mat'));
        
        bins = 10;  
        X = cat(1,X,VI(:)');
        
        for bl=1:bins
            bs = [-3;-0.0885739075035637;-0.0568663768844682;-0.0353849678398447;...
                -0.0177808962429012;-0.00131054146237767;0.0150987043976391;...
                0.0332873622564407;0.0555794681660974;0.0899582511075418;3];
            bands = [bs(bl),bs(bl+1)];
            locs = find(VI>=bands(1)&VI<=bands(2));
           
            HM1{bl} = cat(1,HM1{bl},FD(locs));
        end
        
    end
end

NIH_path = 'E:\Marmoset_NIH_dataset';
Session = dir(fullfile(NIH_path,'*_session_*'));
HM2=cell(10,1);
for idx = 1:39
    
    
    path = fullfile(NIH_path,Session(idx).name);
    for ii=1:99;if ~exist(fullfile(path,'Results',num2str(ii)),'dir');break;end;end
    ScanNum = 1:ii-1;
    
    for nl = ScanNum
        
        fprintf(['marmoset: ',num2str(idx),';scan: ',num2str(nl),'\n'])
        GEPath = fullfile(path,'Results',num2str(nl));
        rp = load(fullfile(GEPath,'rp_smUw2dseq.txt'));
        rp(:,1:3) = rp(:,1:3)/6; rp(:,4:6) = rp(:,4:6)*16;
        rp_ = rp - [rp(1,:);rp(1:end-1,:)];
        FD = sqrt(rp_(:,1).^2+rp_(:,2).^2+rp_(:,3).^2+...
            rp_(:,4).^2+rp_(:,5).^2+rp_(:,6).^2)*1000; % um

        GEPath = fullfile(WholePath,'DiffusionGradient_NIH',Session(idx).name,num2str(nl,'%02d'));
        load(fullfile(GEPath,'VigilanceIndex.mat'));
        

        bins = 10;  
          X = cat(1,X,VI(:)');    
        
        for bl=1:bins
            bs = [-3;-0.0885739075035637;-0.0568663768844682;-0.0353849678398447;...
                -0.0177808962429012;-0.00131054146237767;0.0150987043976391;...
                0.0332873622564407;0.0555794681660974;0.0899582511075418;3];
            bands = [bs(bl),bs(bl+1)];
            locs = find(VI>=bands(1)&VI<=bands(2));
            
            HM2{bl} = cat(1,HM2{bl},FD(locs));
        end        
    end
end


QW=nan(25000*2,10);
for ix=1:10
NUM(ix,1) = numel(HM1{ix});
NUM(ix,2) = numel(HM2{ix});
end



RegBasFunc = [ones(10,1),NUM'];
for iii=1:4
[Beta,~,Residual] = regress(double(nl(:,iii)),RegBasFunc);
data_ready_regress(:,iii) = Residual + Beta(1);

end
