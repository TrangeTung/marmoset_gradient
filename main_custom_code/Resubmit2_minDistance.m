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

CAPpath = fullfile(WholePath,'Figure2');
CAP = spm_read_vols(spm_vol(fullfile(CAPpath,'CAP.nii')));
CAPv = fmask(CAP,lmask); CAPv(isnan(CAPv))=0;
Pupilpath = fullfile(WholePath,'Figure1','Arousal');
P = spm_read_vols(spm_vol(fullfile(Pupilpath,'convCC_Ovl.nii')));
Pv = fmask(P,lmask); Pv(isnan(Pv))=0;
%% dynamic FC estimation
% all VI
VIallM=zeros(512,345);
a_=0;
for idx = 01:62
    path = [ExpTable.Monkey{idx},filesep];
    Session = replace(path,WholePath,'');
    GEEPI = sort([  ExpTable.EPI4D1(idx),ExpTable.EPI4D2(idx),...
        ExpTable.EPI4D3(idx),ExpTable.EPI4D4(idx),...
        ExpTable.EPI4D5(idx),ExpTable.EPI4D6(idx),...
        ExpTable.EPI4D7(idx),ExpTable.EPI4D8(idx)],'ascend');
    GEEPI(isnan(GEEPI)) = [];
    for kk = GEEPI
        fprintf(['session: ',num2str(idx),';folder: ',num2str(kk),'\n']);
        GEPath = fullfile(WholePath,'DiffusionGradient_ION',Session,num2str(kk));
        load(fullfile(GEPath,'VigilanceIndex.mat'));
        a_=a_+1;VIallM(:,a_)=VI;
    end
end
bins=10;VIall = VIallM;
DMall=zeros(bins,bins,size(VIall,2));
for aloop=1:size(VIall,2)
    VI = VIall(:,aloop);
    clear locs*
    for bl=1:bins

        bs = [-3;-0.0885739075035637;-0.0568663768844682;-0.0353849678398447;...
          -0.0177808962429012;-0.00131054146237767;0.0150987043976391;...
          0.0332873622564407;0.0555794681660974;0.0899582511075418;3];
        bands = [bs(bl),bs(bl+1)];
        if bl==1;bands(1)=-3;end
        if bl==bins;bands(2)=3;end
        [locsTime{bl}] = find(VI>=bands(1)&VI<=bands(2));
    end
    
    DM = zeros(bins,bins);
    for xl=1:bins
        for yl=1:bins
            Dis = zeros(1,numel(locsTime{xl}));
            for zl=1:numel(locsTime{xl})
                Dis(zl)=min(abs(locsTime{xl}(zl)-locsTime{yl}));
            end
            DM(xl,yl)=mean(Dis(:));
        end
    end
    DMall(:,:,aloop)=DM;
end

%% 20 
Interp = 5;bins=10;VIall = VIallM;
Cut = 0;
NaNMtx = VIall*0;
NaNMtx(Cut:Interp:size(VIall,2)-Cut,:)=1;
VIall(~NaNMtx) = nan;
DMall_10=zeros(bins,bins,size(VIall,2));
for aloop=1:size(VIall,2)
    VI = VIall(:,aloop);
    if isempty(~isnan(VI));continue;end
    clear locs*
    for bl=1:bins

        bs = [-3;-0.0885739075035637;-0.0568663768844682;-0.0353849678398447;...
          -0.0177808962429012;-0.00131054146237767;0.0150987043976391;...
          0.0332873622564407;0.0555794681660974;0.0899582511075418;3];
        bands = [bs(bl),bs(bl+1)];
        if bl==1;bands(1)=-3;end
        if bl==bins;bands(2)=3;end
        ls = find(VI>=bands(1)&VI<=bands(2));
        if isempty(ls);ls=nan;fprintf([num2str(bl),'\n']);end
         [locsTime{bl}]=ls;
    end
    
    DM = zeros(bins,bins);
    for xl=1:bins
        for yl=1:bins
            Dis = zeros(1,numel(locsTime{xl}));
            for zl=1:numel(locsTime{xl})
                Dis(zl)=min(abs(locsTime{xl}(zl)-locsTime{yl}));
            end
            DM(xl,yl)=mean(Dis(:));
        end
    end
    DMall_10(:,:,aloop)=DM;
end

%% 40 
Interp=10;bins=10;VIall = VIallM;
Cut = 0;
NaNMtx = VIall*0;
NaNMtx(Cut+1:Interp:size(VIall,2)-Cut,:)=1;
VIall(~NaNMtx) = nan;
DMall_20=zeros(bins,bins,size(VIall,2));
for aloop=1:size(VIall,2)
    VI = VIall(:,aloop);
    if isempty(~isnan(VI));break;end
    clear locs*
    for bl=1:bins

        bs = [-3;-0.0885739075035637;-0.0568663768844682;-0.0353849678398447;...
          -0.0177808962429012;-0.00131054146237767;0.0150987043976391;...
          0.0332873622564407;0.0555794681660974;0.0899582511075418;3];
        bands = [bs(bl),bs(bl+1)];
        if bl==1;bands(1)=-3;end
        if bl==bins;bands(2)=3;end
        ls = find(VI>=bands(1)&VI<=bands(2));
        if isempty(ls);ls=nan;fprintf([num2str(bl),'\n']);end
         [locsTime{bl}]=ls;
    end
    
    DM = zeros(bins,bins);
    for xl=1:bins
        for yl=1:bins
            Dis = zeros(1,numel(locsTime{xl}));
            for zl=1:numel(locsTime{xl})
                Dis(zl)=min(abs(locsTime{xl}(zl)-locsTime{yl}));
            end
            DM(xl,yl)=mean(Dis(:));
        end
    end
    DMall_20(:,:,aloop)=DM;
end

F=  figure('Position', [680 190 990 798]);
MyMap=jet(256); %MyMap(1,:)=[100 100 100]/256;
subplot(2,2,1); X1 = nanmean(DMall,3);
X1 = X1/2+X1'/2;
imagesc(2*X1); caxis([0 5]*2);
colormap(MyMap);colorbar;meshgrid on;
axis equal
set(gca,'linewidth',2,'ticklength',[0.04,0.1],'fontsize',15);
subplot(2,2,3); X2 = nanmean(DMall_10,3);
X2 = X2/2+X2'/2;
imagesc(2*X2); caxis([0 20]*2);
colormap(MyMap);colorbar;meshgrid on;
axis equal
set(gca,'linewidth',2,'ticklength',[0.04,0.1],'fontsize',15);
subplot(2,2,4); X3 = nanmean(DMall_20,3);
X3 = X3/2+X3'/2;
imagesc(2*X3); caxis([0 40]*2);
colormap(MyMap);colorbar;meshgrid on;
axis equal
set(gca,'linewidth',2,'ticklength',[0.04,0.1],'fontsize',15);

