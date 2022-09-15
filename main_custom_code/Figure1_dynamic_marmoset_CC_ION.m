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
a
Pupilpath = fullfile(WholePath,'Figure1','Arousal');
P = spm_read_vols(spm_vol(fullfile(Pupilpath,'convCC_Ovl.nii')));
Pv = fmask(P,lmask); Pv(isnan(Pv))=0;
%% dynamic FC estimation
%
NUM=5522;
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
        
        GEPath = fullfile(WholePath,'DiffusionGradient_ION',Session,num2str(kk));
        load(fullfile(GEPath,'rawdata.mat'));
        load(fullfile(GEPath,'VigilanceIndex.mat'));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        VI = VI(randperm(numel(VI)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        V = rawdata.data;
        NANidx = rawdata.nanidx;
        
        RegBasFunc = [ones(numel(VI),1),VI];
        Vr = V*0;
        for iii = 1:size(V,1)
            [Beta,~,Residual] = regress(double(V(iii,:))',RegBasFunc);
            Vr(iii,:) = Residual + Beta(1);
        end
        
        V_ = zeros(NUM,512);        V_(~NANidx,:)=V;
        Vr_ = zeros(NUM,512);       Vr_(~NANidx,:)=Vr;
        clear rawdata V Vr fMRI_4D
        
        %% sort drowsiness phase 10 bin connectivity matrix
        
        bins = 10; 
        DCCa = DCCsimple(V_);
        DCCb = DCCsimple(Vr_);

                     
        Ct=zeros([NUM,NUM,bins]);
        Ctr=zeros([NUM,NUM,bins]);
        CtNUM=zeros(bins,1);
        for bl=1:bins
           
            bs = [-3;-0.0885739075035637;-0.0568663768844682;-0.0353849678398447;...
                -0.0177808962429012;-0.00131054146237767;0.0150987043976391;...
                0.0332873622564407;0.0555794681660974;0.0899582511075418;3];
            bands = [bs(bl),bs(bl+1)];
            locs = find(VI>=bands(1)&VI<=bands(2));
            
           RAM = DCCa(:,:,locs);
           RAMr = DCCb(:,:,locs);
           
           Ct(:,:,bl) = nansum(RAM,3); 
           Ctr(:,:,bl) = nansum(RAMr,3); 
           clear C S RAM W RAMr
           CtNUM(bl) = numel(locs);
           
        end
        clear DCCa DCCb
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Vol(1,1) = struct(  'fname',    fullfile(GEPath,'SW15_VIweighted.nii'),...
                            'dim',double(size(Ct)),'mat',      diag([1;1;1;1]),...
                            'n',        [1,1],   'pinfo',    [1;0;0],...
                            'descrip',  num2str(CtNUM'),      'dt',       [16 0]);
        spm_write_vol(Vol,Ct); 
        Vol(1,1) = struct(  'fname',    fullfile(GEPath,'SW15_VIremoval.nii'),...
                            'dim',double(size(Ctr)),'mat',      diag([1;1;1;1]),...
                            'n',        [1,1],   'pinfo',    [1;0;0],...
                            'descrip',  num2str(CtNUM'),      'dt',       [16 0]);
        spm_write_vol(Vol,Ctr); 
        clear Ct Ctr  Vol V_ Vr_
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
        
end
%}

%% Dynamic gradient estimationa=65;b=86;c=64;d=512; 
AtlasNII = fullfile(codepath,'marmoset_atlas','Trange_atlas_RikenBMA_cortex.nii');
AtlasMask = spm_read_vols_4D(spm_vol(AtlasNII));
%


a=65;b=86;c=64;d=512; % 4D-EPI size
lmask = AtlasMask>20;
Xmask = imresize3(lmask,round([a,b,c]/2),'nearest');
NUM=5522;
%
bins = 10;
for bl=1:bins
    a_=0;NUM=5522;
    Vec=zeros([NUM,NUM,345]);
    wNUM=zeros(1,345);
    
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
        
        a_=a_+1;
        GEPath = fullfile(WholePath,'DiffusionGradient_ION',Session,num2str(kk));
        SWname = fullfile(GEPath,'SW15_VIweighted.nii');
        hdr = spm_vol(SWname);
        FC = spm_slice_vol(hdr,spm_matrix([0 0 bl]),hdr.dim(1:2),0);
        x = str2num(hdr.descrip);
        swNUM = x(bl);
        
        Vec(:,:,a_)=FC;
        wNUM(a_)=swNUM;
        
    end
end

filepath = fullfile(WholePath,'DiffusionGradient_ION','Mean');
load(fullfile(filepath,'Gradients.mat'));
G_Ref = Gd;



FCz = nansum(Vec,3)/sum(wNUM); 
clear Vec
NANidx = isnan(nanmean(FCz,1));
G = GradientMaps('kernel','cs','approach','dm','n_components',20,'align','pa');
data = FCz(~NANidx,~NANidx);
G = G.fit(data,'reference',G_Ref(~NANidx,:));

Gd = G.aligned{1}*50;
Gx = zeros(numel(NANidx),size(Gd,2));
Gx(~NANidx,:) = Gd;
Gd=Gx;


M = funmask(Gd,Xmask);
xhdr=spm_vol(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.nii.gz'));
M_ = zeros([xhdr.dim,size(M,4)]);
for m=1:size(M,4)
    X = M(:,:,:,m);
    X = imresize3(X,xhdr.dim,'linear');
    X = smooth3(X,'box',7);
    M_(:,:,:,m) = X;
end

dest_ = fullfile(WholePath,'DiffusionGradient_ION','Dynamic_VIweighted',num2str(bl,'%02d'));
if ~exist(dest_,'dir');mkdir(dest_);end
save(fullfile(dest_,'G.mat'),'G');


filename = fullfile(dest_,'PG_align.nii');
clear Vol
for k = 1:size(M,4)
    Vol(k,1) = struct(  'fname',    filename,...
                        'dim',      xhdr.dim,...
                        'mat',      xhdr.mat,...
                        'n',        [k,1],...
                        'pinfo',    [1;0;0],...
                        'descrip',  '',...
                        'dt',       [16 0]);
end
spm_write_vol_4D(Vol,M_); 


bar = [-1 -0.0001 0.0001 1]*5;
for k = 1:5
    F = SurfaceMap(fullfile(dest_,'PG_align.nii'),bar/k,codepath,k,'marmoset');
    TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'.tiff']);
    saveas(F,TIFname)
    close all
    X = imread(fullfile(dest_,['Surface_PG_align_',num2str(k),'.tiff']));
    x_ = X(61:400,65:680,:);
    y_ = X(531:870,81:700,:);
    z = cat(2,x_,y_);
    TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'_cut.tiff']);
    imwrite(z,TIFname);
    pause(3);
end

    
end

%}

