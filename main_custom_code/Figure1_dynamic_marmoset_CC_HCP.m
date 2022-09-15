clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
addpath(genpath(codepath));
CortexMaskname = fullfile(WholePath,'DiffusionGradient_HCP',...
    'Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_1mm.nii');
ihdr=spm_vol(CortexMaskname);
Cmask = round(spm_read_vols(ihdr));

ArousalTemplate =  fullfile(WholePath,'DiffusionGradient_HCP',...
    'HUMAN_arousal_template.nii');
ArousalMap = spm_read_vols(spm_vol(ArousalTemplate));
a=91;b=109;c=91;
Av = zeros(1000,1);
for ls = 1:1000
    lm = Cmask==ls;
    lmask = imresize3(lm,[a,b,c],'nearest');
    V = nanmean(fmask(ArousalMap,lmask),1);
    Av(ls,:)=V;
end

filepath = fullfile(WholePath,'DiffusionGradient_HCP');
RESTdir = dir(fullfile(filepath,'*_REST*'));

%{
NUM = 1000;
for idx = 1:numel(RESTdir)
    idx
    path = fullfile(RESTdir(idx).folder,RESTdir(idx).name);
    
    for lo=1:10;if ~exist(fullfile(path,num2str(lo)),'dir');break;end;end
        SS = lo-1;
    
    
    for xl=1:SS
        
        
        tic;
        V_ = csvread(fullfile(path,num2str(xl),'TimeSeries_1000ROIs.csv'));
        data =(V_-mean(V_,2))./std(V_,0,2); 
        data(isnan(data))=0;
        %{
        VI = corr(data,Av);
        z = hilbert(VI);
        radius = angle(z);
        dest = fullfile(path,num2str(xl));
        save(fullfile(dest,'VigilanceIndex.mat'),'VI','radius');        
        %}
        %% sort drowsiness phase 10 bin connectivity matrix
        GEPath = fullfile(path,num2str(xl));
        load(fullfile(GEPath,'VigilanceIndex.mat'));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        VI = VI(randperm(numel(VI)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bins=10;
        Cw = 20; Cs = 1;  Cut = Cw+2*3*Cs;
        [ CtDy ] = tapered_sliding_window(data',Cw, Cs);
        
        
        Ct=zeros([NUM,NUM,bins]);
        CtNUM = zeros(bins,1);
        for bl=1:bins
            bs = [-3;-0.174030372980496;-0.111451106809598;-0.0672818347172635;...
                -0.0306325839124537;0.00286190819950864;0.0360616667538984;...
                0.0711538182233653;0.112381012460658;0.170875221060087;3];
            bands = [bs(bl),bs(bl+1)];
            locs = find(VI>=bands(1)&VI<=bands(2));
            locs(locs<Cut | locs>numel(VI)-Cut+1)=[];
            
            RAM = CtDy(:,:,locs);
            Ct(:,:,bl)=nansum(RAM,3); 
            CtNUM(bl) = numel(locs);
        end
        
        Vol(1,1) = struct(  'fname',    fullfile(GEPath,'DCC_VIweighted_NULL.nii'),...
                            'dim',double(size(Ct)),'mat',      diag([1;1;1;1]),...
                            'n',        [1,1],   'pinfo',    [1;0;0],...
                            'descrip',  num2str(CtNUM'),      'dt',       [16 0]);
        spm_write_vol(Vol,Ct); 
        
        clear Ct CtNUM RAM CtDy V_ data
        
        
        toc;
    end
    
    
end

%}



%% Dynamic gradient estimationa=65;b=86;c=64;d=512; 
%
bins = 10;
for bl=1:bins
    a_=0;NUM=1000;
    Vec=zeros([NUM,NUM,1755]);
    wNUM=zeros(1,1755);
    
for idx = 1:numel(RESTdir)
    idx
    path = fullfile(RESTdir(idx).folder,RESTdir(idx).name);
    
    for lo=1:10;if ~exist(fullfile(path,num2str(lo)),'dir');break;end;end
        SS = lo-1;
    
    
    for xl=1:SS
        
        a_ = a_+1;
        %% sort arousal phase 10 bin connectivity matrix
        GEPath = fullfile(path,num2str(xl));
        SWname = fullfile(GEPath,'DCC_VIweighted_NULL.nii');
        hdr = spm_vol(SWname);
        FC = spm_slice_vol(hdr,spm_matrix([0 0 bl]),hdr.dim(1:2),0);
        x = str2num(hdr.descrip);
        swNUM = x(bl);
        
        Vec(:,:,a_)=FC;
        wNUM(a_)=swNUM;
        
    end
end


FCz = nansum(Vec,3)/sum(wNUM); 


Ref = xlsread(fullfile(WholePath,'DiffusionGradient_HCP','GradientTemplate.xlsx'));
Vec(:,:,a_+1:end)=[];

G = GradientMaps('kernel','cs','approach','dm','align','pa');
G = G.fit(FCz,'reference',Ref);

Gd = G.aligned{1};
Gd(:,2) = -Gd(:,2);
dest = fullfile(filepath,'Dynamic_VIweighted_NULL',num2str(bl,'%02d'));
if ~exist(dest,'dir');mkdir(dest);end
save(fullfile(dest,'G.mat'),'G');

% Gd = (Gd-mean(Gd,1))./std(Gd,0,1);
MY_PG_to_NII_TIFF(Gd*10,dest,'FunctionPG',ihdr,codepath)

    
end






function MY_PG_to_NII_TIFF(Ref_f,root_path,foldername,ref_hdr,codepath)

ihdr = ref_hdr;
Cmask = round(spm_read_vols(ihdr));

RAM=zeros([size(Cmask),10]);
for iii=1:size(Ref_f,2)
    X = zeros(size(Cmask));
    for jjj=1:size(Ref_f,1)
        X(Cmask==jjj+1)=Ref_f(jjj,iii);
    end
    X = smooth3(X,'box',15);
    X = X+flip(X,1);
    RAM(:,:,:,iii)=X;
end
dest = root_path;
dest_ = fullfile(dest,foldername);
if ~exist(dest_,'dir');mkdir(dest_);end
filename = fullfile(dest_,[foldername,'.nii']);

for k = 1:10
    Vol(k,1) = struct(  'fname',    filename,...
                        'dim',      ihdr.dim,    'mat',      ihdr.mat,...
                        'n',        [k,1],       'pinfo',    [1;0;0],...
                        'descrip',  '',          'dt',       [16 0]);
end
spm_write_vol_4D(Vol,RAM);


dest = root_path;
dest_ = fullfile(dest,foldername);
if ~exist(dest_,'dir');mkdir(dest_);end
filename = fullfile(dest_,[foldername,'.nii']);
delete(fullfile(dest_,'fx.func.gii'));
% MA = [20,2,10,2,1];

for k=1:5
    bar = [-1 -0.0001 0.0001 1]*10;%*MA(k);
    F = SurfaceMap(filename,bar,codepath,k,'human');
    TIFname = fullfile(dest_,['PG_align_',num2str(k),'.tiff']);
    saveas(F,TIFname)
    
    
    vx = imread(TIFname);
    a = vx(58:457,84:673,:);
    b = vx(528:927,84:673,:);
    z = cat(2,a,b);
    TIFFname = fullfile(dest_,['PG_align_',num2str(k),'_cut.tiff']);
    imwrite(uint8(z),TIFFname)
    
end
close all
end
