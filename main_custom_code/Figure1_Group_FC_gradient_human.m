% clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
addpath(genpath(codepath));
CortexMaskname = fullfile(WholePath,'DiffusionGradient_HCP',...
    'Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_1mm.nii');
ihdr=spm_vol(CortexMaskname);
Cmask = round(spm_read_vols(ihdr));

%% data package
%{
FILE_PATH = 'F:\HCP_3T\HCP_NII_9110991\';
LIST = loadtxt(fullfile(codepath,'sub_460.txt'));
%
for loop=4:4:numel(LIST)
    FDir = dir(fullfile(FILE_PATH,[num2str(LIST{loop}),'*']));
    fprintf(['ALL: 460; subject ',num2str(loop,'%03d'),'\n']);
    for floop=1:numel(FDir)
        cd(fullfile(FDir(floop).folder,FDir(floop).name));
        GZdir = dir('*.nii.gz');
        for gloop=1:numel(GZdir)
            try
                hdr = spm_vol(GZdir(gloop).name);
                img = spm_read_vols(hdr);
            catch
                break;
            end
            
            a=91;b=109;c=91;d=1200;
            S = zeros(1000,d);
            for ls = 1:1000
                lm = Cmask==ls;
                lmask = imresize3(lm,[a,b,c],'nearest');
                V = nanmean(fmask(img,lmask),1);
                S(ls,:)=V;
            end
            
            dest_ = fullfile(WholePath,'DiffusionGradient_HCP',FDir(floop).name,num2str(gloop));
            if ~exist(dest_,'dir');mkdir(dest_);end
            CSV = fullfile(dest_,'TimeSeries_1000ROIs.csv');
            csvwrite(CSV,S);
            
            clear I S V  img hdr
        end
        
    end
    
    
end
%}


%% HCP Group

a_=0;

filepath = fullfile(WholePath,'DiffusionGradient_HCP');
RESTdir = dir(fullfile(filepath,'*_REST*'));

for idx = 1:50%numel(RESTdir)
    idx
    path = fullfile(RESTdir(idx).folder,RESTdir(idx).name);

    for lo=1:10;if ~exist(fullfile(path,num2str(lo)),'dir');break;end;end
        SS = lo-1;


    for xl=1:SS

        a_ = a_+1;
        %% sort arousal phase 10 bin connectivity matrix
        GEPath = fullfile(path,num2str(xl));
        V_ = csvread(fullfile(path,num2str(xl),'TimeSeries_1000ROIs.csv'));
        data =(V_-mean(V_,2))./std(V_,0,2);
        data(isnan(data))=0;


        FC = corr(data');

        if a_==1;Vec=zeros([size(FC),1755]);end
        Vec(:,:,a_)=FC;

    end
end



Ref = xlsread(fullfile(WholePath,'DiffusionGradient_HCP','GradientTemplate.xlsx'));
Vec(:,:,a_+1:end)=[];
FCz = nanmean(Vec,3);

G = GradientMaps('kernel','cs','approach','dm','align','pa');
G = G.fit(FCz,'reference',Ref);

Gd = G.aligned{1};
dest = fullfile(filepath,'Mean');
Gd = (Gd-mean(Gd,1))./std(Gd,0,1);
MY_PG_to_NII_TIFF(Gd*10,dest,'FunctionPG',ihdr,codepath)
1;
%}



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
