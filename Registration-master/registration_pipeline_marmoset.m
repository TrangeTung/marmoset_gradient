function MY_register_Histology_Marmoset(genename,atlas_file,WholePath)



addpath(genpath('D:\Registration-master'))

% inputs
atlas_file = 'D:\GeneAtlas\S0.nii';
genename='HRH3';

data_base_dir = 'H:\Marmoset_Gene\MarmosetGeneHistology';
files_ = dir(fullfile(data_base_dir,[genename,'*']));
file_dir = files_.name;
% go back to case 1, tro to enable edit mode
input_dir = 'D:\GeneAtlas\INPUT1\';
output_dir = 'D:\GeneAtlas\OUTPUT\';
detailed_output_dir = output_dir;

dx=25;dy=25;dz=200;
create_location_csv_Marmoset(input_dir, genename, dx,dy,dz) % μm

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1 is to create an initial slice alignment
% we use a simple initialization where we locate the center of mass of each
% slice and align by translation
pattern = ['*',genename,'*.jpg'];
find_centers_for_initialization_nissl(input_dir, pattern, detailed_output_dir);
close all;

r = 10; % match to this many neighbors
% this downsampling and iterations is enough for a good initial guess
% not enough for a full accurate reconstruction
downs = [4 2];
niter = 40;
skip_thick = -1; % no thick slices to be skipped
load_initializer = 1;
e = 0.1;
atlas_free_rigid_alignment(input_dir, pattern, detailed_output_dir, r, downs, niter, e, skip_thick, load_initializer)
close all;

load([output_dir,'initializer_A.mat']);

Nissl_dir = dir([input_dir ['*',genename,'*.jpg']]);
for sl = 1:numel(Nissl_dir)
    Is = rgb2gray(imread(fullfile(input_dir,Nissl_dir(sl).name)));
    tform=affine2d(AJ(:,:,sl)');
    IJa = imwarp(Is,tform);
    Iq = IJa(round(size(IJa,1)/2+(1:size(Is,1))-size(Is,1)/2),...
        round(size(IJa,2)/2+(1:size(Is,2))-size(Is,2)/2));
    if sl==1;I=zeros([size(Is),numel(Nissl_dir)]);end
    IX(:,:,sl)=Is;
end
I = IX;
[no,xo]=hist(double(I(:)),256);
ma = xo(no==max(no));
I(629:702,11:112,:)=0;
I(616:714,701:831,:)=0;
I(I==0)=ma;I=uint8(abs(ma-I));
J=histeq(I);

path_out=[output_dir '\' '2dseq.nii'];
mat = [dx/1000,0,0,-size(I,1)/2*dx/1000;
       0,dy/1000,0,-size(I,2)/2*dy/1000;
       0,0,dz/1000,-size(I,3)/2*dz/1000;
       0,0,0,1];
   mat=mat*(10/6);mat(4,4)=1;
Vol(1,1) = struct(  'fname',    path_out,...
                    'dim',      size(I),...
                    'mat',      mat,...
                    'n',        [1,1],...
                    'pinfo',    [1;0;0],...
                    'descrip',  '',...
                    'dt',       [16 0]);
spm_write_vol(Vol,J);

Vol.fname = [output_dir '\' '2dseq_raw.nii'];
spm_write_vol(Vol,I);

source = MY_select_file_for_SPM(output_dir,'^2dseq.nii',[1 Inf]);
ref    = MY_select_file_for_SPM(output_dir,'^S0.nii',[1 Inf]);
all_func = MY_select_file_for_SPM(output_dir,'^2dseq_raw.nii',[1 Inf]);
coreg_mlb = MY_get_default_coreg_batch_struct(ref, source, all_func);
F = spm_figure('GetWin');
spm_jobman('run',coreg_mlb);

% I=spm_read_vols(spm_vol('c2dseq.nii'));
% lmask = ~isnan(I);
% ihdr = spm_vol('S0.nii');
% img = spm_read_vols(ihdr);
% ihdr.fname='cS0.nii';spm_write_vol(ihdr,img.*lmask);
