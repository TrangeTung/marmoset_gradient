function MY_register_Histology_Marmoset(genename,atlas_file,files_dir)


input_dir = fullfile(files_dir,'downs_pic\');
output_dir = fullfile(files_dir,'1st_register\');
mkdir(output_dir);
detailed_output_dir = output_dir;

dx=25;dy=25;dz=250;
create_location_csv_Marmoset(input_dir, genename, dx,dy,dz) % μm

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1 is to create an initial slice alignment
% we use a simple initialization where we locate the center of mass of each
% slice and align by translation
pattern = ['*',genename,'*.jpg'];
% find_centers_for_initialization_nissl(input_dir, pattern, detailed_output_dir);
% close all;

%{
r = 25; % match to this many neighbors
% this downsampling and iterations is enough for a good initial guess
% not enough for a full accurate reconstruction
downs = [4 2];
niter = 40;
skip_thick = -1; % no thick slices to be skipped
load_initializer = 1;
e = 0.1;
atlas_free_rigid_alignment(input_dir, pattern, detailed_output_dir, r, downs, niter, e, skip_thick, load_initializer)
close all;
%}

downs = [6,3,1];
niter = [100,100,5];
affine_for_initial_alignment(atlas_file, input_dir, pattern, detailed_output_dir, downs, niter)
close all;

% 
% 
% 
% copyfile(atlas_file,output_dir);
% 
% initializer_A_URL = [output_dir,'initializer_A.mat'];
% atlas_file = fullfile(output_dir,'Gene_S0.nii');
% MY_reslice_atlas_to_Histology(atlas_file,input_dir,initializer_A_URL,dx,dy,dz)
%     
% Nissl_dir = dir([input_dir ['*',genename,'*.jpg']]);
% for sl = 1:numel(Nissl_dir)
%     Is = imread(fullfile(input_dir,Nissl_dir(sl).name));
%     tform=affine2d(AJ(:,:,sl)');
%     IJa = imwarp(Is,tform);
%     Iq = IJa(round(size(IJa,1)/2+(1:size(Is,1))-size(Is,1)/2),...
%         round(size(IJa,2)/2+(1:size(Is,2))-size(Is,2)/2));
%     if sl==1;I=zeros([size(Is),numel(Nissl_dir)]);end
%     IX(:,:,sl)=Is;
% end
% 
% 
% path_out=[output_dir '\' '2dseq.nii'];
% mat = [dx/1000,0,0,-size(I,1)/2*dx/1000;
%        0,dy/1000,0,-size(I,2)/2*dy/1000;
%        0,0,dz/1000,-size(I,3)/2*dz/1000;
%        0,0,0,1];
%    mat=mat*(10/6);mat(4,4)=1;
% Vol(1,1) = struct(  'fname',    path_out,...
%                     'dim',      size(IX),...
%                     'mat',      mat,...
%                     'n',        [1,1],...
%                     'pinfo',    [1;0;0],...
%                     'descrip',  '',...
%                     'dt',       [16 0]);
% spm_write_vol(Vol,IX);
% 
% 
% copyfile(atlas_file,output_dir);
% ref    = MY_select_file_for_SPM(output_dir,'^2dseq.nii',[1 Inf]);
% source = MY_select_file_for_SPM(output_dir,'^Gene_S0.nii',[1 Inf]);
% coreg_mlb = MY_get_default_coreg_batch_struct(ref, source, source);
% F = spm_figure('GetWin');
% spm_jobman('run',coreg_mlb);

% I=spm_read_vols(spm_vol('c2dseq.nii'));
% lmask = ~isnan(I);
% ihdr = spm_vol('S0.nii');
% img = spm_read_vols(ihdr);
% ihdr.fname='cS0.nii';spm_write_vol(ihdr,img.*lmask);
