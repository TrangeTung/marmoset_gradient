function MY_3D_coreg_slice_to_template(files_dir,codepath)

JPG_dir = dir(fullfile(files_dir,'1st_register',['atlas_*.jpg']));
dx=25;dy=25;dz=250;
for sloop=1:numel(JPG_dir)
    filename = fullfile(JPG_dir(sloop).folder,JPG_dir(sloop).name);
    IMG = imread(filename);
    if sloop==1;I=zeros([size(IMG),numel(JPG_dir)]);end
    I(:,:,sloop)=IMG;
end
path_out=[files_dir '\' 'atlas.nii'];
mat = [0,dy/1000,0,-size(I,2)/2*dy/1000;
    0,0,-dz/1000,+size(I,3)/2*dz/1000;
    -dx/1000,0,0,+size(I,1)/2*dx/1000;
    0,0,0,1];
mat=mat*(6);mat(4,4)=1;
Vol(1,1) = struct(  'fname',    path_out,   'dim',      size(I),...
    'mat',      mat,        'n',        [1,1],...
    'pinfo',    [1;0;0],...
    'descrip',  '',         'dt',       [16 0]);
spm_write_vol(Vol,I/2);
clear Vol


JPG_dir = dir(fullfile(files_dir,'affine_precise','struc',['struc_*.jpg']));
dx=25;dy=25;dz=250;
for sloop=1:numel(JPG_dir)
    filename = fullfile(JPG_dir(sloop).folder,JPG_dir(sloop).name);
    IMG = imread(filename);
    IMG(IMG==0)=nan;
    [no,xo] = hist(double(IMG(:)),256);
    backg = round(xo(no==max(no)));
    IMG(IMG==backg)=nan;
    IMG(isnan(IMG))=backg;
    if sloop==1;I=zeros([size(IMG),numel(JPG_dir)]);end
    I(:,:,sloop)=IMG;
end
path_out=[files_dir '\' 'struc.nii'];
mat = [0,dy/1000,0,-size(I,2)/2*dy/1000;
    0,0,-dz/1000,+size(I,3)/2*dz/1000;
    -dx/1000,0,0,+size(I,1)/2*dx/1000;
    0,0,0,1];
mat=mat*(6);mat(4,4)=1;
Vol(1,1) = struct(  'fname',    path_out,   'dim',      size(I),...
    'mat',      mat,        'n',        [1,1],...
    'pinfo',    [1;0;0],...
    'descrip',  '',         'dt',       [16 0]);
spm_write_vol(Vol,I);
clear Vol



JPG_dir = dir(fullfile(files_dir,'affine_precise','receptor_density',['density_*.jpg']));
dx=25;dy=25;dz=250;
for sloop=1:numel(JPG_dir)
    filename = fullfile(JPG_dir(sloop).folder,JPG_dir(sloop).name);
    IMG = imread(filename);
    if sloop==1;J=zeros([size(IMG),numel(JPG_dir)]);end
    J(:,:,sloop)=IMG;
end
path_out=[files_dir '\' 'density.nii'];
mat = [0,dy/1000,0,-size(J,2)/2*dy/1000;
    0,0,-dz/1000,+size(J,3)/2*dz/1000;
    -dx/1000,0,0,+size(J,1)/2*dx/1000;
    0,0,0,1];
mat=mat*(6);mat(4,4)=1;
Vol(1,1) = struct(  'fname',    path_out,   'dim',      size(I),...
    'mat',      mat,        'n',        [1,1],...
    'pinfo',    [1;0;0],...
    'descrip',  '',         'dt',       [16 0]);
spm_write_vol(Vol,J);
clear Vol


source = MY_select_file_for_SPM(files_dir,'^atlas.nii',[1 Inf]);
other1  = MY_select_file_for_SPM(files_dir,'^density.nii',[1 Inf]);
other2  = MY_select_file_for_SPM(files_dir,'^struc.nii',[1 Inf]);
ref    = MY_select_file_for_SPM(codepath,'^Template_sym_T2_6X.nii',[1 Inf]);
coreg_mlb = MY_get_default_oldnormalize_batch_struct(ref, source, [source;other1;other2]);
% coreg_mlb{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 1;
% coreg_mlb{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 1;
coreg_mlb{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
F = spm_figure('GetWin');
spm_jobman('run',coreg_mlb);


end