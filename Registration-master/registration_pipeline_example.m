% 1. specify filename for vtk atlas
% 2. specify directory for input images at 14.72 um
% they should have a csv file in this directory that indicates their
% geometry, if not, consider running create_location_csv.m first
% this file should have the form
% filename, nx, ny, nz, dx, dy, dz, x0, y0, z0
% 3. specify output directory for ALL transformation parameters, they are
% saved as matlab .mat files
% 4. specify output directory for vtk format transformations
% 5. (optional) specify a file containing segmentations, for QC figures only
%
% NOTE all directories should be strings with a trailing slash

warning('This example will run registration on a compressed dataset at low resolution.  See registration_pipeline_MBA.m for parameters used in our study and brainarchitecture.org for outputs at full resolution.')
tic
close all;
%%
% inputs
seg_file = 'D:\Registration-master\marmoset_atlas_vtk\annotation_50.nii';
atlas_file = 'D:\Registration-master\marmoset_atlas_vtk\Template_sym_S0.nii';
config_file = 'D:\Registration-master\example_nissl_fluoro_config.ini';


% go back to case 1, tro to enable edit mode
input_dir = 'D:\GeneAtlas\INPUT1\';
output_dir = 'D:\GeneAtlas\OUTPUT\';
detailed_output_dir = output_dir;

genename='A2M';
dx=25;dy=25;dz=250;
create_location_csv_Marmoset(input_dir, genename, 25,25,250) % μm

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1 is to create an initial slice alignment
% we use a simple initialization where we locate the center of mass of each
% slice and align by translation
pattern = '*Nissl*.jpg';
find_centers_for_initialization_nissl(input_dir, pattern, detailed_output_dir);
close all;

r = 10; % match to this many neighbors
% this downsampling and iterations is enough for a good initial guess
% not enough for a full accurate reconstruction
downs = [2,1]*3;
niter = 40;
skip_thick = -1; % no thick slices to be skipped
load_initializer = 1;
e = 0.1;
atlas_free_rigid_alignment(input_dir, pattern, detailed_output_dir, r, downs, niter, e, skip_thick, load_initializer)
close all;

% initial affine
downs = [1,2];
niter = 60;
affine_for_initial_alignment(atlas_file, input_dir, pattern, detailed_output_dir, downs, niter)
close all;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2 is to run mapping
% let's start by matching each fluoro slice to its nearest nissl slice
nissl_pattern = pattern;
files = dir([input_dir '*.jpg']);
files = {files.name};
fluoro_pattern = ['*',genename,'*.jpg'];

align_fluoro_to_nissl(input_dir, nissl_pattern, fluoro_pattern, detailed_output_dir);
close all;


% now 3D to 2D transformations for nissl
warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId')
warning('off','MATLAB:griddedInterpolant:MeshgridEval3DWarnId')

% here we will run registration at either only low res
% or only low and medium res
ThreeD_to_2D_registration(atlas_file, input_dir, pattern, config_file, detailed_output_dir)
close all;

% interleave these two transformations in appropriate format
combine_nissl_and_fluoro_transforms(detailed_output_dir)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 3 is to generate standard vtk outputs
apply_deformation({seg_file,atlas_file}, input_dir, detailed_output_dir, output_dir);
close all;

% bregma shift outputs
% convert from origin at the center of the image to the origin at the
% estimated bregma location

toc

