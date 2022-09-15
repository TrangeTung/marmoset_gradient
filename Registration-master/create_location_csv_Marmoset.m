function create_location_csv_Marmoset(directory,genename, dx, dy, dz, outdir,ext)
% because the sorting key may change, this should ONLY be used carefully
% also this file hard codes 15 micron pixels
% note it is actually 0.46*32 = 14.72
% width is 20, but distance from nissl to nissl is 40
% e.g. run like this
% create_location_csv_MBA('/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/PTM830/', 14.72, 14.72, 20)
% here we will identify the nissl images, and the "other" images

% keyboard
if nargin < 7
    ext = 'jpg';
end
if nargin < 6
    outdir = directory; % same input and output
end
disp(['input directory ' directory])

files_g = dir([directory '*',genename,'*.' ext]); % get the gene only

files_ = arrayfun(@(x)strsplit(x.name,'_'),files_g,'UniformOutput',false);
ar = cellfun(@(x)strsplit(x{end},'.'),files_,'UniformOutput',false);
array_g = cellfun(@(x)str2num(x{1}),ar,'UniformOutput',false);

% now we write it out
fid = fopen([outdir 'geometry.csv'],'wt');
% first print the fields
fprintf(fid,'filename, nx, ny, nz, dx, dy, dz, x0, y0, z0\n');



for i=1:numel(array_g)
    fprintf(fid,'%s, ', files_g(i).name);
    info = imfinfo([directory files_g(i).name]);
    fprintf(fid,'%d, ', info.Width);
    fprintf(fid,'%d, ', info.Height);
    fprintf(fid,'%d, ', 1);
    fprintf(fid,'%f, ', dx);
    fprintf(fid,'%f, ', dy);
    fprintf(fid,'%f, ', dz);
    
    if i==1
        z = (1:numel(array_g))*dz;
        z = z - (z(1) + z(end))/2;
    end
    
   
    
    % for offset we will subtract the mean
    x = (0 : info.Width-1)*dx;
    x = x - mean(x);
    fprintf(fid,'%f, ', x(1));
    y = (0 : info.Height-1)*dy;
    y = y - mean(y);
    fprintf(fid,'%f, ', y(1));
    % now for the z offset 
    fprintf(fid,'%f, ', z(i));
    fprintf(fid,'\n');
    
end






% close the file
fclose(fid);