function MY_Precise_affine_slice_registration(files_dir,genename)


% source
source_dir = dir(fullfile(files_dir,'downs_pic',['*',genename,'*.jpg']));
for sloop=1:numel(source_dir)
    filename = fullfile(source_dir(sloop).folder,source_dir(sloop).name);
    IMG = imread(filename);
    if sloop==1;I=zeros([size(IMG),numel(source_dir)]);end
    I(:,:,sloop)=IMG;
end
% ref
ref_dir = dir(fullfile(files_dir,'1st_register',['*','atlas','*.jpg']));
for sloop=1:numel(ref_dir)
    filename = fullfile(ref_dir(sloop).folder,ref_dir(sloop).name);
    IMG = imread(filename);
    if sloop==1;J=zeros([size(IMG),numel(ref_dir)]);end
    J(:,:,sloop)=IMG;
end

other_dir = dir(fullfile(files_dir,'receptor_density',['*',genename,'*.jpg']));
for sloop=1:numel(other_dir)
    filename = fullfile(other_dir(sloop).folder,other_dir(sloop).name);
    IMG = imread(filename);
    if sloop==1;K=zeros([size(IMG),numel(other_dir)]);end
    K(:,:,sloop)=IMG;
end


for sloop = 1:size(I,3)
    moving = I(:,:,sloop)/255;
    fixed = J(:,:,sloop)/255;
    
    %{
    
    [no,xo]=hist(double(moving(:)),256);
    thrs = ceil(xo(no==max(no))*10)/10;
    Q=im2bw(moving,thrs);
    moving(~Q)=0;
    %}
    
        
    [optimizer, metric] = imregconfig('multimodal');
    %movingRegisteredDefault = imregister(moving, fixed, 'affine', optimizer, metric);
    %figure, imshowpair(movingRegisteredDefault, fixed,'montage');
    %optimizer.MaximumIterations = 300;
    tformSimilarity=imregtform(moving,fixed,'affine',optimizer,metric);
    Rfixed = imref2d(size(fixed));
    movingRegisteredRigid = imwarp(moving,tformSimilarity,'OutputView',Rfixed);
    
    movingRegisteredRigid = imhistmatch(movingRegisteredRigid,fixed);
    [D,movingRegisteredRigid] = ...
      imregdemons(movingRegisteredRigid,fixed,1000,'AccumulatedFieldSmoothing',2.0,'PyramidLevels',5);
    
    if sloop==1;Jw=zeros(size(I));end
    Jw(:,:,sloop) = movingRegisteredRigid;
    %subplot(6,10,sloop);imshowpair(fixed,movingRegisteredRigid1)
    
    
    K_movingRegisteredRigid = imwarp(K(:,:,sloop),tformSimilarity,'OutputView',Rfixed);
    K_movingRegisteredRigid = imwarp(K_movingRegisteredRigid,D);
    %figure, imshowpair(movingRegisteredRigid1, fixed);
    if sloop==1;Kw=zeros(size(I));end
    Kw(:,:,sloop) = K_movingRegisteredRigid;
    
end

dest = fullfile(files_dir,'affine_precise','struc');
if ~exist(dest,'dir');mkdir(dest);end
for sloop=1:size(Jw,3)
    filename = fullfile(dest,['struc_slice_',num2str(sloop,'%02d'),'.jpg']);
    IMG=Jw(:,:,sloop)*255;
    imwrite(uint8(IMG),filename);
end

dest = fullfile(files_dir,'affine_precise','receptor_density');
if ~exist(dest,'dir');mkdir(dest);end
for sloop=1:size(Jw,3)
    filename = fullfile(dest,['density_slice_',num2str(sloop,'%02d'),'.jpg']);
    IMG=Kw(:,:,sloop);
    imwrite(uint8(IMG),filename);
end

clear J_ I_
for sloop=1:size(I,3);J_{sloop}=Jw(:,:,sloop);end
tail=ceil(size(I,3)/6)*6-size(I,3);
if tail~=0;for ii=1:tail;J_{sloop+ii}=0*Jw(:,:,sloop);end;end
J__=cell2mat(reshape(J_,6,[]));
for sloop=1:size(I,3);I_{sloop}=J(:,:,sloop);end
if tail~=0;for ii=1:tail;I_{sloop+ii}=0*Jw(:,:,sloop);end;end
I__=cell2mat(reshape(I_,6,[]));
F=figure('Position',[1 41 1920 962]);
imshowpair(I__*255,J__*255)
saveas(F,fullfile(files_dir,'affine_precise','affine.tiff'))
close all

