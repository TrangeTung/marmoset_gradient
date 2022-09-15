function MY_reslice_atlas_to_Histology(atlas_file,input_dir,initializer_A_URL,dx,dy,dz)

JPGdir = dir(fullfile(input_dir,'*.jpg'));
for jloop=1:numel(JPGdir)
    X = imread(fullfile(JPGdir(jloop).folder,JPGdir(jloop).name));
    if jloop==1;J=zeros([size(X),numel(JPGdir)]);end
    J(:,:,jloop)=X;
end

load(initializer_A_URL);
ref = atlas_file;
ihdr=spm_vol(ref);
I = spm_read_vols(ihdr);
[x,y,z] = size(I);
vox = sqrt(sum(ihdr.mat(1:3,1:3).^2));
vx = vox(1)*10^3;
vy = vox(2)*10^3;
vz = vox(3)*10^3;
 
mul = [vx/dx,vy/dy,vz/dz];
Iw = imresize3(I,round(size(I).*mul));

tform = diag([1;1;1;1]);
tform(1:3,1:3)=Ainv(1:3,1:3);
Ia = imwarp(Iw,affine3d(tform));

v = MY_kspace3d(Iw,Ainv);

xstep=round(Ainv(1,4)/dx);
ystep=round(Ainv(2,4)/dy);
zstep=round(Ainv(3,4)/dz);

Iq = Ia((round(size(Ia,1)/2+(1:size(J,1))-size(J,1)/2))+xstep,...
        (round(size(Ia,2)/2+(1:size(J,2))-size(J,2)/2))+ystep,...
        (round(size(Ia,3)/2+(1:size(J,3))-size(J,3)/2))+zstep);

   
    

end