clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));
EyeTrackPath = 'D:\EyeTrack\';

cd(EyeTrackPath);
EyeTrackDir = dir('20*');
%% step1 extract the train images
%{
bla = [];
for iii = 2:numel(EyeTrackDir)
    cd(fullfile(EyeTrackDir(iii).folder,EyeTrackDir(iii).name));
    subDir = dir('*');
    subDir(1:2) = []; a=[];
    for jjj = 1:numel(subDir)
        switch subDir(jjj).isdir
            case true; a(jjj)=1;
            case false; a(jjj)=0;
        end
    end
    AVIPath = [];
    if isempty(find(a,1))
        AVIPath{1} = pwd;
    else
        for jjj = 1:numel(subDir)
            AVIPath{jjj} = fullfile(subDir(jjj).folder,subDir(jjj).name);
        end
    end
    
    for jj = 1:numel(AVIPath)
        cd(AVIPath{jj})

        OpenLevel = xlsread('OpenLevel.xlsx');
        dest = fullfile(AVIPath{jj});

        AVIdir = dir('de_*.avi');
        for idx = 1%:numel(AVIdir)
            EyeLevel = OpenLevel(:,idx);
            EyeLevel(isnan(EyeLevel))=[];
            EyeLevel(end-100:end)=[];
            
            UD(idx) = max(EyeLevel)-min(EyeLevel);
        end
        bar = max(UD)/mean(EyeLevel);
        idx = find(UD==max(UD));
        EyeLevel = OpenLevel(:,idx);
        EyeLevel(isnan(EyeLevel))=[];
        EyeLevel(end-100:end)=[];
        NUM = min(40,max(ceil(bar*100),20));
            
        [no,xo]=hist(EyeLevel,NUM);
        for iter=1:NUM
            r = find(abs(EyeLevel-xo(iter))==min(abs(EyeLevel-xo(iter))));
            train(iter)=min(r);
        end
        
        clear Cls
        fileinfo = VideoReader(AVIdir(idx).name);
        for iter=1:NUM
            Cls(:,:,iter) = rgb2gray(read(fileinfo,train(iter)));
        end
        bla = cat(3,bla,Cls);
        delete('Dictionary.nii');
        Vol = struct(  'fname',    fullfile(['Dictionary.nii']),...
                       'dim',      size(Cls),...
                       'mat',      diag([1 1 1 1]),...
                       'n',        [1,1],...
                       'pinfo',    [1;0;0],...
                       'descrip',  '',...
                       'dt',       [2 0]);
        spm_write_vol_4D(Vol,Cls);
        
    end
end
Vol = struct(  'fname',    fullfile(['D:\EyeTrack\train\DictionaryAll.nii']),...
               'dim',      size(bla),...
               'mat',      diag([1 1 1 1]),...
               'n',        [1,1],...
               'pinfo',    [1;0;0],...
               'descrip',  '',...
               'dt',       [2 0]);
spm_write_vol_4D(Vol,bla);
%}

%{
cd('D:\EyeTrack\train');
IMGall=spm_read_vols(spm_vol('DictionaryAll.nii'));
[a,b,c]=size(IMGall);
IMGVec=reshape(IMGall,[a*b c]);
lidall=spm_read_vols(spm_vol('lid.nii'));
pupilall=spm_read_vols(spm_vol('pupil.nii'));

bla = [];p=0; part=0;
for iii = 1:numel(EyeTrackDir)
    cd(fullfile(EyeTrackDir(iii).folder,EyeTrackDir(iii).name));
    subDir = dir('*');
    subDir(1:2) = []; a=[];
    for jjj = 1:numel(subDir)
        switch subDir(jjj).isdir
            case true; a(jjj)=1;
            case false; a(jjj)=0;
        end
    end
    AVIPath = [];
    if isempty(find(a,1))
        AVIPath{1} = pwd;
    else
        for jjj = 1:numel(subDir)
            AVIPath{jjj} = fullfile(subDir(jjj).folder,subDir(jjj).name);
        end
    end
    
    for jj = 1:numel(AVIPath)
        cd(AVIPath{jj})
        OpenLevel = xlsread('OpenLevel.xlsx');
        dest = fullfile(AVIPath{jj});

        AVIdir = dir('de_*.avi');
        for idx = 1%:numel(AVIdir)
            EyeLevel = OpenLevel(:,idx);
            EyeLevel(isnan(EyeLevel))=[];
            EyeLevel(end-100:end)=[];
            
            UD(idx) = max(EyeLevel)-min(EyeLevel);
        end
        bar = max(UD)/mean(EyeLevel);
        idx = find(UD==max(UD));
        EyeLevel = OpenLevel(:,idx);
        EyeLevel(isnan(EyeLevel))=[];
        EyeLevel(end-100:end)=[];
        NUM = min(40,max(ceil(bar*100),20));
        
        [no,xo]=hist(EyeLevel,NUM);
        for iter=1:NUM
            r = find(abs(EyeLevel-xo(iter))==min(abs(EyeLevel-xo(iter))));
            train(iter)=min(r);
        end
        
        clear IMG lid pupil
        fileinfo = VideoReader(AVIdir(idx).name);
        for iter=1:NUM
            img0 = rgb2gray(read(fileinfo,train(iter)));
            [a,b,c]=size(img0);
            Vec=reshape(img0,[a*b c]);
            t = corr(double(Vec),IMGVec);
            p = find(t==max(t));
            IMG(:,:,iter) = IMGall(:,:,min(p));
            lid(:,:,iter)=lidall(:,:,min(p));
            pupil(:,:,iter)=pupilall(:,:,min(p));        
        end
        
        
        part = max(part,NUM);
        Seg = p+1:p+NUM;
        
        IMG=IMGall(:,:,Seg);
        lid=lidall(:,:,Seg);
        pupil=pupilall(:,:,Seg);

        
        p = p+part;
        
        delete('IMG.nii')
        Vol = struct(  'fname',    fullfile(pwd,['IMG.nii']),...
                       'dim',      size(IMG),...
                       'mat',      diag([1 1 1 1]),...
                       'n',        [1,1],...
                       'pinfo',    [1;0;0],...
                       'descrip',  '',...
                       'dt',       [2 0]);
        spm_write_vol_4D(Vol,IMG);
        Vol = struct(  'fname',    fullfile(['label1.nii']),...
                       'dim',      size(lid),...
                       'mat',      diag([1 1 1 1]),...
                       'n',        [1,1],...
                       'pinfo',    [1;0;0],...
                       'descrip',  '',...
                       'dt',       [2 0]);
        spm_write_vol_4D(Vol,lid);
        Vol = struct(  'fname',    fullfile(['label2.nii']),...
                       'dim',      size(pupil),...
                       'mat',      diag([1 1 1 1]),...
                       'n',        [1,1],...
                       'pinfo',    [1;0;0],...
                       'descrip',  '',...
                       'dt',       [2 0]);
        spm_write_vol_4D(Vol,pupil);
        
    end
    if iii==5
       1; 
    end
end
%}
%% step2 the train dataset
for iii = 2:numel(EyeTrackDir)
    cd(fullfile(EyeTrackDir(iii).folder,EyeTrackDir(iii).name));
    subDir = dir('*');
    subDir(1:2) = []; a=[];
    for jjj = 1:numel(subDir)
        switch subDir(jjj).isdir
            case true; a(jjj)=1;
            case false; a(jjj)=0;
        end
    end
    AVIPath = [];
    if isempty(find(a,1))
        AVIPath{1} = pwd;
    else
        for jjj = 1:numel(subDir)
            AVIPath{jjj} = fullfile(subDir(jjj).folder,subDir(jjj).name);
        end
    end
    
    for jj = 1:numel(AVIPath)
        cd(AVIPath{jj})
        
        cd('D:\EyeTrack\20190113_nih\6165m\')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        OpenLevel = xlsread('OpenLevel.xlsx');
        dest = fullfile(AVIPath{jj});

        AVIdir = dir('de_*.avi');
        for idx = 5%:numel(AVIdir)
            EyeLevel = OpenLevel(:,idx);
            EyeLevel(isnan(EyeLevel))=[];
            EyeLevel(end-100:end)=[];
            
            UD(idx) = max(EyeLevel)-min(EyeLevel);
        end
        bar = max(UD)/mean(EyeLevel);
        idx = find(UD==max(UD));
        EyeLevel = OpenLevel(:,idx);
        EyeLevel(isnan(EyeLevel))=[];
        EyeLevel(end-100:end)=[];
        NUM = min(40,max(ceil(bar*100),20))*2;
            
        [no,xo]=hist(EyeLevel,NUM);
        for iter=1:NUM
            r = find(abs(EyeLevel-xo(iter))==min(abs(EyeLevel-xo(iter))));
            train(iter)=min(r);
        end
        
        clear Cls
        fileinfo = VideoReader(AVIdir(idx).name);
        for iter=1:NUM
            Cls(:,:,iter) = rgb2gray(read(fileinfo,train(iter)));
        end
        delete('Dictionary.nii');
        Vol = struct(  'fname',    fullfile(['IMG.nii']),...
                       'dim',      size(Cls),...
                       'mat',      diag([1 1 1 1]),...
                       'n',        [1,1],...
                       'pinfo',    [1;0;0],...
                       'descrip',  '',...
                       'dt',       [2 0]);
        spm_write_vol_4D(Vol,Cls);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IMG=spm_read_vols(spm_vol('IMG.nii'));
%         rmdir('image','s');mkdir('image');
        for i=1:size(IMG,3)
            Ig = uint8(IMG(:,:,i));
            imwrite(Ig,fullfile(pwd,'image',[num2str(i-1),'.png']));
        end        
        IMG=spm_read_vols(spm_vol('label1.nii'));
%         rmdir('label1','s');mkdir('label1');
        for i=1:size(IMG,3)
            Ig = uint8(IMG(:,:,i)*255);
            imwrite(Ig,fullfile(pwd,'label1',[num2str(i-1),'.png']));
        end        
        IMG=spm_read_vols(spm_vol('label2.nii'));
%         rmdir('label2','s');mkdir('label2');
        for i=1:size(IMG,3)
            Ig = uint8(IMG(:,:,i)*255);
            imwrite(Ig,fullfile(pwd,'label2',[num2str(i-1),'.png']));
        end        
    end
end

% cd('D:\EyeTrack\20191221_nih\61808');
% IMG=spm_read_vols(spm_vol('Dictionary.nii'));
% dest = 'D:\EyeTrack\20191221_nih\61808\image';
% for i=1:size(IMG,3)
%     Ig = uint8(IMG(:,:,i));
%     A = imresize(Ig,[240,320],'bilinear');
%     imwrite(A,fullfile(dest,[num2str(i-1),'.png']));
% end
IMG1=spm_read_vols(spm_vol('lib.nii'));
dest = 'D:\EyeTrack\20191221_nih\61808\image';
for i=1:size(IMG1,3)
    Ig = uint8(IMG1(:,:,i));
    A = imresize(Ig,[240,320],'bilinear');
    imwrite(A,fullfile(dest,[num2str(i-1),'.png']));
end

% cd('D:\EyeTrack\20191221_nih\61808');
% IMG=spm_read_vols(spm_vol('lmask.nii'));
% dest = 'D:\EyeTrack\20191221_nih\61808\label';
% for i=1:size(IMG,3)
%     I = uint8(IMG(:,:,i)*255);
%     A = imresize(I,[240,320],'nearest');
%     imwrite(A,fullfile(dest,[num2str(i-1),'.png']));
% end
IMG1=spm_read_vols(spm_vol('libmask.nii'));
dest = 'D:\EyeTrack\20191221_nih\61808\label';
for i=1:size(IMG1,3)
    Ig = uint8(IMG1(:,:,i)*255);
    A = imresize(Ig,[240,320],'nearest');
    imwrite(A,fullfile(dest,[num2str(i-1),'.png']));
end


cd('D:\EyeTrack\20191221_nih\61808');
fileinfo = VideoReader('de_61808m_run1.avi');
NFrames = fileinfo.NumberOfFrames;
dest = 'D:\EyeTrack\20191221_nih\61808\test';
for i=1:1000%NFrames
    A = rgb2gray(read(fileinfo,i*2));
    imwrite(A,fullfile(dest,[num2str(i-1),'.png']));
end



cd('D:\unet-master\data\membrane\test')
for i=1
   img = imread('999.png');
   lmask = imread('999_predict.png');
   img = imresize(img,[240,320],'bilinear');
   lmask = imresize(lmask,[240,320],'bilinear');
%    lmask = lmask/max(lmask(:));
   lmask(lmask<12.8)=0;
   lmask(lmask>=12.8)=1;
   BW = imfill(lmask,'holes');
   figure;imshow(img+20*BW) 
    
end

%}
