clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));
EyeTrackPath = 'D:\EyeTrack\';

%% Denoise and Down-sample the video
cd(EyeTrackPath);
EyeTrackDir = dir('120*');

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
    
   MY_Denoise_DownSample_AVIVideo(AVIPath);
end

function MY_Denoise_DownSample_AVIVideo(AVIPath)

for jj = 1:numel(AVIPath)
    cd(AVIPath{jj})
    AVIdir = dir('*.avi');
    for idx = 1:numel(AVIdir)
        fileinfo = VideoReader(AVIdir(idx).name);
        NFrames = fileinfo.NumberOfFrames;
        
        
        Id = uint8(zeros([fileinfo.Height,fileinfo.Width,floor(NFrames/6)]));
        
        parfor kk = 1:floor(NFrames/6)
            
            
            I6 = read(fileinfo,[(kk-1)*6+1 kk*6]);
            
            Ia = zeros([fileinfo.Height,fileinfo.Width,6]);
            for fr = 1:6
                I = rgb2gray(I6(:,:,:,fr));
                I = medfilt2(I,[10 5]);
                Ia(:,:,fr) = I;
            end
            
            Id(:,:,kk) = uint8(mean(Ia,3));    
        end
        

        w = VideoWriter(['de_',AVIdir(idx).name]);
        w.FrameRate = 10;
        open(w);
        
        for kk = 1:floor(NFrames/6)
            writeVideo(w, Id(:,:,kk));
        end
        close(w);
        fclose all;

    end
    % delete(fullfile(AVIPath{jj},'run*.avi'))
end
end