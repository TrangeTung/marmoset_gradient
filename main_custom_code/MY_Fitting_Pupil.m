clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));
EyeTrackPath = 'D:\EyeTrack\';

cd(EyeTrackPath);
EyeTrackDir = dir('20*');

%% pupil lid fitting
%{
for iii = 41:numel(EyeTrackDir)
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
        AVIdir = dir('de_*.avi');
        
        Checkdir = dir('Check0*.avi');
        if ~isempty(Checkdir);continue;end
        
        for idx = 1:numel(AVIdir)
            avi_name = AVIdir(idx).name;
            lid_name = ['seg1_',replace(avi_name,'.avi','.mp4')];
            pupil_name = ['seg2_',replace(avi_name,'.avi','.mp4')];
            avi_info = VideoReader(avi_name);
            lid_info = VideoReader(lid_name);
            pupil_info = VideoReader(pupil_name);
            NFrames = avi_info.NumberOfFrames;
           
            cd(AVIPath{jj})
            Vid = uint8(zeros(240,320,3,NFrames));
            clear C_* R_* S_*
            for kk = 1:NFrames
                %% loading
                If = rgb2gray(read(avi_info,kk));
                lmask1 = rgb2gray(read(lid_info,kk));
                lmask1 = lmask1-min(lmask1(:));
                lmask2 = rgb2gray(read(pupil_info,kk));
                lmask2 = lmask2-min(lmask2(:));
                lmask1(lmask1<255*0.2)=0;
                lmask1(lmask1>=255*0.2)=1;
                lmask2(lmask2<255*0.25)=0;
                lmask2(lmask2>=255*0.25)=1;
                
                lmask1 = imfill(lmask1, 'holes');
                lmask2 = imfill(lmask2, 'holes');
                
                if numel(find(lmask1==1))>10
                    imLabel = bwlabel(lmask1);
                    stats = regionprops(imLabel,'Area');
                    [b,index]=sort([stats.Area],'descend');
                    lmask1=ismember(imLabel,index(1));
                end              
                %% eyelid
                if numel(find(lmask1==1))>1000
                    %{
                    BW = edge(lmask1,'canny');
                    [lmaskt,lmaskb] = MY_seg_eyelid_mask(BW);
                    [rt,Ct] = MY_circular_fitting_LMS(lmaskt==1);
                    [rb,Cb] = MY_circular_fitting_LMS(lmaskb==1);
                    Iga = insertShape(zeros(size(If)),'FilledCircle',[Ct(2),Ct(1),rt],...
                        'LineWidth',1,'Color','w','Opacity',1);
                    Igb = insertShape(zeros(size(If)),'FilledCircle',[Cb(2),Cb(1),rb],...
                        'LineWidth',1,'Color','w','Opacity',1);
                    Igab = rgb2gray(Iga)+rgb2gray(Igb);
                    Igab(Igab~=2)=0;
                    EYElid = edge(Igab,'canny');
                    %}
                    Igab = 2*lmask1;
                    EYElid = edge(Igab,'canny');
                    %{
                    if ~isempty(numel(lmask2==1)) 
                        %% pupil
                        pix = floor(sum(Igab(:))/4000);
                        SE = strel('line',pix,90);
                        er = imdilate(EYElid,SE);
                        BW = edge(lmask2,'canny');
                        BW(er==1)=0;
                        [rp,Cp] = MY_circular_fitting_LMS(BW==1);
                        if ~isnan(rp)
                            Igp = insertShape(zeros(size(If)),'FilledCircle',[Cp(2),Cp(1),rp],...
                                'LineWidth',1,'Color','w','Opacity',1);
                            Ip = rgb2gray(Igp);
                            EYEpupil = edge(Ip,'canny');
                        else
                            rp = NaN; Cp = [NaN;NaN]; Ip = zeros(size(If));
                        end
                    else
                        rp = NaN; Cp = [NaN;NaN]; Ip = zeros(size(If));
                    end
                    %}
                    
                    %% pupil area
                    lmask = lmask2.*uint8(Igab);
                    EYEpupil = edge(lmask,'canny');
                    Area = numel(find(lmask>0));
                else
                    rt = NaN; Ct = [NaN;NaN];
                    rb = NaN; Cb = [NaN;NaN];
                    Igab = zeros(size(If));
                    Area = 0;
                    EYEpupil=zeros(240,320);
                    EYElid=zeros(240,320);
                end
                %% Results
                %R_Pupil(kk,:) = rp;      C_Pupil(kk,:) = Cp; 
                %R_Lid(kk,:) = [rt,rb];   C_Lid(kk,:) = [Ct;Cb]; 
                %S_Pupil(kk,:) = sum(Ip(:));
                S_Lid(kk,:) = sum(Igab(:))/2;
                S_Overlay(kk,:) = Area;
                
                Vid(:,:,:,kk) = uint8(repmat(If,[1 1 3]))+...
                    255*uint8(cat(3,zeros(240,320),zeros(240,320),EYElid))+...
                    255*uint8(cat(3,EYEpupil,zeros(240,320),zeros(240,320)));
            end
            %% results
            Excel = fullfile(AVIPath{jj},'FittingResults.xlsx');
            xlswrite(Excel,  S_Lid,     'S_Lid',     [char(double('A')+(idx-1)*1+1-1),'1']);
            xlswrite(Excel,  S_Overlay, 'S_Overlay', [char(double('A')+(idx-1)*1+1-1),'1']);
            %% Video
            w = VideoWriter(['Check0_',avi_name]);
            w.FrameRate = 10;
            open(w);
                        
            for f = 1:floor(NFrames)
                Ig = squeeze(Vid(:,:,:,f));
                Igf = insertText(Ig,[0 0],['Frame ',num2str(f,'%04d')],'FontSize',12,'AnchorPoint','LeftTop');
                writeVideo(w, Igf);
            end
            close(w);
            fclose all;
        end
    end
end

%}

%% correct error fitting (S_Lid > 30000)
%{
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
        AVIdir = dir('de_*.avi');
        
        Checkdir = dir('Check1*.avi');
        if ~isempty(Checkdir);continue;end

        FitRes_Lid     = xlsread('FittingResults.xlsx','S_Lid');
        FitRes_Overlay = xlsread('FittingResults.xlsx','S_Overlay');

        for idx = 1:numel(AVIdir)
            S_Lid = FitRes_Lid(:,idx);
            S_Overlay = FitRes_Overlay(:,idx);
            
            avi_name = AVIdir(idx).name;
            check_name = ['Check0_',avi_name];
            check_info = VideoReader(check_name);
            NFrames = check_info.NumberOfFrames;
            Vid = read(check_info,[1 NFrames]);
            
            S_ERROR = find(S_Lid>2.5*10^4);
            NUM_ERROR = numel(S_ERROR);
            if ~isempty(S_ERROR)
                lid_name = ['seg1_',replace(avi_name,'.avi','.mp4')];
                pupil_name = ['seg2_',replace(avi_name,'.avi','.mp4')];

                avi_info = VideoReader(avi_name);
                lid_info = VideoReader(lid_name);
                pupil_info = VideoReader(pupil_name);
                
                
                for kk = 1:NUM_ERROR
                    frame = S_ERROR(kk);
                    % loading
                    If = rgb2gray(read(avi_info,frame));
                    lmask1 = rgb2gray(read(lid_info,frame));
                    lmask1 = lmask1-min(lmask1(:));
                    lmask2 = rgb2gray(read(pupil_info,frame));
                    lmask2 = lmask2-min(lmask2(:));
                    lmask1(lmask1<255*0.2)=0;
                    lmask1(lmask1>=255*0.2)=1;
                    lmask2(lmask2<255*0.25)=0;
                    lmask2(lmask2>=255*0.25)=1;
                    
                    lmask1 = imfill(lmask1, 'holes');
                    lmask2 = imfill(lmask2, 'holes');
                    lmask2 = lmask2.*lmask1;
                    
                    EYElid=zeros(240,320);
                    EYEpupil=zeros(240,320);
                    if ~isempty(find(lmask1~=0));EYElid=edge(lmask1,'canny');end
                    if ~isempty(find(lmask2~=0));EYEpupil=edge(lmask2,'canny');end
                    
                    Vid_Check1 = uint8(repmat(If,[1 1 3]))+...
                                255*uint8(cat(3,zeros(240,320),zeros(240,320),EYElid))+...
                                255*uint8(cat(3,EYEpupil,zeros(240,320),zeros(240,320)));
                    Igf = insertText(Vid_Check1,[0 0],['Frame ',num2str(frame,'%04d')],...
                                'FontSize',12,'AnchorPoint','LeftTop');
                    Vid(:,:,:,frame) = Igf;
                    
                    S_Lid(frame,:) = numel(find(lmask1~=0));
                    S_Overlay(frame,:) = numel(find(lmask2~=0));              
                end
            end
            % check1_*.avi
            Excel = fullfile(AVIPath{jj},'FittingResults.xlsx');
            xlswrite(Excel,  S_Lid,     'S_Lid1',     [char(double('A')+(idx-1)*1+1-1),'1']);
            xlswrite(Excel,  S_Overlay, 'S_Overlay1', [char(double('A')+(idx-1)*1+1-1),'1']);
            %% Video
            w = VideoWriter(['Check1_',avi_name]);
            w.FrameRate = 10;
            open(w);

            for f = 1:floor(NFrames)
                Ig = squeeze(Vid(:,:,:,f));
                writeVideo(w, Ig);
            end
            close(w);
            fclose all;
        end
    end
end
%}

%% check the S_Lid and S_pupil
%
for iii = 3:numel(EyeTrackDir)
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
        AVIdir = dir('de_*.avi');
        FitRes_Lid     = xlsread('FittingResults.xlsx','S_Lid1');
        FitRes_Overlay = xlsread('FittingResults.xlsx','S_Overlay1');
        
        for kk=1:numel(AVIdir)
            
            F = figure;set(F,'position',[0 0 1900 900]);
            suptitle([AVIPath{jj}]);
            
            S_Lid = FitRes_Lid(:,kk); 
            S_Lid(isnan(S_Lid))=[];
            S_Lid=medfilt1(S_Lid,20);
            S_Overlay = FitRes_Overlay(:,kk); 
            S_Overlay(isnan(S_Overlay))=[];
            S_Overlay=medfilt1(S_Overlay,20);
            % Lid figure
            subplot(2,1,1);
            
            
            
            subplot(2,1,2);
            plot(S_Overlay,'r'); xlim([0 11000]);
            
            
            
            hgexport(figure(F), fullfile(AVIPath{jj},strcat('FittingResults')), hgexport('factorystyle'), 'Format', 'tiff');
            close all
            
        end
        
    end
    
end

function MY_generate_patch_figure(S)
    plot(S,'b','linewidth',2); xlim([0 11000]);%ylim([-1 2000])
    Drowy = S==0;
    S = (S-median(S(S~=0)))/median(S(S~=0));
    Qawake = S<0.5 & S>-0.9;
    Fawake = S>=0.5;
    
    
    
    
end
