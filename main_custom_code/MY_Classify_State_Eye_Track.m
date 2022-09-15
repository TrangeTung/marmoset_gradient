clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));
EyeTrackPath = 'D:\EyeTrack\';

%% Denoise and Down-sample the video
cd(EyeTrackPath);
EyeTrackDir = dir('20*');
%{
for iii = 36%12*3+1:numel(EyeTrackDir)
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
        
        Is = [];
        for idx = numel(AVIdir):-1:1
            fileinfo = VideoReader(AVIdir(idx).name);
            NFrames = fileinfo.NumberOfFrames;
            
            stat = zeros([NFrames,1]);
            % 0 -- close  % 1 -- open but no pupil size
            Detect = vision.CascadeObjectDetector('LeftEye','MergeThreshold',16);
            for kk = 1:NFrames
               If = rgb2gray(read(fileinfo,kk));
               if kk==1; I=uint8(zeros([size(If) NFrames])); end
               bd = step(Detect,If);
               if isempty(bd); stat(kk)=0;else;stat(kk)=1;end
               I(:,:,kk) = If;
            end
            save(fullfile(AVIPath{jj},['stat',num2str(idx),'.txt']),'stat','-ascii');
            
            I(:,:,stat==0)=[];
            if ~isempty(I)
                % gray image clustering
                [a,b,c] = size(I);
                Iv = reshape(I,[a*b c]);
                
                Ks = ceil(sqrt(c));
                if Ks<=50;Ks=min(Ks,c);end
                if Ks~=1
                    [IDX, C] = kmeans(double(Iv'), Ks, 'Distance'  , 'city') ;
                    Cls = reshape(C',[a,b,Ks]);
                else
                    IDX=1;Cls=I;
                end
                save(fullfile(AVIPath{jj},['IDX',num2str(idx),'.txt']),'IDX','-ascii');
                Vol = struct(  'fname',    fullfile(AVIPath{jj},['lib',num2str(idx),'.nii']),...
                               'dim',      [a,b,Ks],...
                               'mat',      diag([1 1 1 1]),...
                               'n',        [1,1],...
                               'pinfo',    [1;0;0],...
                               'descrip',  '',...
                               'dt',       [2 0]);
                spm_write_vol_4D(Vol,Cls);
                
            end
            
            
        end
    end
end



%}


%% pupile radius fitting
%{

for iii = 2% 1:numel(EyeTrackDir)
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
    
    for jj = 1%1:numel(AVIPath)
        cd(AVIPath{jj})
        AVIdir = dir('de_*.avi');
        % load pupil library
        %         libdir = dir('lib*.nii');
        %         lmaskdir = dir('lmask*.nii');
        %         lib=[];lmask=[];
        %         for i=1:numel(libdir)
        %             li = spm_read_vols_4D(spm_vol(libdir(idx).name));
        %             lm = spm_read_vols_4D(spm_vol(lmaskdir(idx).name));
        %             lib=cat(3,lib,li);
        %             lmask=cat(3,lmask,lm);
        %         end
        %         [a,b,c] = size(lib);
        %         libVec = reshape(lib,[a*b c]);
        
        
        for idx = 1%1:numel(AVIdir)
            
            
            lib = spm_read_vols_4D(spm_vol('lib1.nii'));
            lmask = spm_read_vols_4D(spm_vol('lmask1.nii'));
            [a,b,c] = size(lib);
            libVec = reshape(lib,[a*b c]);
            
            
            
            EyeState = load(['stat',num2str(1),'.txt']);
            % 0 -- close  % 1 -- open but no pupil size
            % 2 -- pupil radius
            ClsIDX = load(['IDX',num2str(idx),'.txt']);
            
            fileinfo = VideoReader(AVIdir(idx).name);
            NFrames = fileinfo.NumberOfFrames;
            
            Radius = zeros(NFrames,1);
            Center = zeros(NFrames,2); % 1-- row; 2-- column
            p=0;
            
            for kk = 1:NFrames
                I = read(fileinfo,kk);
                Ig = rgb2gray(I);
                h= fspecial('gaussian',[5 5],10);
                Ig = uint8(filter2(h,Ig));
                
                if kk==1;Igc=uint8(zeros([size(I) NFrames]));end
                Igc(:,:,:,kk) = I;
                
                    %if EyeState(kk)==1
                    if kk==1
                        [a,b,c]=size(Ig);
                        Vec = reshape(Ig,[a*b c]);
                        CC = corr(double(Vec),libVec);
                        if p~=0;CC(p) = CC(p)*1.01;end
                        p = find(CC==max(CC));

                        %q=q+1; p=ClsIDX(q);
                        pmask = lmask(:,:,p);
                    else
                        pmask=prior;
                        
                    end
                    if ~isempty(find(pmask==1))
                        EyeState(kk) = 2;
                        text = {'active awake'};
                        
                        % library mask or -1 frame , better?
                        eyel = (pmask==1);
                        eyer = (pmask==2);
                        
                        SE = strel('line',5,1);
                        eyeld = imdilate(eyel,SE);
                        eyerd = imdilate(eyer,SE);
                        BW = edge(Ig,'canny');
                        eyel = eyeld.*BW;
                        eyer = eyerd.*BW;
                        
                        [r,C] = MY_circular_fitting_LMS(eyel+eyer);
                        
                        if r<47 || r>65 || isnan(r);
                            SE = strel('line',30,1);
                            eyeld = imdilate((pmask==1),SE); eyerd = imdilate((pmask==2),SE);
                            BW = edge(Ig,'canny');
                            eyel = eyeld.*BW; eyer = eyerd.*BW;
                            
                            [r,C] = MY_circular_fitting_LMS(eyel+eyer);
                        end
                        
                        Igf = insertShape(Ig,'circle',[C(2),C(1),r],'LineWidth',1,'Color','w');
                        Igf = insertShape(Igf,'circle',[C(2),C(1),1],'LineWidth',3,'Color','r');
                        figure;imshowpair(Igf,double(Ig)+(eyel+eyer)*100,'montage');
                        
                        Radius(kk) = r; Center(kk,:) = C;
                        Igc(:,:,:,kk) = Igf;
                        
                    else
                        EyeState(kk) = 1;
                        text = {'quiet awake'};
                    end
                %end
                
            end
            
            w = VideoWriter(['Check_',AVIdir(idx).name]);
            w.FrameRate = 10;
            open(w);
            
            for kk = 1:floor(NFrames)
                writeVideo(w, Igc(:,:,:,kk));
            end
            close(w);
            fclose all;
            
            Excel = fullfile(AVIPath{jj},'Radius_Center.xlsx');
            xlswrite(Excel,{'Radius','CenterX','CenterY'},['RUN',num2str(idx)],'A1');
            xlswrite(Excel,[Radius,Center],['RUN',num2str(idx)],'A2');
            clear I*;
            
        end
        
    end
end
%}


%% resort the lib and mask the edge
%{
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
        libdir = dir('lib*.nii');
        if ~isempty(libdir)
        IMG=[];
        for i = 1:numel(libdir)
            ram=spm_read_vols_4D(spm_vol(libdir(i).name));
            IMG=cat(3,IMG,ram);
        end
        [a,b,c]=size(IMG);
        if c>100
            Iv = reshape(IMG,[a*b c]);
            Ks = 100;
            [IDX, C] = kmeans(double(Iv'), Ks, 'Distance'  , 'city') ;
            Cls = reshape(C',[a,b,Ks]);
        else
            Cls = IMG; Ks=c;
        end
        
        Vol = struct(  'fname',    fullfile(AVIPath{jj},['Dictionary.nii']),...
                       'dim',      [a,b,Ks],...
                       'mat',      diag([1 1 1 1]),...
                       'n',        [1,1],...
                       'pinfo',    [1;0;0],...
                       'descrip',  '',...
                       'dt',       [2 0]);
        spm_write_vol_4D(Vol,Cls);
        
        for i = 1:Ks
           I = Cls(:,:,i);
           BW=edge(I,'canny');
           Cls(:,:,i) = Cls(:,:,i)+(BW)*50;
        end
        Vol = struct(  'fname',    fullfile(AVIPath{jj},['DictionaryEdge.nii']),...
                       'dim',      [a,b,Ks],...
                       'mat',      diag([1 1 1 1]),...
                       'n',        [1,1],...
                       'pinfo',    [1;0;0],...
                       'descrip',  '',...
                       'dt',       [2 0]);
        spm_write_vol_4D(Vol,Cls);
        delete('lib*.nii');
        delete('stat*.txt');
        delete('IDX*.txt');
        end
    end
    
end

%% Eyelid fitting
%{
%{
for iii = 2% 1:numel(EyeTrackDir)
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
        
        % load pupil library
        libdir = dir('lib*.nii');
        lmaskdir = dir('lmask*.nii');
%         lib=[];lmask=[];
%         for idx=1:numel(libdir)
%             li = spm_read_vols_4D(spm_vol(libdir(idx).name));
%             lm = spm_read_vols_4D(spm_vol(lmaskdir(idx).name));
%             lib=cat(3,lib,li);
%             lmask=cat(3,lmask,lm);
%         end
%         lmask(lmask==1|lmask==2)=0;
%         lmask(lmask==3)=1;
%         lmask(lmask==4)=2;
        lib = spm_read_vols_4D(spm_vol('libqqq1.nii'));
        lmask = spm_read_vols_4D(spm_vol('lmaskqqq1.nii'));
        [a,b,c]=size(lib);
        libVec = reshape(lib,[a*b,c]);
        
        libp=[];lmaskp=[];
        for idx=1:numel(libdir)
            li = spm_read_vols_4D(spm_vol(libdir(idx).name));
            lm = spm_read_vols_4D(spm_vol(lmaskdir(idx).name));
            libp=cat(3,libp,li);
            lmaskp=cat(3,lmaskp,lm);
        end
        [a,b,c]=size(libp);
        libpVec = reshape(libp,[a*b,c]);
        
        % find the anchor of eyelid
        Im=uint8(mean(libp,3));
        BW=edge(Im,'canny');
        varargin = struct('masklib',lmask,'imglib',lib,'img',Im);
        pmask = MY_eyelid_fitting(1,varargin);
        %[rt,Ct] = MY_circular_fitting_LMS(pmask==1);
        %[rb,Cb] = MY_circular_fitting_LMS(pmask==2);
        varargin = struct('masklib',pmask,'imglib',[],'img',Im);
        [t,b]= MY_eyelid_fitting(2,varargin);
        rt = t(1); Ct = t(2:3);
        rb = b(1); Cb = b(2:3);
        Iga = insertShape(zeros(size(Im)),'FilledCircle',[Ct(2),Ct(1),rt],...
            'LineWidth',1,'Color','w','Opacity',1);
        Igb = insertShape(zeros(size(Im)),'FilledCircle',[Cb(2),Cb(1),rb],...
            'LineWidth',1,'Color','w','Opacity',1);
        Igab = rgb2gray(Iga)+rgb2gray(Igb);
        Igab(Igab~=2)=0;
        po=find(mean(Igab)>0);
        cla = min(po); clb = find(Igab(:,cla)~=0);
        cra = max(po); crb = find(Igab(:,cra)~=0);
        I = insertShape(zeros(size(Im)),'FilledCircle',[mean(cla),mean(clb),2],...
            'LineWidth',1,'Color','w','Opacity',1);
        I = insertShape(I              ,'FilledCircle',[mean(cra),mean(crb),2],...
            'LineWidth',1,'Color','w','Opacity',1);
        anchor = rgb2gray(I);
        
        
        for idx = 5:numel(AVIdir)
            
            fileinfo = VideoReader(AVIdir(idx).name);
            NFrames = fileinfo.NumberOfFrames;
            
            for kk = 1:NFrames
                I = read(fileinfo,kk);
                
                Ig = rgb2gray(I);
                h= fspecial('gaussian',[15 15],20);
                Ig = uint8(filter2(h,Ig));
                
%                 if kk==100 || kk==NFrames
                    varargin = struct('masklib',lmask,'imglib',lib,'img',Ig);
                    pmask = MY_eyelid_fitting(1,varargin);
%                 else
%                     pmask = eyeti + eyebi*2;
%                 end
                
                varargin = struct('masklib',pmask,'imglib',anchor,'img',Ig);
                [t,b]= MY_eyelid_fitting(2,varargin);
                rt = t(1); Ct = t(2:3);
                rb = b(1); Cb = b(2:3);
                
                
                if isnan(rt) || isnan(rb)
                    varargin = struct('masklib',lmask,'imglib',lib,'img',Ig);
                    pmask = MY_eyelid_fitting(1,varargin);
                    varargin = struct('masklib',pmask,'imglib',[],'img',Ig);
                    [t,b]= MY_eyelid_fitting(2,varargin);
                    rt = t(1); Ct = t(2:3);
                    rb = b(1); Cb = b(2:3);
                end
                
                Igf = insertShape(Ig,'circle',[Ct(2),Ct(1),rt],'LineWidth',1,'Color','w');
                Igf = insertShape(Igf,'circle',[Cb(2),Cb(1),rb],'LineWidth',1,'Color','w');
                
                % insertShape eyelids
                Iga = insertShape(zeros(size(Ig)),'FilledCircle',[Ct(2),Ct(1),rt],...
                    'LineWidth',1,'Color','w','Opacity',1);
                Igb = insertShape(zeros(size(Ig)),'FilledCircle',[Cb(2),Cb(1),rb],...
                    'LineWidth',1,'Color','w','Opacity',1);
                Igab = rgb2gray(Iga)+rgb2gray(Igb);
                Igab(Igab~=2)=0;
                BWab = edge(Igab,'canny');
                po=find(mean(Igab)>0);
                BWab(:,1:min(po)+1)=false;
                BWab(:,max(po)-1:end)=false;
                % refresh the eyelid template
                a = insertShape(zeros(size(Ig)),'Circle',[Ct(2),Ct(1),rt],...
                    'LineWidth',1,'Color','w','Opacity',1);
                b = insertShape(zeros(size(Ig)),'Circle',[Cb(2),Cb(1),rb],...
                    'LineWidth',1,'Color','w','Opacity',1);
                eyebi = rgb2gray(a).*BWab; eyebi(eyebi~=0)=1;
                eyeti = rgb2gray(b).*BWab; eyeti(eyeti~=0)=1;
                
                BW = edge(Ig,'canny');
                LightBulb = (Ig>max(Ig(:))*0.7);
                SE = strel('disk',5);
                LB = imdilate(LightBulb,SE);
                SE = strel('disk',2);
                AB = imdilate(BWab,SE);
                BW(Igab~=2 | AB | LB)=0;
                
                %{
                    if kk~=0
                        A = mean(BW);  B = mean(Ig);
                        A(B<max(B)*0.6)=0;
                        [~,LOCS]=findpeaks(A,'NPeaks',2,'MinPeakDistance',80);
                        BW1=false(size(BW));
                        BW1(:,LOCS(1)-10:LOCS(1)+10)= BW(:,LOCS(1)-10:LOCS(1)+10);
                        BW2=BW(:,LOCS(2)-10:LOCS(2)+10);   BW2=bwareaopen(BW2,15);
                        BW1(:,LOCS(2)-10:LOCS(2)+10)=BW2;
                        bw=bwareaopen(BW1,5);
                        [r,C] = MY_circular_fitting_LMS(bw);
                    else
                        SE = strel('line',10,0);
                        BW1 = imdilate(BWprior,SE);
                        bw=BW1.*BW;
                        [r,C] = MY_circular_fitting_LMS(bw);
                        if r<50 || r>70
                            SE = strel('line',20,0);
                            BW1 = imdilate(BWprior,SE);
                            bw=BW1.*BW;
                            [r,C] = MY_circular_fitting_LMS(bw);
                        end
                        if r<50 || r>70
                            A = mean(BW);  B = mean(Ig);
                            A(B<max(B)*0.6)=0;
                            [~,LOCS]=findpeaks(A,'NPeaks',2,'MinPeakDistance',80);
                            BW1=false(size(BW));
                            BW1(:,LOCS(1)-10:LOCS(1)+10)= BW(:,LOCS(1)-10:LOCS(1)+10);
                            BW2=BW(:,LOCS(2)-10:LOCS(2)+10);   BW2=bwareaopen(BW2,15);
                            BW1(:,LOCS(2)-10:LOCS(2)+10)=BW2;
                            bw=bwareaopen(BW1,5);
                            [r,C] = MY_circular_fitting_LMS(bw);
                        end
                    end
                    
                    b = insertShape(Ig,'Circle',[C(2),C(1),r],...
                        'LineWidth',1,'Color','w','Opacity',1);
                    figure;imshowpair(b,double(Ig)+bw*200,'montage');
                    
                    c = insertShape(zeros(size(Ig)),'FilledCircle',[C(2),C(1),r],...
                        'LineWidth',1,'Color','w','Opacity',1);
                    Igabc = Igab+double(rgb2gray(c));
                    Igabc(Igabc~=3)=0;
                    BWprior = edge(Igabc,'canny');
                    Igabd = imerode(Igab,SE);
                    BWprior(Igabd==0)=0;
                %}
                
                [a,b,c]=size(Ig);
                Vec = reshape(Ig,[a*b c]);
                CC = corr(double(Vec),libpVec);
                %if kk~=100 && kk~=NFrames ;CC(p) = CC(p)*1.01;end
                p = find(CC==max(CC));
                pmask = lmaskp(:,:,p);
                eyel = (pmask==1);
                eyer = (pmask==2);
                
                SEd = strel('line',1,90);
                for d=1:3:30
                    SE = strel('line',d,0);
                    eyeld = imdilate(eyel,SE);
                    eyel2 = eyeld.*BW;
%                     eyel2 = imerode(eyel2,SEd);%eyel2 = imdilate(eyel2,SEd);
                    if numel(find(eyel2==1))>=5;eyel=eyel2;break;end
                end
                for d=1:3:50
                    SE = strel('line',d,0);
                    eyerd = imdilate(eyer,SE);
                    eyer2 = eyerd.*BW;
%                     eyer2 = imerode(eyer2,SE);%eyer2 = imdilate(eyer2,SE);
                    if numel(find(eyer2==1))>=5;eyer=eyer2;break;end
                end
                
                
%                 if numel(find(eyel==1))<5;eyel=el;end
%                 if numel(find(eyer==1))<5;eyer=er;end
                
%                 [~,~,para] = hough_circle(BW,0.01,0.1,15,60,0.5);
                
                
                [r,C] = MY_circular_fitting_LMS(eyel+eyer);
                
%                 if r<50 || r>75 || isnan(r);
%                     SE = strel('line',15,0);
%                     eyeld = imdilate((pmask==1),SE); eyerd = imdilate((pmask==2),SE);
%                     eyel = eyeld.*BW; eyer = eyerd.*BW;
%                     SE = strel('line',2,90);
%                     eyel = imerode(eyel,SE);eyel = imdilate(eyel,SE);
%                     eyer = imerode(eyer,SE);eyer = imdilate(eyer,SE);
% %                     if numel(find(eyel==1))<5;eyel=el;end
% %                     if numel(find(eyer==1))<5;eyer=er;end
%                     [r,C] = MY_circular_fitting_LMS(eyel+eyer);
%                 end
                
%                 el=eyel;er=eyer;
%                 if r<45 || r>75 || isnan(r);
%                 Igf = insertShape(Igf,'circle',[C(2),C(1),r],'LineWidth',1,'Color','w');
%                 Igf = insertShape(Igf,'circle',[C(2),C(1),1],'LineWidth',3,'Color','r');     
%                 figure;imshowpair(Igf,double(Ig)+(eyel+eyer)*100,'montage');
%                 end
                
                Radius(kk)=r;
                Center(kk,:)=C;
                Radius_top(kk)=rt;     Center_top(kk,:)=Ct;     
                Radius_bottom(kk)=rb;  Center_bottom(kk,:)=Cb;     
                clear I*
            end
            %}
            %% double check the pupil
            %{
            a = gradient(Radius);
            a(abs(a)>std(a)*2)=nan;
            w = find(isnan(a)|a==0);
            a=w-1;b=w;c=w+1;a(a<1)=[];c(c>NFrames)=[];
            d = union(a,b);wiredFrames=union(c,d);
            for kk = wiredFrames
                I = read(fileinfo,kk);
                Ig = rgb2gray(I);
                h= fspecial('gaussian',[10 10],10);
                Ig = uint8(filter2(h,Ig));
                varargin = struct('masklib',lmask,'imglib',lib,'img',Ig);
                pmask = MY_eyelid_fitting(1,varargin);
                varargin = struct('masklib',pmask,'imglib',[],'img',Ig);
                [t,b]= MY_eyelid_fitting(2,varargin);
                rt = t(1); Ct = t(2:3);
                rb = b(1); Cb = b(2:3);
                
                Iga = insertShape(zeros(size(Ig)),'FilledCircle',[Ct(2),Ct(1),rt],...
                    'LineWidth',1,'Color','w','Opacity',1);
                Igb = insertShape(zeros(size(Ig)),'FilledCircle',[Cb(2),Cb(1),rb],...
                    'LineWidth',1,'Color','w','Opacity',1);
                
                Igab = rgb2gray(Iga)+rgb2gray(Igb);
                Igab(Igab~=2)=0;
                BWab = edge(Igab,'canny');
                
                BW = edge(Ig,'canny');
                LightBulb = (Ig>max(Ig(:))*0.7);
                SE = strel('disk',10);
                LB = imdilate(LightBulb,SE);
                %BWab = edge(Igab,'canny');
                SE = strel('disk',2+3);
                AB = imdilate(BWab,SE);
                BW(Igab~=2 | AB | LB)=0;
                
                
                [a,b,c]=size(Ig);
                Vec = reshape(Ig,[a*b c]);
                CC = corr(double(Vec),libpVec);
                if p~=0;CC(p) = CC(p)*1.01;end
                p = find(CC==max(CC));
                pmask = lmaskp(:,:,p);
                eyel = (pmask==1);
                eyer = (pmask==2);
                
                for d=1:20
                    SE = strel('line',d,0);
                    eyeld = imdilate(eyel,SE);
                    eyel2 = eyeld.*BW;
                    if numel(find(eyel2==1))>=5;eyel=eyel2;break;end
                end
                for d=1:20
                    SE = strel('line',d,0);
                    eyerd = imdilate(eyer,SE);
                    eyer2 = eyerd.*BW;
                    if numel(find(eyer2==1))>=10;eyer=eyer2;break;end
                end
                
                
                [r,C] = MY_circular_fitting_LMS(eyel+eyer);
                
                Radius(kk)=r;
                Center(kk,:)=C;
                Radius_top(kk)=rt;     Center_top(kk,:)=Ct;     
                Radius_bottom(kk)=rb;  Center_bottom(kk,:)=Cb;     
            end
            
            %}

            Radius1 = MY_deNAN_despike_TimeSer(Radius',20);
            Center1(:,1) = MY_deNAN_despike_TimeSer(Center(:,1),10);
            Center1(:,2) = MY_deNAN_despike_TimeSer(Center(:,2),10);
            Radius_top1 = MY_deNAN_despike_TimeSer(Radius_top',50);
            Center_top1(:,1) = MY_deNAN_despike_TimeSer(Center_top(:,1),50);
            Center_top1(:,2) = MY_deNAN_despike_TimeSer(Center_top(:,2),30);
            Radius_bottom1 = MY_deNAN_despike_TimeSer(Radius_bottom',20);
            Center_bottom1(:,1) = MY_deNAN_despike_TimeSer(Center_bottom(:,1),30);
            Center_bottom1(:,2) = MY_deNAN_despike_TimeSer(Center_bottom(:,2),10);
            

            % OutExcel
            Exc = fullfile(AVIPath{jj},'Radius_Center.xlsx');
            xlswrite(Exc,{'R_pupil','C_pupil01','C_pupil02',...
                'R_eyelid_top','C_eyelid_top01','C_eyelid_top02',...
                'R_eyelid_bottom','C_eyelid_bottom01','C_eyelid_bottom02'},...
                ['Video',num2str(idx)],'A1');
            NUM = [Radius1,Center1,Radius_top1,Center_top1,Radius_bottom1,Center_bottom1];
            xlswrite(Exc,NUM,['Video',num2str(idx)],'A2');

            
            w = VideoWriter(['Check0_',AVIdir(idx).name]);
            w.FrameRate = 10;
            open(w);
                        
            for f = 1:floor(NFrames)
                Ig = rgb2gray(read(fileinfo,f));
                Igf = insertShape(Ig,'circle',[Center1(f,2),Center1(f,1),Radius1(f)],'LineWidth',1,'Color','w');
                Igf = insertShape(Igf,'circle',[Center1(f,2),Center1(f,1),1],'LineWidth',3,'Color','r');     
                Igf = insertText(Igf,[0 0],['Frame ',num2str(f,'%04d')],'FontSize',12,'AnchorPoint','LeftTop');
                
                %********************
                Iga = insertShape(zeros(size(Ig)),'FilledCircle',[Center_top1(f,2),Center_top1(f,1),Radius_top1(f)],...
                    'LineWidth',1,'Color','w','Opacity',1);
                Igb = insertShape(zeros(size(Ig)),'FilledCircle',[Center_bottom1(f,2),Center_bottom1(f,1),Radius_bottom1(f)],...
                    'LineWidth',1,'Color','w','Opacity',1);
                Igab = rgb2gray(Iga)+rgb2gray(Igb);
                Igab(Igab~=2)=0;
                BWab = edge(Igab,'canny');
                eye = repmat(BWab+anchor,[1 1 3]);
                eye(:,:,1:2)=0;
                %*********************
                
                Vid = Igf+uint8(eye*200);
                writeVideo(w, Vid);
            end
            close(w);
            fclose all;

            clear I*;
            
        end
        
    end
end
%}











%% denoise pupil fitting
%{
for iii = 2% 1:numel(EyeTrackDir)
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
    
    for jj = 1%1:numel(AVIPath)
        cd(AVIPath{jj})
        AVIdir = dir('de_*.avi');
        
        for idx = 1%1:numel(AVIdir)
            Radius = load(fullfile(AVIPath{jj},['Radius_',num2str(idx),'.txt']));
            Center = load(fullfile(AVIPath{jj},['Center_',num2str(idx),'.txt']));
            
            Radius(Radius==0)=nan;
            Center(Center==0)=nan;
            Radius=smooth(Radius,5);
            Center(:,1)=smooth(Center(:,1),5);
            Center(:,2)=smooth(Center(:,2),5);
            
            fileinfo = VideoReader(AVIdir(idx).name);
            NFrames = fileinfo.NumberOfFrames;
            
            
            for kk = 1:6180%NFrames
                I = read(fileinfo,kk);
                if kk==1;Igc=uint8(zeros([size(I) NFrames]));end
                Igc(:,:,:,kk) = I;
                
                Igf = insertShape(I,'circle',[Center(kk,2),Center(kk,1),Radius(kk)],'LineWidth',1,'Color','w');
                Igf = insertShape(Igf,'circle',[Center(kk,2),Center(kk,1),1],'LineWidth',3,'Color','r');
                Igc(:,:,:,kk) = Igf;
            end
            
            w = VideoWriter(['Check2_',AVIdir(idx).name]);
            w.FrameRate = 10;
            open(w);

            for kk = 1:6180%floor(NFrames)
                writeVideo(w, Igc(:,:,:,kk));
            end
            close(w);
            fclose all;
            save(fullfile(AVIPath{jj},['Radius2_',num2str(idx),'.txt']),'Radius','-ascii');
            save(fullfile(AVIPath{jj},['Center2_',num2str(idx),'.txt']),'Center','-ascii');
            clear I*;
        end
    end
end
%}
%}

%% eye lid open level
for iii = 26:numel(EyeTrackDir)
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

        for idx = 4:numel(AVIdir)
            fileinfo = VideoReader(AVIdir(idx).name);
            NFrames = fileinfo.NumberOfFrames;
            I = read(fileinfo,[1 NFrames]);
            [a,b,c,d] = size(I);
            Iv = reshape(I,[a*b*c d]);
            Ivm = mean(Iv,1);
            
%             [b,a]=butter(3,0.001,'high');
%             OpenLevel=filtfilt(b,a,Ivm);
            
            Excel = 'OpenLevel.xlsx';
            xlswrite(Excel,Ivm(:),'sheet1',[char(double('A')+idx-1),'2']);
        end
        
    end
end