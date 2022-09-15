function MY_downsample_and_density_receptor(pic_dir,genename)

    files_ = dir(fullfile(pic_dir,'RawPic',[genename,'*']));
    bar = waitbar(0,'loading...'); 
    for files_loop = 1:size(files_,1)
        filename = fullfile(files_(files_loop).folder,files_(files_loop).name); 
        
        I = rgb2gray(imread(filename));
        % remove watermarking
        [x,y]=size(I);
        a=1;b=round(y*0.15);a_=round(y*0.8);b_=y;
        c=round(x*0.85);d=x;
        X = I(c:d,a:b);I(c:d,a:b) = median(X(:));
        X = I(c:d,a_:b_);I(c:d,a_:b_) = median(X(:));
        
        Is = imresize(I,0.1);
        BW=edge(I,'approxcanny');
        I_= mat2cell(BW,100*ones(1,size(I,1)/100),100*ones(1,size(I,2)/100));
        BW_ = cell2mat(cellfun( @(x)(mean(x,'all')),I_,'UniformOutput',false));
        BW__ = imresize(BW_,10);
        
        bw = BW__>0.1;
        bw_ = imfill(bw,'holes');
        bw__ = bwareaopen(bw_,round(numel(find(bw_==true))/4));
        
        low_in = mean(double(Is(bw_))) - 2*std(double(Is(bw_)));
        high_in = mean(double(Is(bw_))) + 2*std(double(Is(bw_)));
        J = imadjust(Is,[max(low_in,0) min(high_in,255)]/255,[0 1]);
        
       
        %
        
        I0=I;
        BW = imresize(bw__,10);
        %I(~BW)=0;
        Iall = mean(double(I(BW)));
        I_= mat2cell(I,100*ones(1,size(I,1)/100),100*ones(1,size(I,2)/100));
        C = cell2mat(cellfun( @(x)(mean(x(:))/Iall),I_,'UniformOutput',false));
        C_ = imresize(C,10);
        D = 1-C_; 
        D(D<0)=0;
        
        D = imresize(D,[720,840]);
        If = imresize(255-J,[720,840]);
        
        if files_loop == 1
            Ds = zeros([size(D),size(files_,1)]);
            Ifs = zeros([size(If),size(files_,1)]);
        end
        Ds(:,:,files_loop)=D;
        Ifs(:,:,files_loop)=If;
        
        str=['finish ...',num2str(100*files_loop/size(files_,1),'%0.1f'),'%'];   
        waitbar(files_loop/size(files_,1),bar,str)
        
    end
    close(bar)
    dest = fullfile(pic_dir,'downs_pic');
    if ~exist(dest,'dir');mkdir(dest);end
    for loop = 1:size(Ifs,3);X=uint8(Ifs(:,:,loop));imwrite(X,fullfile(dest,files_(loop).name));end
    
    dest = fullfile(pic_dir,'receptor_density');
    if ~exist(dest,'dir');mkdir(dest);end
    for loop = 1:size(Ds,3);X=uint8(Ds(:,:,loop)*255);imwrite(X,fullfile(dest,files_(loop).name));end

end