clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));



age_band = {[0,36];[36,72];[72,120]};

for ageloop = 1:3
    
    AGEBAND = age_band{ageloop};
    
bins = 10;
for bl=1:bins
    a_=0; CUM_ = 0;
    NUM=5522;
    %Vec_CUM=zeros([NUM,NUM,63]);
    Vec = zeros([NUM,NUM,270]);
    wNUM=zeros(1,345);
    
    
    %% ION 
    ExpRecording = [codepath 'ExpIndex.xlsx'];
    [~,~,CellData] = xlsread(ExpRecording);
    ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
    WholeMaskNII = fullfile(codepath,'marmoset_atlas','Trange_Template_mask_V3.nii');
    lmask = spm_read_vols(spm_vol(WholeMaskNII));


    AtlasNII = fullfile(codepath,'marmoset_atlas','Trange_atlas_RikenBMA_cortex.nii');
    AtlasMask = spm_read_vols_4D(spm_vol(AtlasNII));
    %

    [~,~,CellDATA]= xlsread(fullfile(codepath,'Age_Gender_list.xlsx'),'ION');
    GenderTable = cell2table(CellDATA(2:end,:),'VariableNames',CellDATA(1,:));

    a=65;b=86;c=64;d=512; % 4D-EPI size
    lmask = AtlasMask>20;
    Xmask = imresize3(lmask,round([a,b,c]/2),'nearest');
    NUM=5522;
    %


    for idx = 1:62

        path = [ExpTable.Monkey{idx},filesep];
        Age = GenderTable.age(idx);
        
        if Age<=AGEBAND(2) & Age>AGEBAND(1)
            Session = replace(path,WholePath,'');
            RARE =  ExpTable.RARE(idx);
            SEEPI = sort([ExpTable.SEEPI1(idx),ExpTable.SEEPI2(idx)],'ascend');
            GEEPI = sort([  ExpTable.EPI4D1(idx),ExpTable.EPI4D2(idx),...
                ExpTable.EPI4D3(idx),ExpTable.EPI4D4(idx),...
                ExpTable.EPI4D5(idx),ExpTable.EPI4D6(idx),...
                ExpTable.EPI4D7(idx),ExpTable.EPI4D8(idx)],'ascend');
            GEEPI(isnan(GEEPI)) = [];


            for kk = GEEPI



                GEPath = fullfile(WholePath,'DiffusionGradient_ION',Session,num2str(kk));
                SWname = fullfile(GEPath,'SW15_VIweighted.nii');
                hdr = spm_vol(SWname);
                FC = spm_slice_vol(hdr,spm_matrix([0 0 bl]),hdr.dim(1:2),0);
                x = str2num(hdr.descrip);
                swNUM = x(bl);

                a_=a_+1;
                wNUM(a_) = swNUM;
                Vec(:,:,a_)=FC;
                

                
            end
        end
    end

    %% NIH
    [~,~,CellDATA]= xlsread(fullfile(codepath,'Age_Gender_list.xlsx'),'NIH');
    GenderTable = cell2table(CellDATA(2:end,:),'VariableNames',CellDATA(1,:));
   
    NIH_path = 'E:\Marmoset_NIH_dataset';
    Session = dir(fullfile(NIH_path,'*_session_*'));

    for idx = 1:26


        path = fullfile(NIH_path,Session(idx).name);
        Age = GenderTable.age(idx);
        
        if  Age<=AGEBAND(2) & Age>AGEBAND(1)

            for ii=1:99;if ~exist(fullfile(path,'Results',num2str(ii)),'dir');break;end;end
            ScanNum = 1:ii-1;

            for nl = ScanNum

                fprintf(['marmoset: ',num2str(idx),';scan: ',num2str(nl),'\n'])

                %


                
                GEPath = fullfile(WholePath,'DiffusionGradient_NIH',Session(idx).name,num2str(nl,'%02d'));
                SWname = fullfile(GEPath,'SW15_VIweighted.nii');
                hdr = spm_vol(SWname);
                FC = spm_slice_vol(hdr,spm_matrix([0 0 bl]),hdr.dim(1:2),0);
                x = str2num(hdr.descrip);
                swNUM = x(bl);
                
                a_=a_+1;
                Vec(:,:,a_)=FC;
                wNUM(a_)=swNUM;
                
                

            end
        end
    end
    
    
    %% Summary
    Vec(:,:,a_+1:end)=[];
    wNUM(a_+1:end)=[];
    FCz = nansum(Vec,3)/sum(wNUM); 
    clear Vec Vec_CUM

    load(fullfile(WholePath,'DiffusionGradient_ION','Mean','Gradients.mat'));
    Ref = Gd;

    NANidx = isnan(nanmean(FCz,1));
    G = GradientMaps('kernel','cs','approach','dm','n_components',20,'align','pa');
    G = G.fit(FCz(~NANidx,~NANidx),'reference',Ref);

    Gd = G.aligned{1};
    Gx = zeros(numel(NANidx),size(Gd,2));
    Gx(~NANidx,:) = Gd;

    M = funmask(Gd,Xmask);
    xhdr=spm_vol(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.nii.gz'));
    M_ = zeros([xhdr.dim,size(M,4)]);
    for m=1:size(M,4)
        X = M(:,:,:,m);
        X = imresize3(X,xhdr.dim,'linear');
        X = smooth3(X,'box',7);
        M_(:,:,:,m) = X;
    end

    dest_ = fullfile(WholePath,'DiffusionGradient_merge',['Age_',num2str(ageloop)],'Dynamic_VIweighted',num2str(bl,'%02d'));
    if ~exist(dest_,'dir');mkdir(dest_);end
    save(fullfile(dest_,'G.mat'),'G');


    filename = fullfile(dest_,'PG_align.nii');
    clear Vol
    for k = 1:size(M,4)
        Vol(k,1) = struct(  'fname',    filename,...
                            'dim',      xhdr.dim,...
                            'mat',      xhdr.mat,...
                            'n',        [k,1],...
                            'pinfo',    [1;0;0],...
                            'descrip',  '',...
                            'dt',       [16 0]);
    end
    spm_write_vol_4D(Vol,M_); 

    delete(fullfile(dest_,'fx.func.gii'));
    bar = [-1 -0.0001 0.0001 1]*5;
    for k = 1:5
        F = SurfaceMap(fullfile(dest_,'PG_align.nii'),bar/k,codepath,k,'marmoset');
        TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'.tiff']);
        saveas(F,TIFname)
        close all
        X = imread(fullfile(dest_,['Surface_PG_align_',num2str(k),'.tiff']));
        x_ = X(61:400,65:680,:);
        y_ = X(531:870,81:700,:);
        z = cat(2,x_,y_);
        TIFname = fullfile(dest_,['Surface_PG_align_',num2str(k),'_cut.tiff']);
        imwrite(z,TIFname);
        pause(3);
    end
    

end
end