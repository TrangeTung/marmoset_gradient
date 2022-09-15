clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
addpath(genpath(codepath));

%% topup preprocessing
%{ 
NIH_path = 'H:\Marmoset_connectome_NIH\';
Name_path = fullfile('E:\rsfMRI_marmoset\DiffusionGradient_NIH');
Name = loadtxt(fullfile(Name_path,'list.txt'));

WholePath = 'H:\Marmoset_NIH_dataset';

a_=0;
for idx = 1:26
    
    
    path = fullfile(NIH_path,Name{idx});
    LIST = loadtxt(fullfile(path,'list_fmri'));
    
    
    for ii=1:numel(LIST)
        session = [Name{idx},'_session_',num2str(ii,'%02d')];
        dest = fullfile(WholePath,session,'Results');
        if ~exist(dest,'dir');mkdir(dest);end
        
        a_=a_+1;
        %% GEEPI
        ROOT = dir(fullfile(path,LIST{ii},'BOLD_*_*.nii.gz'));
        for rl=1:numel(ROOT)
            dest_ = fullfile(dest,num2str(rl));
            if ~exist(dest_,'dir');mkdir(dest_);end
            copyfile(fullfile(path,LIST{ii},ROOT(rl).name),dest_);
            cd(dest_); eval(['!rename,',ROOT(rl).name,',2dseq.nii.gz;']);
            
            if contains(ROOT(rl).name,'up');PEdirec=1;end
            if contains(ROOT(rl).name,'down');PEdirec=-1;end
            
            pars = repmat([PEdirec 0 0 1],512,1);
            
            filename = fullfile(dest_,'GE.txt');
            dlmwrite(filename, num2str(pars), 'delimiter', '', 'precision', 6);
        end
        
        %% SEEPI
        ihdr1 = spm_vol(fullfile(path,LIST{ii},'SEEPI_up.nii.gz'));
        img1 = spm_read_vols(ihdr1);
        ihdr2 = spm_vol(fullfile(path,LIST{ii},'SEEPI_down.nii.gz'));
        img2 = spm_read_vols(ihdr2);
        
        ihdr = [ihdr1;ihdr2];
        img = cat(4,img1,img2);
        for ix=1:numel(ihdr)
            ihdr(ix).fname=fullfile(dest,'SEEPI.nii');
            ihdr(ix).n=[ix,1];
        end
        spm_write_vol_4D(ihdr,img);
        
        acqpars_up = repmat([1 0 0 1],8,1);
        acqpars_down = repmat([-1 0 0 1],8,1);
        
        pars = [acqpars_up;acqpars_down];
        filename = fullfile(dest,'SE.txt');
        dlmwrite(filename, num2str(pars), 'delimiter', '', 'precision', 6);
        
        %% T2
        dest_ = fullfile(dest,'T2');
        if ~exist(dest_,'dir');mkdir(dest_);end
        copyfile(fullfile(path,LIST{ii},'InplaneT2.nii.gz'),dest_);
        cd(dest_); eval(['!rename,','InplaneT2.nii.gz',',2dseq.nii.gz;']);
        
        folder_txt = fullfile(WholePath,'folder.txt');
        fid=fopen(folder_txt,'a+');
        fprintf(fid,'%s\n',['Path_',num2str(a_),'="',session,'";']);
        fprintf(fid,'%s\n',['Folder_',num2str(a_),'=(',num2str(1:numel(ROOT)),');']);
        fclose(fid);
    end
    
    
end
%}



%% spatial normalize
NIH_path = 'E:\Marmoset_NIH_dataset';
Session = dir(fullfile(NIH_path,'*_session_*'));

for idx = 1:39
    
    
    path = fullfile(NIH_path,Session(idx).name);
    for ii=1:99;if ~exist(fullfile(path,'Results',num2str(ii)),'dir');break;end;end
    GEEPI = 1:ii-1;
    
    
    fclose('all');
    spm('defaults', 'FMRI');
    set(spm('CreateIntWin','off'),'Visible','on');
    
    for flag_stage = 4:10
        
        if flag_stage == 4
            %% mask all EPI and RARE
            cd([path '\Results']);
            EPI_mask = spm_read_vols(spm_vol('EPI_mask.nii'));
            T2_mask = spm_read_vols(spm_vol('T2_mask.nii'));
            MY_mask_images(path,GEEPI,'Uw2dseq.nii.gz',EPI_mask,'mUw2dseq.nii','EPI');
            MY_mask_images(path,1,'2dseq.nii.gz',T2_mask,'m2dseq.nii','T2');
            
        end
        if flag_stage == 5
            %% Slicetiming
            for kk = GEEPI
                EPI_TR = 2;
                Nslice = 38;
                all_func = MY_find_images_in_all_scans(path,'Results',{kk},'^mUw2dseq','.nii',[1 Inf],'separate_cells');
                slicetiming_mlb = MY_get_default_slicetiming_batch_struct(all_func,Nslice,EPI_TR);
                disp('Start to process Slicetiming!')
                spm_jobman('run',slicetiming_mlb);
            end
        end
        if flag_stage == 6
            %% Realignment
            all_func = MY_find_images_in_all_scans(path,'Results',{GEEPI(:)},'^smUw2dseq','.nii',[1 Inf],'separate_cells');
            realign_mlb = MY_get_default_realign_batch_struct(all_func);
            F = spm_figure('GetWin');
            disp('Start to process realignment !')
            spm_jobman('run',realign_mlb);
            hgexport(figure(F), fullfile([path,'\Results\'],strcat('realign')), hgexport('factorystyle'), 'Format', 'tiff');
            clear realign_mlb all_func;
        end
    
        
        if flag_stage == 7
            %% Head Motion removal
            lmask = spm_read_vols(spm_vol(fullfile(path,'Results','EPI_mask.nii')));
            for kk = GEEPI
                GEPath = fullfile(path,'Results',num2str(kk));
                NIIfile = fullfile(GEPath,'rsmUw2dseq.nii');
                head = spm_vol(NIIfile);
                fMRI_data = spm_read_vols(head);
                data_ready_regress = fmask(fMRI_data,lmask);
                
                disp('Start to process Head Motion removal !')
                rp = load(fullfile(GEPath,'rp_smUw2dseq.txt'));
                PCs = load(fullfile(GEPath,'PCs.txt'));
                [~,drp] = gradient(rp);
                RegBasFunc = [ones([length(rp),1]) rp drp PCs(:,1)];
                
                for iii = 1:size(data_ready_regress,1)
                    [Beta,~,Residual] = regress(double(data_ready_regress(iii,:))',RegBasFunc);
                    data_ready_regress(iii,:) = Residual + Beta(1);
                end
                
                fMRI_data = funmask(data_ready_regress,lmask);
                clear data_ready_regress
                
                %% fMRI after band pass filter
                disp('Start to process band pass filtering !')
                EPI_TR = 2;
                fMRI_data = Smooth_temporal_f(fMRI_data,lmask, EPI_TR);
                for k = 1:size(fMRI_data,4)
                    Vol(k,1) = struct(  'fname',    fullfile(GEPath,'RrsmUw2dseq.nii'),...
                        'dim',      double(head(1).dim),...
                        'mat',      head(1).mat,...
                        'n',        [k,1],...
                        'pinfo',    [1;0;0],...
                        'descrip',  ['clean'],...
                        'dt',       [16 0]);
                end
                spm_write_vol_4D(Vol,fMRI_data);
                clear fMRI_data
            end
        end
       if flag_stage == 8
            %% Func2T2 Coregistration
            for kk = GEEPI
                source{1,1} = [path '\Results\' num2str(kk) '\RrsmUw2dseq.nii,1'];
                ref{1,1} = [path '\Results\T2\m2dseq.nii,1'];
                all_func = MY_find_images_in_all_scans(path,'Results',{kk},'^RrsmUw2dseq','.nii',[1 Inf],'all_mixed');
                Func2T2W_mlb = MY_get_default_coreg_batch_struct(ref, source, all_func);
                disp('Start to process Func2T2W coregistration!');
                F = spm_figure('GetWin');
                spm_jobman('run',Func2T2W_mlb);
                hgexport(figure(F), fullfile([path '\Results\' num2str(kk)], 'coreg'), hgexport('factorystyle'), 'Format', 'tiff');
            end
            clear Func2T2W_mlb other ref source;
        end
        if flag_stage == 9
            %% T22Template coregistration
            ref{1,1} = [codepath '\marmoset_atlas\Trange_Template_sym_T2_6X_V3.nii,1'];
            source{1,1} = [path '\Results\T2\m2dseq.nii,1'];
            all_func = MY_find_images_in_all_scans(path,'Results',{GEEPI(:)},'^cRrsmUw2dseq','.nii',[1 Inf],'all_mixed');
            all_func = [all_func;[path '\Results\T2\m2dseq.nii,1']];
            OldNormalize_mlb = MY_get_default_oldnormalize_batch_struct(ref, source, all_func);
            disp('Start to process OldNormalize!');
            F = spm_figure('GetWin');
            spm_jobman('run',OldNormalize_mlb);
            hgexport(figure(F), fullfile([path '\Results\'], strcat('oldnormalize')), hgexport('factorystyle'), 'Format', 'tiff');
            
            %*************** delete this part after this project***********
            Check_func = char(ref{1,1},[path '\Results\T2\nm2dseq.nii,1']);
            for kk = GEEPI
                delete(fullfile(path,'Results',num2str(kk),'cRrsmUw2dseq.nii'))
                Check_func = char(Check_func,fullfile(path,'Results',num2str(kk),'ncRrsmUw2dseq.nii,1'));
            end
            F = spm_figure('GetWin');
            spm_check_registration(Check_func);
            hgexport(figure(F), fullfile([path '\Results\'], strcat('Coreg_Check')), hgexport('factorystyle'), 'Format', 'tiff');
            %**************************************************************
        end
        if flag_stage == 10
            %% smooth_space
            all_func = MY_find_images_in_all_scans(path,'Results',{GEEPI(:)},'^ncRrsmUw2dseq','.nii',[1 Inf],'all_mixed');
            Smooth_mlb = MY_get_default_smooth_batch_struct(all_func);
            disp('Start to process Smooth!');
            spm_jobman('run',Smooth_mlb);
            clear Smooth_mlb;
            
        end

        
    end
end