%% preparation
% According to your experimental parameters and certain file direction, Change these bellow.
clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));
ExpRecording = [codepath 'ExpIndex.xlsx'];
[~,~,CellData] = xlsread(ExpRecording);
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
WholeMaskNII = fullfile(codepath,'marmoset_atlas','Trange_Template_mask_V3.nii');


%% regressor choices
Reg_choices = {'rp';};

%
for idx = 3:15
    idx
    path = [ExpTable.Monkey{idx},filesep];
    RARE =  ExpTable.RARE(idx);
    SEEPI = sort([ExpTable.SEEPI1(idx),ExpTable.SEEPI2(idx)],'ascend');
    GEEPI = sort([  ExpTable.EPI4D1(idx),ExpTable.EPI4D2(idx),...
        ExpTable.EPI4D3(idx),ExpTable.EPI4D4(idx),...
        ExpTable.EPI4D5(idx),ExpTable.EPI4D6(idx),...
        ExpTable.EPI4D7(idx),ExpTable.EPI4D8(idx)],'ascend');
    GEEPI(isnan(GEEPI)) = [];
        
    fclose('all');
    spm('defaults', 'FMRI');
    set(spm('CreateIntWin','off'),'Visible','on');
    
    for flag_stage = 4:10
        if flag_stage == 1
            %% Bruker2nifiti
            Bruker2nifti_multislice(path,RARE,'marmoset');
            Bruker2nifti_multislice(path,SEEPI,'marmoset');
            Bruker2nifti_multislice(path,GEEPI,'marmoset');
            % *************************************************************
            cd([path,filesep,'Results'])
            eval(['!rename,',num2str(SEEPI(1)),[',SE',num2str(SEEPI(1))]]);
            eval(['!rename,',num2str(SEEPI(2)),[',SE',num2str(SEEPI(2))]]);
            % Send to BaiDuCloud_folder
            BaiDuCloud_Path = fullfile(WholePath,'BaiDuCloud_folder',[path(numel(WholePath)+1:end)]);
            mkdir(BaiDuCloud_Path);
            copyfile(fullfile(path,'Results','T2','2dseq.nii'),BaiDuCloud_Path);
            % *************************************************************
        end
        if flag_stage == 2
            %% GEPI //// Pre-set the FSL topup parameters
            %{
            a=0;b=0;
            for kk = GEEPI
                pars = MY_get_TopUp_pars(path,kk,'');
                NIIfile = fullfile(path,'Results',num2str(kk),'2dseq.nii');
                if pars(1,1) == 1
                    pars_1 = pars(1:8,:);a = 1;
                    head_1 = spm_vol(NIIfile); Img_1 = spm_read_vols(head_1);
                else
                    pars_0 = pars(1:8,:);b = 1;
                    head_0 = spm_vol(NIIfile); Img_0 = spm_read_vols(head_0);
                end
                if (a==1 && b==1);break; end
            end
            pars = [pars_1;pars_0];
            Img(:,:,:,1:8) = Img_1(:,:,:,1:8);
            Img(:,:,:,9:16) = Img_0(:,:,:,1:8);
            
            for k = 1:size(Img,4)
                Vol(k,1) = struct(  'fname',    fullfile(path,'Results','GEEPI.nii'),...
                    'dim',      double(head_1(1).dim),...
                    'mat',      head_1(1).mat,...
                    'n',        [k,1],...
                    'pinfo',    [1;0;0],...
                    'descrip',  ['GEEPI'],...
                    'dt',       [4 0]);
            end
            spm_write_vol_4D(Vol,Img);
            
            
            filename = fullfile(path,'Results','GE.txt');
            dlmwrite(filename, num2str(pars), 'delimiter', '', 'precision', 6);
            
            
            FSL_Path = fullfile(WholePath,'FSL_folder',[path(numel(WholePath)+1:end)]);
            mkdir(FSL_Path);
            copyfile(fullfile(path,'Results','GEEPI.nii'),FSL_Path);
            copyfile(fullfile(path,'Results','GE.txt'),FSL_Path);
            %}
            
            %% SEPI //// Pre-set the FSL topup parameters
            %
            a=0;b=0;
            for kk = SEEPI
                pars = MY_get_TopUp_pars(path,kk,'SE');
                NIIfile = fullfile(path,'Results',['SE',num2str(kk)],'2dseq.nii');
                if pars(1,1) == 1
                    pars_1 = pars(1:8,:);a = 1;
                    head_1 = spm_vol(NIIfile); Img_1 = spm_read_vols(head_1);
                else
                    pars_0 = pars(1:8,:);b = 1;
                    head_0 = spm_vol(NIIfile); Img_0 = spm_read_vols(head_0);
                end
                if (a==1 && b==1);break; end
            end
            pars = [pars_1;pars_0];
            Img(:,:,:,1:8) = Img_1(:,:,:,1:8);
            Img(:,:,:,9:16) = Img_0(:,:,:,1:8);
            
            for k = 1:size(Img,4)
                Vol(k,1) = struct(  'fname',    fullfile(path,'Results','SEEPI.nii'),...
                    'dim',      double(head_1(1).dim),...
                    'mat',      head_1(1).mat,...
                    'n',        [k,1],...
                    'pinfo',    [1;0;0],...
                    'descrip',  ['SEEPI'],...
                    'dt',       [4 0]);
            end
            spm_write_vol_4D(Vol,Img);
            
            
            filename = fullfile(path,'Results','SE.txt');
            dlmwrite(filename, num2str(pars), 'delimiter', '', 'precision', 6);
            
            
            FSL_Path = fullfile(WholePath,'FSL_folder',[path(numel(WholePath)+1:end)]);
            mkdir(FSL_Path);
            copyfile(fullfile(path,'Results','SEEPI.nii'),FSL_Path);
            copyfile(fullfile(path,'Results','SE.txt'),FSL_Path);
            %}
            
            %% GEEPI
            for kk = GEEPI
                pars = MY_get_TopUp_pars(path,kk,'');
                filename = fullfile(path,'Results',num2str(kk),'acqparsGE.txt');
                dlmwrite(filename, num2str(pars), 'delimiter', '', 'precision', 6);
                
                copyfile(fullfile(path,'Results',num2str(kk),'2dseq.nii'),FSL_Path);
                copyfile(fullfile(path,'Results',num2str(kk),'acqparsGE.txt'),FSL_Path);
                cd(FSL_Path);
                eval(['!rename,2dseq.nii,2dseq_',num2str(kk),'.nii']);
                eval(['!rename,acqparsGE.txt,acqparsGE_',num2str(kk),'.txt']);
            end
        end
        if flag_stage == 3
            %% Unzip .gz file and rename Uw2dseq.nii
            FSL_TopUp_Path = fullfile(WholePath,'FSL_folder',[path(numel(WholePath)+1:end)]);
            Unzip_path = fullfile(path,'Results');
            % SEEPI
            gunzip(fullfile(FSL_TopUp_Path,'Trange_unwrapped_images.nii.gz'),Unzip_path);
            cd(Unzip_path);
            eval(['!rename,Trange_unwrapped_images.nii,Uw_SEEPI.nii']);
            % GEEPI
            for kk = GEEPI
                Unzip_path = fullfile(path,'Results',num2str(kk));
                % delete(fullfile(Unzip_path,'Uw2dseq.nii'));
                gunzip(fullfile(FSL_TopUp_Path,['Trange_2dseq_',num2str(kk),'.nii.gz']),Unzip_path);
                cd(Unzip_path);
                eval(['!rename,',['Trange_2dseq_',num2str(kk),'.nii'],',Uw2dseq.nii']);
            end
        end
        if flag_stage == 4
            %% mask all EPI and RARE
            cd([path 'Results']);
            EPI_mask = spm_read_vols(spm_vol('EPI_mask.nii'));
            T2_mask = spm_read_vols(spm_vol('T2_mask.nii'));
            MY_mask_images(path,GEEPI,'Uw2dseq.nii',EPI_mask,'mUw2dseq.nii','EPI');
            MY_mask_images(path,RARE,'2dseq.nii',T2_mask,'m2dseq.nii','T2');
            
        end
        if flag_stage == 5
            %% Slicetiming
            for kk = GEEPI
                Segments = MY_search_bruker_method('Segments',kk,path);
                EPI_TR = MY_search_bruker_method('EPI_TR',kk,path)/1000*Segments;
                Nslice = MY_search_bruker_method('Nslice',kk,path);
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
            hgexport(figure(F), fullfile([path,'Results\'],strcat('realign')), hgexport('factorystyle'), 'Format', 'tiff');
            clear realign_mlb all_func;
        end
        if flag_stage == 7
            %% RETROICOR RETROICOR_ prefix
            BW = spm_read_vols_4D(spm_vol(fullfile(path,'Results','EPI_mask.nii')));
            for kk = GEEPI
                disp('Start to process RETROICOR !')
                % RETROICOR
                GEPath = fullfile(path,'Results',num2str(kk));
                NIIfilename = fullfile(GEPath,'rsmUw2dseq.nii');
                TR = MY_search_bruker_method('EPI_TR',kk,path);
                Phyfilename = fullfile(path,'Physiology',num2str(kk),[num2str(kk),'.lvm']);
                MY_RETROICOR_correction(NIIfilename,Phyfilename,BW,TR);

            end
            
        end
        
        if flag_stage == 8
            %% Head Motion removal
            lmask = spm_read_vols(spm_vol(fullfile(path,'Results','EPI_mask.nii')));
            for kk = GEEPI
                GEPath = fullfile(path,'Results',num2str(kk));
                NIIfile = fullfile(GEPath,'GrsmUw2dseq.nii');
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
                Segments = MY_search_bruker_method('Segments',kk,path);
                EPI_TR = MY_search_bruker_method('EPI_TR',kk,path)/1000*Segments;
                fMRI_data = Smooth_temporal_f(fMRI_data,lmask, EPI_TR);
                for k = 1:size(fMRI_data,4)
                    Vol(k,1) = struct(  'fname',    fullfile(GEPath,'RGrsmUw2dseq.nii'),...
                        'dim',      double(head(1).dim),...
                        'mat',      head(1).mat,...
                        'n',        [k,1],...
                        'pinfo',    [1;0;0],...
                        'descrip',  ['clean'],...
                        'dt',       [4 0]);
                end
                spm_write_vol_4D(Vol,fMRI_data);
                clear fMRI_data
            end
        end
        
        if flag_stage == 9
            %% double slice number to get off the wrapped situation
            MY_Double_Slice_with_Empty([path,'Results'],'m2dseq.nii','T2','rare');
            MY_Double_Slice_with_Empty([path,'Results'],'RGrsmUw2dseq.nii',GEEPI,'epi');
            %% T22Template coregistration
            ref{1,1} = [codepath 'marmoset_atlas\Trange_Template_sym_T2_6X_V3.nii,1'];
            source{1,1} = [path 'Results\T2\Sm2dseq.nii,1'];
            all_func = MY_find_images_in_all_scans(path,'Results',{GEEPI(:)},'^SRGrsmUw2dseq','.nii',[1 Inf],'all_mixed');
            all_func = [all_func;[path 'Results\T2\Sm2dseq.nii,1']];
            OldNormalize_mlb = MY_get_default_oldnormalize_batch_struct(ref, source, all_func);
            disp('Start to process OldNormalize!');
            F = spm_figure('GetWin');
            spm_jobman('run',OldNormalize_mlb);
            hgexport(figure(F), fullfile([path 'Results\'], strcat('oldnormalize')), hgexport('factorystyle'), 'Format', 'tiff');
            
            %*************** delete this part after this project***********
            Check_func = char(ref{1,1},[path 'Results\T2\nSm2dseq.nii,1']);
            for kk = GEEPI
                delete(fullfile(path,'Results',num2str(kk),'RGrsmUw2dseq.nii'))
                Check_func = char(Check_func,fullfile(path,'Results',num2str(kk),'nSRGrsmUw2dseq.nii,1'));
            end
            F = spm_figure('GetWin');
            spm_check_registration(Check_func);
            hgexport(figure(F), fullfile([path 'Results\'], strcat('Coreg_Check')), hgexport('factorystyle'), 'Format', 'tiff');
            %**************************************************************
        end
        if flag_stage == 10
            %% smooth_space
            all_func = MY_find_images_in_all_scans(path,'Results',{GEEPI(:)},'^nSRGrsmUw2dseq','.nii',[1 Inf],'all_mixed');
            Smooth_mlb = MY_get_default_smooth_batch_struct(all_func);
            disp('Start to process Smooth!');
            spm_jobman('run',Smooth_mlb);
            clear Smooth_mlb;
            
        end

        
    end
end