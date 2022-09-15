Eg_path = 'D:\rsfMRI_marmoset\20181103_091548_20181103_NIH_marmoset_restfmri_6143m_1_1\Results\';
cd(Eg_path);

%% SEEPI
Raw_Img = spm_read_vols(spm_vol('GEEPI.nii'));
PhaseEncoding1 = mean(Raw_Img(:,:,:,1:8),4);
PhaseEncoding2 = mean(Raw_Img(:,:,:,9:16),4);
Raw_Gray2D_1 = MY_reshape_3D_to_2D(PhaseEncoding1);
Raw_Gray2D_2 = MY_reshape_3D_to_2D(PhaseEncoding2);

UnWrap_Img = spm_read_vols(spm_vol('Uw_GEEPI.nii'));
PhaseEncoding1 = mean(UnWrap_Img(:,:,:,1:8),4);
PhaseEncoding2 = mean(UnWrap_Img(:,:,:,9:16),4);
UnWrap_Gray2D_1 = MY_reshape_3D_to_2D(PhaseEncoding1);
UnWrap_Gray2D_2 = MY_reshape_3D_to_2D(PhaseEncoding2);

q = [Raw_Gray2D_1,UnWrap_Gray2D_1;Raw_Gray2D_2,UnWrap_Gray2D_2];
[QsizeX,QsizeY] = size(q);
q(round((QsizeX/2-1:QsizeX/2+1)),:) = max(q(:));
q(:,round((QsizeY/2-1:QsizeY/2+1))) = max(q(:));
p = [Raw_Gray2D_1-Raw_Gray2D_2,UnWrap_Gray2D_1-UnWrap_Gray2D_2];
p(:,round((QsizeY/2-1:QsizeY/2+1))) = max(p(:));

F = figure;set(F,'Position',[0,0,1200,600]*1.5);
subplot(6,1,1:4);imshow(q,[]);colorbar;
title('GEEPI Raw Images (Left) and Unwrapped Images (Right)','Fontsize',16)
subplot(6,1,5:6);imshow(p,[]);colorbar;caxis([-2000,2000]);
title('GEEPI Gray Value Differences of Raw Images (Left) and Unwrapped Images (Right)','Fontsize',16)
% print(F,fullfile(WholePath,'TopUpComparison_SEEPI.tiff'),'-dtiff','-r600');

%% GEEPI
PhaseEncoding1 = mean(spm_read_vols(spm_vol(fullfile(Eg_path,'17','2dseq.nii,1'))),4);
PhaseEncoding2 = mean(spm_read_vols(spm_vol(fullfile(Eg_path,'18','2dseq.nii,1'))),4);
Raw_Gray2D_1 = MY_reshape_3D_to_2D(PhaseEncoding1);
Raw_Gray2D_2 = MY_reshape_3D_to_2D(PhaseEncoding2);

PhaseEncoding1 = mean(spm_read_vols(spm_vol(fullfile(Eg_path,'17','Uw2dseq.nii'))),4);
PhaseEncoding2 = mean(spm_read_vols(spm_vol(fullfile(Eg_path,'18','Uw2dseq.nii'))),4);
UnWrap_Gray2D_1 = MY_reshape_3D_to_2D(PhaseEncoding1);
UnWrap_Gray2D_2 = MY_reshape_3D_to_2D(PhaseEncoding2);

q = [Raw_Gray2D_1,UnWrap_Gray2D_1;Raw_Gray2D_2,UnWrap_Gray2D_2];
[QsizeX,QsizeY] = size(q);
q(round((QsizeX/2-1:QsizeX/2+1)),:) = max(q(:));
q(:,round((QsizeY/2-1:QsizeY/2+1))) = max(q(:));
p = [Raw_Gray2D_1-Raw_Gray2D_2,UnWrap_Gray2D_1-UnWrap_Gray2D_2];
p(:,round((QsizeY/2-1:QsizeY/2+1))) = max(p(:));

F = figure;set(F,'Position',[0,0,1200,600]*1.5);
subplot(6,1,1:4);imshow(q,[]);colorbar;
title('GEEPI Raw Images (Left) and Unwrapped Images (Right)','Fontsize',16)
subplot(6,1,5:6);imshow(p,[]);colorbar;caxis([-3000,3000]);
title('GEEPI Gray Value Differences of Raw Images (Left) and Unwrapped Images (Right)','Fontsize',16)
print(F,fullfile(WholePath,'TopUpComparison_GEEPI.tiff'),'-dtiff','-r600');



function Gray2D = MY_reshape_3D_to_2D(Gray3D)

pump = 3;
blank = 1;
[SizeX,SizeY,SizeZ] = size(Gray3D);
Gray2D = ones([SizeY*pump,SizeX*(SizeZ+blank)/pump])*min(Gray3D(:));
for j = 1:SizeZ
    Block_Row = ceil(j/((SizeZ+blank)/pump));
    Block_Column = (j - (SizeZ+blank)/pump*(Block_Row-1));
   
    Gray2D(((Block_Row-1)*SizeY+1):Block_Row*SizeY,((Block_Column-1)*SizeX+1):Block_Column*SizeX) = Gray3D(:,:,j)' ;
end
end


