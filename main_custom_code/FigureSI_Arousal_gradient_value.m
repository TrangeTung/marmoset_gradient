
clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
cd(codepath);
addpath(genpath(codepath));


ICApath = fullfile(WholePath,'SubMit','ICA_final','processed');
filename = fullfile(ICApath,'normICAv1_all_vtsplit_premotor.nii.gz');
ihdr = spm_vol(filename);
img = spm_read_vols(ihdr);
ICAmask = img>0;
ICAname = {'ventralsoma';'dorsalsoma';'frontalpole';'parahip';'OPFC';...
    'AUD';'frontalpaietal';'DMN';'visualMT';'visual0';...
    'visual1';'visual2';'visual3';'mACC';'premotor'};
ICAcolor = {[061,139,198];[080,185,231];[235,051,066];[212,179,130];[207,062,079];...
    [062,181,174];[236,091,083];[245,173,051];[051,134,093];[225,230,051];...
    [171,194,119];[078,189,097];[057,166,097];[233,051,153];[067,084,160]};

% gradient
bins = 10; V_ION=[];
for bl=1:bins
    dest = fullfile(WholePath,'DiffusionGradient_ION','Dynamic_VIweighted',num2str(bl,'%02d'));
    filename = fullfile(dest,'PG_align.nii');
    ihdr = spm_vol(filename);
    PG = spm_read_vols(ihdr);
    [a,b,c,d] = size(PG);
    
    for i=1:numel(ICAname)
        X = ICAmask(:,:,:,i);
        Xmask = imresize3(X,[a,b,c],'nearest');
        V_ION(i,bl,:) = nanmean(fmask(PG,Xmask),1); 
    end
end

bins = 10; V_NIH=[];
for bl=1:bins
    dest = fullfile(WholePath,'DiffusionGradient_NIH','Dynamic_VIweighted',num2str(bl,'%02d'));
    filename = fullfile(dest,'PG_align.nii');
    ihdr = spm_vol(filename);
    PG = spm_read_vols(ihdr);
    [a,b,c,d] = size(PG);
    
    for i=1:numel(ICAname)
        X = ICAmask(:,:,:,i);
        Xmask = imresize3(X,[a,b,c],'nearest');
        V_NIH(i,bl,:) = nanmean(fmask(PG,Xmask),1); 
    end
end

V = (V_ION+V_NIH)/2; 
F = figure('Position', [680 120 984 858]);
for grd = 1:4
    subplot(2,2,grd);
%    if grd==1;Reorder = [12,11,10,09,13,08,04,07,05,03,14,06,15,01,02];end 
%    if grd==2;Reorder = [01,06,11,12,14,02,15,04,05,13,03,09,10,08,09];end 
%    if grd==3;Reorder = [12,11,10,09,13,08,04,07,05,03,14,06,15,01,02];end 
%    if grd==4;Reorder = [12,11,10,09,13,08,04,07,05,03,14,06,15,01,02];end 
%     VX = squeeze(V(:,:,grd));
%     ram = VX(:,5)+VX(:,6);
%     maxmin = [find(ram==max(ram)),find(ram==min(ram))];
    for ic=1:15
        plot(squeeze(V(ic,:,grd)),'Color',ICAcolor{ic}/255,'linewidth',1.5);
        hold on;
    end
    ylabel('Gradients');
    set(gca,'box','off','linewidth',1.5,'fontsize',15,'tickdir','out')
end






