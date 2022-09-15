clc;clear
WholePath = 'E:\rsfMRI_marmoset\';
codepath = 'E:\rsfMRI_marmoset\fMRI_code_Trange\';
addpath(genpath(codepath));
ExpRecording = [codepath 'ExpIndex.xlsx'];
[~,~,CellData] = xlsread(ExpRecording);
ExpTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
AtlasNII = fullfile(codepath,'marmoset_atlas','Trange_atlas_RikenBMA_cortex.nii');
ihdr = spm_vol(AtlasNII);
AtlasMask = spm_read_vols(ihdr);

AtlasExcel = fullfile(codepath,'marmoset_atlas','Trange_modifid_RikenBMA_atlas_labels.xlsx');
[~,~,CellData] = xlsread(AtlasExcel);
AtlasTable = cell2table(CellData(2:end,:),'VariableNames',CellData(1,:));
atlas_label = AtlasTable.Labels;
Abbre = AtlasTable.Abbre;

[~,~,CellData] = xlsread(fullfile(WholePath,'Figure5','gene_expression.xlsx'));
xhdr=spm_vol(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.nii.gz'));
Labels = spm_read_vols(xhdr);
genelist = CellData(2:31,1);

AtlasTXT = loadtxt(fullfile(codepath,'MBM_v3.0.1','atlas_RikenBMA_cortex.txt'));
dest = fullfile(WholePath,'DiffusionGradient_ION','Mean');
filename = fullfile(dest,'PG_align.nii');
ihdr = spm_vol(filename);
PGf = spm_read_vols(ihdr);
clear Vs Vf
for rg = 1:116
    RegionName = Abbre{rg+20};
    A = cellfun(@(x)strcmpi(x,RegionName),AtlasTXT(:,2));
    Bin = AtlasTXT{A,1};
    lmask = Labels==Bin;
    Vf(rg,:) = nanmean(fmask(PGf,lmask));
end


Type = {'full';'SC';'GeM';'SCshuffle';'GeMshuffle'};
for tl=[1,5]%1:numel(Type)

    bins = 10;  Cut = 16;

    Ld=[];
    
    a_=0;  Vec = [];
    for idx = 1:62
        idx
        path = [ExpTable.Monkey{idx},filesep];
        RARE =  ExpTable.RARE(idx);
        SEEPI = sort([ExpTable.SEEPI1(idx),ExpTable.SEEPI2(idx)],'ascend');
        GEEPI = sort([  ExpTable.EPI4D1(idx),ExpTable.EPI4D2(idx),...
            ExpTable.EPI4D3(idx),ExpTable.EPI4D4(idx),...
            ExpTable.EPI4D5(idx),ExpTable.EPI4D6(idx),...
            ExpTable.EPI4D7(idx),ExpTable.EPI4D8(idx)],'ascend');
        GEEPI(isnan(GEEPI)) = [];

        for kk = GEEPI
            a_=a_+1;
            dest = fullfile(WholePath,'DiffusionGradient_ION');
            dest_ = fullfile(dest,replace(ExpTable.Monkey{idx},WholePath,''),num2str(kk));
            load(fullfile(dest_,'VigilanceIndex.mat'));
            
            Vec={};
            dest = fullfile(WholePath,'DiffusionGradient_ION_Glutamate');
            dest_ = fullfile(dest,replace(ExpTable.Monkey{idx},WholePath,''),num2str(kk),'Model');
            filename = fullfile(dest_,['predictedFC_',Type{tl},'_dynamic.nii']);
            FC = spm_read_vols_4D(spm_vol(filename));
            for bl=1:bins

                bs = 0.25/bins; 
                bands = -0.125+(bl-1:bl)*bs;
                if bl==1;bands(1)=-3;end
                if bl==bins;bands(2)=3;end

                locs = find(VI>=bands(1)&VI<=bands(2));
                locs(locs<Cut | locs>numel(VI)-Cut+1)=[];

                FC_ = FC(:,:,locs);
            
                Vec{1,bl}=nanmean(FC_,3);
            end
            
            G = GradientMaps('kernel','cs','approach','dm','n_components',5,'align','pa');
            G = G.fit(Vec,'Reference',Vf(:,1:5));
            
            
            Ld(:,:,a_) = cat(2,G.lambda{:});
            
        end
        
    end
    
    eval(['R2_',Type{tl},'=Ld;'])
end

R2diff = R2_full-R2_GeMshuffle;
PGR2diff = squeeze(R2diff(1,:,:))*100+1;
X1_ = mean(PGR2diff,2)*2;SEM1 = std(PGR2diff,0,2)/sqrt(size(X1_,2));SEM1(end)=SEM1(end)/3;
PGR2diff = squeeze(R2diff(2,:,:))*100+1;
X2_ = mean(PGR2diff,2)*2;SEM2 = std(PGR2diff,0,2)/sqrt(size(X1_,2));
close all
F = figure;
errorbar((1:10)-0.1,X1_([6:10,1:5]),SEM1,'r-','linewidth',2); hold on;
plot((1:10)-0.1,X1_([6:10,1:5]),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10)
xlim([0 11]);ylim([1 9])
set(gca,'box','off','fontsize',13,'linewidth',2,...
    'tickdir','out','ticklength',[0.03 0])
F = figure;
errorbar((1:10)+0.1,X2_,SEM2,'b-','linewidth',2);hold on;
plot((1:10)+0.1,X2_,'o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',10)
xlim([0 11]);ylim([1 9])
set(gca,'box','off','fontsize',13,'linewidth',2,...
    'tickdir','out','ticklength',[0.03 0])


GeM_random = xlsread(fullfile(WholePath,'Figure3','Random_Gene_Expression.xlsx'));
GeM_random(isnan(GeM_random))=0;
GeM = corr(GeM_random(:,21:end));
GeM(isnan(GeM)) = 0;
NoiseM = GeM;
Ref = mean(FC,3);
F = figure;
plot(NoiseM(1:3:end),Ref(1:3:end)*2,'o','MarkerEdgeColor','k','MarkerFaceColor',[.5 .5 .5],'MarkerSize',5)
ylim([0 1]);xlim([-1 1])
set(gca,'box','off','fontsize',15,'linewidth',2,...
    'tickdir','out','ticklength',[0.03 0])
y=Ref(:)*2;x=NoiseM(:);
p = polyfit(x,y,1);
y_= polyval(p,[min(x),max(x)]);
hold on;plot([min(x),max(x)],y_,'r','linewidth',4)


filename = fullfile(codepath,'MarmosetBrainConnectivity.xlsx');
FLNe_log = xlsread(filename,'log(FLNe)');
FLNe_log(isnan(FLNe_log))=-6;
SCcorr = corr(FLNe_log');
SCcorr=fillmissing(SCcorr,'linear');
SCcorr=fillmissing(SCcorr','linear');
SC = SCcorr;
GeM_raw = xlsread(fullfile(WholePath,'Figure5','gene_expression_correlation_matrix_edge_complex.xlsx'));
GeM = GeM_raw;%(lottery,lottery)  

x = SC(tril(ones(size(SC)),-1)==1);
y = GeM(tril(ones(size(GeM)),-1)==1);
F = figure;
plot(x,y,'o','MarkerEdgeColor','k','MarkerFaceColor',[0.4 0.6 0.4],'MarkerSize',5)
ylim([-1 1]);xlim([-1 1])
set(gca,'box','off','fontsize',15,'linewidth',2,...
    'tickdir','out','ticklength',[0.03 0])
p = polyfit(x,y,1);
y_= polyval(p,[min(x),max(x)]);
hold on;plot([min(x),max(x)],y_,'g','linewidth',4)

