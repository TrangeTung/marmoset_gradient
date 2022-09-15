close all
clear X
Type = {'full';'SC';'GeM';'SCshuffle';'GeMshuffle'};

for tl=1:5%:numel(Type)
    
    bins = 10;  Cut = 16;
    
    for bl=1:bins
        
        dest = fullfile(WholePath,'DiffusionGradient_ION','DynamicModel',Type{tl},num2str(bl,'%02d'));
        load(fullfile(dest,'G.mat'));
        X(:,bl,tl)=G.lambda{1};
    end
end
% full lambda
F=figure('Position', [680 678 339 300]);
plot(squeeze(X(1,:,1)),'r');hold on;
plot(squeeze(X(2,:,1)),'b')
set(gca,'xtick','','linewidth',1.5,'box','off');

% SCshuffle lambda
F=figure('Position', [680 678 339 300]);
plot(squeeze(X(1,:,4)),'r');hold on;
plot(squeeze(X(2,:,4)),'b')
set(gca,'xtick','','linewidth',1.5,'box','off');

% GeMshuffle lambda
F=figure('Position', [680 678 339 300]);
plot(squeeze(X(1,:,5)),'r');hold on;
plot(squeeze(X(2,:,5)),'b')
set(gca,'xtick','','linewidth',1.5,'box','off');


%% Arousal related R2
Type = {'full';'SC';'GeM';'SCshuffle';'GeMshuffle'};
for tl=4:5
    
    bins = 10;  Cut = 16;
    
    for bl=1:bins
        
        dest = fullfile(WholePath,'DiffusionGradient_ION','DynamicModel',Type{tl},num2str(bl,'%02d'));
        data = csvread(fullfile(dest,'Connectome.csv'));
        FCz = data;
        if tl>=4
            dest_ = fullfile(WholePath,'DiffusionGradient_ION','DynamicModel',Type{1},num2str(bl,'%02d'));
            data_ = csvread(fullfile(dest_,'Connectome.csv'));
            FCz = data_-data;
        end
        
        Vec = data_(:)';
        VecP = FCz(:)';
        mu = mean(VecP,2).*mean(Vec,2)-mean(Vec.*VecP,2);
        md = mean(VecP,2).^2 - mean(VecP.^2,2);
        b = mean(Vec,2) - mu./md.*mean(VecP,2);
        Vecf = mu./md.*VecP+b;
        SSR = sum((Vecf-mean(Vec,2)).^2,2);
        SST = sum((Vec-mean(Vec,2)).^2,2);
        R2(bl,tl-3) = SSR./SST;

        
    end
    
end
F = figure;plot(R2)
ylabel('R2');
set(gca,'linewidth',2,'fontsize',15,'box','off','xtick','');

%% topographic flow
ICApath = fullfile(WholePath,'SubMit','ICA_final','processed');
filename = fullfile(ICApath,'normICAv1_all_vtsplit_premotor.nii.gz');
ihdr = spm_vol(filename);
img = spm_read_vols(ihdr);
ICAmask = img>0;
ICAname = {'ventralsoma';'dorsalsoma';'frontalpole';'parahip';'OPFC';...
    'AUD';'frontalpaietal';'DMN';'visualMT';'visual0';'visual1';...
    'visual2';'visual3';'mACC';'premotor'};

% gradient
Type = {'full';'SC';'GeM';'SCshuffle';'GeMshuffle'};
for tl=[1 4 5]

bins = 10; V_ION=[];
for bl=1:bins
    dest = fullfile(WholePath,'DiffusionGradient_ION','DynamicModel',Type{tl},num2str(bl,'%02d'));
    filename = fullfile(dest,'PG.nii');
    ihdr = spm_vol(filename);
    PG = spm_read_vols(ihdr);
    [a,b,c,d] = size(PG);
    
    for i=1:numel(ICAname)
        X = ICAmask(:,:,:,i);
        Xmask = imresize3(X,[a,b,c],'nearest');
        V_ION(i,bl,:) = nanmean(fmask(PG,Xmask),1); 
    end
end

Map = [0:255;zeros(1,256);255:-1:0]';
MyMap = 1/255*interp1(1:256,Map,10:25:256);
F = figure( 'Position', [680 607 409 371]);
for iy = 1:9
    for ix=1:numel(ICAname)
       X = V_ION(ix,iy:iy+1,2);
       Y = V_ION(ix,iy:iy+1,1);
       
       quiver(X(1),Y(1),X(2)-X(1),Y(2)-Y(1),'color',MyMap(iy,:),...
           'linewidth',1,'MaxHeadSize',1);
       hold on;
       if iy==5
          text(X(1),Y(1),ICAname{ix},'HorizontalAlignment','center') ;
       end
       
    end 
end
xlabel('Gradient 2');ylabel('Gradient 1');
set(gca,'linewidth',2,'fontsize',15,'box','off');

end