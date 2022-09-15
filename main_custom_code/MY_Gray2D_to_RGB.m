function RGBmap = MY_Gray2D_to_RGB(Gray2D,lmask,lb,ub)


defaulRGBmap = parula(64);%

RGBmap = ones([size(Gray2D),3]);

map = (Gray2D-lb)/(ub-lb);
map(map<=0) = 10^-9;
map(map>=1) = 1;

Mapnan = isnan(map);
map(Mapnan) = 10^-9;

RGBmap(:,:,1) = reshape(defaulRGBmap(ceil(map*64),1)*255,size(map));
RGBmap(:,:,2) = reshape(defaulRGBmap(ceil(map*64),2)*255,size(map));
RGBmap(:,:,3) = reshape(defaulRGBmap(ceil(map*64),3)*255,size(map));

[a,b,c]=size(RGBmap);
Vmap = reshape(RGBmap,[a*b c]);
Vmap(lmask==0,1) = 0;
Vmap(lmask==0,2) = 0;
Vmap(lmask==0,3) = 0;

RGBmap = uint8(reshape(Vmap,[a,b,c]));
RGBmap = flip(flip(RGBmap,1),2);

% bar=zeros([a,round(a/2),3]);
% cb = flip(permute(repmat(defaulRGBmap,[1 1 10]),[1 3 2]),1);
% [x,y,~]=size(cb);
% cbr = imresize3(cb*255,[round(a*0.8),round(a*0.8/x*y),3]);
% [x,y,~]=size(cbr);
% bar(round(a*0.1)+(1:x),round(a/10)+(1:y),:)=cbr;
% 
% I=insertText(bar,[round(a/10)+y,round(a*0.1)],num2str(ub,'%1.2f'),...
%     'FontSize', 15,'TextColor',[255 255 255],'BoxColor',[0 0 0],'BoxOpacity',0,'AnchorPoint','LeftTop');
% I=insertText(I,[round(a/10)+y,round(a*0.1)+x],num2str(lb,'%1.2f'),...
%     'FontSize', 15,'TextColor',[255 255 255],'BoxColor',[0 0 0],'BoxOpacity',0,'AnchorPoint','LeftBottom');
% 
% RGBmap = cat(2,RGBmap,I);

end