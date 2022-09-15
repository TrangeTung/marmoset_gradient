
function cRGB = MY_display_ColorMap_RGB(cdata,bar_value)


%% threshold map
%{
%{
Num=numel(cdata);
bar = sort(bar_value,'ascend');
negmin = bar(1);
negmax = bar(2);
posmin = bar(3);
posmax = bar(4);

cdata(cdata>posmax)=posmax;
cdata(cdata<negmin)=negmin;

posBW = cdata>posmin;
negBW = cdata<negmax;
%}


cdata(posBW)=cdata(posBW)/posmax;
cdata(negBW)=cdata(negBW)/negmin*(-1);
% white-->red // white-->blue
%{
cRGB = ones(Num,3)*0.8;
cRGB(negBW,1)=1-abs(cdata(negBW));
cRGB(negBW,2)=1-abs(cdata(negBW));
cRGB(negBW,3)=1;
cRGB(posBW,1)=1;
cRGB(posBW,2)=1-cdata(posBW);
cRGB(posBW,3)=1-cdata(posBW);
zeoBW = (cdata<posmin & cdata>negmax)|isnan(cdata);
cRGB(zeoBW,1)=1;
cRGB(zeoBW,2)=1;
cRGB(zeoBW,3)=1;
%}

% red-->yellow // blue-->green
%
cRGB = ones(Num,3)*0.8;
cRGB(negBW,1)=0;
cRGB(negBW,2)=abs(cdata(negBW));
cRGB(negBW,3)=1-abs(cdata(negBW));
cRGB(posBW,1)=1;
cRGB(posBW,2)=cdata(posBW);
cRGB(posBW,3)=0;
%}


%}

%% no threshold map
%
defaultMap = autumn(200);
% defaultMap = [flipud(hot(200));];%repmat([0.5,0.5,0.5],10,1)
% defaultMap = parula(200);%
% map0=[26+(255-26)/255*(1:255)',136+(255-136)/255*(1:255)',159+(255-159)/255*(1:255)'];
% map1=[242+(255-242)/255*(255:-1:1)',154+(255-154)/255*(255:-1:1)',0+(255-0)/255*(255:-1:1)'];
% Map = [map0;map1]/255;
% defaultMap = interp1(1:size(Map,1),Map,(1:200)/200*size(Map,1),'pchip');
% HCP colorbar
Map = [74,45,123;69,48,123;69,50,125;69,55,126;62,62,129;56,65,132;52,68,134;
    54,71,136;45,77,136;36,80,137;28,83,137;24,89,138;9,93,139;14,95,139;
    3,100,137;2,104,137;17,107,138;10,112,139;15,116,139;0,118,138;
    22,121,139;16,125,140;21,127,139;17,132,140;20,133,136;10,139,137;
    17,141,137;12,146,135;15,147,134;18,150,132;29,154,132;18,159,126;
    27,162,127;15,165,126;34,168,124;33,170,121;33,172,118;37,176,115;
    39,179,111;47,180,107;41,181,103;57,180,97;64,182,96;62,183,93;
    68,183,87;74,182,82;84,184,74;98,188,69;112,188,67;130,195,65;
    143,199,62;151,203,59;168,207,56;178,210,53;191,214,50;204,218,46;
    216,220,41;228,221,41;239,222,63;245,197,95;248,157,44;242,120,30;
    228,85,39;190,62,39]/255;
defaultMap = interp1(1:size(Map,1),Map,(1:200)/200*size(Map,1),'pchip');
% defaultMap = [repmat([.8 .8 .8],100,1);[ones(100,1),1-(0.01:0.01:1)',1-(0.01:0.01:1)']];

red = [ones(100,1),1-(0.01:0.01:1)',1-(0.01:0.01:1)'];
blue = [(0.01:0.01:1)',(0.01:0.01:1)',ones(100,1)];
% defaultMap = ([blue;red]);

map_nothre_reshape = cdata;
nothrebar = [min(bar_value(:)) max(bar_value(:))];
tmap = ones([size(map_nothre_reshape,1),3]);
map_normalize = (map_nothre_reshape-nothrebar(1))/(nothrebar(2)-nothrebar(1));
map_normalize(map_normalize<=0) = 0.000000001;
map_normalize(map_normalize>=1) = 1;
map_normalize1=map_normalize;
map_normalize(isnan(map_normalize)) = 0.0000000001;
tmap(:,1) = reshape(defaultMap(ceil(map_normalize*200),1),size(map_normalize));
tmap(:,2) = reshape(defaultMap(ceil(map_normalize*200),2),size(map_normalize));
tmap(:,3) = reshape(defaultMap(ceil(map_normalize*200),3),size(map_normalize));
tmap(isnan(map_normalize1),:) = repmat([0.75,0.75,0.75],numel(find(isnan(map_normalize1))),1);
cRGB = tmap;
%}

end
