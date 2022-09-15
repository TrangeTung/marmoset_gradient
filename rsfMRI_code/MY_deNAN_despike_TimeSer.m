function DeSer = MY_deNAN_despike_TimeSer(Ser,thr)

LOCS = (isnan(Ser));
SE = strel('line',round(thr/2),0);
LOCS = imdilate(LOCS,SE);
LOCS = imerode(LOCS,SE);
LOCSd = LOCS-[LOCS(2:end);LOCS(1)];
LOCSds = find(LOCSd==-1);
LOCSde = find(LOCSd==+1);
for k = 1:numel(LOCSds)
    LOCSa = Ser(LOCSds(k)-1);
    LOCSb = Ser(LOCSde(k)+1);
    Seg = LOCSds(k):LOCSde(k);
    X = (Seg-Seg(1)+1)/max((Seg-Seg(1)+1));
    RR = interp1(0:1,[LOCSa,LOCSb],X);
    Ser(Seg) = RR;
end

b = envelope(gradient(Ser));
LOCS = (b>thr);
SE = strel('line',round(thr/2),0);
LOCS = imdilate(LOCS,SE);
LOCS = imerode(LOCS,SE);
LOCSd = LOCS-[LOCS(2:end);LOCS(1)];
LOCSds = find(LOCSd==-1);
LOCSde = find(LOCSd==+1);
if ~isempty(LOCSds)
if LOCSds(1)>LOCSde(1);LOCSde(1)=[];end
if LOCSds(1)==1;LOCSde(1)=[];LOCSds(1)=[];end
for k = 1:numel(LOCSde)
    LOCSa = Ser(LOCSds(k)-1);
    LOCSb = Ser(LOCSde(k)+1);
    Seg = LOCSds(k):LOCSde(k);
    X = (Seg-Seg(1)+1)/max((Seg-Seg(1)+1));
    RR = interp1(0:1,[LOCSa,LOCSb],X);
    Ser(Seg) = RR;
end
end
DeSer=medfilt1(Ser,5);
end