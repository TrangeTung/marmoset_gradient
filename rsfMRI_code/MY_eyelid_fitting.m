function varargout= MY_eyelid_fitting(flag,varargin)

switch flag
    case 1
        [lmask,lib,Ig] = parse_inputs(varargin);
        % template based eyelid capture
        [a,b,c] = size(lib);
        libVec = reshape(lib,[a*b c]);
        [a,b,c]=size(Ig);
        Vec = reshape(Ig,[a*b c]);
        
        for i=1:size(lib,3)
            CC(i)=MY_two_IMAGE_nmi(Ig,lib(:,:,i));
        end
        %CC = corr(double(Vec),libVec);
        p = find(CC==max(CC));
        pmask = lmask(:,:,p);
%         moving = lib(:,:,p);
%         fixed = Ig;
%         [optimizer, metric] = imregconfig('monomodal');
%         tformSimilarity = imregtform(moving,fixed,'similarity',optimizer,metric);
%         Rfixed = imref2d(size(fixed));
%         pmask = imwarp(pmask,tformSimilarity,'OutputView',Rfixed);
        
        varargout={pmask};
        
    case 2
        
        [pmask,anchor,Ig] = parse_inputs(varargin);
        % circular fitting (top and bottom eyelids)
        eyet = (pmask==1);
        eyeb = (pmask==2);
        SE = strel('disk',10);
        eyetd = imdilate(eyet,SE);
        eyebd = imdilate(eyeb,SE);
        
        if isempty(anchor);anchor=zeros(size(pmask));end
        
        BW = edge(Ig,'canny');
        eyet = eyetd.*BW; eyetl=bwareaopen(eyet,15)+anchor;
        eyeb = eyebd.*BW; eyebl=bwareaopen(eyeb,15)+anchor;
        
        [rt,Ct] = MY_circular_fitting_LMS(eyetl);
        [rb,Cb] = MY_circular_fitting_LMS(eyebl);
        
        varargout= {[rt;Ct],[rb;Cb]};
        
        
end

end


function  [lmask,lib,Ig] = parse_inputs(varargin)

in= varargin{1};
lmask = in{1}.masklib;
lib =  in{1}.imglib;
Ig = in{1}.img;

end