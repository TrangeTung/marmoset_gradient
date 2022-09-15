function MY_Double_Slice_with_Empty(Root_path,NIIfile,folder,file_type)

for kk = folder
    switch lower(file_type)
        case {'epi'}
            subpath = num2str(kk);
        otherwise
            subpath = 'T2';
    end
    head = spm_vol(fullfile(Root_path,subpath,NIIfile));
    Img = spm_read_vols(head);
    newImg = zeros([size(Img,1) size(Img,2) size(Img,3)*2 size(Img,4)]);
    % new Img
    newImg(:,:,size(Img,3)/2+1:size(Img,3)/2+size(Img,3),:) = Img;
    
    or_mat =    head(1).mat;
    or_mat(4,3) = or_mat(4,3)*2;

    for k = 1:size(Img,4)
        Vol(k,1) = struct(  'fname',    fullfile(Root_path,subpath,['S' NIIfile]),...
            'dim',      [size(Img,1) size(Img,2) size(Img,3)*2],...
            'mat',      or_mat,...
            'n',        [k,1],...
            'pinfo',    [1;0;0],...
            'descrip',  'doubleSlice',...
            'dt',       [16 0]);
    end
    spm_write_vol_4D(Vol,newImg);
    
    
end


end