function Gray2D = MY_reshape_3D_to_2D(Gray3D,rowNum)



[X,Y,Z] = size(Gray3D);
blank = ceil(Z/rowNum)*rowNum - Z;

% Gray3D = Gray3D - min(Gray3D(:));
% Gray3D = Gray3D/max(Gray3D(:))*256;
Gray2D = ones([Y*rowNum,X*(Z+blank)/rowNum])*mean(Gray3D(:))*nan;
for j = 1:Z
    Block_Row = ceil(j/((Z+blank)/rowNum));
    Block_Column = (j - ((Z+blank)/rowNum)*(Block_Row-1));
   
    Gray2D(((Block_Row-1)*Y+1):Block_Row*Y,((Block_Column-1)*X+1):Block_Column*X) = Gray3D(:,:,j)' ;
end
end

