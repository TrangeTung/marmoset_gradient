function [hough_space,hough_circle,para] = hough_circle(BW,step_r,step_angle,r_min,r_max,p)  
  
  
circleParaXYR=[];  
para=[];  
  
[m,n] = size(BW);  
size_r = round((r_max-r_min)/step_r)+1;
size_angle = round(2*pi/step_angle);  
  
hough_space = zeros(m,n,size_r);  
  
[rows,cols] = find(BW);
ecount = size(rows);
  
% Hough
  
% a = x-r*cos(angle)  
% b = y-r*sin(angle)  
  
for i=1:ecount  
    for r=1:size_r 
        for k=1:size_angle 
            a = round(rows(i)-(r_min+(r-1)*step_r)*cos(k*step_angle));  
            b = round(cols(i)-(r_min+(r-1)*step_r)*sin(k*step_angle));  
            if(a>0&a<=m&b>0&b<=n)  
            hough_space(a,b,r) = hough_space(a,b,r)+1;
            end  
        end  
    end  
end  
  
  
max_para = max(max(max(hough_space)));
index = find(hough_space>=max_para*p);
length = size(index);
hough_circle = false(m,n);  
%hough_circle = zeros(m,n);  
for i=1:ecount  
    for k=1:length  
        par3 = floor(index(k)/(m*n))+1;  
        par2 = floor((index(k)-(par3-1)*(m*n))/m)+1;  
        par1 = index(k)-(par3-1)*(m*n)-(par2-1)*m;  
        if((rows(i)-par1)^2+(cols(i)-par2)^2<(r_min+(par3-1)*step_r)^2+5&...  
                (rows(i)-par1)^2+(cols(i)-par2)^2>(r_min+(par3-1)*step_r)^2-5)  
              hough_circle(rows(i),cols(i)) = true;  
        end  
    end  
end                 
  
for k=1:length  
    par3 = floor(index(k)/(m*n))+1;
    par2 = floor((index(k)-(par3-1)*(m*n))/m)+1;  
    par1 = index(k)-(par3-1)*(m*n)-(par2-1)*m;  
    circleParaXYR = [circleParaXYR;par1,par2,par3];  
    hough_circle(par1,par2)= true; 
    %fprintf(1,'test1:Center %d %d \n',par1,par2);  
end  
  
while size(circleParaXYR,1) >= 1  
    num=1;  
    XYR=[];  
    temp1=circleParaXYR(1,1);  
    temp2=circleParaXYR(1,2);  
    temp3=circleParaXYR(1,3);  
    c1=temp1;  
    c2=temp2;  
    c3=temp3;  
    temp3= r_min+(temp3-1)*step_r;  
   if size(circleParaXYR,1)>1       
     for k=2:size(circleParaXYR,1)  
      if (circleParaXYR(k,1)-temp1)^2+(circleParaXYR(k,2)-temp2)^2 > temp3^2  
         XYR=[XYR;circleParaXYR(k,1),circleParaXYR(k,2),circleParaXYR(k,3)]; 
      else    
      c1=c1+circleParaXYR(k,1);  
      c2=c2+circleParaXYR(k,2);  
      c3=c3+circleParaXYR(k,3);  
      num=num+1;  
      end   
    end  
   end   
      %fprintf(1,'sum %d %d radius %d\n',c1,c2,r_min+(c3-1)*step_r);  
      c1=round(c1/num);  
      c2=round(c2/num);  
      c3=round(c3/num);  
      c3=r_min+(c3-1)*step_r;  
      %fprintf(1,'num=%d\n',num)  
      %fprintf(1,'Center %d %d radius %d\n',c1,c2,c3);     
      para=[para;c1,c2,c3]; 
      circleParaXYR=XYR;  
end  