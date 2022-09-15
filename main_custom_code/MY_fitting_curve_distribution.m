function [x_,y_,Rsquare]=MY_fitting_curve_distribution(x,y,type)

EX = sum(x.*y/sum(y));
DX = sum((x.^2).*(y/sum(y)))-EX^2;


switch type
    case 'gaussian'
        
        gaussEqn = 'c*exp(-((x-a)/b).^2)';

        p0(1) = EX;   p0(2) = sqrt(DX);   
        p0(3) = 1;%max(y(:))/max(normpdf(x,p0(1),p0(2)));
        
        [f1,goodness] = fit(x',y',gaussEqn,'Start', p0);
        a=f1.a;b=f1.b;c=f1.c;
        
        x_ = x;
        y_ = c*exp(-((x_-a)/b).^2);
        
        Rsquare = goodness.rsquare;
    case 'gamma'
        GamEqn = 'c*gampdf(x,a,b)';
        
        p0(1) = EX*EX/DX;   p0(2) = 1/(EX/DX);%std(y.*x)   
        p0(3) = max(y(:))/max(gampdf(x,p0(1),p0(2)));
        if p0(3)<10^-2;p0(3)=1;end
        
        try
            [f2,goodness] = fit(x',y',GamEqn,'Start', p0);
            a=f2.a;b=f2.b;c=f2.c;
        catch
            p1=lsqcurvefit(@(p,x)p(3)*gampdf(x,p(1),p(2)),p0,x(:),y(:));
            a=p1(1);b=p1(2);c=p1(3);
            x_ = x;
            y_ = c*gampdf(x_,a,b);
            [~,goodness] = fit(y',y_','a*x+b','Start', [1,0]);
        end
        
        x_ = x;
        y_ = c*gampdf(x_,a,b);
        
        
        Rsquare = goodness.rsquare;
        
        
    case 'exponential'
        ExpEqn ='b*exp(a*x)';
        
        ys = smooth(x,y,0.1,'loess');
        ym = max(ys); xm = x(ys==ym);
        
        p0(1) = (log(ym)-log(y(end)))/(xm-x(end));
        p0(2) = ym/exp(p0(1)*xm);
        
        [f3,goodness] = fit(x',y',ExpEqn,'Start', p0);
        a=f3.a;b=f3.b;
        
        x_ = x;
        y_ = b*exp(a*x_);
        
        Rsquare = goodness.rsquare;

end




    
end

