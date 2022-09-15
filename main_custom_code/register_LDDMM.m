%% Register two images
%  Changed: Dec 31st, 2011
%
function [Mp,sx,sy,sz,vx,vy,vz] = register_LDDMM(F,M,opt)
    if nargin<3;  opt = struct();  end;
    if ~isfield(opt,'sigma_fluid');      opt.sigma_fluid     = 10;              end;
    if ~isfield(opt,'sigma_diffusion');  opt.sigma_diffusion = 10;              end;
    if ~isfield(opt,'sigma_i');          opt.sigma_i         = .1;              end;
    if ~isfield(opt,'sigma_x');          opt.sigma_x         = 1.0;              end;
    if ~isfield(opt,'niter');            opt.niter           = 250;              end;
    if ~isfield(opt,'vx');               opt.vx              = zeros(size(M));   end;
    if ~isfield(opt,'vy');               opt.vy              = zeros(size(M));   end;
    if ~isfield(opt,'vz');               opt.vz              = zeros(size(M));   end;
    if ~isfield(opt,'imagepad');         opt.imagepad        = 1.2;              end;
    if ~isfield(opt,'stop_criterium');   opt.stop_criterium  = 1e-4;             end;
    if ~isfield(opt,'do_display');       opt.do_display      = 1;                end;
    if ~isfield(opt,'do_plotenergy');    opt.do_plotenergy   = 1;                end;
    if ~isfield(opt,'do_avi');           opt.do_avi          = 0;                end;
    
    if opt.do_avi
        if ~isfield(opt, 'aviobj')
            opt.aviobj = avifile('spectral-demons.avi','compression','None', 'fps',10);
            opt.do_closeavi = 1; % create and close avi file here
        end
        if ~isfield(opt, 'do_closeavi'); opt.do_closeavi = 0; end;
        global aviobj;
        aviobj = opt.aviobj;
    end;
    %% Pad image
    [F,lim] = imagepad(F,opt.imagepad);
    [M,lim] = imagepad(M,opt.imagepad);
    
    %% T is the deformation from F to M
    if ~isempty(opt.vx) && ~isempty(opt.vy) && ~isempty(opt.vz)
        vx = imagepad(opt.vx,opt.imagepad);
        vy = imagepad(opt.vy,opt.imagepad);
        vz = imagepad(opt.vz,opt.imagepad);
    end
    e  = zeros(1,opt.niter);
    e_min = 1e+100;      % Minimal energy
    
    %% Iterate update fields
    for iter=1:opt.niter
        % Find update
        [ux,uy,uz] = findupdate(F,M,vx,vy,vz,opt.sigma_i,opt.sigma_x);
        % Regularize update
        ux    = imgaussian(ux,opt.sigma_fluid);
        uy    = imgaussian(uy,opt.sigma_fluid);
        uz    = imgaussian(uz,opt.sigma_fluid);
        % Compute step (max half a pixel)
        step  = opt.sigma_x;
        
        % Update correspondence (demons) - additive
        vx = vx + step*ux;
        vy = vy + step*uy;
        vz = vz + step*uz;
        % Update correspondence (demons) - composition
        %[vx,vy,vz] = compose(vx,vy,vz,step*ux,step*uy,step*uz);
        
        % Regularize deformation
        vx = imgaussian(vx,opt.sigma_diffusion);
        vy = imgaussian(vy,opt.sigma_diffusion);
        vz = imgaussian(vz,opt.sigma_diffusion);
        % Get Transformation
        [sx,sy,sz] = expfield(vx,vy,vz);  % deformation field
        
        % Compute energy
        e(iter) = energy(F,M,sx,sy,sz,opt.sigma_i,opt.sigma_x);
        disp(['Iteration: ' num2str(iter) ' - ' 'energy: ' num2str(e(iter))]);
        if e(iter)<e_min
            sx_min = sx; sy_min = sy; sz_min = sz;
            vx_min = vx; vy_min = vy; vz_min = vz;
            e_min  = e(iter);
        end
        
        % Stop criterium
        if iter>1 && abs(e(iter) - e(max(1,iter-5))) < e(1)*opt.stop_criterium
            break;
        end
        if opt.do_display
            % display deformation
            subplot(2,4,7); showvector(ux,uy,uz,4,0,lim); title('Update');
            subplot(2,4,8); showgrid  (sx,sy,sz,4,lim);   title('Transformation');
            drawnow;
            
            % Display registration
            Mp     = iminterpolate(M,sx,sy,sz);
            diff   = (F-Mp).^2;
            showimage(F,'Fixed', M,'Moving', Mp,'Warped', diff,'Diff', 'lim',lim,'nbrows',2,'caxis',[0 256]); drawnow;
            % Plot energy
            if opt.do_plotenergy
                subplot(2,2,3)
                hold on;
                plot(1:iter,e(1:iter),'r-'); xlim([0 opt.niter]);
                xlabel('Iteration'); ylabel('Energy');
                hold off;
                drawnow
            end
        end
        
        if opt.do_avi; aviobj = addframe(aviobj,getframe(gcf)); end;
        
    end
    
    %% Get Best Transformation
    vx = vx_min;  vy = vy_min;  vz = vz_min;
    sx = sx_min;  sy = sy_min;  sz = sz_min;
    %% Transform moving image
    Mp = iminterpolate(M,sx,sy,sz);
    
    %% Unpad image
    Mp = Mp(lim(1):lim(2),lim(3):lim(4),lim(5):lim(6));
    vx = vx(lim(1):lim(2),lim(3):lim(4),lim(5):lim(6));
    vy = vy(lim(1):lim(2),lim(3):lim(4),lim(5):lim(6));
    vz = vz(lim(1):lim(2),lim(3):lim(4),lim(5):lim(6));
    sx = sx(lim(1):lim(2),lim(3):lim(4),lim(5):lim(6));
    sy = sy(lim(1):lim(2),lim(3):lim(4),lim(5):lim(6));
    sz = sz(lim(1):lim(2),lim(3):lim(4),lim(5):lim(6));
    if opt.do_avi && opt.do_closeavi
        aviobj = close(aviobj);
    end
    
end



%% Display vector field
%  Changed: Dec 11th, 2011
%
function showvector(ux,uy,uz,downsample,scale,lim)
    if nargin<3; downsample = 1; end;
    
    sizex = size(ux,1);
    sizey = size(ux,2);
    sizez = size(ux,3);
    
    % Use mid slice as image
    slice = ceil(sizez/2);
    ux = ux(:,:,slice);
    uy = uy(:,:,slice);
    
    ux  = ux(1:downsample:end, 1:downsample:end);
    uy  = uy(1:downsample:end, 1:downsample:end);
    
    if nargin<5; scale = 3;                     end; % Scale vector to show small ones
    if nargin<6; lim   = [0 sizex-1 0 sizey-1 0 sizez-1]; end; % Display whole image
    [y,x] = ndgrid((1:downsample:sizex)+downsample/2, (1:downsample:sizey)+downsample/2); % coordinate image
    quiver(x,y,ux,uy,scale);                  % show vectors
    daspect([1 1 1]);
    axis([lim(3)-0.5 lim(4)+0.5 lim(1)-0.5 lim(2)+0.5]);      % which vector to show
    set(gca,'YDir','reverse');
    
end

%% Display several images on the same figure
%  Changed: Dec 9th, 2011
%
function showimage(varargin)
    % Check parameters
    nb_args   = size(varargin,2);
    nb_images = nb_args;
    nb_cols   = 0;
    nb_rows   = 1;
    row       = 1;
    crange    = [0 1]; % default image intensities
    
    for i=1:nb_args
        if ischar(varargin{i})
            if isequal(varargin{i},'lim')
                lim       = varargin{i+1};
                nb_images = nb_images-2;
            elseif isequal(varargin{i},'nbcols')
                nb_cols   = varargin{i+1};
                nb_images = nb_images-2;
            elseif isequal(varargin{i},'nbrows')
                nb_rows   = varargin{i+1};
                nb_images = nb_images-2;
            elseif isequal(varargin{i},'row')
                row       = varargin{i+1};
                if row>nb_rows; nb_rows = row; end;
                nb_images = nb_images-2;
            elseif isequal(varargin{i},'caxis')
                crange    = varargin{i+1};
                nb_images = nb_images-2;
            else
                nb_images = nb_images-1;
            end
        end
    end
    
    if nb_cols==0; nb_cols = nb_images; end;
    
    % Display images
    iter_image = 1;
    for iter_arg=1:nb_args
        if ~ischar(varargin{iter_arg})
            I = varargin{iter_arg};
            
            % Use mid slice as image
            slice = ceil(size(I,3)/2);
            I = I(:,:,slice);
            subplot(nb_rows,nb_cols,(row-1)*nb_cols + iter_image);
            imagesc(I,crange);
            daspect([1 1 1]);
            if exist('lim'); axis([lim(3)-0.5 lim(4)+0.5 lim(1)-0.5 lim(2)+0.5]); end
            axis off;
            if iter_arg+1<=nb_args && ischar(varargin{iter_arg+1})
                title(varargin{iter_arg+1});
            end
            iter_image = iter_image+1;
        end
        if iter_image>nb_images
            break;
        end
    end
        
end
%% Display vector field
%  Changed: Dec 11th, 2011
%
function showgrid(ux,uy,uz,downsample,lim)
    if nargin<3; downsample = 1; end;
    
    sizex = size(ux,1);
    sizey = size(ux,2);
    sizez = size(ux,3);
    
    % Use mid slice as image
    slice = ceil(sizez/2);
    ux = ux(:,:,slice);
    uy = uy(:,:,slice);
    
    ux  = ux(1:downsample:end, 1:downsample:end);
    uy  = uy(1:downsample:end, 1:downsample:end);
    
    if nargin<5; scale = 3;                     end; % Scale vector to show small ones
    if nargin<6; lim   = [0 sizex-1 0 sizey-1 0 sizez-1]; end; % Display whole image
    [y,x] = ndgrid((1:downsample:sizex)+downsample/2, (1:downsample:sizey)+downsample/2); % coordinate image
    z = zeros(size(x));
    mesh(x+ux,y+uy,z); view(2);
    daspect([1 1 1]);
    axis([lim(3) lim(4) lim(1) lim(2)] + .5 + [downsample 0 downsample 0]/2);      % which vector to show
    set(gca,'YDir','reverse');
    
end
% J = jacobian(S)
%   Determinant of Jacobian of a displacement field
%
% Herve Lombaert, Jan. 8th, 2013
%
function det_J = jacobian(sx,sy,sz)
    % Gradients
    [gx_y,gx_x,gx_z] = gradient(sx);
    [gy_y,gy_x,gy_z] = gradient(sy);
    [gz_y,gz_x,gz_z] = gradient(sz);
    
    % Add identity
    gx_x = gx_x + 1;
    gy_y = gy_y + 1;
    gz_z = gz_z + 1;
    
    % Determinant
    det_J = gx_x.*gy_y.*gz_z + ...
            gy_x.*gz_y.*gx_z + ...
            gz_x.*gx_y.*gy_z - ...
            gz_x.*gy_y.*gx_z - ...
            gy_x.*gx_y.*gz_z - ...
            gx_x.*gz_y.*gy_z;
end
%% Interpolate image
%  Changed: Dec 31st, 2011
%
% In Matlab, (x,y) are:
%
%      j               x
%   o------>        o------>
%   |               |
% i |  . I(i,j)   y |  . I(x,y)
%   |               |
%   v               v
%
function I = iminterpolate(I,sx,sy,sz,mode)
    if size(size(I),2)==4; I = iminterpolate_multichannel(I,sx,sy,sz,mode); return; end;
    
    if nargin<5; mode = 'linear'; end;
    
    % Find update points on moving image
    nx = size(I,1); ny = size(I,2); nz = size(I,3);
    [y,x,z] = ndgrid(1:nx, 1:ny, 1:nz); % coordinate image
    x_prime = x + sx; % updated x values (1st dim, rows)
    y_prime = y + sy; % updated y values (2nd dim, cols)
    z_prime = z + sz; % updated z values (3rd dim, slices)
    
    % Interpolate updated image
    I = interpn(y,x,z,I,y_prime,x_prime,z_prime,mode,0); % moving image intensities at updated points
    
end
%% Exponentiate vector field
%  Changed: Dec 31st, 2011
%
function [vx,vy,vz] = expfield(vx, vy, vz)
    % Find n, scaling parameter
    normv2 = vx.^2 + vy.^2 + vz.^2;
    m = sqrt(max(normv2(:)));
    n = ceil(log2(m/0.5)); % n big enough so max(v * 2^-n) < 0.5 pixel)
    n = max(n,0);          % avoid null values
    
    % Scale it (so it's close to 0)
    vx = vx * 2^-n;
    vy = vy * 2^-n;
    vz = vz * 2^-n;
    % square it n times
    for i=1:n
        [vx,vy,vz] = compose(vx,vy,vz, vx,vy,vz);
    end
end
%% Compose two vector fields
%  Changed: Dec 31st, 2011
%
function [vx,vy,vz] = compose(ax,ay,az, bx,by,bz)
    % Pad Image
    imagepadscale = 1.2; % just to get too many outside values
    [ax,lim] = imagepad(ax,imagepadscale);
    [ay,lim] = imagepad(ay,imagepadscale);
    [az,lim] = imagepad(az,imagepadscale);
    [bx,lim] = imagepad(bx,imagepadscale);
    [by,lim] = imagepad(by,imagepadscale);
    [bz,lim] = imagepad(bz,imagepadscale);
    % Coordinates
    nx = size(ax,1);
    ny = size(ax,2);
    nz = size(ax,3);
    
    [y,x,z] = ndgrid(1:nx, 1:ny, 1:nz); % coordinate image
    % Where points are going
    xp  = iminterpolate(x+ax, bx,by,bz);
    yp  = iminterpolate(y+ay, bx,by,bz);
    zp  = iminterpolate(z+az, bx,by,bz);
    
    % Update field
    vx = xp - x;
    vy = yp - y;
    vz = zp - z;
    
    % Zero vectors going outside the image
    zr  = (xp==0 & yp==0 & zp==0);
    vx(zr) = 0;
    vy(zr) = 0;
    vz(zr) = 0;
    % Unpad image
    vx = vx(lim(1):lim(2),lim(3):lim(4),lim(5):lim(6));
    vy = vy(lim(1):lim(2),lim(3):lim(4),lim(5):lim(6));
    vz = vz(lim(1):lim(2),lim(3):lim(4),lim(5):lim(6));
    
end

%% Get energy
function e = energy(F,M,sx,sy,sz,sigma_i,sigma_x)
    % Intensity difference
    Mp     = iminterpolate(M,sx,sy,sz);
    diff2  = (F-Mp).^2;
    area   = numel(M);
    % Transformation Gradient
    jac = jacobian(sx,sy,sz);
    
    % Three energy components
    e_sim  = sum(diff2(:)) / area;
    e_reg  = sum(jac(:).^2) / area;
    
    % Total energy
    e      = e_sim + (sigma_i^2/sigma_x^2) * e_reg;
end
%% Pad image
%  Changed: Dec 31st, 2011
%
function [I,lim] = imagepad(I,scale)
    if nargin<2; scale = 2; end; % default, pad image twice as big
    
    if size(size(I),2)>3; [I,lim]=imagepad_multichannel(I,scale); return; end; % check if I is multichannel
    
    Ip  = zeros(ceil(size(I)*scale));
    lim = bsxfun(@plus, floor(size(I)*(scale-1)/2), [[1 1 1]; size(I)]); % image limits
    Ip(lim(1):lim(2),lim(3):lim(4),lim(5):lim(6)) = I;                   % pad image
    I = Ip;
end
%% Pad image
%  Changed: Dec 31st, 2011
%
function [I,lim] = imagepad_multichannel(I,scale)
    if nargin<2; scale = 2; end; % default, pad image twice as big
    nchannels = size(I,4);
    for i=1:nchannels
        Ip(:,:,:,i) = imagepad(I(:,:,:,i),scale);
    end
        
    imagesize = [size(I,1) size(I,2) size(I,3)];
    lim = bsxfun(@plus, floor(imagesize*(scale-1)/2), [[1 1 1]; imagesize]); % image limits
    I = Ip;
end

%% Find update between two images
%  Changed: Dec 31st, 2011
%
function [ux,uy,uz] = findupdate(F,M,vx,vy,vz,sigma_i,sigma_x)
    % Get Transformation
    [sx,sy,sz] = expfield(vx,vy,vz);
    % Interpolate updated image
    M_prime = iminterpolate(M,sx,sy,sz);     % intensities at updated points
    area    = size(M,1)*size(M,2)*size(M,3); % area of moving image
    
    % image difference
    diff = F - M_prime;
    
    % fixed image gradient
    [gx_f,gy_f,gz_f] = gradient(F);          % image gradient
    normg2_f  = gx_f.^2 + gy_f.^2 + gz_f.^2; % squared norm of gradient
    
    % moving image gradient
    [gx,gy,gz] = gradient(M_prime);          % image gradient
    normg2  = gx.^2 + gy.^2 + gz.^2;         % squared norm of gradient
    
    % update is Idiff / (||J||^2+(Idiff^2)/sigma_x^2) J, with Idiff = F(x)-M(x+s), and J = Grad(M(x+s));
    scale = diff ./ (normg2 + diff.^2*sigma_i^2/sigma_x^2);
    scale(normg2==0) = 0;
    scale(diff  ==0) = 0;
    scale = scale .* ((scale>=0) + (scale<0) .* sign(gx.*gx_f + gy.*gy_f + gz.*gz_f)); % avoid collapsing gradients (change sign if moving goes backward)
    ux = gx .* scale;
    uy = gy .* scale;
    uz = gz .* scale;
    
    % Zero non overlapping areas
    %ux(F==0)       = 0; uy(F==0)       = 0;
    %ux(M_prime==0) = 0; uy(M_prime==0) = 0;
end

%% Apply gaussian filter to image
%  Changed: Dec 31st, 2011
%
function I = imgaussian(I,sigma)
    if sigma==0; return; end; % no smoothing
    
    if size(size(I))==4; I = imgaussian_multichannel(I,sigma); return; end;
    
    % Create Gaussian kernel
    radius   = ceil(3*sigma);
    [y,x,z]  = ndgrid(-radius:radius,-radius:radius,-radius:radius); % kernel coordinates
    h        = exp(-(x.^2 + y.^2 + z.^2)/(2*sigma^2));
    h        = h / sum(h(:));
    
    % Filter image
    I = imfilter(I,h);
end
%% Apply gaussian filter to image
%  Changed: Dec 31st, 2011
%
function I = imgaussian_multichannel(I,sigma)
    if sigma==0; return; end; % no smoothing
    
    nchannels = size(I,4);
    for i=1:nchannels
        Ip(:,:,:,i) = imgaussian(I(:,:,:,i),sigma);
    end
    
    I = Ip;
end

