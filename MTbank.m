function [X, params] = MTbank(Xt,Xr,maxspeed)

    w = single(sqrt(size(Xt,2)));
    ns = round((w/2)^2);
    Xt=single(Xt);
    Xr=single(Xr);
    
    if nargin < 3
        maxspeed = single(max(Xr(:)));
    end
    maxspeed=single(maxspeed);
    
    
    
%     if nargin < 4
%         sigmas = 1.8;
%     end
        
   % bwspeed = 1;
    %bwtheta = 2.5;
   bwtheta = 1;%single([.4,.72,2]);
   bwspeed = 1;%single([0.5,1,2,4]);
   sigmas = 2;%2;%single([0.1, 0.5, 1, 2]);
    
    [xi yi] = meshgrid(1:w,1:w);
    
    xr = [];
    yr = [];
    
    iter = 1;

    
 X = zeros(size(Xt,1),3,8,ns,length(bwtheta),length(bwspeed),length(sigmas),'single');
 
%     filts=zeros(w^2,(w/2)^2,length(sigmas),'single');
filts=zeros(w^2,ns,length(sigmas),'single');

        
   for jj = 1:w/2
       for ii = 1:w/2
                   
            x0 = single((ii-1)*2+1.5);
            y0 = single((jj-1)*2+1.5);

            for yy=1:length(sigmas)
            filt = exp(-((xi-x0).^2+(yi-y0).^2)/2/sigmas(yy)^2);
            filt = filt(:);
            filts(:,iter,yy) = filt;
            end
            xr = [xr;x0];
            yr = [yr;y0];
            iter = iter+1;
        end
    end
    
    filts = filts/36;
   %% 
    spds =single([0.8,2,3]);
    themaxes=zeros(length(spds),length(bwspeed),'single');
for nn=1:length(bwspeed)
    for ll = 1:3
        prefspeed = spds(ll)/3*maxspeed;
        Xrs = 0:.01:maxspeed;
        themaxes(ll,nn) = single(max(2*(exp(-bwspeed(nn)^-2/2*((Xrs-prefspeed).^2))-exp(-bwspeed(nn)^-2/2*((Xrs+prefspeed).^2)))));
    end
end
    
    thresh = @(x) x.*(x>0);
    
    xs = zeros(3,8,ns,length(bwtheta),length(bwspeed),length(sigmas),'single');
    ys = zeros(3,8,ns,length(bwtheta),length(bwspeed),length(sigmas),'single');
    thetas = zeros(3,8,ns,length(bwtheta),length(bwspeed),length(sigmas),'single');
    speeds = zeros(3,8,ns,length(bwtheta),length(bwspeed),length(sigmas),'single');
    bwthetas = zeros(3,8,ns,length(bwtheta),length(bwspeed),length(sigmas),'single');
    bwspeeds = zeros(3,8,ns,length(bwtheta),length(bwspeed),length(sigmas),'single');
    sigmass = zeros(3,8,ns,length(bwtheta),length(bwspeed),length(sigmas),'single');
    for hh=1:length(sigmas)
        for uu=1:length(bwspeed)
            for oo=1:length(bwtheta)
                for kk = 1:8
                    preftheta = single((kk-1)/4*pi);
                    
                    X2 = 1/sinh(bwtheta(oo))*(exp(bwtheta(oo)*cos(Xt-preftheta))-1);
                    
                    for ll = 1:3
                        prefspeed = spds(ll)/3*maxspeed;
                        
                        X1 = 2*(exp(-bwspeed(uu)^-2/2*((Xr-prefspeed).^2))-exp(-bwspeed(uu)^-2/2*((Xr+prefspeed).^2)));
                        
                        Xf = X1.*X2/themaxes(ll,uu);
                        
%                         X(:,ll,kk,:,oo,uu,hh) = thresh(Xf*squeeze(filts(:,:,hh)));
                        X(:,ll,kk,:,oo,uu,hh) = (Xf*squeeze(filts(:,:,hh)));
                        
                        xs(ll,kk,:,oo,uu,hh) = xr;
                        ys(ll,kk,:,oo,uu,hh) = yr;
                        thetas(ll,kk,:,oo,uu,hh) = preftheta;
                        speeds(ll,kk,:,oo,uu,hh) = prefspeed;
                        bwthetas(ll,kk,:,oo,uu,hh)= bwtheta(oo);
                        bwspeeds(ll,kk,:,oo,uu,hh)= bwspeed(uu);
                        sigmass(ll,kk,:,oo,uu,hh)= sigmas(hh);
                    end
                end
            end
        end
    end
    params.xs = xs(:);
    params.ys = ys(:);
    params.thetas = thetas(:);
    params.speeds = speeds(:);
    params.sigmass = sigmass(:);
    params.bwspeeds = bwspeeds(:);
    params.bwthetas = bwthetas(:);
    
%     X = reshape(X,size(X,1),ns*8*3*length(bwtheta)*length(bwspeed)*length(sigmas));
    %X=single(X);