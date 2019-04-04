function [] = polarmosaic(y,rg,r1,r2)
    precision = 8;
    if nargin < 2
        themin = min(y);
        themax = max(y);
    else
        themin = rg(1);
        themax = rg(2);
    end
    %y = (y - themin)/(themax-themin);
    y = max(min(y,themax),themin);
    
    x = ((0:length(y))-.5)/length(y)*2*pi;
    
    xs = zeros(4*precision,length(y));
    ys = zeros(4*precision,length(y));
    for ii = 1:length(y)
        %get wedge coordinates
        [xc,yc] = wedgecoords(x(ii),x(ii+1),r1,r2,precision);
        xs(:,ii) = xc;
        ys(:,ii) = yc;
    end
            
    %Plot as a patch
    %size(y)
    h = fill(xs,ys,y','Clipping','off','LineStyle','none');
    %get(h,'Faces')
    
    thesel = (1:precision)';
    
    xint = xs(thesel ,:);
    xint = xint(:);
    yint = ys(thesel ,:);
    yint = yint(:);
    
    h = line([xint;xint(1)],[yint;yint(1)]);
    set(h,'Color',[0.9,0.9,0.9]*.9,'LineWidth',1);
    
    thesel = 2*precision + (precision:-1:1)';
    
    
    xint = xs(thesel,:);
    xint = xint(:);
    yint = ys(thesel ,:);
    yint = yint(:);
    
    h = line([xint;xint(1)],[yint;yint(1)]);
    set(h,'Color',[0.9,0.9,0.9]*.9,'LineWidth',1,'Clipping','off');
    set(gca,'CLim',[themin,themax],'Clipping','off');
end

function [xs,ys] = wedgecoords(t1,t2,r1,r2,precision)
   
    ts = [t1,t2,t2,t1,t1];
    rs = [r1,r1,r2,r2,r1];
    
    xs = [];
    ys = [];
    
    for ii = 1:4
        tb = ts(ii+1);
        ta = ts(ii);
        rb = rs(ii+1);
        ra = rs(ii);
        dr = rb-ra;
        dt = tb-ta;
        for jj = (1:precision)/precision
            rn = ra+jj*dr;
            tn = ta+jj*dt;
            xs = [xs;rn*cos(tn)];
            ys = [ys,rn*sin(tn)];
        end
    end
end