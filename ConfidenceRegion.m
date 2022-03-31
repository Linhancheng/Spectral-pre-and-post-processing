function ConfidenceRegion(xdat,alpha,distribution)
%   ConfidenceRegion(xdat,alpha,distribution)
%   xdat£ºp*N or N*p£¬p = 1,2 or 3
%   alpha£º[0,1]£¬default: 0.05
%   distribution£º'norm'or 'experience'
%   CopyRight£ºxiezhh
%   date£º2011.4.14
%
%   example1£ºx = normrnd(10,4,100,1);
%             ConfidenceRegion(x)
%             ConfidenceRegion(x,'exp')
%   example2£ºx = mvnrnd([1;2],[1 4;4 25],100);
%             ConfidenceRegion(x)
%             ConfidenceRegion(x,'exp')
%   example3£ºx = mvnrnd([1;2;3],[3 0 0;0 5 -1;0 -1 1],100);
%             ConfidenceRegion(x)
%             ConfidenceRegion(x,'exp')
% 
if nargin == 1
    distribution = 'norm';
    alpha = 0.05;
elseif nargin == 2
    if ischar(alpha)
        distribution = alpha;
        alpha = 0.05;
    else
        distribution = 'norm';
    end
end

if ~isscalar(alpha) || alpha>=1 || alpha<=0
    error('')
end
if ~strncmpi(distribution,'norm',3) && ~strncmpi(distribution,'experience',3)
    error('')
end

[m,n] = size(xdat);
p = min(m,n);  
if ~ismember(p,[1,2,3])
    error('')
end

if m < n
    xdat = xdat';
end
xm = mean(xdat); 
n = max(m,n);  

switch p
    case 1    
        xstd = std(xdat);
        if strncmpi(distribution,'norm',3)
            lo = xm - xstd*norminv(1-alpha/2,0,1); 
            up = xm + xstd*norminv(1-alpha/2,0,1);
            
        else
            lo = prctile(xdat,100*alpha/2); 
            up = prctile(xdat,100*(1-alpha/2)); 
            
        end
        
        xin = xdat(xdat>=lo & xdat<=up);
        xid = find(xdat>=lo & xdat<=up);
        plot(xid,xin,'.')
        hold on
        xout = xdat(xdat<lo | xdat>up);
        xid = find(xdat<lo | xdat>up);
        plot(xid,xout,'r+')
        h = refline([0,lo]);
        set(h,'color','k','linewidth',2)
        h = refline([0,up]);
        set(h,'color','k','linewidth',2)
        xr = xlim;
        yr = ylim;
        text(xr(1)+range(xr)/20,lo-range(yr)/20,'lower limit',...
            'color','g','FontSize',15,'FontWeight','bold')
        text(xr(1)+range(xr)/20,up+range(yr)/20,'upper limit',...
            'color','g','FontSize',15,'FontWeight','bold')
        title(TitleText)
        hold off
    case 2    
        x = xdat(:,1);
        y = xdat(:,2);
        s = inv(cov(xdat));  
        xd = xdat-repmat(xm,[n,1]);
        rd = sum(xd*s.*xd,2);
        if strncmpi(distribution,'norm',3)
            r = chi2inv(1-alpha,p);
            %r = p*(n-1)*finv(1-alpha,p,n-p)/(n-p)/n;
            TitleText = ' ';
        else
            r = prctile(rd,100*(1-alpha));
            TitleText = ' ';
        end
        plot(x(rd<=r),y(rd<=r),'.','MarkerSize',16) 
        hold on
        plot(x(rd>r),y(rd>r),'r+','MarkerSize',10)  
        plot(xm(1),xm(2),'k+'); 
        h = ellipsefig(xm,s,r,1); 
        xlabel('X')
        ylabel('Y')
        title(TitleText)
        hold off;
    case 3    
        x = xdat(:,1);
        y = xdat(:,2);
        z = xdat(:,3);
        s = inv(cov(xdat));  
        xd = xdat-repmat(xm,[n,1]);
        rd = sum(xd*s.*xd,2);
        if strncmpi(distribution,'norm',3)
            r = chi2inv(1-alpha,p);
            %r = p*(n-1)*finv(1-alpha,p,n-p)/(n-p)/n;
            TitleText = ' ';
        else
            r = prctile(rd,100*(1-alpha));
            TitleText = ' ';
        end
        plot3(x(rd<=r),y(rd<=r),z(rd<=r),'.','MarkerSize',16)  
        hold on
        plot3(x(rd>r),y(rd>r),z(rd>r),'r+','MarkerSize',10)  
        plot3(xm(1),xm(2),xm(3),'k+');  
        h = ellipsefig(xm,s,r,2);  
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        title(TitleText)
        hidden off;
        hold off;
end

function  h = ellipsefig(xc,P,r,tag)

[V, D] = eig(P);
if tag == 1
    aa = sqrt(r/D(1));
    bb = sqrt(r/D(4));
    t = linspace(0, 2*pi, 60);
    xy = V*[aa*cos(t);bb*sin(t)];  
    h =plot(xy(1,:)+xc(1),xy(2,:)+xc(2), 'k', 'linewidth', 2);
else
    aa = sqrt(r/D(1,1));
    bb = sqrt(r/D(2,2));
    cc = sqrt(r/D(3,3));
    [u,v] = meshgrid(linspace(-pi,pi,30),linspace(0,2*pi,30));
    x = aa*cos(u).*cos(v);
    y = bb*cos(u).*sin(v);
    z = cc*sin(u);
    xyz = V*[x(:)';y(:)';z(:)'];  
    x = reshape(xyz(1,:),size(x))+xc(1);
    y = reshape(xyz(2,:),size(y))+xc(2);
    z = reshape(xyz(3,:),size(z))+xc(3);
    h = mesh(x,y,z);  
end