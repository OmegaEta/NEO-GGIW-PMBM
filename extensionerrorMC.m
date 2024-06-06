function [extensionerror]=extensionerrorMC(x,X,y,Y)
% Return the sum of non overlapping areas;
% x,X are the estimated center and covariance matrix of single target
% y,Y are the true center and covariance matrix of single target
% this function use Monte Carlo
% all rights reserved

%num of shape
Nx=size(x,2);
Ny=size(y,2);

%parameter
Npi=20;%3.14/20
Npoint=1e5;%iteration number of MC
N=2;

%preprocess
mu=[x y];
sigma=cat(3,X,Y);
xx=[];yy=[];
for i=1:Nx+Ny
    sqrtP=N*sqrtm_2by2(sigma(:,:,i));
    theta=0:pi/Npi:2*pi;
    
    xy = sqrtP*[cos(theta); sin(theta)];
    xy = xy + [mu(1,i)*ones(1,length(theta)) ; mu(2,i)*ones(1,length(theta))];
    
    xx = [xx;xy(1,:)';NaN];
    yy = [yy;xy(2,:)';NaN];
end

%range of MC
x_range=[min(xx)-1 max(xx)+1];
y_range=[min(yy)-1 max(yy)+1];
xlen=x_range(2)-x_range(1);
ylen=y_range(2)-y_range(1);

%MC
W=([xlen 0;0 ylen]*rand(2,Npoint))+repmat([x_range(1);y_range(1)],[1,Npoint]);
if_inside=zeros(2,Npoint);
for i=1:Nx
    %1 if in estimated target; 0 if not
    temp=W-repmat(x(:,i),[1 Npoint]);
    [R,S]=eig(X(:,:,i));
    S=sqrtm_2by2(S);iR=inv(R);iS=inv(S);
    temp=1/N*iS*iR*temp;
    if_inside(1,:)=((temp(1,:).^2+temp(2,:).^2)<=1)|if_inside(1,:);
end
for i=1:Ny
    %1 if in true target; 0 if not
    temp=W-repmat(y(:,i),[1 Npoint]);
    [R,S]=eig(Y(:,:,i));
    S=sqrtm_2by2(S);iR=inv(R);iS=inv(S);
    temp=1/N*iS*iR*temp;
    if_inside(2,:)=((temp(1,:).^2+temp(2,:).^2)<=1)|if_inside(2,:);
end
MC_result=xor(if_inside(1,:),if_inside(2,:));

extensionerror=sum(MC_result)/Npoint*xlen*ylen;

% draw
% hold on;
% plot(xx,yy,'linewidth',1);
% scatter(W(1,MC_result),W(2,MC_result),'.b');
end