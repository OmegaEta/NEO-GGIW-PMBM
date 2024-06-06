function [x,y] = Sigmacircle(cx,cy,P,N,format)

if nargin < 4
    N=1;
    format = '-b';
elseif nargin < 5
    format = '-b';
end
N2 = 20;

num = size(cx,1);
x=[];
y=[];

for i = 1 : num
    if size(P(:,:,i),1) == 2
        sqrtP = N*sqrtm_2by2(P(:,:,i));
    else
        sqrtP = N*sqrtm(P(:,:,i));
    end

    phi = 0:pi/N2:2*pi;

    xy = sqrtP*[cos(phi); sin(phi)];
    xy = xy + [cx(i)*ones(1,length(phi)) ; cy(i)*ones(1,length(phi))];

    x = [x;xy(1,:)';NaN];
    y = [y;xy(2,:)';NaN];
    
%     %plot(x,y,'k-','linewidth',1);
%     plot(x,y,'-','color',[0.7 0.7 0.7],'linewidth',1);
%     fill(xy(1,:),xy(2,:),[0.9 0.9 0.9]);
    
  
    
    

end
    %     c = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'w'];
    %    红   绿  蓝 青绿 洋红 黄  黑  白
    hold on;
    plot(x,y,format,'linewidth',1);
    %     pbaspect([1 1 1]);
    %     axis([-200 200 -200 200])
end


% function [x y] = Sigmacircle(cx,cy,P,N,N2)
% 
% if nargin < 4 %输入参数的个数
%     N=1;
%     N2 = 20;
% elseif nargin < 5
%     N2 = 20;
% end
% 
% if size(P,1) == 2
%     sqrtP = N*sqrtm_2by2(P);
% else
%     sqrtP = N*sqrtm(P);
% end
% 
% phi = 0:pi/N2:2*pi;
% xy = sqrtP*[cos(phi); sin(phi)];
% xy = xy + [cx*ones(1,length(phi)) ; cy*ones(1,length(phi))];
% 
% x = xy(1,:)';
% y = xy(2,:)';
% plot(x,y,'k-','linewidth',2);