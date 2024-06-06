function GGIW_=initShape(GGIW,W,num)
% NEOpartition
J=NEOpartition(W,num,30);
%
shape=repmat(struct('a',0,'b',0,'m',zeros(4,1),'P',zeros(4,4),'v',0,'V',[]),[num,1]);
% init shape
for i=1:num
    if isempty(J{i,1})||size(J{i,1},2)<5
        if isempty(J{i,1})
            z_bar = mean(W,2);
        else
            z_bar = mean(W(:,J{i,1}),2);
        end
        
        shape(i).a=100;
        shape(i).b=10;
        shape(i).m=[z_bar-GGIW.m(1:2,1);0;0];
        shape(i).P=diag([4,4,4,4]);
        shape(i).v=22.5;
        shape(i).V=32.5*eye(2);
    else
        card_W = size(W(:,J{i,1}),2);
        z_bar = mean(W(:,J{i,1}),2);
        temp = (W(:,J{i,1}) - z_bar);
        Z = temp*temp';
    
        shape(i).a=card_W;
        shape(i).b=1;
        shape(i).m=[z_bar-GGIW.m(1:2,1);0;0];
        shape(i).P=diag([2,2,2,2]);
        shape(i).v=card_W;
        shape(i).V=Z;
        
        if shape(i).v<7
            shape(i).v=7;
        end
    end
end
% return GGIW with shapes
GGIW_=GGIW;
GGIW_.shape=shape;
end