TITLE={'gospa','ExtensionError','LocationError','MissedError','FalseError'};
figure(1)
% preprocess estimates
chosenMC=1;
cx=[];cy=[];P=zeros(2,2,0);
for i=1:length(estimates{chosenMC,1})
    num_target=length(estimates{chosenMC,1}{i}.g);
    for j=1:num_target
        cx=[cx;(estimates{chosenMC,1}{i}.xshape(1,:,j)+estimates{chosenMC,1}{i}.x(1,j))'];
        cy=[cy;(estimates{chosenMC,1}{i}.xshape(2,:,j)+estimates{chosenMC,1}{i}.x(2,j))'];
        P=cat(3,P,estimates{chosenMC,1}{i}.Xshape(:,:,:,j));
    end
end
% preprocess estimates_original
cx_original=[];cy_original=[];P_original=zeros(2,2,0);
for i=1:length(estimates_original{chosenMC,1})
    cx_original=[cx_original;estimates_original{chosenMC,1}{i}.x(1,:)'];
    cy_original=[cy_original;estimates_original{chosenMC,1}{i}.x(2,:)'];
    P_original=cat(3,P_original,estimates_original{chosenMC,1}{i}.X);
end
% cx_original=[];cy_original=[];P_original=zeros(2,2,0);
% for i=1:length(estimates{chosenMC,1})
%     cx_original=[cx_original;estimates{chosenMC,1}{i}.x(1,:)'];
%     cy_original=[cy_original;estimates{chosenMC,1}{i}.x(2,:)'];
%     P_original=cat(3,P_original,estimates{chosenMC,1}{i}.X);
% end
% preprocess Z
Z=[];
for i=1:length(Scenario.Z{chosenMC,1})
    Z=[Z Scenario.Z{chosenMC,1}{i,1}];
end
%draw
Sigmacircle(cx,cy,P,2,'-r');
Sigmacircle(cx_original,cy_original,P_original,2,'--b');
scatter(Z(1,:),Z(2,:),'.k');
xlabel('X (m)');
ylabel('Y (m)');
box on;
axis([-200 200 -200 200]);
for i=1:5
figure(i+1);
hold on;
box on;
if i==3
    plot(1:K,1/2*GOSPA02(:,i),'-r');
else
    plot(1:K,GOSPA02(:,i),'-r');
end
plot(1:K,GOSPA02_original(:,i),'--b');
xlabel('Time (s)');
ylabel('value');
legend('NEO-GGIW-PMBM','GGIW-PMBM');
title(TITLE(i))
end
figure(7);
hold on;
plot(1:K,simulation_time02,'-or')
plot(1:K,simulation_time_original02,'--*b');
xlabel('Time');
ylabel('Time (s)');
legend('NEO-GGIW-PMBM','GGIW-PMBM');
title('time cost')