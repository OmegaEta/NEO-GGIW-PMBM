function P=NEOpartition(W,num,N2)
%meaningless situation
if size(W,2)<num||size(W,2)<3
    P=cell(num,1);
    return;
end
% test code
% clc;clear;close all;
% dbstop if error
% P1=[6.79 -3.8278;-3.8278 2.3700];
% P2=[6.79 3.8278;3.8278 2.3700];
% r= ellipse_uniform_sampling(cat(3,P1,P2),[0 0;0 0],poissrnd(2*12));
% % P1=[5.8481 -1.5;-1.5 0.6519];
% % P2=[5.8481 1.5;1.5 0.6519];
% % r= ellipse_uniform_sampling(cat(3,P1,P2),[0 1;0 -1],poissrnd(2*12));
% 
% hold on;
% axis equal;
% 
% W=r;
% num=2;

N=3;
if nargin<3
    N2=10;
end

sigma=cov(W(1,:),W(2,:));
[~,S]=eig(sigma);
S=eye(2).*max(max(S));

mu=mean(W,2);

sqrtP=N*sqrtm_2by2(S);

phi = 0:pi/N2:2*pi-pi/N2;

xy = sqrtP*[cos(phi); sin(phi)];
xy = xy + [mu(1)*ones(1,length(phi)) ; mu(2)*ones(1,length(phi))];

x = xy(1,:);
y = xy(2,:);

% figure(1);
% plot([x x(1)],[y y(1)],'linewidth',1);
% scatter(W(1,:),W(2,:),'.');

pervector=zeros(size(phi,2),size(W,2));
for i=1:size(phi,2)
    for j=1:size(W,2)
        theta=acos((W(1,j)-xy(1,i))/norm(W(:,j)-xy(:,i)));
        if(W(2,j)-xy(2,i)<0)
            theta=2*pi-theta;
        end
        pervector(i,j)=theta;
    end
end

% figure(2);
% hold on;
% for j=1:size(W,2)
%     scatter(phi,pervector(:,j)','.');
% end

precision=0.005;
% precision=0.001;
gaussian_x=0:precision:(2-precision)*pi;
gaussian_kernel_mat=zeros(size(phi,2),size(gaussian_x,2));
gaussian_kernel_mat_component=zeros(size(phi,2),size(gaussian_x,2),size(W,2));
for i=1:size(phi,2)
    for j=1:size(W,2)
        gaussian_kernel_mat_component(i,:,j)=normpdf(gaussian_x,pervector(i,j),pi/40);
    end
end
gaussian_kernel_mat=sum(gaussian_kernel_mat_component,3);

% figure(3);
% hold on;
% for i=1:size(phi,2)
% plot(gaussian_x,gaussian_kernel_mat(i,:))
% end

del=zeros(1,size(W,2));
del_list=false(1,size(W,2));
prob_line=zeros(3,num);
for i=1:num
    [phii,vectorr]=find(gaussian_kernel_mat==max(max(gaussian_kernel_mat)),1);
    prob_line(:,i)=[xy(:,phii);gaussian_x(1,vectorr)];
    k=tan(prob_line(3,i));
    for j=1:size(W,2)
        if del_list(1,j)==0
            del(1,j)=(W(2,j)-prob_line(2,i)-k*(W(1,j)-prob_line(1,i)))^2/(1+k^2);
        else
            del(1,j)=inf;
        end
    end
    [~,idx]=sort(del);
    del_num=floor(size(W,2)/num);
    del_list(1,idx(1:del_num))=ones(1,del_num);

    gaussian_kernel_mat=sum(gaussian_kernel_mat_component(:,:,~del_list),3);
end

% color=['r--','b--'];
% for i=1:num
%     figure(1);
%     hold on;
%     k=tan(prob_line(3,i));
%     prob_x=-7:0.1:7;
%     prob_y=k*(prob_x-prob_line(1,i))+prob_line(2,i);
%     plot(prob_x,prob_y,color(i));
% end

%partitions and measurements for every partition cell 
P=cell(num,1);
meas_cell=cell(num,1);

%init partitions
%compute the distance of every measurement and every line
for i=1:size(prob_line,2)
    k=tan(prob_line(3,i));
    temp=W-repmat(prob_line(1:2,i),[1 size(W,2)]);
    square_distance(i,:)=((temp(2,:)-k*temp(1,:)).^2)/(k^2+1);
end
%assign each measurement to the nearest line
for i=1:size(W,2)
    [~,idx]=min(square_distance(:,i));
    meas_cell{idx,1}=[meas_cell{idx,1} i];
end
%init Lp
Lp=0;
meas_cell_list{1,1}=meas_cell;
gaussian_kernel_mat_now=zeros(size(phi,2),size(gaussian_x,2),num);
for i=1:num
    gaussian_kernel_mat_now(:,:,i)=sum(gaussian_kernel_mat_component(:,:,meas_cell{i,1}),3);
    Lp=Lp+max(max(gaussian_kernel_mat_now(:,:,i)));
end

%stochastic optimization
max_iterator=200;
num_repetition=0;
max_repetition=size(W,2);
P_set=cell(max_iterator,1);
L_set=zeros(max_iterator,1);
P_set{1,1}=meas_cell;
L_set(1,1)=Lp;
for iterator=1:max_iterator
    %randomly choose a measurement
    mea_selected=randi(size(W),1);
    %the cell index of chosen measurement and the
    %measurement cell after removing the chosen measurement
    Loc = cellfun(@(x) x==mea_selected,meas_cell,'Un',0);
    cell_idx = find(cellfun(@(x) any(x(:)),Loc));
    mea_after_move = meas_cell{cell_idx}(meas_cell{cell_idx}~=mea_selected);

    gaussian_kernel_mat_after_remove=gaussian_kernel_mat_now(:,:,cell_idx)-gaussian_kernel_mat_component(:,:,mea_selected);
    Lp_after_remove=max(max(gaussian_kernel_mat_after_remove));
    
    %compute Wp
    Wp=zeros(num+num,1);
    for i=1:length(meas_cell)
        if i==cell_idx
            Wp(i,1)=0;
        else
            temp=meas_cell;
            temp{i,1}=sort([meas_cell{i,1} mea_selected]);
            temp{cell_idx,1}=mea_after_move;
            gaussian_kernel_mat_after_add=zeros(size(phi,2),size(gaussian_x,2),length(meas_cell));
            if ~any(cellfun(@(x) isequal(x,temp),meas_cell_list))
                meas_cell_list{end+1,1}=temp;
                gaussian_kernel_mat_after_add(:,:,i)=gaussian_kernel_mat_now(:,:,i)+gaussian_kernel_mat_component(:,:,mea_selected);
                Lp_after_add=max(max(gaussian_kernel_mat_after_add(:,:,i)));
                Wp(i,1)=Lp_after_remove+Lp_after_add-Lp;
            else
                Wp(i,1)=-inf;
            end
        end
    end
    
    for i=1:length(meas_cell)
        if i==cell_idx
            Wp(i+num,1)=-inf;
        else
            temp=[meas_cell{cell_idx,1} meas_cell{i,1}];
            [IDX,~]=kmeanspp(W(:,[meas_cell{cell_idx,1} meas_cell{i,1}]),2);
            cell1{i,1}=temp(IDX==1);
            cell2{i,1}=temp(IDX==2);
            if isequal(cell1{i,1},meas_cell{cell_idx,1})||isequal(cell1{i,1},meas_cell{i,1})
                Wp(i+num,1)=-inf;
            else
                temp1=meas_cell;
                temp2=meas_cell;
                temp1{i,1}=cell1{i,1};
                temp1{cell_idx,1}=cell2{i,1};
                temp2{i,1}=cell2{i,1};
                temp2{cell_idx,1}=cell1{i,1};
                if ~any(cellfun(@(x) isequal(x,temp1)|isequal(x,temp2),meas_cell_list))
                    meas_cell_list{end+1,1}=temp1;
                    meas_cell_list{end+1,1}=temp2;
                    gaussian_kernel_mat_after_kmeans1=zeros(size(phi,2),size(gaussian_x,2),length(meas_cell));
                    gaussian_kernel_mat_after_kmeans2=zeros(size(phi,2),size(gaussian_x,2),length(meas_cell));

                    gaussian_kernel_mat_after_kmeans1(:,:,i)=sum(gaussian_kernel_mat_component(:,:,cell1{i,1}),3);
                    gaussian_kernel_mat_after_kmeans2(:,:,i)=sum(gaussian_kernel_mat_component(:,:,cell2{i,1}),3);
                    Lp_after_Kmeans1=max(max(gaussian_kernel_mat_after_kmeans1(:,:,i)));
                    Lp_after_Kmeans2=max(max(gaussian_kernel_mat_after_kmeans2(:,:,i)));
                    Wp(i+num,1)=Lp_after_Kmeans1+Lp_after_Kmeans2-Lp;
                else
                    Wp(i+num,1)=-inf;
                end
            end
        end
    end
    
    [wp,~] = normalizeLogWeights(Wp);
%     I = find(cumsum(exp(wp))>rand(),1); %sampling
    [~,I]=max(wp);
    if I==cell_idx
        num_repetition = num_repetition + 1;
        if (num_repetition > max_repetition)
            P_set{iterator+1,1}=meas_cell;
            L_set(iterator+1,1)=Lp;
            break;
        end
    elseif I<=num
        Lp=Wp(I,1)+Lp;
        meas_cell{I,1}=[meas_cell{I,1} mea_selected];
        meas_cell{cell_idx,1}=mea_after_move;
        gaussian_kernel_mat_now(:,:,I)=gaussian_kernel_mat_after_add(:,:,I);
        gaussian_kernel_mat_now(:,:,cell_idx)=gaussian_kernel_mat_after_remove;
    elseif I<=2*num
        Lp=Wp(I,1)+Lp;
        meas_cell{I-num,1}=cell1{I-num,1};
        meas_cell{cell_idx,1}=cell2{I-num,1};
        gaussian_kernel_mat_now(:,:,I-num)=gaussian_kernel_mat_after_kmeans1(:,:,I-num);
        gaussian_kernel_mat_now(:,:,cell_idx)=gaussian_kernel_mat_after_kmeans2(:,:,I-num);
    end
        
    if I~=cell_idx
    num_repetition = 0;
    end
    
    P_set{iterator+1,1}=meas_cell;
    L_set(iterator+1,1)=Lp;
end

% [~,idx]=max(L_set);
% P=P_set{idx,1};
P=meas_cell;


prob_line_after=zeros(3,num);
for i=1:num
%     temp=zeros(size(phi,2),size(gaussian_x,2));
%     for j=1:size(phi,2)
%         temp(j,:)=sum(gaussian_kernel_mat_component(P{i,1},:,j),1);
%     end
    temp=sum(gaussian_kernel_mat_component(:,:,P{i,1}),3);
    [phii,vectorr]=find(temp==max(max(temp)),1);
    prob_line_after(:,i)=[xy(:,phii);gaussian_x(1,vectorr)];
%     figure(1);
%     hold on;
%     k=tan(prob_line_after(3,i));
%     prob_x=-7:0.1:7;
%     prob_y=k*(prob_x-prob_line_after(1,i))+prob_line_after(2,i);
%     plot(prob_x,prob_y,'b');
%     if i==1
%         scatter(W(1,P{i,1}),W(2,P{i,1}),'bx');
%     else
%         scatter(W(1,P{i,1}),W(2,P{i,1}),'ro');
%     end
end

% figure(4);
% plot(1:length(L_set),L_set);

% least squares method for optimization
for i=1:size(W,2)
    distance=Inf(size(prob_line_after,2),1);
    for j=1:size(prob_line_after,2)
        k=tan(prob_line_after(3,j));
        distance(j)=((W(2,i)-prob_line_after(2,j)-k*(W(1,i)-prob_line_after(1,j))).^2)/(k^2+1);
    end
    [~,idx]=min(distance);
    Loc = cellfun(@(x) x==i,P,'Un',0);
    cell_idx = find(cellfun(@(x) any(x(:)),Loc));
    if idx~=cell_idx&&length(P{cell_idx,1})>=3
        P{cell_idx,1}(P{cell_idx,1}==i)=[];
        P{idx,1}=sort([P{idx,1} i]);
    end
end

% figure(1);
% hold on;
% for i=1:num
%         if i==1
%         scatter(W(1,P{i,1}),W(2,P{i,1}),'bh');
%     else
%         scatter(W(1,P{i,1}),W(2,P{i,1}),'rh');
%     end
% end

end