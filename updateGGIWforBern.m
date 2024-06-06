function [GGIW,lik,lik_rough] = updateGGIWforBern(GGIW,W,model)

d = 2;

card_W = size(W,2);

GGIW_.a = GGIW(end).a + card_W;
GGIW_.b = GGIW(end).b + 1;

z_bar = mean(W,2);
epsilon = z_bar - model.measmodel.h(GGIW(end).m);
H = model.measmodel.H(GGIW(end).m);

X_hat = GGIW(end).V/(GGIW(end).v - 2*d - 2);
X_hat = (X_hat + X_hat')/2;

S = H*GGIW(end).P*H' + X_hat/card_W;
S = (S + S')/2;

Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';

K = GGIW(end).P*H'*iS;

GGIW_.m = GGIW(end).m + K*epsilon;
GGIW_.P = GGIW(end).P - K*H*GGIW(end).P;

temp = (W - z_bar);
Z = temp*temp';

X_sqrt = sqrtm_2by2(X_hat);
S_sqrt_inv = sqrtm(iS);
N = X_sqrt*S_sqrt_inv*(epsilon*epsilon')*S_sqrt_inv'*X_sqrt';

GGIW_.v = GGIW(end).v + card_W;
GGIW_.V = GGIW(end).V + N + Z;

lik_rough = (GGIW(end).v-d-1)/2*log(det2(GGIW(end).V)) - (GGIW_.v-d-1)/2*log(det2(GGIW_.V))...
        + gamma2ln((GGIW_.v-d-1)/2) - gamma2ln((GGIW(end).v-d-1)/2)...
        + log(det2(X_hat))/2 - log(det_S)/2 + gammaln(GGIW_.a)...
        - gammaln(GGIW(end).a) + GGIW(end).a*log(GGIW(end).b) - GGIW_.a*log(GGIW_.b)...
        - (card_W*log(pi)+log(card_W))*d/2;

% predict orientation
GGIW(end).shape=predictShapeRotate(GGIW(end).m(3:4),GGIW_.m(3:4),GGIW(end).shape);

% partition based on closest distance
Z=W-repmat(model.measmodel.h(GGIW_.m),[1,card_W]);
W=W-repmat(z_bar,[1,card_W]);
J=cell(length(GGIW(end).shape),1);
aver_m=[0;0];
for i=1:length(GGIW(end).shape)
    aver_m=aver_m+GGIW(end).shape(i).m(1:2,1);
end
aver_m=aver_m/length(GGIW(end).shape);

% precompute subobjects to reduce computational complexity
S=zeros(2,2,length(GGIW(end).shape));
for i=1:length(GGIW(end).shape)
    H = model.measmodel.H(GGIW(end).shape(i).m);
    Stemp = GGIW(end).shape(i).V/(GGIW(end).shape(i).v-2*d-2) + H*GGIW(end).shape(i).P*H' + model.measmodel.R;
    S(:,:,i) = (Stemp + Stemp')/2;
end

for i=1:length(GGIW(end).shape)
    nu=W-(GGIW(end).shape(i).m(1:2,1)-aver_m);
    dist(i,:)= sum((inv(chol(S(:,:,i)))'*nu).^2);
end
for i=1:card_W
    idx=find(dist(:,i)==min(dist(:,i)));
    I=randperm(length(idx),1);
    J{idx(I),1}=[J{idx(I),1} i];
end


lik=0;

for i=1:length(GGIW(end).shape)
    if isempty(J{i,1})
        [GGIW_.shape(i,1),lik_]=misdetectionGGIW(GGIW(end).shape(i));
    else
        [GGIW_.shape(i,1),lik_]=updateGGIW(GGIW(end).shape(i),Z(:,J{i,1}),model);
    end
    lik=lik+lik_;
end

GGIW(end) = GGIW_;
end

