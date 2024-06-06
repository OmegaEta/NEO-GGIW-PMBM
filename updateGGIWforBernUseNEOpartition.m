function [GGIW,lik,lik_rough] = updateGGIWforBernUseNEOpartition(GGIW,W,model)

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
% adjust the shape because the shapes are relative to the rough ellipse
GGIW(end).shape=predictShapeRotate(GGIW(end).m(3:4),GGIW_.m(3:4),GGIW(end).shape);

% NEOpartition
Z=W-repmat(model.measmodel.h(GGIW_.m),[1,card_W]);
J=NEOpartition(Z,2);

% update shape and compute lik

lik1=0;
lik2=0;
for i=1:length(GGIW(end).shape)
    if isempty(J{i,1})
        [GGIW1.shape(i,1),lik1_]=misdetectionGGIW(GGIW(end).shape(i));
    else
        [GGIW1.shape(i,1),lik1_]=updateGGIW(GGIW(end).shape(i),Z(:,J{i,1}),model);
    end
    lik1=lik1+lik1_;
end
for i=1:length(GGIW(end).shape)
    if isempty(J{3-i,1})
        [GGIW2.shape(i,1),lik2_]=misdetectionGGIW(GGIW(end).shape(i));
    else
        [GGIW2.shape(i,1),lik2_]=updateGGIW(GGIW(end).shape(i),Z(:,J{3-i,1}),model);
    end
    lik2=lik2+lik2_;
end
if lik1>lik2
    for i=1:length(GGIW(end).shape)
        GGIW_.shape(i,1)=GGIW1.shape(i,1);
    end
    lik=lik1;
else
    for i=1:length(GGIW(end).shape)
        GGIW_.shape(i,1)=GGIW2.shape(i,1);
    end
    lik=lik2;
end

GGIW(end) = GGIW_;
end

