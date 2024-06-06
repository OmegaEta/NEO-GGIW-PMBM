function [GGIW,lik]=misdetectionGGIW(GGIW)

lik=(GGIW.b/(GGIW.b+1))^GGIW.a;
lik=log(lik);

GGIW.b=GGIW.b+1;
end