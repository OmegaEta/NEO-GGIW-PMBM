function [Bern,lik] = detectionBernfinal(Bern,C,model)

[Bern1.GGIW,lik1,lik1_rough] = updateGGIWforBern(Bern.GGIW,C,model);

% compare to NEOpartition
if size(C,2)>=3
    [Bern2.GGIW,lik2,lik2_rough] = updateGGIWforBernUseNEOpartition(Bern.GGIW,C,model);
else
    lik2=-inf;
end
if lik1>lik2
    Bern.GGIW=Bern1.GGIW;
    lik=lik1_rough;
else
    Bern.GGIW=Bern2.GGIW;
    lik=lik2_rough;
end

lik = lik + log(Bern.r) + log(model.Pd) + log(Bern.w_death(end));
Bern.r = 1;

Bern.t_death = Bern.t_death(end);
Bern.w_death = 1;

end
