function track = pseudoTrack(GGIW)

wp=repmat(1/length(GGIW),[length(GGIW) 1]);

[w_hat,Bern.GGIW] = GGIW_merge_wrap(wp,GGIW);

Bern.r = 1;
lik = log(w_hat);

Bern.w_death = 1;

Bern.t_death=1;

track = struct('Bern',Bern,'lik',lik,'assocHistory',[]);
end

