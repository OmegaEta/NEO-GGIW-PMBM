function [partitions,lik,new_merge_table] = ObjectsSO(MBM,j,PPP,W,gating_matrix_d,gating_matrix_u,model)

m = size(W,2);

%Pre-store the results of detectionPPP update to avoid duplicate computation
%%%%%%%%%%%%%

Lstore = -inf;
Mstore = false(1,m);

%%%%%%%%%%%%%

track_indices = find(MBM.table(j,:)>0);     %extract the jth global hypothesis
n_tt = length(track_indices);               %number of tracks in the jth global hypothesis
gating_matrix_dj = zeros(m,n_tt);           %extract gating matrix for each track

merge_table=MBM.merge_table{j};         %extract the jth table of merged tracks
if size(merge_table,2)<n_tt
    merge_table(1,n_tt)=0;%ensure the validity of the length of merge_table 
end
prob_merge_mat=zeros(n_tt,n_tt);           %matrix for tracks which may merge with each other

n_merged_cell=sum(merge_table>0);       %number of the tracks which are merged
n_merge_cell=sum(unique(merge_table)~=0);%number of the tracks which are generated through merging processing
n_tt_now=n_tt-n_merged_cell+n_merge_cell;%number of tracks ObjectsSO will truly process in the following

idx_tracks2n_tt=cell(n_tt_now,1);%the table to restore the new tracks' index to the actual tracks' index 
idx_n_tt2tracks=zeros(1,n_tt);%the table to get the new tracks' index by the actual tracks' index 

%get tracks¡¢gating_matrix_dj¡¢idx_tracks2n_tt¡¢idx_n_tt2tracks
idx=1;
for i = 1:n_tt
    tracks_truth(i,1) = MBM.track{track_indices(i)}(MBM.table(j,track_indices(i)));
    if merge_table(i)==0
        gating_matrix_dj(:,idx) = gating_matrix_d{track_indices(i)}(:,MBM.table(j,track_indices(i)));
        tracks(idx,1) = MBM.track{track_indices(i)}(MBM.table(j,track_indices(i)));
        idx_tracks2n_tt{idx,1}=[i];
        idx_n_tt2tracks(1,i)=idx;
        idx=idx+1;
    end
end
for i=1:n_merge_cell
    GGIWs=[];
    for k=find(merge_table==i)
        gating_matrix_dj(:,idx) = gating_matrix_dj(:,idx)|gating_matrix_d{track_indices(k)}(:,MBM.table(j,track_indices(k)));
        GGIWs=[GGIWs;MBM.track{track_indices(k)}(MBM.table(j,track_indices(k))).Bern.GGIW(end)];
    end
    tracks(idx,1) = pseudoTrack(GGIWs);
    idx_tracks2n_tt{idx,1}=find(merge_table==i);
    idx_n_tt2tracks(1,merge_table==i)=idx;
    idx=idx+1;
end

%initialize prob_merge_mat by merge_table
if any(merge_table)
    for i=1:n_merge_cell
        for ii=find(merge_table==i)
            for jj=find(merge_table==i)
                if ii~=jj
                    prob_merge_mat(ii,jj)=1;
                end
            end
        end
    end
end

%Pre-computing likelihood of new measurement being a cell of its own
Lc = zeros(m,1);
for i = 1:m
    if any(gating_matrix_u(i,:))
        [~,Lc(i)] = detectionPPP(PPP.w(gating_matrix_u(i,:)),...
            PPP.GGIW(gating_matrix_u(i,:)),W(:,i),model);
    else
        Lc(i) = log(model.lambda_fa);
    end
end

%Pre-computing misdetection likelihood
Lmiss = zeros(n_tt_now,1);
for i = 1:n_tt_now
    [~,Lmiss(i)] = misdetectionBern(tracks(i,1).Bern,model);
end

%Maximum likelihood initialisation
nu = length(PPP.w);

meas_cell = cell(n_tt_now+nu,1);
ini_likelihood = -inf(n_tt_now+1,m);
for i = 1:n_tt_now
    ghat = tracks(i,1).Bern.GGIW(end).a/tracks(i,1).Bern.GGIW(end).b;
    zhat = model.measmodel.h(tracks(i,1).Bern.GGIW(end).m);
    Shat = model.measmodel.H(tracks(i,1).Bern.GGIW(end).m)*tracks(i,1).Bern.GGIW(end).P*model.measmodel.H(tracks(i,1).Bern.GGIW(end).m)' ...
        + model.measmodel.R + tracks(i,1).Bern.GGIW(end).V/(tracks(i,1).Bern.GGIW(end).v-2*2-2);
    Shat = (Shat + Shat')/2;
    log_pD_ghat_r = log(model.Pd) + log(ghat) + log(tracks(i,1).Bern.r);
    for j = 1:m
        if gating_matrix_dj(j,i)
            ini_likelihood(i,j) = log_pD_ghat_r + log_mvnpdf(W(:,j),zhat,Shat);
        end
    end
end
ini_likelihood(n_tt_now+1,:) = Lc';

l = 1;
for j = 1:m
    [~,idx] = max(ini_likelihood(:,j));
    if idx <= n_tt_now
        meas_cell{idx} = [meas_cell{idx} j];
    else
        idx_u = find(gating_matrix_u(j,:)==true,1);
        if ~isempty(idx_u)
            meas_cell{n_tt_now+idx_u,1} = [meas_cell{n_tt_now+idx_u,1} j];
        else
            meas_cell{n_tt_now+nu+l,1} = j;
            l = l + 1;
        end
    end
end

idx = ~cellfun('isempty',meas_cell);
idx(1:n_tt_now) = true;
meas_cell =  meas_cell(idx);

N = length(meas_cell);
Lp = zeros(N,1);
for i = 1:n_tt_now
    if isempty(meas_cell{i})
        Lp(i) = Lmiss(i);
    else
        [~,Lp(i)] = detectionBern(tracks(i,1).Bern,W(:,meas_cell{i}),model);
    end
end
for i = n_tt_now+1:N
    in_gate = sum(gating_matrix_u(meas_cell{i},:),1)>=1;
    if any(in_gate)
        [~,Lp(i)] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,meas_cell{i}),model);
    else
        Lp(i) = log(model.lambda_fa);
    end
end


T = model.num_iterations*(m+n_tt_now);     %number of iterations
partitions = cell(T,1);
partitions{1} = meas_cell;
lik = zeros(T,1);
lik(1) = sum(Lp);

max_repetition = max(20,ceil(model.max_repetition*(m+n_tt_now)/2));
num_repetition = 0;
%Stochastic optimisation
for t = 1:T
%     num_iteration=num_iteration+1;
    %randomly select a measurement
    mea_selected = randi(m,1);
    %find the corresponding cell index
    Loc = cellfun(@(x) x==mea_selected,meas_cell,'Un',0);
    cell_idx = find(cellfun(@(x) any(x(:)),Loc));
    N = length(meas_cell);
    
    mea_after_move = meas_cell{cell_idx}(meas_cell{cell_idx}~=mea_selected);
    Wp = -inf(2*N+4,1);
    Lp_original = Lp;
    meas_cell_original = meas_cell;
    
    %selected measurement moved to an existing cell
    Lp_move1 = -inf(N,1);
    if cell_idx <=n_tt_now %selected measurement from a target cell
        if isempty(mea_after_move)
            Lp_move2 = Lmiss(cell_idx);
        else
            [~,Lp_move2] = detectionBern(tracks(cell_idx).Bern,W(:,mea_after_move),model);
        end
    else %selected measurement from a new cell
        if isempty(mea_after_move)
            Lp_move2 = 0;
        else
            if length(mea_after_move) == 1
                Lp_move2 = Lc(mea_after_move);
            else
                in_gate = sum(gating_matrix_u(mea_after_move,:),1)>=1;
                if any(in_gate)
                    mstore = ismember(1:m,mea_after_move);
                    idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                    if ~any(idx)
                        [~,Lp_move2] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,mea_after_move),model);
                        Mstore = [Mstore;mstore];
                        Lstore = [Lstore;Lp_move2];
                    else
                        Lp_move2 = Lstore(idx);
                    end
                else
                    Lp_move2 = -inf;
                end
            end
        end
    end
    for i = 1:N
        if i==cell_idx
            Wp(i) = 0;
        else
            if i <= n_tt_now
                if gating_matrix_dj(mea_selected,i)
                    [~,Lp_move1(i)] = detectionBern(tracks(i).Bern,W(:,[meas_cell{i} mea_selected]),model);
                end
            else
                in_gate = sum(gating_matrix_u([meas_cell{i} mea_selected],:),1)>=1;
                if any(in_gate)
                    mstore = ismember(1:m,[meas_cell{i} mea_selected]);
                    idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                    if ~any(idx)
                        [~,Lp_move1(i)] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,[meas_cell{i} mea_selected]),model);
                        Mstore = [Mstore;mstore];
                        Lstore = [Lstore;Lp_move1(i)];
                    else
                        Lp_move1(i) = Lstore(idx);
                    end
                end
            end
            Wp(i) = Lp_move1(i) + Lp_move2 - Lp(cell_idx) - Lp(i);
        end
    end
    
    %selected measurement moved to a new cell
    if cell_idx <=n_tt_now
        if size(idx_tracks2n_tt{cell_idx,1},2)==1
            Wp(N+1) = Lp_move2 + Lc(mea_selected) - Lp(cell_idx);
        else
            Wp(N+1) = -inf;
        end
    elseif cell_idx > n_tt_now && ~isempty(mea_after_move)
        Wp(N+1) = Lp_move2 + Lc(mea_selected) - Lp(cell_idx);
    end
    
    if length(meas_cell{cell_idx}) > 1
        %selected cell merged with an existing cell
        Lp_merge = -inf(N,1);
        if cell_idx <=n_tt_now %selected cell is a target cell
            %selected cell becomes a new cell
            in_gate = sum(gating_matrix_u(meas_cell{cell_idx},:),1)>=1;
            if any(in_gate)
                mstore = ismember(1:m,meas_cell{cell_idx});
                idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                if ~any(idx)
                    [~,Lp_new] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,meas_cell{cell_idx}),model);
                    Mstore = [Mstore;mstore];
                    Lstore = [Lstore;Lp_new];
                else
                    Lp_new = Lstore(idx);
                end
            else
                Lp_new = -inf;
            end
            if size(idx_tracks2n_tt{cell_idx,1},2)==1
                Wp(N+2) = Lmiss(cell_idx) + Lp_new - Lp(cell_idx);
            else
                Wp(N+2) = -inf;
            end
            
            for i = 1:n_tt_now %merged with a target cell
                if i~=cell_idx
                    if any(gating_matrix_dj([meas_cell{cell_idx}],i))
                        [~,Lp_merge(i)] = detectionBern(tracks(i).Bern,W(:,[meas_cell{cell_idx} meas_cell{i}]),model);
                    end
                    Wp(i+N+2) = Lp_merge(i) + Lmiss(cell_idx) - Lp(cell_idx) - Lp(i);
                end
            end
            for i = n_tt_now+1:N %merged with a new cell
                in_gate = sum(gating_matrix_u([meas_cell{cell_idx} meas_cell{i}],:),1)>=1;
                if any(in_gate)
                    mstore = ismember(1:m,[meas_cell{cell_idx} meas_cell{i}]);
                    idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                    if ~any(idx)
                        [~,Lp_merge(i)] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,[meas_cell{cell_idx} meas_cell{i}]),model);
                        Mstore = [Mstore;mstore];
                        Lstore = [Lstore;Lp_merge(i)];
                    else
                        Lp_merge(i) = Lstore(idx);
                    end
                end
                Wp(i+N+2) = Lp_merge(i) + Lmiss(cell_idx) - Lp(cell_idx) - Lp(i);
            end
        else %selected cell is a new cell
            
            for i = 1:n_tt_now %merged with a target cell
                if any(gating_matrix_dj([meas_cell{cell_idx}],i))
                    [~,Lp_merge(i)] = detectionBern(tracks(i).Bern,W(:,[meas_cell{cell_idx} meas_cell{i}]),model);
                end
                Wp(i+N+2) = Lp_merge(i) - Lp(cell_idx) - Lp(i);
            end
            for i = n_tt_now+1:N %merged with a new cell
                if i~=cell_idx
                    in_gate = sum(gating_matrix_u([meas_cell{cell_idx} meas_cell{i}],:),1)>=1;
                    if any(in_gate)
                        mstore = ismember(1:m,[meas_cell{cell_idx} meas_cell{i}]);
                        idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                        if ~any(idx)
                            [~,Lp_merge(i)] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,[meas_cell{cell_idx} meas_cell{i}]),model);
                            Mstore = [Mstore;mstore];
                            Lstore = [Lstore;Lp_merge(i)];
                        else
                            Lp_merge(i) = Lstore(idx);
                        end
                    end
                    Wp(i+N+2) = Lp_merge(i) - Lp(cell_idx) - Lp(i);
                end
            end
        end
        
        %selected cell splited into two parts
        if cell_idx>n_tt_now||size(idx_tracks2n_tt{cell_idx,1},2)==1
            if length(meas_cell{cell_idx}) > 2
                [IDX,~] = kmeanspp(W(:,meas_cell{cell_idx}),2);
                Lp_split1 = -inf(2,1);
                Lp_split2 = -inf(2,1);
                if cell_idx <=n_tt_now %selected cell is a target cell
                    if length(meas_cell{cell_idx}(IDX==2)) > 1
                        [~,Lp_split1(1)] = detectionBern(tracks(cell_idx).Bern,W(:,meas_cell{cell_idx}(IDX==1)),model);
                        in_gate = sum(gating_matrix_u(meas_cell{cell_idx}(IDX==2),:),1)>=1;
                        if any(in_gate)
                            mstore = ismember(1:m,meas_cell{cell_idx}(IDX==2));
                            idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                            if ~any(idx)
                                [~,Lp_split2(1)] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,meas_cell{cell_idx}(IDX==2)),model);
                                Mstore = [Mstore;mstore];
                                Lstore = [Lstore;Lp_split2(1)];
                            else
                                Lp_split2(1) = Lstore(idx);
                            end
                        end
                    end
                    Wp(2*N+3) = Lp_split1(1) + Lp_split2(1) - Lp(cell_idx);
                    
                    if length(meas_cell{cell_idx}(IDX==1)) > 1
                        [~,Lp_split1(2)] = detectionBern(tracks(cell_idx).Bern,W(:,meas_cell{cell_idx}(IDX==2)),model);
                        in_gate = sum(gating_matrix_u(meas_cell{cell_idx}(IDX==1),:),1)>=1;
                        if any(in_gate)
                            mstore = ismember(1:m,meas_cell{cell_idx}(IDX==1));
                            idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                            if ~any(idx)
                                [~,Lp_split2(2)] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,meas_cell{cell_idx}(IDX==1)),model);
                                Mstore = [Mstore;mstore];
                                Lstore = [Lstore;Lp_split2(2)];
                            else
                                Lp_split2(2) = Lstore(idx);
                            end
                        end
                    end
                    Wp(2*N+4) = Lp_split1(2) + Lp_split2(2) - Lp(cell_idx);
                    
                else %selected cell is a new cell
                    if length(meas_cell{cell_idx}(IDX==1)) > 1
                        in_gate = sum(gating_matrix_u(meas_cell{cell_idx}(IDX==1),:),1)>=1;
                        if any(in_gate)
                            mstore = ismember(1:m,meas_cell{cell_idx}(IDX==1));
                            idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                            if ~any(idx)
                                [~,Lp_split1(1)] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,meas_cell{cell_idx}(IDX==1)),model);
                                Mstore = [Mstore;mstore];
                                Lstore = [Lstore;Lp_split1(1)];
                            else
                                Lp_split1(1) = Lstore(idx);
                            end
                        end
                    end
                    if length(meas_cell{cell_idx}(IDX==2)) > 1
                        in_gate = sum(gating_matrix_u(meas_cell{cell_idx}(IDX==2),:),1)>=1;
                        if any(in_gate)
                            mstore = ismember(1:m,meas_cell{cell_idx}(IDX==2));
                            idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                            if ~any(idx)
                                [~,Lp_split2(1)] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,meas_cell{cell_idx}(IDX==2)),model);
                                Mstore = [Mstore;mstore];
                                Lstore = [Lstore;Lp_split2(1)];
                            else
                                Lp_split2(1) = Lstore(idx);
                            end
                        end
                    end
                    Wp(2*N+3) = Lp_split1(1) + Lp_split2(1) - Lp(cell_idx);
                end
            end
        end
    end
    
    [Wp,~] = normalizeLogWeights(Wp);
    I = find(cumsum(exp(Wp))>rand(),1); %sampling
    if I==cell_idx
        Lp = Lp_original;
        meas_cell = meas_cell_original;
        num_repetition = num_repetition + 1;
        if (num_repetition > max_repetition)
            break;
        end
    elseif I <= N %move selected measurement to an existing cell
        meas_cell{I} = [meas_cell{I} mea_selected];
        Lp(I) = Lp_move1(I);
        meas_cell{cell_idx} = mea_after_move;
        if isempty(meas_cell{cell_idx})
            if cell_idx <=n_tt_now
                Lp(cell_idx) = Lmiss(cell_idx);
            else
                meas_cell(cell_idx) = [];
                Lp(cell_idx) = [];
            end
        else
            Lp(cell_idx) = Lp_move2;
        end
%         num_movetoexistingcell=num_movetoexistingcell+1;
    elseif I==N+1 %move selected measurement to a new cell
        meas_cell{N+1,1} = mea_selected;
        Lp(N+1,1) = Lc(mea_selected);
        meas_cell{cell_idx} = mea_after_move;
        Lp(cell_idx) = Lp_move2;
%         num_movetonewcell=num_movetonewcell+1;
    elseif I==N+2 %selected target cell becomes a new cell
        meas_cell{N+1,1} = meas_cell{cell_idx};
        Lp(N+1,1) = Lp_new;
        meas_cell{cell_idx} = [];
        Lp(cell_idx) = Lmiss(cell_idx);
%         num_cellbecomenewcell=num_cellbecomenewcell+1;
    elseif I<=2*N+2 %merge selected cell with an existing cell
        meas_cell{I-N-2} = [meas_cell{cell_idx} meas_cell{I-N-2}];
        Lp(I-N-2) = Lp_merge(I-N-2);
        if cell_idx <= n_tt_now
            meas_cell{cell_idx} = [];
            Lp(cell_idx) = Lmiss(cell_idx);
        else
            meas_cell(cell_idx) = [];
            Lp(cell_idx) = [];
        end
        if I-N-2<=n_tt_now&&cell_idx<=n_tt_now% update prob_merge_mat
            temp=[idx_tracks2n_tt{I-N-2,1} idx_tracks2n_tt{cell_idx,1}];
            for ii=temp
                for jj=temp
                    if ii~=jj
                        prob_merge_mat(ii,jj)=1;
                    end
                end
            end
        end
%         num_mergewithexistingcell=num_mergewithexistingcell+1;
    elseif I==2*N+3 %split selected cell
        meas_cell{N+1,1} = meas_cell{cell_idx}(IDX==2);
        Lp(N+1,1) = Lp_split2(1);
        meas_cell{cell_idx} = meas_cell{cell_idx}(IDX==1);
        Lp(cell_idx) = Lp_split1(1);
%         num_split1=num_split1+1;
    elseif I==2*N+4 %split selected cell
        meas_cell{N+1,1} = meas_cell{cell_idx}(IDX==1);
        Lp(N+1,1) = Lp_split2(2);
        meas_cell{cell_idx} = meas_cell{cell_idx}(IDX==2);
        Lp(cell_idx) = Lp_split1(2);
%         num_split2=num_split2+1;
    end
    
    if I~=cell_idx
        num_repetition = 0;
    end
    
    partitions{t,1} = cellfun(@(x) sort(x),meas_cell,'Un',0);
    lik(t) = sum(Lp);
end

idx_empty = cellfun('isempty',partitions);
lik = lik(~idx_empty);
partitions = partitions(~idx_empty);

[lik,idx] = unique(lik);
partitions = partitions(idx);

%merge in chaotic era to simplify data association
%restore partitions
if any(merge_table)
    for i=1:length(partitions)
        restore_partition=cell(length(partitions{i,1})+n_merged_cell-n_merge_cell,1);
        for j=1:length(restore_partition)
            if j<=n_tt
                restore_partition{j,1}=partitions{i,1}{idx_n_tt2tracks(1,j),1};
            else
                restore_partition{j,1}=partitions{i,1}{j-n_tt+n_tt_now,1};
            end
        end
        partitions{i,1}=restore_partition;
    end
end

% try to split merged tracks
if any(merge_table)
    for i=1:max(merge_table)
        if length(partitions{1,1}{idx_tracks2n_tt{i+n_tt-n_merged_cell,1}(1,1),1})>2
            [IDX,C] = kmeanspp(W(:,partitions{1,1}{idx_tracks2n_tt{i+n_tt-n_merged_cell,1}(1,1),1}),2);
            min_dist=inf;
            for k=idx_tracks2n_tt{i+n_tt-n_merged_cell,1}
                Shat= tracks_truth(k,1).Bern.GGIW(end).V/(tracks_truth(k,1).Bern.GGIW(end).v-2*2-2);
                Shat = (Shat + Shat')/2;
                [~,D]=eig(Shat);
                dist=max(max(D));
                if dist<min_dist
                    min_dist=dist;
                end
            end
            if norm(C(:,1)-C(:,2))>2*min_dist
                lik_split=-inf(length(idx_tracks2n_tt{i+n_tt-n_merged_cell,1}),1);
                idx_maxlik=zeros(length(idx_tracks2n_tt{i+n_tt-n_merged_cell,1}),1);
                W_idx=cell(2,1);
                W_idx{1,1}=partitions{1,1}{idx_tracks2n_tt{i+n_tt-n_merged_cell,1}(1,1),1}(IDX==1);
                W_idx{2,1}=partitions{1,1}{idx_tracks2n_tt{i+n_tt-n_merged_cell,1}(1,1),1}(IDX==2);
                idx=1;
                for k=idx_tracks2n_tt{i+n_tt-n_merged_cell,1}
                    [~,lik_split1]=updateGGIWforBernUseNEOpartition(tracks_truth(k,1).Bern.GGIW,W(:,W_idx{1,1}),model);
                    [~,lik_split2]=updateGGIWforBernUseNEOpartition(tracks_truth(k,1).Bern.GGIW,W(:,W_idx{2,1}),model);
                    [lik_split(idx) idx_maxlik(idx)]=max([lik_split1 lik_split2]);
                    idx=idx+1;
                end
                [~,idx]=max(lik_split);
                
                for j=1:length(partitions)
                    for k=idx_tracks2n_tt{i+n_tt-n_merged_cell,1}
                        if k==idx_tracks2n_tt{i+n_tt-n_merged_cell,1}(1,idx)
                            partitions{j,1}{k,1}=W_idx{idx_maxlik(idx),1};
                        else
                            partitions{j,1}{k,1}=W_idx{3-idx_maxlik(idx),1};
                        end
                    end
                end
                
                prob_merge_mat(idx_tracks2n_tt{i+n_tt-n_merged_cell,1}(1,idx),:)=zeros(1,n_tt);
                prob_merge_mat(:,idx_tracks2n_tt{i+n_tt-n_merged_cell,1}(1,idx))=zeros(n_tt,1);
            end
        end
    end
end

new_merge_table=zeros(1,n_tt);
%repartition
if any(any(prob_merge_mat))
    %compute merge matrix with warshall
    merge_mat=warshall(prob_merge_mat);
    %compute new merge_table
    idx_now=1;
    for i=1:n_tt
        if ~any(merge_mat(i,:))||new_merge_table(1,i)>0
            continue;
        else
            new_merge_table(1,merge_mat(i,:))=idx_now;
            idx_now=idx_now+1;
        end
    end
    %repartition after merging cell
    partition_mat=zeros(length(lik),m);                 %turn partitions to matrix in order to use unique function
    merge_cell_w=[];                        %measurements to partition
    for i=1:length(lik)
        for j=1:length(partitions{i,1})
            %compute partition_mat (convert to matrix)
            partition_mat(i,partitions{i,1}{j,1})=j;
            %compute merge_cell_w (collect all measurements)
            if j<=n_tt&&new_merge_table(1,j)>0
                merge_cell_w=[merge_cell_w partitions{i,1}{j,1}];
            end
        end
    end
    %ensure measurement unique
    merge_cell_w=unique(merge_cell_w);
    %compute quantity of shape
    shape_num=0;
    shape_track_table=[];
    shape_track_th_table=[];
    for i=1:n_tt
        if new_merge_table(1,i)>0
            shape_num=shape_num+length(tracks_truth(i,1).Bern.GGIW(end).shape);
            %table between shape and tracks
            shape_track_table=[shape_track_table repmat(i,[1 length(tracks_truth(i,1).Bern.GGIW(end).shape)])];
            shape_track_th_table=[shape_track_th_table 1:length(tracks_truth(i,1).Bern.GGIW(end).shape)];
        end
    end
    %repartition
    pseudo_likelihood=-inf(length(merge_cell_w),shape_num);
    for k = 1:shape_num
        i=shape_track_table(k);
        shape_th=shape_track_th_table(k);
        ghat = tracks_truth(i,1).Bern.GGIW(end).shape(shape_th).a/tracks_truth(i,1).Bern.GGIW(end).shape(shape_th).b;
        zhat = model.measmodel.h(tracks_truth(i,1).Bern.GGIW(end).m)+model.measmodel.h(tracks_truth(i,1).Bern.GGIW(end).shape(shape_th).m);
        Shat = model.measmodel.H(tracks_truth(i,1).Bern.GGIW(end).shape(shape_th).m)*tracks_truth(i,1).Bern.GGIW(end).shape(shape_th).P*model.measmodel.H(tracks_truth(i,1).Bern.GGIW(end).shape(shape_th).m)' ...
            + model.measmodel.R + tracks_truth(i,1).Bern.GGIW(end).shape(shape_th).V/(tracks_truth(i,1).Bern.GGIW(end).shape(shape_th).v-2*2-2);
        Shat = (Shat + Shat')/2;
        log_pD_ghat_r = log(model.Pd) + log(ghat) + log(tracks_truth(i,1).Bern.r);
        for j = 1:length(merge_cell_w)
            pseudo_likelihood(j,k) = log_pD_ghat_r + log_mvnpdf(W(:,merge_cell_w(j)),zhat,Shat);
        end
    end
    
    [~,idx]=max(pseudo_likelihood,[],2);
    %generate pseudo_partition and partitions
    pseudo_partition=cell(n_tt,1);
    for i=1:shape_num
        pseudo_partition{shape_track_table(i),1}=[pseudo_partition{shape_track_table(i),1} merge_cell_w(idx==i)];
    end
    %restore partitions
    for i=1:n_tt
        if new_merge_table(1,i)>0
            %set 0 to measurements which relate to merging
            partition_mat(partition_mat==i)=0;
        end
    end
    %unique
    [partition_mat,~,~]=unique(partition_mat,'rows');
    %restore partition_mat to partitions
    partitions=cell(size(partition_mat,1),1);
    %get partitions by processing every row of partition_mat
    for i=1:size(partition_mat,1)
        k=max(partition_mat(i,:));
        if k>n_tt
            partitions{i,1}=cell(k,1);
        else
            partitions{i,1}=cell(n_tt,1);
        end
        %
        for j=1:max(k,n_tt)
            if j>n_tt||new_merge_table(1,j)==0
                partitions{i,1}{j,1}=find(partition_mat(i,:)==j);
            else
                partitions{i,1}{j,1}=pseudo_partition{j,1};
            end
        end
    end
    %compute lik
    lik=zeros(size(partition_mat,1),1);
    for i=1:length(lik)
        for j=1:length(partitions{i,1})
            if j<=n_tt
                if ~isempty(partitions{i,1}{j,1})
                    [~,lik_] = detectionBern(tracks_truth(j).Bern,W(:,partitions{i,1}{j,1}),model);
                else
                    [~,lik_] = misdetectionBern(tracks_truth(j).Bern,model);
                end
            else
                if ~isempty(partitions{i,1}{j,1})
                    in_gate = sum(gating_matrix_u(partitions{i,1}{j,1},:),1)>=1;
                    if any(in_gate)
                        [~,lik_] = detectionPPP(PPP.w(in_gate),PPP.GGIW(in_gate),W(:,partitions{i,1}{j,1}),model);
                    else
                        lik_=-100;
                        %lik_=-inf£¬bug if all lik=-inf
                        continue;
                    end
                end
            end
            lik(i)=lik(i)+lik_;
        end
    end
end

end

function R=warshall(A)
n=size(A,1);
R=logical(A);
for i=1:n
    for j=1:n
        if R(j,i)==1
            for k=1:n
                R(j,k)=logical(R(j,k)+R(i,k));
            end
        end
    end
end
end
