load('ParallelTarget.mat')
load('ParallelTargetGroundTruth_original.mat')
% load('ParallelTargetVXdistance10T30Poisson8s.mat')
% load('ParallelTargetGroundTruth_originalVXdistance10T30s.mat')

K = length(Scenario.Z{1});
model_original.K = K;

% Effective window length for the gamma prediction
w_e_gamma = 20;
% Effective window length for the extent prediction
w_e_extent = 15;

model_original.tao = 1/(log(w_e_extent)-log(w_e_extent-1));
model_original.eta = 1/(1-1/w_e_gamma);
model_original.Ts = 1;   %sampling interval
sigma_v = 0.1;  %standard deviation of motion noise
sigma_r = 0.1;  %standard deviation of measurement noise
model_original.motionmodel = motionmodel.cvmodel(model_original.Ts,sigma_v);
model_original.measmodel = measmodel.cvmeasmodel(sigma_r);

%re-construct the data structure

% % generate tracks (ground truth)
% X = cell(K,1);
% E = cell(K,1);
% N = zeros(K,1);
% groundTruth = cell(K,1);
% % generate tracks (ground truth)
% for targetnum = 1:length(targetTracks)
%     for k = targetTracks(targetnum).birthTime:targetTracks(targetnum).deathTime
%         targetstate = targetTracks(targetnum).x(1:model.motionmodel.d,k-targetTracks(targetnum).birthTime+1);
%         targetextent = targetTracks(targetnum).X(:,:,k-targetTracks(targetnum).birthTime+1);
%         X{k} = [X{k} targetstate];
%         E{k} = cat(3,E{k},targetextent);
%         N(k) = N(k) + 1;
%     end
% end
% 
% for k = 1:K
%     groundTruth{k}.x = X{k};
%     groundTruth{k}.X = E{k};
% end

%target existence probability
model_original.Ps = 0.99;
%target detection probability
model_original.Pd = Scenario.detection_prob;

%range of the surveillance area
range_c = [-1 1;-1 1]*200;
%Poisson false alarm (clutter) rate
lambda_c = Scenario.false_alarm_rate;
%Poisson clutter intensity
model_original.lambda_fa = lambda_c/prod(range_c(:,2)-range_c(:,1));

% target initial state
nbirths = 2;
xstart = zeros(model_original.motionmodel.d,nbirths);

xstart(:,1) = [-110 75 0 0];
xstart(:,2) = [-110 -75 0 0];

% xstart(:,1) = [-118 66 0 0];
% xstart(:,2) = [-118 -66 0 0];
%Birth model
d = 2;
model_original.birth.w = 0.5*ones(nbirths,1);
model_original.birth.GGIW = repmat(struct('a',100,'b',10,'m',[],'P',diag([10,10,10,10]),'v',22.5,'V',72.5*eye(2)),[nbirths,1]);
for i = 1:nbirths
    model_original.birth.GGIW(i).m = xstart(:,i);
end

% Gating parameters
Pg = 0.999;
model_original.gamma= chi2inv(Pg,model_original.measmodel.d);
model_original.Qd = 1 - model_original.Pd*Pg;

% Thresholds
model_original.threshold_r = 1e-2;   %existence probability of Bernoulli component
model_original.threshold_u = 1e-2;   %weight of mixture component in PPP
model_original.threshold_w = 0.1;   %1e-2 weight of global hypothesis (multi-Bernoulli)
model_original.threshold_s = 1e-4;   %weight of the trajectory is still alive
model_original.recycle = 1e-1;       %recycling threshold
model_original.merge = 4;            %merge threshold used to merge similar GGIWs
model_original.M = 100;              %cap of number of MBM components in PMBM
model_original.num_iterations = 3;   %controls the number of iterations used in SO
model_original.max_repetition = 1;   %controls the number of iterations used in SO

%extract target state from Bernoulli components with existence probability
%larger than this threshold 
model_original.exist_r = 0.5;        

