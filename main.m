clear;clc;close all;
dbstop if error
%Choose a scenario: Scenario 1: 27 targets born at four different
%locations; Scenario 2: targets move in proximity.
scenario = 2;

%Parameter setting
if scenario == 1
    modelparas1;
    modelparas1_original;
elseif scenario==2
    modelparas2;
    modelparas2_original;
else
    modelparas3;
    modelparas3_original;
end

%Number of Monte Carlo Simulations
numMC = length(Scenario.Z);
numMC = 5;%for test
%Parameters used in GOSPA metric
c = 20;
p = 1;

%Number of time steps
K = model.K;

GOSPA = zeros(K,5,numMC);
trajectoryEstimates = cell(numMC,1);
simulation_time = zeros(K,1);

GOSPA_original = zeros(K,5,numMC);
trajectoryEstimates_original = cell(numMC,1);
simulation_time_original = zeros(K,1);

estimates = cell(numMC,1);
estimates_original = cell(numMC,1);

for t = 1:numMC
    Z = Scenario.Z{t};
    % Initialisation
    PPP.w = log(model.birth.w);
    PPP.GGIW = model.birth.GGIW;
    MBM.w = [];     % Global hypotheses weights
    MBM.track = {}; % Locl hypotheses trees
    MBM.table = []; % Global hypotheses look-up table
    PPP_original.w = log(model_original.birth.w);
    PPP_original.GGIW = model_original.birth.GGIW;
    MBM_original.w = [];     % Global hypotheses weights
    MBM_original.track = {}; % Locl hypotheses trees
    MBM_original.table = []; % Global hypotheses look-up table
    
    estimates{t,1} = cell(K,1);
    estimates_original{t,1} = cell(K,1);
    
    for k = 1:K
        [t k]
        tic
        %Update step
        [PPP,MBM] = updatePMBM(PPP,MBM,Z{k},k,model);
        %Extract estimates (both estimate of the current time and the
        %estimate of the full trajectory) 
        [estimates{t,1}{k},trajectoryEstimates{t}{k}] = estimator(MBM,model);
        %Prediction Step
        if k < K
            [PPP,MBM] = predictPMBM(PPP,MBM,model);
        end
        simulation_time(k) = simulation_time(k)+toc;
    end

    %Evaluate filtering performance using GOSPA
    for k=1:K
        GOSPA(k,:,t) = GOSPAmetric(estimates{t,1}{k},groundTruth{k},c,p);
    end
    
    %original PMBM
    for k = 1:K
        [t k]
        tic
        %Update step
        [PPP_original,MBM_original] = updatePMBM_original(PPP_original,MBM_original,Z{k},k,model_original);
        %Extract estimates (both estimate of the current time and the
        %estimate of the full trajectory) 
        [estimates_original{t,1}{k},trajectoryEstimates_original{t}{k}] = estimator_original(MBM_original,model_original);
        %Prediction Step
        if k < K
            [PPP_original,MBM_original] = predictPMBM_original(PPP_original,MBM_original,model_original);
        end
        simulation_time_original(k) = simulation_time_original(k)+toc;
    end
    
    %Evaluate filtering performance using GOSPA
    for k=1:K
        GOSPA_original(k,:,t) = GOSPAmetric_original(estimates_original{t,1}{k},groundTruth_original{k},c,p);
    end
end
simulation_time02=simulation_time/K;
simulation_time_original02=simulation_time_original/K;
GOSPA02 = sum(GOSPA,3) / numMC;
GOSPA02_original = sum(GOSPA_original,3) / numMC;
displayresult;