% GENERATES ROC AND PR CURVES USING THRESHOLD AVERAGING.
% Inputs: 
% Generator = name of the script in the Data_Generators folder which will
%             produce simulated data.
% its_max = number of ROC and PR curves to average over (>=2).
% t = sampling time points 1*(n+1) row vector.
% snr_cell = signal-to-noise ratio for cellular dynamics.
% snr_meas = signal-to-noise ration for measurement.
% mod_S = number of cells which are averaged to produce data.
% Outputs:
% data = a typical simulated datasets.
% genes = names of the model variables.
% TPR = true positive rates.
% FPR = false positive rates.
% P = precisions.
% R = recalls.
% schemes = names of inference schemes which are assessed.
% table = details of the individual inference schemes.
% © Chris Oates 2011.

function [data,genes,TPR,FPR,AUR,P,R,schemes,table] = curves(Generator,its_max,t,snr_cell,snr_meas,mod_S)

% search Methods folder for inference schemes to assess
cd('../Schemes')
files = dir('*.m');
M = length(files);
cd('../Data Generators')
gen = str2func(Generator);

% initialise storage
[~,A,~] = gen(t,snr_cell,snr_meas,mod_S);
D = length(A);
TPR = zeros(M,1+D^2); 
FPR = zeros(M,1+D^2); 
AUR = zeros(its_max,M);
P = zeros(M,1+D^2); 
R = zeros(M,1+D^2); 
schemes = cell(M,1);
table = zeros(M,4);
cd('../Schemes')

for its = 1:its_max

    % simulate the experiment
    cd('../Data Generators')
    [y,A,genes] = gen(t,snr_cell,snr_meas,mod_S);
    cd('../Schemes')
    if its == 1
        data = y;
    end
    
    % assess the inference schemes
    for m = 1:M

        % apply scheme m
        method = files(m).name;
        method(end-1:end) = [];
        schemes(m) = cellstr(method);
        f_method = str2func(method);
        [A_est,C] = f_method(y,t);
        table(m,:) = C;
        
        % assess the performance of scheme m
        cd('../Assessment')
        [tpr,fpr,aur,p,r] = curve(A_est,A);
        
        % record performance
        TPR(m,:) = TPR(m,:) + tpr/its_max;
        FPR(m,:) = FPR(m,:) + fpr/its_max;
        AUR(its,m) = aur;
        P(m,:) = P(m,:) + p/its_max;
        R(m,:) = R(m,:) + r/its_max;
        cd('../Schemes')

    end 
end
cd('../Assessment')
end

% GENERATE SINGLE ROC AND PR CURVES
% Inputs:
% A = predicted adjacency matrix.
% A_true = true adjacency matrix.
% Outputs:
% tpr = true positive rates.
% fpr = false positive rates.
% aur = area under ROC curve.
% p = precisions.
% r = responses.
% © Chris Oates 2011.

function [tpr,fpr,aur,p,r] = curve(A,A_true)

true_network = ceil(abs(A_true));
D = length(A_true);

% scale A for thresholding
A = abs(A);
A = reshape(A,1,D^2);
[~,order] = sort(A);
for i = 1:D^2
    A(order(i)) = (2*i-1)/(2*D^2);
end
A = reshape(A,D,D);

% initialise storage
tp = zeros(1,D^2+1); fp = zeros(1,D^2+1); tn = zeros(1,D^2+1); fn = zeros(1,D^2+1);

% perform the computation
for k = 0:D^2
    
    % compute the network estimator
    threshold = k / D^2;
    inferred_network = im2bw(A,threshold);
    
    % assess the predictions
    true_positives = true_network & inferred_network;
    false_positives = inferred_network - true_positives;
    true_negatives = (~true_network) & (~inferred_network);
    false_negatives = (~inferred_network) - true_negatives;
    
    % do not assess self-loops
    tp(k+1) = sum_offdiagonal(true_positives);
    fp(k+1) = sum_offdiagonal(false_positives);
    tn(k+1) = sum_offdiagonal(true_negatives);
    fn(k+1) = sum_offdiagonal(false_negatives);
    
end

% compute ROC curve
fpr = fp ./ (tn + fp); tpr = tp ./ (tp + fn);
aur = -trapz(fpr,tpr);

% compute PR curve
p = tp ./ (tp + fp); r = tp ./ (tp + fn);
end

function z = sum_offdiagonal(M)
    z = sum(sum(M)) - sum(diag(M));
end

