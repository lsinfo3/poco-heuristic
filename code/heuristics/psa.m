function [M, statsM] = psa(T0, rho, m, s, distanceMatrix, k, seed)
%PSA Pareto Simulated Annealing
%   T0:     starting temperature (~50)
%   rho:    cooling rate (~.9)
%   m:      number of iterations per temperature level
%   s:      number of elements in the set of generating solutions S
%   k:      number of controllers per placement
%   seed:   seed for the RNG

DEBUG = false;

t0 = tic;

rng(seed);

% logger for debugging
if DEBUG
    f = fopen('C:/Users/stas/git/poco/code/psa/probabilities.csv', 'w');
end

nStats = 6;
alpha = 1.05;
diameter = max(distanceMatrix(:));

% n: number of nodes
n = size(distanceMatrix, 1);
T = T0;

% Random placements
S = generateS(n, s, k);
[avgLatencyN2C, maxLatencyN2C, controllerlessNodes, controllerImbalance, avgLatencyC2C, maxLatencyC2C] = ...
    calculateMetrics(distanceMatrix, S);
D = zeros(size(S, 1));

statsS = [avgLatencyN2C' / diameter, maxLatencyN2C' / diameter, controllerlessNodes', controllerImbalance' / n, avgLatencyC2C' / diameter, maxLatencyC2C' / diameter];
% M is the set of Pareto optima from S
M = S(paretoGroup(statsS), :);

Lambda = zeros(s, nStats);

itercount = 0;

if DEBUG
    fprintf('[%s] ### Starting PSA ###\n', datestr(now, ' YYYY-mm-dd HH:MM:SS '))
end

while T > 1
    
    % Permute rowwise (allows more efficient implementation of the next step).
    S = shuffleRows(S);
    
    % For each placement x in S, draw a random neighbor y by replacing the
    % location of one controller location with a random location.
    % draw up to k/2 new neighbors, this increases diversity while in higher temperature regions
    Y = drawNeighbors5(S, n, ceil(T*k / (2* T0)));
    % Calculate the stats for the new placements.
    [avgLatencyN2C, maxLatencyN2C, controllerlessNodes, controllerImbalance, avgLatencyC2C, maxLatencyC2C] = ...
        calculateMetrics(distanceMatrix, [Y; M]);
    statsYM = [avgLatencyN2C' / diameter, maxLatencyN2C' / diameter, controllerlessNodes', controllerImbalance' / n, avgLatencyC2C' / diameter, maxLatencyC2C' / diameter];
    statsY = statsYM(1:size(Y, 1), :);
    % Update M.
    SYM = [S; Y; M];
    optima = paretoGroup([statsS; statsYM]);
    M = SYM(optima, :);
    
    if (itercount == 0)
        
        % (Initially random) objective weights for each solution (normalized rowwise).
        Lambda = rand(s, nStats);
        Lambda = bsxfun(@times, Lambda, 1./(sum(Lambda, 2)));
        
    else
        
        % For each solution s in S, find a solution closest to it.
        
        % D contains the Lambda weighted, source dependent distance between elements in S.
        [p,q] = meshgrid(1:size(statsS,1), 1:size(statsS,1));
        pairs = [p(:) q(:)];

        dists = abs(statsS(pairs(:,1),:)-statsS(pairs(:,2),:));
        D(sub2ind(size(D), pairs(:, 1), pairs(:, 2))) = sum(Lambda(pairs(:, 1), :) .* dists((pairs(:, 1) - 1) * size(S, 1) + pairs(:, 2), :), 2);

        D = D + diag(Inf(size(diag(D))));
        
        % sidx == indices of x' in S (x' are the elements in S closest to x)
        [~, sidx] = min(D, [], 2);
        
        % combinations of x and x' as indices of S
        combs = [(1:size(S,1))', sidx];
        
        % update Lambda
        % \lambda^{x_i}_j * \alpha, if (f_j(x_i) - f_j(x'_i)) \leq 0
        % \lambda^{x_i}_j / \alpha, else
        A = ones(size(Lambda)) * alpha;
        A(logical(dists((combs(:, 1) - 1) * size(S, 1) + combs(:, 2), :) > 0)) = 1 / alpha;
        Lambda = A .* Lambda;
        
    end
    
    % update S
    % if x dominates y, reject y (keep x) (I) FIXME actually, there's still a probability to accept y
    % if y dominates x, accept y (replace x with y) (II)
    % else, accept y with probability P (exponentiated weighted sum, temperature dependent) (III)
    YdomX = find(sum(statsY <= statsS) == nStats);
    nondom = setdiff(1:s, [YdomX]);
    % (II)
    S(YdomX, :) = Y(YdomX, :);
    
    % (III)
    P = min(1, exp(sum(Lambda(nondom, :) .* (statsS(nondom, :) - statsY(nondom, :)), 2)));
    % alternative P
%     P = min(1, exp(max(Lambda(nondom, :) .* (statsS(nondom, :) - statsY(nondom, :)), [], 2)));

    acc = logical(rand(length(nondom), 1) < P);
    S(acc, :) = Y(acc, :);
    
    [avgLatencyN2C, maxLatencyN2C, controllerlessNodes, controllerImbalance, avgLatencyC2C, maxLatencyC2C] = ...
    calculateMetrics(distanceMatrix, S);
    statsS = [avgLatencyN2C' / diameter, maxLatencyN2C' / diameter, controllerlessNodes', controllerImbalance' / n, avgLatencyC2C' / diameter, maxLatencyC2C' / diameter];
    
    itercount = itercount + 1;
    
    % update T
    if mod(itercount, m) == 0
%         fprintf('[%s] Changing temperature from %.2f to %.2f (step %03d / %03d)\n', datestr(now, ' YYYY-mm-dd HH:MM:SS '), ...
%                     T, T * rho, ceil(-log2(T0)/log2(rho)) - ceil(-log2(T)/log2(rho)), ceil(-log2(T0)/log2(rho)))
        T = T * rho;
    end
    
end

    [avgLatencyN2C, maxLatencyN2C, controllerlessNodes, controllerImbalance, avgLatencyC2C, maxLatencyC2C] = ...
        calculateMetrics(distanceMatrix, M);
%    statsM = [avgLatencyN2C' / diameter, maxLatencyN2C' / diameter, controllerlessNodes', controllerImbalance' / n, avgLatencyC2C' / diameter, maxLatencyC2C' / diameter];
    statsM.maxNumberOfControllerlessNodes=controllerlessNodes;
    statsM.avgLatencyN2C=avgLatencyN2C./diameter;
    statsM.maxLatencyN2C=maxLatencyN2C./diameter;
    statsM.avgLatencyC2C=avgLatencyC2C./diameter;
    statsM.maxLatencyC2C=maxLatencyC2C./diameter;
    statsM.controllerImbalance=controllerImbalance./n;    
    if DEBUG
        fclose(f);
        fprintf('[%s] ### PSA FINISHED IN %.2f SECONDS ###\n', datestr(now, ' YYYY-mm-dd HH:MM:SS '), toc(t0))
    end

end


function S = generateS(n, s, k)
%GENERATES Generate the set of generating solutions S (i.e., controller placements).
%   n:  network size
%   s:  output size, |S|
%   k:  number of controllers

S = zeros(s + 1, k);
ii = 1;
nuniq = size(unique(S, 'rows'), 1);
while nuniq ~= s + 1
    S(ii, :) = randsample(n, k);
    nuniq2 = size(unique(S, 'rows'), 1);
    if nuniq2 > nuniq
        ii = ii + 1;
        nuniq = nuniq + 1;
    end
end

S = S(1:s, :);

end

function B = shuffleRows(A)

[nRows, nCols] = size(A);
[vartilde, idx] = sort(rand(nRows, nCols), 2);
% convert column indices into linear indices
idx = (idx - 1) * nRows + ndgrid(1:nRows, 1:nCols);
% rearrange
B = A;
B(:) = B(idx);

end

function X = drawNeighbors(S, n)

X = S;
k = size(S, 2);

for i = 1:(size(S, 1))
    
    r = S(i, :);
    id = randi(k);
    sub = randi(n);
    while ismember(sub, r)
        sub = randi(n);
    end
    r(id) = sub;
    X(i, :) = r;
    
end

end

function X = drawNeighbors5(S, n, nNeighbors)

k = size(S, 2);

X = S;

for ii = 1:size(X,1)
    % draw k+1 candidates (this guarantees that at least one of them is new as randsample's default is w/o replacement)
    c = randsample(n, k + nNeighbors);
    % choose nNeighbors that are not part of the current placement
    % => the new placement has still k distinct controllers and does not equal the original one
    X(ii, 1:nNeighbors) = c(find(~ismember(c, X(ii, :)), nNeighbors, 'first'));
end

end

