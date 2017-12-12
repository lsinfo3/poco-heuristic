function [centers, assignment] = capacitatedKMedoids(distanceMatrix, k, rho, seed)

rng(seed, 'twister');

% fprintf('%d ', seed)

% If a new set of centers doesn't improve costs by at least this amount, the algorithm stops.
TAU = eps;

n = size(distanceMatrix, 1);

% Get initial controller locations via k-Medoids
[vartilde, centers] = kmedioids(distanceMatrix, k);
% Order required for later
centers = sort(centers);

% Generate bipartite graph (bpg) by using the distances between the centers and all nodes
ctrlIDs = reshape(repmat(1:k, ceil(n/k) + rho, 1), (ceil(n/k) + rho) * k, 1);
bpg = distanceMatrix(centers, :);
bpg = bpg(ctrlIDs, :);

% assignment(i) = val => node i connects to controller ctrlIDs(val)
% WARNING: For a nonsquare bpg, assignment is the index vector of the
% larger dimension assigned to the smaller dimension.
% => add dummy row to enforce dimension relationship
[assignment, costs] = lapjv([bpg; Inf(1, n)]);
% centerIDs(i) = nodeID of the center node i is assigned to
centerIDs = centers(ctrlIDs(assignment));

% Cell array with nodeIDs of nodes in each cluster, i.e., clusters{centers(i)} = vector of nodes in cluster i
clusters = accumarray(centerIDs', 1:n, [], @(x) {x});
% clusters{i} = vector of nodes in cluster i
clusters = clusters(~cellfun('isempty', clusters));
% For each cluster, find the cluster member which minimizes the sum of distances within the cluster
getCenters = @(cluster) min(sum(distanceMatrix(cluster, cluster), 1));
% newCosts are the costs per cluster if the new centers are used
% clusters{i}(newCenterIdxs(i)) = nodeID of new center for cluster i
[newCosts, newCenterIdxs] = cellfun(getCenters, clusters);
newCosts = sum(newCosts);
% Get abovementioned nodeIDs
newCenters = cellfun(@(cl, id) cl(id), clusters, num2cell(newCenterIdxs));

%fprintf('Centers: [ %s ] -> [ %s ]\n', sprintf('%d\t', sort(centers)), sprintf('%d\t', sort(newCenters)))
%fprintf('costs -> newCosts: %f -> %f\n', costs, newCosts)

ctr = 0;

while (costs - newCosts) > TAU
    centers = sort(newCenters');
    % generate new bpg, reassign n2c, recalculate costs / centers
    bpg = distanceMatrix(centers, :);
    bpg = bpg(ctrlIDs, :);
    % backup assignment in case the new assignment produces higher costs
    tmpAssignment = assignment;
    [assignment, costs] = lapjv([bpg; Inf(1, n)]);
    % No need to continue if the new assignment produces higher costs
    if costs > newCosts
        assignment = tmpAssignment;
        break;
    end
    centerIDs = centers(ctrlIDs(assignment));
    clusters = accumarray(centerIDs', 1:n, [], @(x) {x});
    clusters = clusters(~cellfun('isempty', clusters));
    getCenters = @(cluster) min(sum(distanceMatrix(cluster, cluster), 1));
    [newCosts, newCenterIdxs] = cellfun(getCenters, clusters);
    newCosts = sum(newCosts);
    % Get abovementioned nodeIDs
    newCenters = cellfun(@(cl, id) cl(id), clusters, num2cell(newCenterIdxs));
    ctr = ctr + 1;
%    fprintf('Centers: [ %s ] -> [ %s ]\n', sprintf('%d\t', sort(centers)), sprintf('%d\t', sort(newCenters)))
%    fprintf('Iteration %03d - newCosts: %f ##\n', ctr, newCosts);
end

centers = newCenters;

assignment = ctrlIDs(assignment);

end
