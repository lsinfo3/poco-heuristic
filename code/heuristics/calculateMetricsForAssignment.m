% CALCULATEMETRICS Calculates all placement metrics for given distance 
% matrix, placements, and node weights.
%
%   [U,V,W,X,Y,Z] = CALCULATEMETRICS(A,B) calculates and returns 
%   placement metrics for given distance matrix A and placements B, where B
%   is a matrix containing in each column one placement (consisting of the
%   node IDs given in the rows of this column). 
%   U is a vector containing the average latency between the nodes and the 
%   controllers, V is a vector containing the maximum latency between the 
%   nodes and the controllers, W is a vector containing the maximum number 
%   of controller-less nodes, X is a vector containing the imbalance of the
%   controllers, y is a vector containing the average latency between the 
%   controllers and Z is a vector containing the maximum latency between 
%   the controllers.
%
%   [U,V,W,X,Y,Z] = CALCULATEMETRICS(A,B,C) calculates and returns 
%   placement metrics for given distance matrix A and placements B as well 
%   as node weights C, where C is a row vector containing for each node 
%   a node weight value. 
%
%   For example use cases, see also PLOTEXAMPLE.

%   Copyright 2012-2013 David Hock, Stefan Geißler, Fabian Helmschrott,
%                       Steffen Gebert
%                       Chair of Communication Networks, Uni Würzburg   

function  [avgLatencyN2C,maxLatencyN2C,controllerlessNodes,controllerImbalance,avgLatencyC2C,maxLatencyC2C]=calculateMetricsForAssignment(distanceMatrix,placementArray,assignment,nodeWeights)
% Check Matlab-version
matlabVersion=str2num(regexprep(version('-release'),'[^0-9]',''));

if matlabVersion < 2013
    [uniquePlacementArray,tildevar,backIdx]=unique(placementArray,'rows');
else
    [uniquePlacementArray,tildevar,backIdx]=unique(placementArray,'rows','legacy');
end
%fprintf('%d;%d\n',length(placementArray),length(uniquePlacementArray));
clear tildevar;

% Create temporary matrix T
% tempidx == assignment, temp == n2c distance given the assignment
% distMat(1:n, placementArray(tempidx));
% assignment/tempidx format = column vector of length n with elements \in 1:k
tempidx = assignment;
n = size(distanceMatrix, 1);
nPlacements = size(placementArray, 1);
% equivalent of tempidx but with nodeIDs rather than controllerIDs
tidx2 = placementArray(sub2ind(size(placementArray), repmat(1:nPlacements, n, 1), tempidx));
tidx2 = reshape(tidx2, n, nPlacements);
temp = distanceMatrix(sub2ind(size(distanceMatrix), repmat((1:n)', 1, nPlacements), tidx2));
%[temp,tempidx]=min(reshape(distanceMatrix(:,uniquePlacementArray),size(distanceMatrix,1),size(uniquePlacementArray,1),size(uniquePlacementArray,2)),[],3);

% Calculate load balancing

% Indices of all controllers with at least one node assigned
if matlabVersion < 2013
    usedidx=unique(tempidx);
else
    usedidx=unique(tempidx, 'legacy');
end

% Simple case - only one controller - imbalance is 0
if length(usedidx)<=1        
        controllerImbalance=zeros(1,size(tempidx,2));
else   
    
    if ~exist('nodeWeights','var')
        tempidx4=hist(tempidx,usedidx);
    else
        % If node weights are defined, they are included in calculation of
        % imbalance
        tempidx4=accumarray([reshape(repmat(1:size(tempidx,2),size(tempidx,1),1),numel(tempidx),1) reshape(tempidx,numel(tempidx),1)],repmat(nodeWeights,1,size(tempidx,2)),[], @sum)';
    end
    
    % Transpose necessary to repair one-column vector
    if size(tempidx4,1)==1 && size(tempidx4,2)>1
        tempidx4=tempidx4';
    end
    
    % Calculate imbalance as max minus min 
    controllerImbalance=max(tempidx4(end-length(usedidx)+1:end,:))-min(tempidx4(end-length(usedidx)+1:end,:));
end

% Number of controller-less nodes is calculated with or without node
% weights
if exist('nodeWeights','var')
    controllerlessNodes=sum((temp==Inf).*repmat(nodeWeights',1,size(temp,2)));
else
    controllerlessNodes=sum(temp==Inf);
end

% Maximum and average node-to-controller latency are calculated
temp(temp==Inf)=nan;
avgLatencyN2C=nanmean(temp);
maxLatencyN2C=nanmax(temp);

% Maximum and average inter-controller latency are calculated - matrix needs to be of square size
if diff(size(distanceMatrix))==0
    %Controller to Controller - max
    temp=reshape(distanceMatrix(:,uniquePlacementArray),size(uniquePlacementArray,1)*size(distanceMatrix,1),size(uniquePlacementArray,2));
    maxtemp=reshape(nanmax(temp,[],2),size(distanceMatrix,1),size(uniquePlacementArray,1));
    maxtemp(maxtemp==Inf)=nan;
    maxLatencyC2C=nanmax(maxtemp(sub2ind(size(maxtemp),uniquePlacementArray',repmat(1:size(uniquePlacementArray,1),size(uniquePlacementArray,2),1))),[],1);

    % Controller to Controller - avg
    % avoid considering c2c distances where c1 == c2 FIXME
    temp(temp==0) = nan;
    meantemp=reshape(nanmean(temp,2),size(distanceMatrix,1),size(uniquePlacementArray,1));
    meantemp(meantemp==Inf)=nan;
    avgLatencyC2C=nanmean(meantemp(sub2ind(size(meantemp),uniquePlacementArray',repmat(1:size(uniquePlacementArray,1),size(uniquePlacementArray,2),1))),1);
else
    maxLatencyC2C=nan(size(maxLatencyN2C));
    avgLatencyC2C=nan(size(maxLatencyN2C));
end

diameter = max(distanceMatrix(:));

avgLatencyN2C=avgLatencyN2C(backIdx)/diameter;
maxLatencyN2C=maxLatencyN2C(backIdx)/diameter;
controllerlessNodes=controllerlessNodes(backIdx);
controllerImbalance=controllerImbalance(backIdx)/n;
avgLatencyC2C=avgLatencyC2C(backIdx)/diameter;
maxLatencyC2C=maxLatencyC2C(backIdx)/diameter;

return
