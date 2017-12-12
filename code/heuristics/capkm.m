function [centers, assignments] = capkm(distanceMatrix, k, rhoRange, nReps, seed)
    
    n = size(distanceMatrix, 1);
    centers = zeros(k, nReps * length(rhoRange));
    assignments = zeros(n, nReps * length(rhoRange));
    
    for r = 1:length(rhoRange)
        for i = 1:nReps
            
            repID = (r - 1) * nReps + i;
            [center, assignment] = capacitatedKMedoids(distanceMatrix, k, rhoRange(r), seed + repID);
            centers(:, repID) = center;
            assignments(:, repID) = assignment;
            
        end
    end
    
end
