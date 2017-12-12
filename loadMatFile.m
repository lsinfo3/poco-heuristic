    function [topology,coordinates,tm,nodenames,distanceMatrix]=loadMatFile(filename)
        if exist(filename,'file')
            load(filename);
        end
        if ~exist('topology','var') || ~exist('coordinates','var')
            errordlg('The selected file is not a topology file.','File Error');
            topology=[];
            coordinates=[];
            tm=[];
            distanceMatrix=[];
        end
        if ~exist('tm','var')
            tm=[];
        end
        if ~exist('nodenames','var')
            nodenames=[];
        end
        if ~exist('distanceMatrix','var')
           distanceMatrix=allToAllShortestPathMatrix(topology); 
        end
    end