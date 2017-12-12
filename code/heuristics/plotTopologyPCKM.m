% PLOTTOPOLOGY Plots a given topology to visualize placements and failure scenarios.
%
%   PLOTTOPOLOGYPCKM(A,B,C) plots the topology A with the coordinates B of the 
%   nodes and controller ids C.
%
%   The parameters A,B,C can be followed by parameter/value pairs to 
%   specify additional properties of the plot.
%
%   The following optional parameters are possible. Default values are
%   shown in brackets '{}'.
%     - 'nodeWeigths',{ones(1,size(topology,1))}: using this parameter, a 
%       node weight vector can be passed to PLOTTOPOLOGY. If not set, 
%       uniform node weights are used.
%     - 'ShowNodeToControllerLatency',{'on'}|'off': indicates whether the
%       node to controller latency should be displayed in the plot.
%     - 'ShowControllerImbalance',{'on'}|'off': indicates whether the
%       controller imbalance should be displayed in the plot.
%     - 'ShowControllerToControllerLatency',{'on'}|'off': indicates whether
%       the controller to controller latency should be displayed in the 
%       plot.
%     - 'ShowIds',{'on'}|'off': indicates whether the node ids should be 
%       displayed in the plot.
%     - 'ShowControllerlessHeatmap',{'on'}|'off': indicates whether the 
%       risk of controller-less nodes should be displayed in the plot.
%     - 'ShowNodeWeights',{'on'}|'off': indicates whether the node weights 
%       should be displayed in the plot.
%     - 'Parent',{figure()}: defines the parent of the plot. The value can
%       be either a axes or figure handle. If no value is provided, a new
%       figure is opened.
%     - 'Export',{'off'}|'pdf'|'png'|'jpg': defines as which format the 
%       plot should be exported. If set to 'off', the plot is not exported.
%     - 'Filename','filename': the name the exported file should have.
%         Default filename is a combination of the topology and the used
%         parameters for plotting.
%         Example:
%           V34_E84_(4-11-20)_()_()_()_(010100).pdf
%           For 34 Vertices, 84 Edges, controllers at node 4, 11 and 20,
%           and the following parameters for plotting:
%             - ShowNodeToControllerLatency: off
%             - ShowControllerImbalance: on
%             - ShowControllerToControllerLatency: off
%             - ShowIds: on
%             - ShowControllerlessHeatmap: off
%             - ShowNodeWeigths: off
%     - 'ReferenceDiameter',{nan}: sets the diameter for the topology. If 
%       set to nan, the longest of all shortest path is taken as diameter.
%     - 'Colors',{hsv(5)}: sets a different colormap for the plot.
%     - 'Markers',{'odsv^'}: changes the shapes of the controllers and assigned
%         nodes.
%     - 'DistanceMatrix',{[]}: sets the distance matrix for the topology.
%
%     Furthermore, the PLOTTOPOLOGY command supports the Figure Property
%     'Position'.
%
%   For example use cases, see also PLOTEXAMPLE.

%   Copyright 2012-2014 David Hock, Stefan Geißler, Fabian Helmschrott,
%                       Steffen Gebert
%                       Chair of Communication Networks, Uni Würzburg   
%
function plotTopologyPCKM(topology,coordinates,controllerPlaces,assignments,currPlacement,varargin)

p = inputParser;    % Parser to parse input arguments

p.addRequired('topology', @(x)isnumeric(x) && diff(size(x))==0); % symmetric numerical topology matrix
p.addRequired('coordinates', @(x)isnumeric(x) && size(x,1)==size(topology,1) && size(x,2)==2);
p.addParamValue('nodeWeights', ones(1,size(topology,1)) ,@(x)isnumeric(x) && size(x,1)==size(topology,1) && size(x,2)==1);

% Plot Parameters
p.addParamValue('ShowNodeToControllerLatency', 'off',  @(x)strcmpi(x,'off')||strcmpi(x,'on'));
p.addParamValue('ShowControllerImbalance', 'off',  @(x)strcmpi(x,'off')||strcmpi(x,'on'));
p.addParamValue('ShowControllerToControllerLatency', 'off',  @(x)strcmpi(x,'off')||strcmpi(x,'on'));
p.addParamValue('ShowIds', 'off',  @(x)strcmpi(x,'off')||strcmpi(x,'on'));
p.addParamValue('ShowControllerlessHeatmap', 'off',  @(x)strcmpi(x,'off')||strcmpi(x,'on'));
p.addParamValue('ShowNodeWeights', 'off',  @(x)strcmpi(x,'off')||strcmpi(x,'on'));
p.addParamValue('Parent', [],  @(x)strcmpi(get(x,'type'),'axes')||strcmpi(get(x,'type'),'figure'));
p.addParamValue('Position', [0 0 800 600],  @(x)isnumeric(x) && length(x)==4);
p.addParamValue('Export', 'off',  @(x)strcmpi(x,'off')||strcmpi(x,'pdf')||strcmpi(x,'png')||strcmpi(x,'jpg'));


% Topology / Placement Parameters
p.addRequired('controllerPlaces',  @isnumeric);
p.addParamValue('ReferenceDiameter', nan,  @isnumeric);
p.addParamValue('DistanceMatrix', [], @isnumeric);
p.addRequired('assignments', @(x)isnumeric(x));
p.addRequired('currPlacement', @isnumeric);


% Visualization Parameters
p.addParamValue('Colors', [hsv(5);hsv(min(1,length(controllerPlaces)-5))],  @isnumeric);
p.addParamValue('Markers', 'odsv^pdshov^x+*',  @ischar);

p.parse(topology, coordinates, controllerPlaces, assignments, currPlacement,varargin{:});
% p.Results

plotids=strcmpi(p.Results.ShowIds,'on');  % boolean to show ids
plotlatencyN2C=strcmpi(p.Results.ShowNodeToControllerLatency,'on'); % boolean to color according to node-to-controller distances
plotbal=strcmpi(p.Results.ShowControllerImbalance,'on'); % boolean to color and shape according to balancing
plotlatencyC2C=strcmpi(p.Results.ShowControllerToControllerLatency,'on'); % boolean to color according to inter-controller distances
plottm=strcmpi(p.Results.ShowNodeWeights,'on'); % boolean if tm should be used for node sizes
plotheatmap=strcmpi(p.Results.ShowControllerlessHeatmap,'on'); % boolean if heatmap should be shown
plotexport=~strcmpi(p.Results.Export,'off'); % boolean if plot is used for display only

% Check Matlab-version
matlabVersion=str2num(regexprep(version('-release'),'[^0-9]','')) ;

if isempty(p.Results.Parent)
    figure();
else
    if strcmpi(get(p.Results.Parent,'type'),'figure')        
        figure(p.Results.Parent);
    else
        axes(p.Results.Parent);
    end
end
set(gca,'YDir','normal');
gplot(topology~=inf,coordinates,'-k');
hold on
% gray links
[rlink,clink]=find(topology~=inf);
beginx=coordinates(rlink,1);
beginy=coordinates(rlink,2);
endx=coordinates(clink,1);
endy=coordinates(clink,2);
for m=1:length(beginx)
    plot([beginx(m) endx(m)],[beginy(m) endy(m)],'Color',[0.6 0.6 0.6],'LineWidth',2, 'LineSmoothing','on');
end

if isempty(p.Results.DistanceMatrix)
    distanceMatrix = allToAllShortestPathMatrix(topology);
else
    distanceMatrix = p.Results.DistanceMatrix;
end

if isnan(p.Results.ReferenceDiameter)
    distanceMatrix(distanceMatrix==inf)=nan;
    mydiameter=nanmax(nanmax(distanceMatrix));
else
    mydiameter=p.Results.ReferenceDiameter;
end

mycolors=p.Results.Colors;

mymarkers=p.Results.Markers;

tm=p.Results.nodeWeights;

if ~plottm
    tm=ones(size(topology,1),1);
end

if plotheatmap
    controllersfailed=[];
    nodesfailed=[];
    linksfailed=[];
end

if length(mycolors)<length(controllerPlaces)
    mycolors=repmat(hsv(10),ceil(length(controllerPlaces)/10),1);
end

% to provide from errors due to unexpected input
if isempty(controllerPlaces)
    plotlatencyN2C=0; % boolean to color according to node-to-controller distances
    plotbal=0; % boolean to color and shape according to balancing
    plotlatencyC2C=0; % boolean to color according to inter-controller distances
    plotheatmap=0; % boolean if heatmap should be shown
end

% Plot nodes
% heatmap
maxfail=2;
if plotheatmap
    heatmap=zeros(1,size(topology,1));
    for i=1:maxfail
        failurepatterns=combnk(1:size(topology,1),i);
        for j=1:size(failurepatterns,1) %
            nodesfailed=failurepatterns(j,:);
            topologyreduced=topology;
            % disable nodes
            topologyreduced(nodesfailed,:)=Inf;
            topologyreduced(:,nodesfailed)=Inf;
            modifiedDistanceMatrix=allToAllShortestPathMatrix(topologyreduced);
            for m=1:length(nodesfailed)
                modifiedDistanceMatrix(nodesfailed(m),nodesfailed(m))=Inf;
            end
            accessnodes=[];
            for m=1:size(modifiedDistanceMatrix,1)
                if length(find(modifiedDistanceMatrix(m,:)~=inf))==1
                    modifiedDistanceMatrix(m,:)=Inf;
                    modifiedDistanceMatrix(:,m)=Inf;
                    accessnodes=[accessnodes m];
                end
            end
            temp=min(modifiedDistanceMatrix(:,controllerPlaces),[],2);
            if matlabVersion < 2013
                uncovered=setdiff(find(temp==Inf),[accessnodes nodesfailed]);
            else
                uncovered=setdiff(find(temp==Inf),[accessnodes nodesfailed],'legacy');
            end
            heatmap(uncovered)=heatmap(uncovered)+1; % add 1 to all controller-less nodes in this scenario
        end
    end
    uncovered=[];
    nodesfailed=[];
end

for i = 1:length(controllerPlaces)
    modifiedDistanceMatrix=distanceMatrix;
    notAssignedToThis=find((assignments(:,currPlacement)==i)==0);
    modifiedDistanceMatrix(:,notAssignedToThis)=Inf;
    modifiedDistanceMatrix(notAssignedToThis,:)=Inf;
    [temp,temp2]=min(modifiedDistanceMatrix(:,controllerPlaces(i)),[],2);
    b=find(temp~=inf);
% colored links
    if plotbal || plotlatencyN2C
        [rlink,clink]=find(modifiedDistanceMatrix~=inf);
        beginx=coordinates(rlink,1);
        beginy=coordinates(rlink,2);
        endx=coordinates(clink,1);
        endy=coordinates(clink,2);    
        for m=1:length(beginx)
            if (temp2(rlink(m))==temp2(clink(m)) && (abs(temp(rlink(m))-(temp(clink(m))-topology(rlink(m),clink(m))))<1e-3 || abs(temp(clink(m))-(temp(rlink(m))-topology(rlink(m),clink(m))))<1e-3))
                if plotbal
                    plot([beginx(m) endx(m)],[beginy(m) endy(m)],'Color',mycolors(i,:),'LineWidth',3, 'LineSmoothing','on');
                else
                    plotFadingLink([beginx(m) endx(m)],[beginy(m) endy(m)],getgreenyellowred(temp(rlink(m))/mydiameter),getgreenyellowred(temp(clink(m))/mydiameter));
                end
            end
        end
    end


    for m=b'
            v='o';
            if plotheatmap
                if max(heatmap)>0
                    c=getgreenyellowred(log(1+heatmap(m))/log(5));
                else
                    c=getgreenyellowred(0);
                end
            elseif plotlatencyN2C
                c=getgreenyellowred(temp(m)/mydiameter);
            elseif plotbal
                c=mycolors(assignments(m,currPlacement),:);
            else
                c=[0.6 0.6 0.6];
            end
            if plotbal
                v=mymarkers(1+mod(assignments(m,currPlacement)-1,length(mymarkers)));
                cborder=darken(mycolors(assignments(m,currPlacement),:));
            else
                cborder=darken(c);
            end
            placex=coordinates(m,1);
            placey=coordinates(m,2);
            plot3(placex,placey,1,v,'MarkerFaceColor',c,'MarkerSize',20*tm(m),'Color',cborder,'LineWidth',3, 'LineSmoothing','on');        
    end
end

% Plot controllers
modifiedDistanceMatrix=distanceMatrix;
temp=modifiedDistanceMatrix(controllerPlaces,controllerPlaces);
temp(temp==Inf)=nan;
temp=nanmax(temp,[],2);
for m=1:length(controllerPlaces(:))
    placex=coordinates(controllerPlaces(m),1);
    placey=coordinates(controllerPlaces(m),2);
    if plotlatencyC2C
        c=getgreenyellowred(temp(m)/mydiameter);
    elseif plotlatencyN2C
        c='g';
    elseif plotbal
        c=mycolors(m,:);
    else
        c=[0.6 0.6 0.6];
    end
    v='o';
    if plotbal
        v=mymarkers(1+mod(temp2(m)-1,length(mymarkers)));
        cborder=darken(mycolors(m,:));
    else
        cborder=darken(c);
    end
    plot3(placex,placey,1,v,'MarkerFaceColor',c,'MarkerSize',20*tm(controllerPlaces(m)),'Color',cborder,'LineWidth',3, 'LineSmoothing','on');
    plot3(placex,placey,1,v,'MarkerFaceColor','none','MarkerSize',9+20*tm(controllerPlaces(m)),'Color',cborder,'LineWidth',3, 'LineSmoothing','on');
    plot3(placex,placey,1,v,'MarkerFaceColor','none','MarkerSize',18+20*tm(controllerPlaces(m)),'Color',cborder,'LineWidth',3, 'LineSmoothing','on');
end


if plotids
    for i=1:size(coordinates,1)
        placex=coordinates(i,1);
        placey=coordinates(i,2);
        text('Position',[placex,placey,1],'String',num2str(i),'FontSize',8,'HorizontalAlignment','center');
    end
end

% Adapt visible range and remove axes
set(gca,'Visible','off');
xlim([min(coordinates(:,1))-0.3 max(coordinates(:,1))+0.3]);
ylim([min(coordinates(:,2))-1 max(coordinates(:,2))+1.5]);


set(gcf,'Position',p.Results.Position, 'PaperPositionMode','auto');

if plotexport==1 % Export is activated, save as pdf file
    originalPaperPosition=get(gcf,'PaperPosition');
    if strcmpi(p.Results.Export,'pdf')
        set(gcf,'PaperSize',originalPaperPosition(3:4),'PaperPosition',[0 0 originalPaperPosition(3:4)]);
    end
    plotpattern(1)=plotlatencyN2C;
    plotpattern(2)=plotbal;
    plotpattern(3)=plotlatencyC2C;
    plotpattern(4)=plotids;
    plotpattern(5)=plotheatmap;
    plotpattern(6)=plottm;
    
    if isempty(p.Results.FileName)
        numNodes=size(topology,1);
        numLinks=length(find(topology~=inf & topology>0));
        filename=[sprintf('V%d_E%d_(%s)_(%s)',numNodes,numLinks,...
            regexprep(num2str(controllerPlaces),'\s+','-'),...
            regexprep(num2str(plotpattern),'\s','')), ['.' p.Results.Export]];
    else
        filename=p.Results.FileName;
        if ~strcmpi(filename(end-length(p.Results.Export)+1:end),p.Results.Export)
            filename=[filename '.' p.Results.Export];
        end
    end
    
    saveas(gcf,filename);
end

% getgreenyellowred Creates a colormap from green to dark red depending on
% the value x in [0,1]. Values x<0 are set to green, values x>1 to dark
% red.
    function c=getgreenyellowred(x)
        myyellow=0.5;
        myred=1;
        mycolormap=[];
        
        % Blue RGB component is always zero
        mycolormap(1:200,3)=zeros(1,200);
        
        % Green to Yellow
        mycolormap(1:ceil(200*myyellow),1)=(1:ceil(200*myyellow))/ceil(200*myyellow);
        mycolormap(1:ceil(200*myyellow),2)=ones(1,ceil(200*myyellow));
        
        % Yellow to Red
        mycolormap((ceil(200*myyellow)+1):ceil(200*myred),1)=ones(1,length((ceil(200*myyellow)+1):ceil(200*myred)));
        mycolormap((ceil(200*myyellow)+1):ceil(200*myred),2)=1-(1:length((ceil(200*myyellow)+1):ceil(200*myred)))/length((ceil(200*myyellow)+1):ceil(200*myred));
        
        % Red to Dark Red
        mycolormap((ceil(200*myred)+1):200,1)=1-(1:length((ceil(200*myred)+1):200))/(1.5*length((ceil(200*myred)+1):200));
        mycolormap((ceil(200*myred)+1):200,2)=0;
        
        c=mycolormap(max(1,min(200,ceil(x*200))),:);
    end

end
