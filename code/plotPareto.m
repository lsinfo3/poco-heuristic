% PLOTPARETO Displays the solution space for two given metrics including a set of Pareto-optimal values.
%
%   X = PLOTPARETO(A,B) displays the solution space including a set of
%   Pareto-optimal values for two given metrics A and B, where A and B are
%   two vectors of equal length. The returned array X contains a vector
%   with the ids of all Pareto-optimal values.
%
%   The parameters A and B can be followed by parameter/value pairs to
%   specify additional properties of the plot.
%
%   The following optional parameters are possible. Default values are
%   shown in brackets '{}'.
%     - 'showMeanValues',{'on'}|'off': indicates whether mean values of
%       metrics A and B should be shown
%     - 'ShowSelectedidx',{0}: indicates a particular placement id that
%       should be highlighted in the graph. If set to 0, no particular
%       placement is highlighted.
%     - 'Parent',{figure()}: defines the parent of the plot. The value can
%       be either a axes or figure handle. If no value is provided, a new
%       figure is opened.
%     - 'Export',{'off'}|'pdf'|'png'|'jpg': defines as which format the
%       plot should be exported. If set to 'off', the plot is not exported.
%     - 'Filename','filename': the name the exported file should have.
%         Default filename is 'paretoPlot.EXTENSION' where EXTENSION
%         corresponds to the extension provided in the 'Export' parameter.
%
%     Furthermore, the PLOTPARETO command supports following Figure and
%     Axes Properties that can be directly passed as parameter/value pairs:
%     'XLabel','YLabel','XLim','Ylim','Position'.
%
%   For example use cases, see also PLOTEXAMPLE.

%   Copyright 2012-2013 David Hock, Stefan Geißler, Fabian Helmschrott,
%                       Steffen Gebert
%                       Chair of Communication Networks, Uni Würzburg

function [paretoidx]=plotPareto(xin,yin,varargin)

p = inputParser;    % Parser to parse input arguments

p.addRequired('xin', @isnumeric);
p.addRequired('yin', @(x)isnumeric(x) && length(xin)==length(yin));

% Plot Parameters
p.addParamValue('ShowMeanValues', 'on',  @(x)strcmpi(x,'off')||strcmpi(x,'on'));
p.addParamValue('ShowSelectedidx', 0,  @(x)isnumeric(x) && x>0 && x<=length(xin));
p.addParamValue('Parent', [],  @(x)strcmpi(get(x,'type'),'axes')||strcmpi(get(x,'type'),'figure'));
p.addParamValue('Export', 'off',  @(x)strcmpi(x,'off')||strcmpi(x,'pdf')||strcmpi(x,'png')||strcmpi(x,'jpg'));
p.addParamValue('FileName', '',  @ischar);
p.addParamValue('XLabel', '',  @ischar);
p.addParamValue('YLabel', '',  @ischar);
p.addParamValue('XLim', [],  @(x)isnumeric(x));
p.addParamValue('YLim', [],  @(x)isnumeric(x));
p.addParamValue('Position', [0 0 800 600],  @(x)isnumeric(x) && length(x)==4);
p.addParamValue('dim3in', [],  @(x)isnumeric(x) && length(x)==length(xin));
p.addParamValue('dim4in', [],  @(x)isnumeric(x) && length(x)==length(xin));
p.addParamValue('Dim3Label', '',  @ischar);
p.addParamValue('Dim4Label', '',  @ischar);
p.addParamValue('Dim3Lim', [],  @(x)isnumeric(x));
p.addParamValue('Dim4Lim', [],  @(x)isnumeric(x));

p.parse(xin, yin, varargin{:});
mycolors=hsv(5);
% Check Matlab-version
matlabVersion=str2num(regexprep(version('-release'),'[^0-9]','')) ;

if isempty(p.Results.Parent)
    figure();
else
    if strcmpi(get(p.Results.Parent,'type'),'figure')
        figure(p.Results.Parent);
        clf;
    else
        axes(p.Results.Parent);
        cla;
    end
end

box off;

selectedidx=p.Results.ShowSelectedidx;

meanflag=strcmpi(p.Results.ShowMeanValues,'on');
plotexport=~strcmpi(p.Results.Export,'off'); % boolean if plot is used for display only

% Plot light blue dots indicating the entire search space
% do not plot several identical points to reduce the total number of points

dim3in=p.Results.dim3in;
dim4in=p.Results.dim4in;

indys=[];
if ~isempty(p.Results.XLim)
    indys=[indys find(xin<p.Results.XLim)];
end
if ~isempty(p.Results.YLim)
    indys=[indys find(yin<p.Results.YLim)];
end
if ~isempty(p.Results.Dim3Lim)
    indys=[indys find(dim3in<p.Results.Dim3Lim)];
end
if ~isempty(p.Results.Dim4Lim)
    indys=[indys find(dim4in<p.Results.Dim4Lim)];
end
if matlabVersion < 2013
    indys=unique(indys);
    mythin=unique([xin;yin;dim3in;dim4in]','rows');
else
    indys=unique(indys,'legacy');
    mythin=unique([xin;yin;dim3in;dim4in]','rows', 'legacy');
end
xinthin=mythin(:,1)';
yinthin=mythin(:,2)';
dim3inthin=mythin(:,3)';
dim4inthin=mythin(:,4)';
plot(xinthin,yinthin,'.','MarkerSize',5,'Color',[0.7 0.7 0.7]);
hold on

% Get mean
meanx=mean(xin);
meany=mean(yin);

if isempty(indys)
    indys=1:length(xin);
end

% Find pareto values
if isempty(dim3in) && isempty(dim4in)
    [x,y]=paretobest2fastlogic(xin(indys),yin(indys));
else
    [xp4,yp4,dim3,dim4]=paretobest4fastlogic(xin(indys),yin(indys),dim3in(indys),dim4in(indys));
    [x,y]=paretobest2fastlogic(xin(indys),yin(indys));
%     dim3idx=nan(size(x));
%     for i=1:length(x)
%         dim3idx(i)=dim3(find(xp4==x(i) & yp4==y(i)));
%     end
% 
%     [x2,dim32]=paretobest2fastlogic(xin,dim3in);
%     y2=nan(size(x2));
%     for i=1:length(x2)
%         y2(i)=yp4(find(xp4==x2(i) & dim3==dim32(i)));
%     end
% 
%     [y3,dim33]=paretobest2fastlogic(yin,dim3in);
%     x3=nan(size(y3));
%     for i=1:length(y3)
%         x3(i)=xp4(find(yp4==y3(i) & dim3==dim33(i)));
%     end
end

paretoidx=zeros(1,length(x));
for i=1:length(x)
    paretoidx(i)=find(xin==x(i) & yin==y(i),1,'last');
end

xlims2=get(gca,'xlim');
ylims2=get(gca,'ylim');

if xlims2(2)>max(xinthin) || xlims2(2)<max(xinthin)
    xlims2(2)=max(xinthin)*1.1;
end

if ylims2(2)>max(yinthin) || ylims2(2)<max(yinthin)
    ylims2(2)=max(yinthin)*1.1;
end

if xlims2(1)<0
    xlims2(1)=0;
else
    xlims2(1)=min(xinthin)*0.9;
end

if ylims2(1)<0
    ylims2(1)=0;
else
    ylims2(1)=min(yinthin)*0.9;
end

if xlims2(1)>min(xinthin)
    xlims2(1)=min(xinthin)*0.9;
end

if ylims2(1)>min(yinthin)
    ylims2(1)=min(yinthin)*0.9;
end
% if meanflag
%     textXaxis = text('Position',[xlims(2),meany],'String','Mean','HorizontalAlignment','right','VerticalAlignment','top');
%     textYaxis = text('Position',[meanx,ylims(2)],'String','Mean','HorizontalAlignment','right','VerticalAlignment','top','Rotation',90);
%     plot(xlims,repmat(meany,2,1),'k--','LineWidth',2);
%     plot(repmat(meanx,2,1),ylims,'k--','LineWidth',2);
% end

c=[0.9 0.9 0.9];

% [ux,tildevar,uxidx]=unique(xp4);
% mygray=gray(length(ux));
% for i=1:length(ux)
%     idx=find(xp4==ux(i))
%     [tildevar,idxs]=sortrows([yp4(idx)',dim3(idx)'])
%     plot3(ux(i)*ones(size(idxs)),yp4(idx(idxs)),dim3(idx(idxs)),'-','LineWidth',2,'MarkerSize',10,'Color',[mygray(i,1) 0 0],'MarkerFaceColor',c);
% end
% 
% [uy,tildevar,uyidx]=unique(yp4);
% mygray=gray(length(uy));
% for i=1:length(uy)
%     idx=find(yp4==uy(i))
%     [tildevar,idxs]=sortrows([xp4(idx)',dim3(idx)'])
%     plot3(xp4(idx(idxs)),uy(i)*ones(size(idxs)),dim3(idx(idxs)),'-','LineWidth',2,'MarkerSize',10,'Color',[0 mygray(i,1) 0],'MarkerFaceColor',c);
% end
% 
% [udim3,tildevar,udim3idx]=unique(dim3);
% mygray=gray(length(udim3));
% for i=1:length(udim3)
%     idx=find(dim3==udim3(i))
%     [tildevar,idxs]=sortrows([xp4(idx)',yp4(idx)'])
%     plot3(xp4(idx(idxs)),yp4(idx(idxs)),udim3(i)*ones(size(idxs)),'-','LineWidth',2,'MarkerSize',10,'Color',[0 0 mygray(i,1)],'MarkerFaceColor',c);
% end
% if ~(isempty(dim3in) && isempty(dim4in))
plot(x,y,'-','LineWidth',2,'MarkerSize',10,'Color',darken(c),'MarkerFaceColor',c);
% plot3(x2,y2,dim32,'-','LineWidth',2,'MarkerSize',10,'Color',darken(c),'MarkerFaceColor',c);
% plot3(x3,y3,dim33,'-','LineWidth',2,'MarkerSize',10,'Color',darken(c),'MarkerFaceColor',c);
for i=1:length(dim3);
    if ~strcmp(cellstr(p.Results.Dim3Label),cellstr(p.Results.Dim4Label))
        col='g';
        if ((max(dim4)-min(dim4))~=0)
            col=getgreenyellowred((dim4(i)-min(dim4))/(max(dim4)-min(dim4))*1.5);
        end
    else
        col=[0.9 0.9 0.9];
    end
    if dim3(i)==0 || length(dim3)==1
        size=10;
    else
        size=5+10*(dim3(i)-min(dim3))/(max(dim3)-min(dim3));
    end
    plot(xp4(i),yp4(i),'s','MarkerSize',size,'Color',darken(col),'MarkerFaceColor',col);
end
% else
%     plot(x,y,'o-','LineWidth',2,'MarkerSize',10,'Color',darken(c),'MarkerFaceColor',c);
% end
if selectedidx>0
    plot(xin(selectedidx),yin(selectedidx),'x','MarkerSize',10,'MarkerFaceColor',mycolors(4,:),'Color',darken(mycolors(4,:)),'LineWidth',2);
end

%wlims=[min(dim4) max(dim4)/1.5];
wlims=[min(dim4) max(dim4)/1.2];
xlims=[min(dim3) max(dim3)];
%ylims=get(gca,'xlim');
%zlims=get(gca,'ylim');
ylims=xlims2;
zlims=ylims2;

set(gca, 'Layer','top')
box off;
xlabel(p.Results.XLabel);
ylabel(p.Results.YLabel);
if ~isempty(p.Results.XLim)
    xlim([xlims2(1) p.Results.XLim]);
    ylims=[xlims2(1) p.Results.XLim];
else
    xlim(xlims2);
end
if ~isempty(p.Results.YLim)
    ylim([ylims2(1) p.Results.YLim]);
    zlims=[ylims2(1) p.Results.YLim];
else
    ylim(ylims2);
end
rectangle('Position',[ylims(1)+0.01*diff(ylims),zlims(2)-0.26*diff(zlims),(22/26-0.01)*diff(ylims),0.24*diff(zlims)],'EdgeColor','k','FaceColor','w')
if ~strcmp(cellstr(p.Results.Dim3Label),cellstr(p.Results.Dim4Label))
    text(ylims(1)+0.10*diff(ylims),zlims(2)-0.2*diff(zlims),sprintf('%.4f  ',wlims(1)),'HorizontalAlignment','right','VerticalAlignment','middle')
    plot(ylims(1)+0.10*diff(ylims),zlims(2)-0.2*diff(zlims),'s','MarkerFaceColor','g','Color','k','MarkerSize',10)
    text(ylims(1)+0.22*diff(ylims),zlims(2)-0.2*diff(zlims),sprintf('%.4f  ',sum(wlims)/2),'HorizontalAlignment','right','VerticalAlignment','middle')
    plot(ylims(1)+0.22*diff(ylims),zlims(2)-0.2*diff(zlims),'s','MarkerFaceColor','y','Color','k','MarkerSize',10)
    text(ylims(1)+0.34*diff(ylims),zlims(2)-0.2*diff(zlims),sprintf('>%.4f  ',wlims(2)),'HorizontalAlignment','right','VerticalAlignment','middle')
    plot(ylims(1)+0.34*diff(ylims),zlims(2)-0.2*diff(zlims),'s','MarkerFaceColor','r','Color','k','MarkerSize',10)
    text(ylims(1)+0.40*diff(ylims),zlims(2)-0.2*diff(zlims),[p.Results.Dim4Label],'HorizontalAlignment','left','VerticalAlignment','middle')
end
% text(ylims(1)+0.08*diff(ylims),zlims(2)-0.13*diff(zlims),'\rightarrow','HorizontalAlignment','center','VerticalAlignment','middle')
text(ylims(1)+0.10*diff(ylims),zlims(2)-0.08*diff(zlims),sprintf('%.4f  ',xlims(1)),'HorizontalAlignment','right','VerticalAlignment','middle')
plot(ylims(1)+0.10*diff(ylims),zlims(2)-0.08*diff(zlims),'s','MarkerFaceColor','w','Color','k','MarkerSize',5)
text(ylims(1)+0.22*diff(ylims),zlims(2)-0.08*diff(zlims),sprintf('%.4f  ',sum(xlims)/2),'HorizontalAlignment','right','VerticalAlignment','middle')
plot(ylims(1)+0.22*diff(ylims),zlims(2)-0.08*diff(zlims),'s','MarkerFaceColor','w','Color','k','MarkerSize',10)
text(ylims(1)+0.34*diff(ylims),zlims(2)-0.08*diff(zlims),sprintf('>%.4f  ',xlims(2)),'HorizontalAlignment','right','VerticalAlignment','middle')
plot(ylims(1)+0.345*diff(ylims),zlims(2)-0.08*diff(zlims),'s','MarkerFaceColor','w','Color','k','MarkerSize',15)
% text(ylims(1)+0.08*diff(ylims),zlims(2)-0.07*diff(zlims),'\rightarrow','HorizontalAlignment','center','VerticalAlignment','middle')
text(ylims(1)+0.40*diff(ylims),zlims(2)-0.08*diff(zlims),[p.Results.Dim3Label],'HorizontalAlignment','left','VerticalAlignment','middle')

%set(gca,'Position',[0.08 0.25 0.9 0.7])

if plotexport==1 % Export is activated, save as pdf file
    originalPaperPosition=get(gcf,'PaperPosition');
    set(gcf,'Position',p.Results.Position,'PaperSize',originalPaperPosition(3:4),'PaperPosition',[0 0 originalPaperPosition(3:4)]);

    if isempty(p.Results.FileName)
        filename=['paretoPlot.' p.Results.Export];
    else
        filename=p.Results.FileName;
        if ~strcmpi(filename(end-length(p.Results.Export)+1:end),p.Results.Export)
            filename=[filename '.' p.Results.Export];
        end
    end

    saveas(gcf,filename);
end
end