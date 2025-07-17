function SLIDE_plotGBS_GBO_2D(input)

%% Unpack Data

IPL_min = input.userSettings.IPL_min;
ebsd = input.ebsd;
grains_SLIDE = input.grains;
incList = input.incList;
grainsPlot = grains_SLIDE.inc_1;

%% Extract the necessary data from each deformation increment

IncrementStruct = struct;

for i = incList

    % extract the 2D sliding magnitude
    IncrementStruct.(['inc_' num2str(i)]).slidingMagnitude2D = grains_SLIDE.(['inc_' num2str(i)]).prop.slidingMagnitude2D_Total*1000;
    
    % extract the sliding/opening identification
    openingIdentification = grains_SLIDE.(['inc_' num2str(i)]).prop.OpeningIdentification_Total;
    openingIdentification(isnan(IncrementStruct.(['inc_' num2str(i)]).slidingMagnitude2D)) = nan;
    IncrementStruct.(['inc_' num2str(i)]).openingIdentification = openingIdentification;

    % load the images
    IncrementStruct.(['inc_' num2str(i)]).IBSE = imread(['Data\ZnCoating_IBSE_Increment' num2str(i) '.tif']);

end

%% Plot Settings and Preparation

fontsizelabel = 14;

% General Plot Settings
fontsize = 14;  %fontsize used for plotting
backgroundcolor = [0.75, 0.75, 0.75]; %define the background color of the plot [gray]

% Calculate a box to show that data close to edge of ROI is removed
ebsdRange = [min(ebsd.x(:)) max(ebsd.x(:)) min(ebsd.y(:)) max(ebsd.y(:))];
ebsdRange = ebsdRange + [1 -1 1 -1]*IPL_min/1000;
ebsdRangeRect = [ebsdRange(1) ebsdRange(3) ebsdRange(2)-ebsdRange(1) ebsdRange(4)-ebsdRange(3)];

%% Plot the opening identification

close all
h = newMtexFigure;

for i = 1:5

    % plot opening IBSE image
    hold on
    imshow(IncrementStruct.(['inc_' num2str(i)]).IBSE);
    set(gca,'Color','none')
    colormap(h.currentAxes,'gray')
    set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',2)
    mtexTitle(['IBSE Increment ' num2str(i)],'Fontsize',fontsize,'interpreter','latex')
    set(gca,'Color',backgroundcolor)
    h.children(i).YLim = [1 4096];
    h.children(i).XLim = [1 4096];
    nextAxis

end

for i = 1:5

    % make 1 pixel opening, otherwise the plot makes all sliding red
    IncrementStruct.(['inc_' num2str(i)]).openingIdentification(1) = 1;

    % plot opening identification
    hold on
    plot(grainsPlot.boundary,IncrementStruct.(['inc_' num2str(i)]).openingIdentification,'LineWidth',5,'micronbar','off');
    set(gca,'Color','none')
    colormap(h.currentAxes,[0 1 0; 1 0 0])
    plot(grainsPlot.boundary,'lineWidth',1,'lineColor',[0 0 0])
    rectangle('Position',ebsdRangeRect,'LineWidth',2,'EdgeColor','r')
    set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',2)
    mtexTitle('GBS/GBO Identification','Fontsize',fontsize,'interpreter','latex')
    set(gca,'Color',backgroundcolor)
    rectangle('Position',ebsdRangeRect,'LineWidth',2,'EdgeColor','r')

    if i~=5
        nextAxis
    else
        cbinc = colorbar('location','eastoutside');
    end

end

%% Positioning the figure

screensize = get(0,'screensize');
h.parent.Position = get(0,'ScreenSize');
drawnow

% Position
spacingx = 10;
spacingy = 50;

offsety = 125;
offsetx = 50;

xlength = screensize(4)/3.5;


% position maps
for i = 1:5
    h.children(i).Position = [offsetx+(xlength+spacingx)*(i-1) offsety+spacingy+xlength xlength xlength];
end
for i = 6:10
    h.children(i).Position = [offsetx+(xlength+spacingx)*(i-6) offsety xlength xlength];
end

drawnow

% Colorbars
baroffset = 5;
barinsertion = 20;
cbSmall = xlength-2*barinsertion;

% color bar for sliding/opening
cbinc.Units = 'pixels';
cbinc.Location = 'southoutside';
cbinc.Position([1 2 3 4]) = [offsetx+barinsertion offsety-20-baroffset cbSmall 20];
cbinc.FontSize = fontsizelabel;
cbinc.Ticks = [0.25 0.75];
cbinc.TickLabels = {'Sliding','Opening'};


end

