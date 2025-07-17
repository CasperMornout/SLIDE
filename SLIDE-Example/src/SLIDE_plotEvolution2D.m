function SLIDE_plotEvolution2D(input)

%% Unpack Data

IPL_min = input.userSettings.IPL_min;
ebsd = input.ebsd;
grains_SLIDE = input.grains;
incList = input.incList;
grainsPlot = grains_SLIDE.inc_1;


%% Extract the necessary data from each deformation increment

IncrementStruct = struct;

for i = incList

    % In the SLIDE code (SLIDE_Solve_Minimization), the displacements are
    % filtered and therefore the strain field is also slightly filtered.
    % Therefore, we take a filtered strain field here to plot.
    IncrementStruct.(['inc_' num2str(i)]).E_eff = ebsd.prop.(['inc_' num2str(i) '_Eequi_Filter']);

    % extract the 2D sliding magnitude
    IncrementStruct.(['inc_' num2str(i)]).slidingMagnitude2D = grains_SLIDE.(['inc_' num2str(i)]).prop.slidingMagnitude2D_Total*1000;
    
end

%% Plot Settings and Preparation

fontsizelabel = 14;

% General Plot Settings
fontsize = 14;  %fontsize used for plotting
backgroundcolor = [0.75, 0.75, 0.75]; %define the background color of the plot [gray]

% color limits for each axis
clim_equi = [0.02 1];
clim_int = [10 800]; % used to be [50 800] 
clim_misori =[0 90];

% Calculate a box to show that data close to edge of ROI is removed
ebsdRange = [min(ebsd.x(:)) max(ebsd.x(:)) min(ebsd.y(:)) max(ebsd.y(:))];
ebsdRange = ebsdRange + [1 -1 1 -1]*IPL_min/1000;
ebsdRangeRect = [ebsdRange(1) ebsdRange(3) ebsdRange(2)-ebsdRange(1) ebsdRange(4)-ebsdRange(3)];

%% Plotting the evolution of strain and sliding

close all
figure('WindowState','maximized')
h = newMtexFigure('ncols',5,'nrows',3);

for i = incList

    % plot strain fields
    if i == 1
    plot(ebsd,IncrementStruct.(['inc_' num2str(i)]).E_eff)
    else
    plot(ebsd,IncrementStruct.(['inc_' num2str(i)]).E_eff,'micronbar','off')
    end
    hold on
    h.currentAxes.CLim = clim_equi;
    colormap(h.currentAxes,'viridis')
    set(gca,'Color',backgroundcolor)
    set(gca,'colorscale','log')
    set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',2)
    hold on
    micronbar1(5,1,[0 0 0],14,'k')
    mtexTitle(['$\hat{E}_{eff}$ (Increment ' num2str(i) ') [-]'],'Fontsize',fontsize,'interpreter','latex')

    if i == 5
        cbStrain = colorbar('location','eastoutside');
    end

        nextAxis

end

for i = 1:5

    % plot sliding magnitude in 2D
    hold on
    plot(grainsPlot.boundary,IncrementStruct.(['inc_' num2str(i)]).slidingMagnitude2D,'LineWidth',5,'micronbar','off');
    set(gca,'Color','none')
    set(gca,'ColorScale','log')
    set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',2)
    colormap(h.currentAxes,'jet')
    hold on
    clim(clim_int)
    mtexTitle('$\bar{\gamma}^{s/o}_{2D}$ [nm]','Fontsize',fontsize,'interpreter','latex')
    set(gca,'Color',backgroundcolor)
    plot(grainsPlot.boundary,'lineWidth',1,'lineColor',[0 0 0])
    rectangle('Position',ebsdRangeRect,'LineWidth',2,'EdgeColor','r')

    if i ~= 5
        nextAxis
    else
        cbSliding = colorbar('location','eastoutside');
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

% color bar for strain
cbStrain.Units = 'pixels';
cbStrain.Position([1 2 3 4]) = [offsetx+(xlength+spacingx)*5+baroffset offsety+xlength+spacingy+barinsertion 20 cbSmall];
cbStrain.FontSize = fontsizelabel; 
cbStrain.Ticks = [0.02 0.5 1];
cbStrain.TickLabels = {'0.02' '0.5' '1'};
cbStrain.Label.String = '';

% color bar for sliding
cbSliding.Units = 'pixels';
cbSliding.Position([1 2 3 4]) = [offsetx+(xlength+spacingx)*5+baroffset offsety+barinsertion 20 cbSmall];
cbSliding.FontSize = fontsizelabel;
cbSliding.Ticks = [10 50 200  800];
cbSliding.TickLabels = {'10' '50' '200' '800'};
cbSliding.Label.String = '';

end

