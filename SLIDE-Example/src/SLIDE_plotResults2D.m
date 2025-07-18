function SLIDE_plotResults2D(input)

%% Unpack Data

IPL_min = input.userSettings.IPL_min;
ebsd = input.ebsd;
grains = input.grains;
inc_chosen = input.userSettings.currentInc;

%% Define sliding magnitudes to plot

slidingMagnitude2d = grains.prop.slidingMagnitude2D_Total*1000;

%% Define an ebsd strain field to plot

E_eff_Filtered = ebsd.prop.(['inc_' num2str(inc_chosen) '_Eequi']);

%% Plot Settings and Preparation

fontsizelabel = 14;

% General Plot Settings
fontsize = 22;  %fontsize used for plotting
backgroundcolor = [0.75, 0.75, 0.75]; %define the background color of the plot [gray]

% color limits for each axis
clim_equi = [0.02 1];
clim_int = [50 800];
clim_misori =[0 90];

% Calculate a box to show that data close to edge of ROI is removed
ebsdRange = [min(ebsd.x(:)) max(ebsd.x(:)) min(ebsd.y(:)) max(ebsd.y(:))];
ebsdRange = ebsdRange + [1 -1 1 -1]*IPL_min/1000;
ebsdRangeRect = [ebsdRange(1) ebsdRange(3) ebsdRange(2)-ebsdRange(1) ebsdRange(4)-ebsdRange(3)];

%% Plotting

close all
figure('WindowState','maximized')
h = newMtexFigure('ncols',4);

% plot strain fields
plot(ebsd,E_eff_Filtered)
hold on
h.currentAxes.CLim = clim_equi;
colormap(h.currentAxes,'viridis')
cbstrain = colorbar('location','southoutside');
set(gca,'Color',backgroundcolor)
set(gca,'colorscale','log')
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',2)
hold on
micronbar1(5,1,[0 0 0],14,'k ')
mtexTitle('$\hat{E}_{eff}$ ($\varepsilon_g \ = \ 5.88 % $) [-]','Fontsize',fontsize,'interpreter','latex')

nextAxis

% plotting the sliding magnitude in 2D
hold on
plot(grains.boundary,slidingMagnitude2d,'LineWidth',5,'micronbar','off');
set(gca,'Color','none')
set(gca,'ColorScale','log')
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',2)
colormap(h.currentAxes,'jet')
hold on
clim(clim_int)
mtexTitle('$\bar{\gamma}^{s/o}_{2D}$ [nm]','Fontsize',fontsize,'interpreter','latex')
set(gca,'Color',backgroundcolor)
plot(grains.boundary,'lineWidth',1,'lineColor',[0 0 0])
rectangle('Position',ebsdRangeRect,'LineWidth',2,'EdgeColor','r')
cbsliding = colorbar('location','southoutside');

nextAxis

% plotting the opening and sliding distinction
openingIdentification = grains.prop.OpeningIdentification_Total;
openingIdentification(isnan(slidingMagnitude2d)) = nan;

hold on
plot(grains.boundary,openingIdentification,'LineWidth',5,'micronbar','off');
set(gca,'Color','none')
colormap(h.currentAxes,[0 1 0; 1 0 0])
plot(grains.boundary,'lineWidth',1,'lineColor',[0 0 0])
rectangle('Position',ebsdRangeRect,'LineWidth',2,'EdgeColor','r')
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',2)
mtexTitle('GBS/GBO Identification','Fontsize',fontsize-4,'interpreter','latex')
set(gca,'Color',backgroundcolor)
rectangle('Position',ebsdRangeRect,'LineWidth',2,'EdgeColor','r')
cbopening = colorbar('location','southoutside');

nextAxis

% plotting the theta angle
hold on
plot(grains.boundary,grains.prop.theta_s_Total,'LineWidth',5,'micronbar','off');
set(gca,'Color','none')
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',2)
colormap(h.currentAxes,'jet')
hold on
cbangle = colorbar('location','southoutside');
clim(clim_misori)
mtexTitle('$\bar{\theta}^{s/o}$ [$^{\circ}$]','Fontsize',fontsize,'interpreter', 'latex')
set(gca,'Color',backgroundcolor)
plot(grains.boundary,'lineWidth',1,'lineColor',[0 0 0])
rectangle('Position',ebsdRangeRect,'LineWidth',2,'EdgeColor','r')


%% Positioning the figure

screensize = get(0,'screensize');
h.parent.Position = get(0,'ScreenSize');
drawnow

% Position
spacingx = 10;
spacingy = 0;

offsety = 0;
offsetx = 50;

xlength = screensize(4)/3;

% position maps
h.children(1).Position = [offsetx offsety+spacingy+xlength xlength xlength];
h.children(2).Position = [offsetx+(xlength+spacingx)*1 offsety+spacingy+xlength xlength xlength];
h.children(3).Position = [offsetx+(xlength+spacingx)*2 offsety+spacingy+xlength xlength xlength];
h.children(4).Position = [offsetx+(xlength+spacingx)*3 offsety+spacingy+xlength xlength xlength];

drawnow

% Colorbars
baroffset = 5;
barinsertion = 20;
barinsertion2 = 100;
cbSmall = xlength-2*barinsertion;
cbBig = xlength*2+spacingx - 2*barinsertion2;

% color bar for strain
cbstrain.Units = 'pixels';
cbstrain.Position([1 2 3 4]) = [offsetx+barinsertion offsety+xlength+spacingy-baroffset-20 cbSmall 20];
cbstrain.FontSize = fontsizelabel;
cbstrain.Ticks = [0.02 0.5 1];
cbstrain.TickLabels = {'0.02' '0.5' '1'};
cbstrain.Label.String = '';

% color bar for sliding
cbsliding.Units = 'pixels';
cbsliding.Position([1 2 3 4]) = [offsetx+barinsertion+xlength+spacingx offsety+xlength+spacingy-baroffset-20 cbSmall 20];
cbsliding.FontSize = fontsizelabel;
cbsliding.Ticks = [50 100 200 400 800];
cbsliding.TickLabels = {'50' '100' '200' '400' '800'};
cbsliding.Label.String = '';

% color bar for sliding/opening
cbopening.Units = 'pixels';
cbopening.Location = 'southoutside';
cbopening.Position([1 2 3 4]) = [offsetx+barinsertion+xlength*2+spacingx*2 offsety+xlength+spacingy-baroffset-20 cbSmall 20];
cbopening.FontSize = fontsizelabel;
cbopening.Ticks = [0.25 0.75];
cbopening.TickLabels = {'Sliding','Opening'};

% color bar for theta
cbangle.Units = 'pixels';
cbangle.Position([1 2 3 4]) = [offsetx+barinsertion+xlength*3+spacingx*3 offsety+xlength+spacingy-baroffset-20 cbSmall 20];
cbangle.FontSize = fontsizelabel;
cbangle.Ticks = [0 45 90];
cbangle.TickLabels = {'0' '45' '90'};
cbangle.Label.String = '';




end

