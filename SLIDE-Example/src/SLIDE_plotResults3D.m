function SLIDE_plotResults3D(input)

%% Unpack Data

IPL_min = input.userSettings.IPL_min;
ebsd = input.ebsd;
grains = input.grains;
inc_chosen = input.userSettings.currentInc;

%% Define sliding magnitudes to plot

slidingMagnitude2d = grains.prop.slidingMagnitude2D_Total*1000;
slidingMagnitude3d = grains.prop.slidingMagnitude3D_Total*1000;
dZTotal_Absolute = abs(grains.prop.dZ_Total*1000);
beta_Absolute = abs(grains.prop.beta_s_Total);

%% Define a (filtered) ebsd strain field to plot

% In the SLIDE code (SLIDE_Solve_Minimization), the displacements are
% filtered and therefore the strain field is also slightly filtered.
% Therefore, we take a filtered strain field here to plot.
E_eff_Filtered = ebsd.prop.(['inc_' num2str(inc_chosen) '_Eequi_Filter']);

%% Plot Settings and Preparation

fontsizelabel = 14;

% General Plot Settings
fontsize = 22;  %fontsize used for plotting
backgroundcolor = [0.75, 0.75, 0.75]; %define the background color of the plot [gray]

% color limits for each axis
clim_equi = [0.02 1];
clim_int = [50 800];
clim_misori =[0 90];

% define a lavender colormap for the density plot
cmapLavender = flip([linspace(75,235,25)' linspace(0,235,25)' linspace(130,255,25)']/255);

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

nextAxis

% plotting the sliding magnitude Z
hold on

plot(grains.boundary,dZTotal_Absolute,'LineWidth',5,'micronbar','off');
hold on
rectangle('Position',ebsdRangeRect,'LineWidth',2,'EdgeColor','r')
clim([0 0.3])
set(gca,'ColorScale','log')
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',2)
colormap(h.currentAxes,'jet')
hold on
clim(clim_int)
mtexTitle('$d\bar{Z}$ [nm]','Fontsize',fontsize,'interpreter','latex')
set(gca,'Color',backgroundcolor)
plot(grains.boundary,'lineWidth',1,'lineColor',[0 0 0])
rectangle('Position',ebsdRangeRect,'LineWidth',2,'EdgeColor','r')

nextAxis

% plotting the 3D sliding magnitude
hold on
plot(grains.boundary,slidingMagnitude3d,'LineWidth',5,'micronbar','off');
set(gca,'ColorScale','log')
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',2)
colormap(h.currentAxes,'jet')
hold on
cbsliding = colorbar('location','southoutside');
clim(clim_int)
mtexTitle('$\bar{\gamma}^{s/o}$ [nm]','Fontsize',fontsize,'interpreter','latex')
set(gca,'Color',backgroundcolor)
plot(grains.boundary,'lineWidth',1,'lineColor',[0 0 0])
rectangle('Position',ebsdRangeRect,'LineWidth',2,'EdgeColor','r')

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

nextAxis

% plotting the beta angle
hold on
plot(grains.boundary,beta_Absolute,'LineWidth',5,'micronbar','off');
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',2)
colormap(h.currentAxes,'jet')
hold on
clim(clim_misori)
mtexTitle('$\bar{\beta}^{s/o}$ [$^{\circ}$]','Fontsize',fontsize,'interpreter', 'latex')
set(gca,'Color',backgroundcolor)
plot(grains.boundary,'lineWidth',1,'lineColor',[0 0 0])
rectangle('Position',ebsdRangeRect,'LineWidth',2,'EdgeColor','r')


nextAxis

% plotting the sliding angle distribution
theta_s_List = grains.prop.theta_s_Total;
alpha_t_List = grains.prop.alpha_t_Total;

theta_s_List = theta_s_List(openingIdentification == 0);
alpha_t_List = alpha_t_List(openingIdentification == 0);

beta_s_List = beta_Absolute(openingIdentification == 0);

dscatter(alpha_t_List,theta_s_List,'plottype','contour','smoothing',2)
colormap(gca,cmapLavender)
ylim([0 90])
xlim([0 90])
axis square
set(gca, 'xtick', [], 'ytick', [],'linewidth',2)

nextAxis
scatter(alpha_t_List,theta_s_List,6,beta_s_List,'filled')
xlabel('$\bar{\alpha}^{t}$ [$^\circ$]','FontSize',18,'Interpreter','latex')
ylabel('$\bar{\theta}^{s}$ [$^\circ$]','FontSize',18,'Interpreter','latex')
axis square
ylim([0 90])
xlim([0 90])
clim([0 90])
set(gca,'Color','none')
colormap(gca,'jet')
cbbeta = colorbar('location','southoutside');
cbbeta.Label.String = 'Out-of-Plane Sliding Angle \beta [\circ]';
set(gca,'xaxisLocation','top')
xticks([0 30 60 90])
yticks([0 30 60 90])
title('Sliding Angles - GBS Only')

%% Positioning the figure

screensize = get(0,'screensize');
h.parent.Position = get(0,'ScreenSize');
drawnow

% Position
spacingx = 10;
spacingy = 100;

offsety = 125;
offsetx = 50;

xlength = screensize(4)/3;

inset9 = 23;


% position maps
h.children(1).Position = [offsetx offsety+spacingy+xlength xlength xlength];
h.children(2).Position = [offsetx+(xlength+spacingx)*1 offsety+spacingy+xlength xlength xlength];
h.children(3).Position = [offsetx+(xlength+spacingx)*2 offsety+spacingy+xlength xlength xlength];
h.children(4).Position = [offsetx+(xlength+spacingx)*3 offsety+spacingy+xlength xlength xlength];

h.children(5).Position = [offsetx offsety xlength xlength];
h.children(6).Position = [offsetx+(xlength+spacingx)*1 offsety xlength xlength];
h.children(7).Position = [offsetx+(xlength+spacingx)*2 offsety xlength xlength];

h.children(8).Position = [offsetx+(xlength+spacingx)*3+inset9 offsety xlength-inset9 xlength-inset9];
h.children(9).Position = [offsetx+(xlength+spacingx)*3+inset9 offsety xlength-inset9 xlength-inset9];

drawnow

h.children(9).YLabel.Position(1) = 0;
h.children(9).XLabel.Position(2) = 90;
h.children(9).YLabel.FontSize = 20;
h.children(9).XLabel.FontSize = 20;

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
cbsliding.Position([1 2 3 4]) = [offsetx+xlength*1.5+spacingx*2+barinsertion2 offsety+xlength+spacingy-baroffset-20 cbBig 20];
cbsliding.FontSize = fontsizelabel;
cbsliding.Ticks = [50 100 200 400 800];
cbsliding.TickLabels = {'50' '100' '200' '400' '800'};
cbsliding.Label.String = '';

% color bar for sliding/opening
cbopening.Units = 'pixels';
cbopening.Location = 'southoutside';
cbopening.Position([1 2 3 4]) = [offsetx+barinsertion offsety-20-baroffset cbSmall 20];
cbopening.FontSize = fontsizelabel;
cbopening.Ticks = [0.25 0.75];
cbopening.TickLabels = {'Sliding','Opening'};

% color bar for theta
cbangle.Units = 'pixels';
cbangle.Position([1 2 3 4]) = [offsetx+xlength+spacingx+barinsertion2 offsety-baroffset-20 cbBig 20];
cbangle.FontSize = fontsizelabel;
cbangle.Ticks = [0 45 90];
cbangle.TickLabels = {'0' '45' '90'};
cbangle.Label.String = '';

% color bar for residual
cbbeta.Units = 'pixels';
cbbeta.Position([1 2 3 4]) = [offsetx+xlength*3+spacingx*3+barinsertion+inset9 offsety-baroffset-20 cbSmall-inset9 20];
cbbeta.FontSize = fontsizelabel;
cbbeta.Ticks = [0 30 60 90];
cbbeta.TickLabels = {'0' '30' '60' '90'};
cbbeta.Label.String = '$\bar{\beta}^{s}$ [$^\circ$]';
cbbeta.Label.FontSize = 20;
cbbeta.Label.Position(2) = -0.2;
cbbeta.Label.Interpreter = 'latex';



end

