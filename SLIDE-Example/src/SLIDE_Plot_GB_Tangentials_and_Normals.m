function [] = SLIDE_Plot_GB_Tangentials_and_Normals(boundary_info,grains,ebsd)

% Final Verification plot
% Plots one figure for each grain
ebsd_temp = ebsd('indexed'); 

for grain_loop = 1:length(grains)

    close all
    figure('WindowState','maximized')

    % grab the X and Y coordinates of the normal vectors for each ebsd point
    xvectors = boundary_info.(['grainentry' num2str(grain_loop)]).n_pix_vec.x;
    yvectors = boundary_info.(['grainentry' num2str(grain_loop)]).n_pix_vec.y;
    
    % create the corresponding normal vector and normalize it
    veccombi = vector3d([xvectors yvectors zeros(size(xvectors))]').normalize;

    % plot all the datapoints and their vectors
    scatter(ebsd_temp.x(ebsd_temp.grainId == grain_loop),ebsd_temp.y(ebsd_temp.grainId == grain_loop))
    hold on
    plot(boundary_info.(['grainentry' num2str(grain_loop)]).xs,boundary_info.(['grainentry' num2str(grain_loop)]).ys,'LineWidth',3)
    grid on
    quiver(ebsd_temp.x,ebsd_temp.y,veccombi.x',veccombi.y')
    daspect([1 1 1])
    scatter(grains(grain_loop).boundary.x,grains(grain_loop).boundary.y,'k')
    legend('EBSD Points','Smoothened Boundary','Normal Vector Approximation','Raw Boundary Data')
    set(gca, 'YDir','reverse')

    title(['Grain ' num2str(grain_loop)])

    drawnow

    saveas(gcf,['Figures/GB_Tangentials_and_Normals/Grain' num2str(grain_loop) '.png'])

end

close all

end

