function SLIDE_Struct = SLIDE_Datapointwise(ebsd,grains,boundary_info)

    fprintf(['Performing datapoint-wise SLIDE for ' num2str(length(grains)) ' Grains \n'])
    
    % Perform SLIDE for each grain
    for grain_loop = 1:length(grains)
    
        % Set the relevant set of datapointwise GB tangentials and normals
        ebsd.prop.normalvectors = boundary_info.(['grainentry' num2str(grain_loop)]).n_pix_vec.normalize;
        ebsd.prop.directionvectors = boundary_info.(['grainentry' num2str(grain_loop)]).b_pix_vec.normalize;
    
        % Perform SLIDE
        ebsdOUT = SLIDE_Solve_Minimization(ebsd.gridify);
    
        % Store the relevant output in a struct for each grain
        SLIDE_Struct.(['grainentry' num2str(grain_loop)]).ebsdID = ebsdOUT;
        SLIDE_Struct.(['grainentry' num2str(grain_loop)]).grain  = grains(grain_loop);
        
    end
    
    fprintf(['Completed datapoint-wise SLIDE for ' num2str(length(grains)) ' Grains \n'])

end

