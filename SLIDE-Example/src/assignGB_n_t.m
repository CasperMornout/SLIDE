function [boundary_info] = assignGB_n_t(ebsd,grains,padding)

fprintf(['Assigning GB Tangentials and Normals to EBSD datapoints for ' num2str(length(grains)) ' Grains\n'])

% Initialize a struct containing all GB information
boundary_info = struct;

% get the x and y pixels of the EBSD dataset
ebsd = ebsd.gridify;
xpixels = ebsd.x;
ypixels = ebsd.y;

% Loop over Grains
for grain_loop = 1:length(grains)

    % Extract the smoothened grain and some information on it
    grain_temp = grains(grain_loop);

    %%%% --- 1) Extract information about the GB tangential and normal --- %%%%

    % Get grain midpoints (discrete coordinates that described the GB location)
    xs = grain_temp.boundary.midPoint(:,1); ys = grain_temp.boundary.midPoint(:,2);

    % Calculate the gradient in x and y direction of the smoothened
    % boundary and put the results in a vector
    gradx = gradient(xs); grady = gradient(ys);
    smooth_coords = vector3d(xs,ys,zeros(size(xs)));

    % Define the direction vectors of each boundary points
    boundary_directions = vector3d(gradx,grady,zeros(size(gradx)));

    % Enforce the x direction of the GB tangential to always be positive
    negative_x_dirs = boundary_directions.x < 0;
    boundary_directions.x(negative_x_dirs) = - boundary_directions.x(negative_x_dirs);
    boundary_directions.y(negative_x_dirs) = - boundary_directions.y(negative_x_dirs);

    % Calculate the angle between tangential and x-direction
    alpha_t = angle(boundary_directions,xvector)/degree;

    % Define the corresponding normal vectors of each boundary points by rotating 90 degrees over the Z axis.
    rot1 = rotation.byAxisAngle(vector3d(0,0,1),pi/2);
    boundary_normals = rotate(boundary_directions,rot1);

    % Enforce the x direction of the GB normal to always be positive
    negative_x_dirs = boundary_normals.x < 0;
    boundary_normals.x(negative_x_dirs) = - boundary_normals.x(negative_x_dirs);
    boundary_normals.y(negative_x_dirs) = - boundary_normals.y(negative_x_dirs);

    %%%% --- 2) Apply the GB tangential and normal to the surrounding pixels --- %%%%

    % Define structs to store data for each grain separately
    boundary_info.(['grainentry' num2str(grain_loop)]).n_pix.x = nan(size(ebsd('indexed')));
    boundary_info.(['grainentry' num2str(grain_loop)]).n_pix.y = nan(size(ebsd('indexed')));
    boundary_info.(['grainentry' num2str(grain_loop)]).b_pix.x = nan(size(ebsd('indexed')));
    boundary_info.(['grainentry' num2str(grain_loop)]).b_pix.y = nan(size(ebsd('indexed')));

    % define a range of pixels close to the GB
    graintemp_range = [min(grain_temp.x) max(grain_temp.x) min(grain_temp.y) max(grain_temp.y)] + padding/1000*[-1 1 -1 1];

    % select all x/y pixels from the ebsd dataset that are either inside or close to the grain
    pixelsCloseToGrain = ones(size(xpixels));
    pixelsCloseToGrain(xpixels<graintemp_range(1)) = 0;
    pixelsCloseToGrain(xpixels>graintemp_range(2)) = 0;
    pixelsCloseToGrain(ypixels<graintemp_range(3)) = 0;
    pixelsCloseToGrain(ypixels>graintemp_range(4)) = 0;

    % extract the relevant x and y coordinates, and their indices in the ebsd data
    xpixels_temp = xpixels(pixelsCloseToGrain == 1);
    ypixels_temp = ypixels(pixelsCloseToGrain == 1);
    pixels_temp_idx = find(pixelsCloseToGrain == 1);

    % get the closest distance from each ebsd pixel to the GB
    % also get the pixel of the GB to which the distance is close
    [closest_distance, closest_pixel] = min(idpm([xpixels_temp ypixels_temp],[smooth_coords.x smooth_coords.y]),[],2);

    % store in a matrix and remove pixels further than the GB padding
    dis_pix_matrix = [closest_distance closest_pixel pixels_temp_idx];
    dis_pix_matrix(dis_pix_matrix(:,1) > padding/1000,:) = [];

    % get the GB normal/tangential for the required pixels
    bn_x_sel = boundary_normals(dis_pix_matrix(:,2)).x;
    bn_y_sel = boundary_normals(dis_pix_matrix(:,2)).y;
    bd_x_sel = boundary_directions(dis_pix_matrix(:,2)).x;
    bd_y_sel = boundary_directions(dis_pix_matrix(:,2)).y;

    % Store the information on datapointwise GB tangentials and normals in
    % a struct, for each grain separately
    boundary_info.(['grainentry' num2str(grain_temp.id)]).n_pix.x(dis_pix_matrix(:,3)) = bn_x_sel;
    boundary_info.(['grainentry' num2str(grain_temp.id)]).n_pix.y(dis_pix_matrix(:,3)) = bn_y_sel;
    boundary_info.(['grainentry' num2str(grain_temp.id)]).b_pix.x(dis_pix_matrix(:,3)) = bd_x_sel;
    boundary_info.(['grainentry' num2str(grain_temp.id)]).b_pix.y(dis_pix_matrix(:,3)) = bd_y_sel;
    % Store the smoothened GB coordinates, as well as normal and tangential
    boundary_info.(['grainentry' num2str(grain_temp.id)]).xs = xs;
    boundary_info.(['grainentry' num2str(grain_temp.id)]).ys = ys;
    boundary_info.(['grainentry' num2str(grain_temp.id)]).n = boundary_normals;
    boundary_info.(['grainentry' num2str(grain_temp.id)]).b = boundary_directions;
    % Store the angle with the direction of tension
    boundary_info.(['grainentry' num2str(grain_temp.id)]).alpha_t_List = alpha_t;

    % also store the datapointwise GB normals and tangentials as 2D vectors
    boundary_info.(['grainentry' num2str(grain_loop)]).n_pix_vec = ...
        vector3d(boundary_info.(['grainentry' num2str(grain_loop)]).n_pix.x, ...
        boundary_info.(['grainentry' num2str(grain_loop)]).n_pix.y, ...
        zeros(size(boundary_info.(['grainentry' num2str(grain_loop)]).n_pix.x)));

    boundary_info.(['grainentry' num2str(grain_loop)]).b_pix_vec = ...
        vector3d(boundary_info.(['grainentry' num2str(grain_loop)]).b_pix.x, ...
        boundary_info.(['grainentry' num2str(grain_loop)]).b_pix.y, ...
        zeros(size(boundary_info.(['grainentry' num2str(grain_loop)]).b_pix.x)));

end

fprintf(['Assigned GB Tangentials and Normals to EBSD datapoints for ' num2str(length(grains)) ' Grains\n'])


end

