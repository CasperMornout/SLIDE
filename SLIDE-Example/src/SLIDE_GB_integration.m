function SLIDE_Struct_new = SLIDE_GB_integration(input)

% unpack the input
IPL_min = input.userSettings.IPL_min;
residualTH = input.userSettings.residualTH;
strainTH = input.userSettings.strainTH;
distPureSliding = input.userSettings.distPureSliding;
grains = input.grains;
SLIDE_Struct = input.SLIDE_Struct;
boundary_info = input.boundary_info;

% unpack the ebsd data
ebsd = input.ebsd.gridify;
xpixels = ebsd.x;
ypixels = ebsd.y;

% Determine a feasible step size for integration based on EBSD pixelsize
Nstep = ceil(sqrt(2)*2*(IPL_min*10^-3/max(unique(gradient(ebsd.prop.x)))));

fprintf(['Performing GB-integration for ' num2str(length(grains)) ' Grains - Increment ' num2str(input.userSettings.currentInc) '\n'])

for grain_loop =  1:length(grains)

    % Obtain the relevant output data from the SLIDE struct
    ebsdOUT = SLIDE_Struct.(['grainentry' num2str(grain_loop)]).ebsdID;

    %%%% 1) Filtering Output Data %%%%

    % Filtering the SLIDE output on strain and residual
    residualFraction = ebsdOUT.prop.residualFraction;
    strainfieldThreshold = ebsdOUT.prop.Eeff;

    filterList = logical(size(residualFraction));
    filterList(residualFraction > residualTH) = 1; % Filtering 1 - residual
    filterList(strainfieldThreshold < strainTH) = 1; % Filtering 2 - strain minimum

    % Take the filtered sliding vectors and make them nan
    s_exp = ebsdOUT.prop.s_exp;
    s_exp(filterList,:) = nan;
    ebsdOUT.prop.s_exp = s_exp;

    %%%% 2) GB-Integration %%%%

    % obtain sliding direction components in matrix form
    slipdirX = reshape(ebsdOUT.prop.s_exp(:,1),size(ebsdOUT));
    slipdirY = reshape(ebsdOUT.prop.s_exp(:,2),size(ebsdOUT));

    % For interpolation, the nans must be set to 0
    slipdirXinterp = slipdirX; slipdirXinterp(isnan(slipdirXinterp)) = 0;
    slipdirYinterp = slipdirY; slipdirYinterp(isnan(slipdirYinterp)) = 0;

    % create gridded interpolant objects for both sliding directions, this
    % will allow fast interpolation later
    F_slipdirX = griddedInterpolant(xpixels',ypixels',slipdirXinterp','cubic');
    F_slipdirY = griddedInterpolant(xpixels',ypixels',slipdirYinterp','cubic');

    % get a list of the boundary coordinates + normal/direction vectors
    bInfoTemp = boundary_info.(['grainentry' num2str(grain_loop)]);
    boundaryCoords = [bInfoTemp.xs bInfoTemp.ys];
    boundaryVecN = normalize(bInfoTemp.n);
    boundaryVecB = normalize(bInfoTemp.b);

    % Initialize arrays
    sExpIntTotal = cell(length(boundaryCoords),1);  % sliding vector (GB-integrated)
    sMagIntTotal = zeros(length(boundaryCoords),1); % sliding magnitude (GB-integrated)
    theta_s_List = zeros(length(boundaryCoords),1); % misorientation to GB direction (Theta_s)
    alpha_s_List = zeros(length(boundaryCoords),1); % misorientation to direction of tension (Alpha_s)

    % loop over every GB coordinate to perform GB-integration
    for GBpix = 1:length(boundaryCoords)

        % define the start and end of the integration line
        lineBegin = boundaryCoords(GBpix,:) - [boundaryVecN(GBpix).x boundaryVecN(GBpix).y]*(IPL_min*10^-3/2);
        lineEnd = boundaryCoords(GBpix,:) + [boundaryVecN(GBpix).x boundaryVecN(GBpix).y]*(IPL_min*10^-3/2);

        % create line segments
        lineCoordX = linspace(lineBegin(1),lineEnd(1),Nstep+2);
        lineCoordY = linspace(lineBegin(2),lineEnd(2),Nstep+2);

        % make the data nan if the pixel is too close to the edge of the
        % EBSD dataset because interpolation would not be possible
        if boundaryCoords(GBpix,1)>max(xpixels(:))-IPL_min*10^-3 || ...
                boundaryCoords(GBpix,1)<min(xpixels(:))+IPL_min*10^-3 || ...
                boundaryCoords(GBpix,2)>max(ypixels(:))-IPL_min*10^-3 || ...
                boundaryCoords(GBpix,2)<min(ypixels(:))+IPL_min*10^-3

            theta_s = nan;
            alpha_s = nan;
            sMagInt = nan;
            sVecInt = vector3d(nan,nan,nan);

        else

            % calculate distance along the integration path
            distCoord = [0 cumsum((diff(lineCoordY).^2 + diff(lineCoordX).^2).^0.5)];

            % interpolate the sliding vector (x and y components) on the
            % line segment
            sVecInterp_x = F_slipdirX(lineCoordX',lineCoordY');
            sVecInterp_y = F_slipdirY(lineCoordX',lineCoordY');

            % interpolate the x and y values of the vector over the length
            % of the line segment
            sVecInterp_x_Int = trapz(distCoord,sVecInterp_x);
            sVecInterp_y_Int = trapz(distCoord,sVecInterp_y);

            % create a total (integrated) sliding vector and magnitude
            sVecInt = vector3d(sVecInterp_x_Int,sVecInterp_y_Int,0);
            sMagInt = norm(sVecInt);

            % Calculate angles between sliding vector (from SLIDE)
            % and the GB tangential (from EBSD) [Theta]
            theta_s0 = angle(sVecInt,boundaryVecB(GBpix))/degree;
            theta_s180 = theta_s0 - 180;
            theta_s_Tot = [theta_s0 theta_s180];
            theta_s = min(abs(theta_s_Tot),[],2);

            alpha_s0 = angle(sVecInt,xvector)/degree;
            alpha_s180 = alpha_s0 - 180;
            alpha_s_Tot = [alpha_s0 alpha_s180];
            alpha_s = min(abs(alpha_s_Tot),[],2);

            % Calculate a number of elements that should have some result
            stepSize = diff(distCoord(1:2));
            numElPureSliding = floor(distPureSliding/stepSize);

            % If the number of SLIDE-identified pixels along the line is
            % too low, we remove the data points
            if numel(sVecInterp_x(sVecInterp_x ~= 0)) < numElPureSliding
                theta_s = nan;
                sMagInt = nan;
                sVecInt = vector3d(nan,nan,nan);
            end

        end

        % save some relevant variables
        sExpIntTotal{GBpix} = sVecInt;
        sMagIntTotal(GBpix) = sMagInt;
        theta_s_List(GBpix) = theta_s;
        alpha_s_List(GBpix) = alpha_s;


    end

    % store the relevant output
    SLIDE_Struct.(['grainentry' num2str(grain_loop)]).sVecInt = sExpIntTotal;
    SLIDE_Struct.(['grainentry' num2str(grain_loop)]).sMagInt = sMagIntTotal;
    SLIDE_Struct.(['grainentry' num2str(grain_loop)]).theta_s_List = theta_s_List;
    SLIDE_Struct.(['grainentry' num2str(grain_loop)]).alpha_s_List = alpha_s_List;

end

SLIDE_Struct_new = SLIDE_Struct;

fprintf(['Completed GB-integration for ' num2str(length(grains)) ' Grains - Increment ' num2str(input.userSettings.currentInc) '\n'])


end

