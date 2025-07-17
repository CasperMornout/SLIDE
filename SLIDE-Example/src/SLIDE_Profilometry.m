function SLIDE_Struct_new = SLIDE_Profilometry(input)

% unpack variables
grains = input.grains;
SLIDE_Struct = input.SLIDE_Struct;
boundary_info = input.boundary_info;
IPL_min = input.userSettings.IPL_min;
Confocal_X = input.Confocal_X;
Confocal_Y = input.Confocal_Y;
Confocal_Z = input.Confocal_Z;
IPL_min_CONF = input.userSettings.IPL_min_CONF;

fprintf(['Performing Profilometry GB-integration for ' num2str(length(grains)) ' Grains - Increment ' num2str(input.userSettings.currentInc) '\n'])

% unpack the ebsd data
ebsd = input.ebsd.gridify;
xpixels = ebsd.x;
ypixels = ebsd.y;

% set a step size for interpolation confocal data. 
% Step size is hardcoded here
NstepConfocal = 100;

% create gridded interpolant for forward-deformed alignment of confocal
% data (for each EBSD pixel, these fields shown the corresponding X and Y
% location in the Confocal dataset)
CONF_Xmapping = griddedInterpolant(xpixels',ypixels',input.xCONF');
CONF_Ymapping = griddedInterpolant(xpixels',ypixels',input.yCONF');

% create a gridded interpolant for Confocal Z Data (to do interpolation)
F_CONF = griddedInterpolant(Confocal_X',Confocal_Y',Confocal_Z','cubic');

% loop over each grain
for grain_loop = 1:length(grains)

    %%%% 0) Gather the relevant data %%%%

    % get a list of the boundary coordinates + normal/direction vectors
    bInfoTemp = boundary_info.(['grainentry' num2str(grain_loop)]);
    boundaryCoords = [bInfoTemp.xs bInfoTemp.ys];
    boundaryVecN = normalize(bInfoTemp.n);

    % take out the SLIDE 2D vectors for this particular grain
    sExpIntTotal = SLIDE_Struct.(['grainentry' num2str(grain_loop)]).sVecInt;

    % Create cells and lists to store the results
    dZList = zeros(length(boundaryCoords),1);
    betaList = zeros(length(boundaryCoords),1);
    slidingMagnitude3D = zeros(length(boundaryCoords),1);

    %%%% 1) Forward deform the GB %%%%

    % Forward deform the GB into the forward-deformed Confocal configuration
    boundaryCoordsFWDComplete = [CONF_Xmapping(boundaryCoords) CONF_Ymapping(boundaryCoords)];

    %%%% 2) Loop over each GB coord to calculate line profiles on Confocal Data %%%%
    for GBpix = 1:length(boundaryCoords)

        % if in-plane sliding is nonexistent, no 3D vector can be
        % constructed, and the profilometry analysis is skipped
        if isnan(sExpIntTotal{GBpix}.x)

            % fill everything in as nan
            dZList(GBpix) = nan;
            betaList(GBpix) = nan;
            slidingMagnitude3D(GBpix) = nan;
            sExpIntTotal{GBpix}.z = nan;

            continue

        end

        % make the data nan if the pixel is too close to the edge of the ROI
        if boundaryCoords(GBpix,1)>max(xpixels(:))-IPL_min*10^-3 || ...
                boundaryCoords(GBpix,1)<min(xpixels(:))+IPL_min*10^-3 || ...
                boundaryCoords(GBpix,2)>max(ypixels(:))-IPL_min*10^-3 || ...
                boundaryCoords(GBpix,2)<min(ypixels(:))+IPL_min*10^-3

            % fill everything in as nan
            dZList(GBpix) = nan;
            betaList(GBpix) = nan;
            slidingMagnitude3D(GBpix) = nan;
            sExpIntTotal{GBpix}.z = nan;

        else

            % forward deform the GB center coordinate
            boundaryCoordsFWD = boundaryCoordsFWDComplete(GBpix,:);

            % define the start and end of the integration line
            lineBegin = boundaryCoordsFWD - [boundaryVecN(GBpix).x boundaryVecN(GBpix).y]*((IPL_min_CONF+400)*10^-3/2);
            lineEnd = boundaryCoordsFWD + [boundaryVecN(GBpix).x boundaryVecN(GBpix).y]*((IPL_min_CONF+400)*10^-3/2);

            % create line segments
            lineCoordX_FWD2 = linspace(lineBegin(1),lineEnd(1),NstepConfocal+2);
            lineCoordY_FWD2 = linspace(lineBegin(2),lineEnd(2),NstepConfocal+2);

            % interpolate CONFOCAL data on the line segment
            CONF_line = F_CONF(lineCoordX_FWD2',lineCoordY_FWD2');
            % calculate the distance along the line
            distX = cumsum(abs(diff(lineCoordX_FWD2)));
            distY = cumsum(abs(diff(lineCoordY_FWD2)));
            distCoordC = [0 (distX.^2+distY.^2).^0.5];

            % get the first and last bit (0.2 micron)
            lengthCuts = 0.2;
            CONF_begin_pix = distCoordC<lengthCuts;
            CONF_end_pix = distCoordC>max(distCoordC)-lengthCuts;
            CONF_lineBegin = mean(CONF_line(CONF_begin_pix));
            CONF_lineEnd = mean(CONF_line(CONF_end_pix));

            % make sure that the difference is taken from left to right
            if lineCoordX_FWD2(1) < lineCoordX_FWD2(end)
                dZ = CONF_lineEnd - CONF_lineBegin;
            elseif lineCoordX_FWD2(1) >= lineCoordX_FWD2(end)
                dZ = CONF_lineBegin - CONF_lineEnd;
            end

            % save the heigt step and also construct 3d Sliding vector
            dZList(GBpix) = dZ;
            sExpIntTotal{GBpix}.z = dZ;
            slidingMagnitude3D(GBpix) = norm(sExpIntTotal{GBpix});

            % calculate the angle between vector and Z axis (beta)
            betaList(GBpix) = 90 - angle(sExpIntTotal{GBpix},zvector)/degree;

        end
    end

    SLIDE_Struct.(['grainentry' num2str(grain_loop)]).dZList = dZList;
    SLIDE_Struct.(['grainentry' num2str(grain_loop)]).sExpInt3D = sExpIntTotal;
    SLIDE_Struct.(['grainentry' num2str(grain_loop)]).beta_s_List = betaList;
    SLIDE_Struct.(['grainentry' num2str(grain_loop)]).sMagInt3D = slidingMagnitude3D;

end

fprintf(['Completed Profilometry GB-integration for ' num2str(length(grains)) ' Grains - Increment ' num2str(input.userSettings.currentInc) '\n'])

SLIDE_Struct_new = SLIDE_Struct;

end

