function grains = combineSLIDEresults(input)

% unpack variables
grains = input.grains;
boundary_info = input.boundary_info;
SLIDE_Struct = input.SLIDE_Struct;

%%%% Step 1 - initalize some arrays to store the relevant data %%%%

% GB-Integration of displacements
slidingMagnitude2D_Total = nan(size(grains.boundary));
slidingVector2D_Total = cell(size(grains.boundary));
theta_s_Total = nan(size(grains.boundary));
alpha_s_Total = nan(size(grains.boundary));
alpha_t_Total = nan(size(grains.boundary));

% IBSE Data
OpeningIdentification_Total = nan(size(grains.boundary));

% Confocal Data
slidingMagnitude3D_Total = nan(size(grains.boundary));
dZ_Total = nan(size(grains.boundary));
slidingVector3D_Total = cell(size(grains.boundary));
beta_s_Total = nan(size(grains.boundary));

% get all the midpoints of every grain
gbCoordsAll = grains.boundary.midPoint;

caughtOpening = 0;
caught3D = 0;

%%%% Step 2 - Loop over all grains to insert the data %%%%
for i = 1:length(grains)

    % find the SLIDE information for the grain
    SLIDElocal = SLIDE_Struct.(['grainentry' num2str(i)]);

    % take the midpoints of the relevant grain from the larger list
    gbInfo = boundary_info.(['grainentry' num2str(i)]);
    gbCoordsResult = [gbInfo.xs gbInfo.ys];

    % get the closest matching coordinates
    [~, sm_pixel] = min(idpm(gbCoordsResult,gbCoordsAll),[],2);

    % insert the sliding magnitude / vector / misori in the larger grain2d variable
    slidingMagnitude2D_Total(sm_pixel) = SLIDElocal.sMagInt;
    slidingVector2D_Total(sm_pixel) = SLIDElocal.sVecInt;
    theta_s_Total(sm_pixel) = SLIDElocal.theta_s_List;
    alpha_s_Total(sm_pixel) = SLIDElocal.alpha_s_List;
    alpha_t_Total(sm_pixel) = gbInfo.alpha_t_List;

    try
        OpeningIdentification_Total(sm_pixel) = SLIDElocal.OpeningIdentification;
    catch
        if caughtOpening == 0
            warning('No Opening Data Found. Did you run without using IBSE Opening Identification?')
            caughtOpening = 1;
        end
    end

    if input.CONFdata
        slidingMagnitude3D_Total(sm_pixel) = SLIDElocal.sMagInt3D;
        dZ_Total(sm_pixel) = SLIDElocal.dZList;
        slidingVector3D_Total(sm_pixel) = SLIDElocal.sExpInt3D;
        beta_s_Total(sm_pixel) = SLIDElocal.beta_s_List;
    else
        if caught3D == 0
            warning('No Profilometry Data Found. Did you run 2D SLIDE?')
            caught3D = 1;
        end

    end

end

%%%% Step 3 - Append all data to the smooth GB %%%%
grains.prop.slidingMagnitude2D_Total = slidingMagnitude2D_Total;
grains.prop.slidingVector2D_Total = slidingVector2D_Total;
grains.prop.theta_s_Total = theta_s_Total;
grains.prop.alpha_s_Total = alpha_s_Total;
grains.prop.alpha_t_Total = alpha_t_Total;
try
    grains.prop.OpeningIdentification_Total = OpeningIdentification_Total;
end
if input.CONFdata
    grains.prop.slidingMagnitude3D_Total = slidingMagnitude3D_Total;
    grains.prop.dZ_Total = dZ_Total;
    grains.prop.slidingVector3D_Total = slidingVector3D_Total;
    grains.prop.beta_s_Total = beta_s_Total;
end

end

