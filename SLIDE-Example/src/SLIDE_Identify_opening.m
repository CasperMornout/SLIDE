function SLIDE_Struct_new = SLIDE_Identify_opening(input)

% unpack variables
inc_chosen = input.userSettings.currentInc;
grains = input.grains;
SLIDE_Struct = input.SLIDE_Struct;
boundary_info = input.boundary_info;
IPL_min = input.userSettings.IPL_min;
IPL_min_IBSE = input.userSettings.IPL_min_IBSE;
openingThreshold = input.userSettings.openingThreshold;
numPixUnder = input.userSettings.numPixUnder;

fprintf(['Performing IBSE Sliding/Opening Analysis for ' num2str(length(grains)) ' Grains - Increment ' num2str(input.userSettings.currentInc) '\n'])

% unpack the ebsd data
ebsd = input.ebsd.gridify;
xpixels = ebsd.x;
ypixels = ebsd.y;

% Read the IBSE image and normalize it so values are between 0 and 1
IBSE_img = imread(input.filepath);
IBSE_img = double(IBSE_img)/255;

% create X/Y meshgrids based on the size of the EBSD image
IBSEPsize = input.psize;
[IBSE_X,IBSE_Y] = meshgrid(0:IBSEPsize:IBSEPsize*(size(IBSE_img,2)-1), 0:IBSEPsize:IBSEPsize*(size(IBSE_img,1)-1));

% create a gridded interpolant for IBSE image (to do interpolation)
F_IBSE_Img = griddedInterpolant(IBSE_X',IBSE_Y',IBSE_img','cubic');

% Read the X/Y maps for IBSE Forward-Deformed Alignment
% These fields show, for each EBSD pixel, the corresponding X and Y
% location in the IBSE dataset
IBSE_Xmapping = griddedInterpolant(xpixels',ypixels',input.xIBSE');
IBSE_Ymapping = griddedInterpolant(xpixels',ypixels',input.yIBSE');

% Define a step size for GB-integration along IBSE data
NstepIBSE = ceil(IPL_min_IBSE*10^-3/IBSEPsize);

% Set the number of lines with which we interpolate over the IBSE data
% this means that we collect data for a region with the
% following dimensions: IPL_min_IBSE x (IBSEPsize * length(distlist))
nrSegments = input.userSettings.nrSegments;
distList = -(nrSegments-1)/2:(nrSegments-1)/2;

% loop over each grain
for grain_loop = 1:length(grains)

    %%%% 0) Gather the relevant data %%%%

    % get a list of the boundary coordinates + normal/direction vectors
    bInfoTemp = boundary_info.(['grainentry' num2str(grain_loop)]);
    boundaryCoords = [bInfoTemp.xs bInfoTemp.ys];
    boundaryVecN = normalize(bInfoTemp.n);

    % Create cells and lists to store the results
    IBSE_lineCoordCell = cell(length(boundaryCoords),2);
    IBSE_lineMatrixCell = cell(length(boundaryCoords),1);

    %%%% 1) Forward deform the GB %%%%

    % Forward deform the GB into the forward-deformed IBSE configuration
    boundaryCoordsFWDComplete = [IBSE_Xmapping(boundaryCoords) IBSE_Ymapping(boundaryCoords)];

    %%%% 2) Loop over each GB coord to perform GB-integration on IBSE Data %%%%
    for GBpix = 1:length(boundaryCoords)

        lineCoordX_FWD2_Cell = cell(1,length(distList));
        lineCoordY_FWD2_Cell = cell(1,length(distList));
        IBSE_LineMatrix = nan(length(distList),NstepIBSE+2);

        % make the data nan if the pixel is too close to the edge of the ROI
        if boundaryCoords(GBpix,1)>max(xpixels(:))-IPL_min*10^-3 || ...
                boundaryCoords(GBpix,1)<min(xpixels(:))+IPL_min*10^-3 || ...
                boundaryCoords(GBpix,2)>max(ypixels(:))-IPL_min*10^-3 || ...
                boundaryCoords(GBpix,2)<min(ypixels(:))+IPL_min*10^-3

            % do nothing

        else

            % forward deform the GB center coordinate
            boundaryCoordsFWD = boundaryCoordsFWDComplete(GBpix,:);

            lineBegin = boundaryCoordsFWD - [boundaryVecN(GBpix).x boundaryVecN(GBpix).y]*(IPL_min_IBSE*10^-3/2);
            lineEnd = boundaryCoordsFWD + [boundaryVecN(GBpix).x boundaryVecN(GBpix).y]*(IPL_min_IBSE*10^-3/2);

            % create line segments
            lineCoordX_FWD2 = linspace(lineBegin(1),lineEnd(1),NstepIBSE+2);
            lineCoordY_FWD2 = linspace(lineBegin(2),lineEnd(2),NstepIBSE+2);

            % find the tangential in the deformed configuration
            lineNormal = normalize(vector3d(lineCoordY_FWD2(end)-lineCoordY_FWD2(1),lineCoordX_FWD2(end)-lineCoordX_FWD2(1),0));

            for dist = 1:length(distList)
                lineCoordX_FWD2_Cell{dist} = lineCoordX_FWD2 + lineNormal.x*distList(dist)*IBSEPsize;
                lineCoordY_FWD2_Cell{dist} = lineCoordY_FWD2 - lineNormal.y*distList(dist)*IBSEPsize;
                IBSE_LineMatrix(dist,:)  = F_IBSE_Img(lineCoordX_FWD2_Cell{dist}',lineCoordY_FWD2_Cell{dist}');
            end

        end

        % Save some variables for plotting purposes
        IBSE_lineCoordCell{GBpix,1} = lineCoordX_FWD2_Cell;
        IBSE_lineCoordCell{GBpix,2} = lineCoordY_FWD2_Cell;
        IBSE_lineMatrixCell{GBpix} = IBSE_LineMatrix;

    end

    % store the relevant output
    SLIDE_Struct.(['grainentry' num2str(grain_loop)]).IBSE_lineCoordCell = IBSE_lineCoordCell;
    SLIDE_Struct.(['grainentry' num2str(grain_loop)]).IBSE_lineMatrix = IBSE_lineMatrixCell;

end


%%%% 4) Identify Opening based on brightness threshold %%%%

% loop over each grain
for grain_loop = 1:length(grains)

    % get a list of the boundary coordinates + normal/direction vectors
    bInfoTemp = boundary_info.(['grainentry' num2str(grain_loop)]);
    boundaryCoords = [bInfoTemp.xs bInfoTemp.ys];

    % obtain the IBSE Lines for all GB segments
    IBSE_lineMatrix = SLIDE_Struct.(['grainentry' num2str(grain_loop)]).IBSE_lineMatrix;

    % Create a matrix to store the opening identification
    OpeningIdentification = zeros(length(boundaryCoords),1);

    % loop to take the mean of the line matrix
    for GBpix = 1:length(boundaryCoords)

        IBSE_lineMatrix_GBsegment = IBSE_lineMatrix{GBpix};
        IBSE_mean = mean(IBSE_lineMatrix_GBsegment);
        IBSE_thresholded = double(IBSE_mean < openingThreshold);

        reshapedList=reshape(find(diff([0 IBSE_thresholded 0])~=0),2,[]);
        [pixOpening,~]=max(diff(reshapedList));

        if pixOpening >= numPixUnder
            OpeningIdentification(GBpix) = 1;
        end

    end

    SLIDE_Struct.(['grainentry' num2str(grain_loop)]).OpeningIdentification = OpeningIdentification;

end

SLIDE_Struct_new = SLIDE_Struct;

fprintf(['Completed IBSE Sliding/Opening Analysis for ' num2str(length(grains)) ' Grains - Increment ' num2str(input.userSettings.currentInc) '\n'])

end

