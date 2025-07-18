function ebsdOUT = SLIDE_Solve_Minimization(ebsdIN)

%%%% Performs SLIDE %%%%
% collect U and V and check for NaN values
data.U = ebsdIN.prop.U.x;
data.V = ebsdIN.prop.U.y;
nanU = isnan(ebsdIN.prop.U.x);
data.U(nanU) = NaN;
data.V(nanU) = NaN;

% create dummy EBSD to contain the output
ebsdOUT = dummyEBSDSimple(ebsdIN.orientations(1),ebsdIN.prop.x,ebsdIN.prop.y);

% calculate the pixelsize for gradient computation
pixelsize = max(unique(gradient(ebsdOUT.prop.x)));

% calculate displacement gradient tensor components
[Dxx, Dxy] = gradient(data.U,pixelsize,pixelsize);
[Dyx, Dyy] = gradient(data.V,pixelsize,pixelsize);

% calculate effective shear strain
Eeff = calcEffectiveE(Dxx,Dxy,Dyx,Dyy);

% store useful fields in ebsdOUT
ebsdOUT.prop.U = data.U;
ebsdOUT.prop.V = data.V;
ebsdOUT.prop.Hxx = Dxx;
ebsdOUT.prop.Hxy = Dxy;
ebsdOUT.prop.Hyx = Dyx;
ebsdOUT.prop.Hyy = Dyy;
ebsdOUT.prop.Eeff = Eeff;

% Store displacement gradient tensor components in a matrix
D = zeros(2,2,length(Dxx(:)));
D(1,1,:) = Dxx(:);
D(1,2,:) = Dxy(:);
D(2,1,:) = Dyx(:);
D(2,2,:) = Dyy(:);

% preparation of arrays
residual_fraction_list = nan(length(Dxx(:)),1); % normalized residual norm R
s_exp_list = nan(length(Dxx(:)),2); % sliding direction s
n_exp_list = ebsdIN('indexed').prop.normalvectors; % GB normal

% define the relevant pixels that need to be analysed 
% any pixel that is has a non-nan GB normal needs to be analyzed
analysisvectors_rel_idx = find(~isnan(n_exp_list.x));

% loop over each pixel that required analysis
for i = analysisvectors_rel_idx'

    % construct n_exp and D as matrices
    n_exp = n_exp_list(i);
    n_exp = n_exp.matrix;
    Di = D(:,:,i);

    % calculate experimental sliding direction by performing a
    % least-squares optimization
    s_exp_i = Di / n_exp(1:2)';
    s_exp_list(i,:) = s_exp_i;

    % calculate the normalized residual norm R
    residual_fraction_list(i) = norm(Di-s_exp_i*n_exp(1:2)')/norm(Di);

end

% assign all outcomes to the output EBSD
ebsdOUT.prop.residualFraction = residual_fraction_list;
ebsdOUT.prop.s_exp = s_exp_list;

end


