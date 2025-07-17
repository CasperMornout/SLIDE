function [ data ] = filterDisplacements(U,V,options)
% filter displacement data to reproduce VIC-2D strain filter smoothing
%   'strain window size', specifically filter size, defined in 'data points' (i.e. in terms of the
%   stepsize used for DIC) This definition is similar to defintion of strain window in VIC-2D 
%   input

% created by Tijmen Vermeij & Johan Hoefnagels

%%% create filter
if options.filter
    % OPTION I: Provide the standard deviation of gaussian filter (in data points)
    if options.option == 1
        filt_std = options.filt_std;
        EffectiveWindowSize = 2*sqrt(-2*log(1-options.cutofffraction)*filt_std^2);
        if options.cutofffraction < 1
            window = ceil(EffectiveWindowSize)+1-rem(ceil(EffectiveWindowSize),2); 
            imageFilter=fspecial('gaussian',window ,filt_std);
            imageFilter(imageFilter<(1-options.cutofffraction)*max(imageFilter(:)))=0;
        else
    %         window = 15*ceil(EffectiveWindowSize)+1-rem(ceil(EffectiveWindowSize),2); 
            window = ceil(15*filt_std) + 1 - rem(ceil(15*filt_std),2);
            imageFilter=fspecial('gaussian',window ,filt_std);
        end
    end

    % OPTION 2: Provide the effective strain window size (in data points)
    if options.option == 2
        EffectiveWindowSize = options.EffectiveWindowSize; 
        filt_std = 0.5*sqrt(-0.5*EffectiveWindowSize^2/log(1-options.cutofffraction));
        window = ceil(EffectiveWindowSize)+1-rem(ceil(EffectiveWindowSize),2); 
        imageFilter=fspecial('gaussian',window ,filt_std);
        imageFilter(imageFilter<(1-options.cutofffraction)*max(imageFilter(:)))=0;
    end
else
    imageFilter = 0;
end


%%% perform filtering
  


% blur displacement fields (use nanconv script which can handle nan values well)
if options.filter
    Uf = nanconv(U,imageFilter, 'nonanout','edge');
    Vf = nanconv(V,imageFilter, 'nonanout','edge');
else
    Uf = U;
    Vf = V;
end

data.imageFilter = imageFilter;
data.U = Uf;
data.V = Vf;



end

