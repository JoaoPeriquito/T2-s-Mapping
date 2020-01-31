% Pixelwise T2(*) fitting using MATLAB fitting toolbox 
% authors: Ludger Starke and Joao Periquito
% INPUT:
%   data           - T2(*) images with different T2(*) weighting
%   TEs            - a vector with the echo times for each image

% OUTPUT:
%   T2smap         - T2(*) map
%   aMap           - M0 map
%   bMap           - R2(*) map
%   r2sMap         - goodness of the fit (r squared)

%eg. [T2map, aMap, bMap, r2sMap] = T2analysis(T2_data,TE_T2); or
%    [T2map, aMap, bMap, r2sMap] = T2analysis(T2s_data,TE_T2s);


function [T2map, aMap, bMap, r2Map] = T2analysis(data, TEs)

Max = max(data(:));
data = data/Max;

[TEs, I] = sort(TEs);
data = data(:,:,I);

dim = size(data);

if numel(dim) ~= 3
    
    error('Wrong dimension\n')

end


%% Model

T2model = fittype( 'a*exp(-x*b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions(T2model);
opts.Algorithm = 'Trust-Region';
% opts.Algorithm = 'Levenberg-Marquardt';
opts.Lower = [0, 0];
opts.Upper = [100*max(data(:)), 10]; 
opts.TolFun = 10^(-4);
opts.TolX = 10^(-4);

aMap = zeros(dim(1), dim(2));
bMap = zeros(dim(1), dim(2));
r2Map = zeros(dim(1), dim(2));

counter = 1;
fprintf('\n')
%% Looping through pixels
for ii = 1:dim(1)
    
    if (ii/dim(1) > counter/20)
        fprintf('%d%% ', counter*5)
        counter = counter + 1;
    end

    for jj = 1:dim(2)
        
%       fprintf('%d - %d : ', ii, jj)

%% Starting values calculation        
        SI = squeeze(data(ii, jj, :));
        startA = SI(1);
        
        [~, minIndex] = min(abs(SI/max(SI) - 0.63));
        startB = 1/TEs(minIndex);
%% Fitting        
        opts.StartPoint = [startA, startB];
        [fitResult, goF] = fit(TEs(:), SI(:), T2model, opts);
%% Results        
        aMap(ii, jj) = fitResult.a;
        bMap(ii, jj) = fitResult.b;
        r2Map(ii, jj) = goF.rsquare;
       
    end
end

fprintf('100%%\n')

T2map = bMap.^(-1);
aMap = aMap*Max;