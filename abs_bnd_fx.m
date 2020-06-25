function [alh] = abs_bnd_fx(wl,cp)
% Function to estimate ALH using hyperspectral particulate beam-attenuation
% (cp) data. This code linearly interpolates the (cp) spectra onto a
% standard wavelength vector, and derives cp residuals (cpr) by subtracting
% a power-law fit. The cpr vectors are projected onto the first three 
% eigenvectors derived by eigendecomposition of a subsampled portion of the
% Tara Oceans Expedition cp dataset. Predictors based on the positive and 
% negative values of the projection are applied in a multiple linear
% regression to estimate chlorophyll a absortion line height (alh).
% 
% ** This code is a beta version and has not been thoroughly tested.
% 
% ** Regional differences may exist in the relationships that were used to 
%    construct this code. Correspondence describing your dataset and this
%    code's performance with your data is extremely appreciated.
%    (correspondance address provided below).
% 
% ** Script does not remove filter artifacts. See Chase et al. 2013, and
%    visit http://misclab.umeoce.maine.edu/software.php for corresponding
%    code to remove filter artifacts associated with WetLabs ac-s.
% 
% ** This code is valid for cp spectra measured with a WetLabs ac-s and 
%    which correspond to a consistent nominal wavelength set.
% 
% INPUTS:
%       [wl]:   1 x N vector [nm];      wavelength
%       [cp]:   M x N matrix [/m];      particulate beam attenuation
% 
% OUTPUTS:
%       [alh]:  M x 1 vector [/m];      absorption line height
%       
%       
%  Comments to help improve or optimize this code, as well as feedback on
%  performance of this code (using any relevant dataset) would be extremely
%  welcome and may be directed with subject line "Absorption Band Effects"
%  to <hhouskee@ucsc.edu>.
% 
% 
%  Best wishes,
% 
%  Henry Houskeeper
%  hhouskee@ucsc.edu
%  Ocean Sciences Department
%  University of California, Santa Cruz
% 
%  VERSION 1; last updated 2020 June 25

fprintf('\n%s\n%s\n%s\n%s\n\n',...
    '*******************************************************', ...
    '        Warning, this code is a beta version', ...
    '    and has not been tested on multiple datasets.', ...
    '*******************************************************')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                  Eigenvectors (derived June 25 2020)                    %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

PSI = [-0.1460362 -0.0994965 -0.2034638
-0.2195018 -0.1299693 -0.330608
-0.3226411 -0.1744869 -0.3635149
-0.4138938 -0.1852109 -0.206234
-0.4669648 -0.120849 0.02655046
-0.4640913 0.002939042 0.2185824
-0.3787169 0.1896811 0.3257801
-0.2471053 0.3438464 0.2603245
-0.1228482 0.4334627 0.1421026
-0.03431109 0.4477047 -0.04362244
0.006117889 0.4070195 -0.2417868
-0.0008747435 0.3367003 -0.389474
-0.02924991 0.268842 -0.4690348];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                  Coefficients (derived June 25 2020)                    %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

f0 = 0.0002657;
f1 = 0.7338176;
f2 = -0.2612777;
f3 = 0.6789917;
f4 = -0.8034422;
f5 = -0.7563318;
f6 = 2.2358288;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%               Standard wavelength vector (June 25 2020)                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

WL = [399.1	403.7	408.2	412.5	416.9	421.7	426.9	431.5	435.7...
    440.5	444.9	450.0	454.8	459.4	464.0	468.7	473.4	478.3...
    483.8   488.1	492.7	497.2	501.9	506.8	511.6	516.6	521.6...
    526.3   531.0	535.3	539.9	544.4	549.2	553.8	558.4	562.8...
    567.3   571.2	575.4	579.5	583.9	588.3	592.7	597.2	601.8...
    606.5   611.2	616.0	620.7	625.2	629.7	634.2	638.7	643.3...
    647.8   652.7	657.3	661.9	666.5	670.9	675.7	679.9	684.2...
    688.4   692.3	696.5	700.2	704.1	707.9	711.5	715.2	718.9];

[mwl,nwl] = size(wl);
if mwl > nwl
    wl = wl';
    [mwl,nwl] = size(wl);
end

[M,N] = size(cp);
if M == N
    fprintf('Warning, square cp matrix, we assume spectra are entered as rows.\n')
elseif M == nwl
    cp = cp';
    [M,N] = size(cp);
end

cp_stand = interp1(wl,cp',WL,'linear')';

cpr = NaN(M,N);
flag = NaN(M,1);
alh = NaN(M,1);

ifit = find( (WL >= 500 & WL <= 640) | (WL >=710 & WL <= 720));
isub = find( WL >= 647.8 & WL <= 700.2 );

% % setting options for fminsearch
opts = optimset('fminsearch');      
opts = optimset(opts,'MaxIter',4000); 
opts = optimset(opts,'MaxFunEvals',2000);
opts = optimset(opts,'TolFun',1e-9);
opts = optimset(opts,'Display','off');
b0=[0.1 1];

for k = 1:M

    i = find(isnan(cp_stand(k,ifit)) == 0);
    [B,FVAL,flag(k,1)] = fminsearch(@least_squares_AD,b0,opts,cp_stand(k,ifit(i)),WL(ifit(i)));
    cpr(k,:) = cp_stand(k,:) - ( B(1)*(532./WL).^B(2) );
    
    if flag(k,1) == 0
        fprintf('Warning: Iterations exceeded. Observation %i set to NaN.\n',k)
    elseif flag(k,1) == -1
        fprintf('Warning: Residual not derived. Observation %i set to NaN.\n',k)
    end
    
end

PCTS = cpr(:,isub) * PSI;

for k = 1:M

    if PCTS(k,1) >= 0
        P(1,k) = PCTS(k,1); P(2,k) = 0;
    else
        P(1,k) = 0; P(2,k) = abs(PCTS(k,1));
    end
  
    if PCTS(k,2) >= 0
        P(3,k) = PCTS(k,2); P(4,k) = 0;
    else
        P(3,k) = 0; P(4,k) = abs(PCTS(k,2));
    end
    
    if PCTS(k,3) >= 0
        P(5,k) = PCTS(k,3); P(6,k) = 0;
    else
        P(5,k) = 0; P(6,k) = abs(PCTS(k,3));
    end
    
    alh(k,1) = f0 + (f1 .* P(1,k)) + (f2 .* P(2,k)) + (f3 .* P(3,k))...
         + (f4 .* P(4,k)) + (f5 .* P(5,k)) + (f6 .* P(6,k));
     
    if alh(k) <= 0
        alh(k) = NaN;
        fprintf('Warning: Negative alh derived. Observation %i set to NaN.\n',k)
    end
    
end

alh(flag<1) = NaN;

return

function y = least_squares_AD(x0,spec,l);
    y=sum((spec-x0(1).*(532./l).^x0(2)).^2);
return

