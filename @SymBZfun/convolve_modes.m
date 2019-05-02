function convS = convolve_modes(obj,qh,qk,ql,en,omega,S,varargin)
inpForm.fname  = {'resfun' 'pars'};
inpForm.defval = {'gauss'   1    };
inpForm.size   = {[1 -2]   [1 -1]};
inpForm.soft   = {false    false};

warnState = warning('off','sw_readparam:UnreadInput');
kwds = sw_readparam(inpForm, varargin{:});

pars = kwds.pars;
% Determine which resolution function to use.
scale_factor = 1;
if ischar(kwds.resfun)
    switch lower(kwds.resfun)
        case {'g','gauss','gaussian'}
            fwhm = pars(1);
            if numel(pars) > 1; scale_factor = pars(2); end
            resfun = @(Emat,center)gauss_internal(Emat,center,@(x)(fwhm+x*0));
        case {'l','lor','lorentz','lorentzian'}
            resfun = @(emat, cen) lorz_internal(emat, cen, pars);
            if numel(pars) > 1; scale_factor = pars(2); end
        case {'v','voi','voigt'}
            resfun = @(emat, cen) voigt_internal(emat, cen, pars);
            if numel(pars) > 2; scale_factor = pars(3); end
        case {'s','sho','simple harmonic oscillator'}
            convS = sho_internal(en, omega, S, pars);
            return
        otherwise
            error('convolve_modes:UnknowResFun','Unknown resolution function %s', kwds.resfun);
    end
else
    fwhm = kwds.resfun;
    if numel(pars) > 0; scale_factor = pars(1); end
    resfun = @(Emat,center)gauss_internal(Emat,center,@(x)(fwhm+x*0));
end

res = zeros(size(omega,1),1);
for i=1:size(omega,2)
    res = res + S(:,i).*resfun( en, omega(:,i));
end
convS = scale_factor * sum( res, 2); % sum over modes

end



%--------------------------------------------------------------------------------------------------
function G = gauss_internal(Emat,center,FWHMfun)
% 1D Gauss function
%
% G = gauss_internal(Emat,center,FWHMfun)
%
% Input:
%
% Emat      Matrix with energy bin values, dimensions of [nQ nE].
% center    Column vector of the center of the Gaussians with nQ elements.
% FWHMfun   Function handle: dE = resfun(E), works on vectors.
%
% Output:
%
% G         Output matrix, with Gaussians, dimensions are [nQ nE].

% resolution for every omega values
sig = repmat(FWHMfun(center)/sqrt(log(256)),[1 size(Emat,2)]);
% intensities
G = exp(-(bsxfun(@minus,center,Emat)).^2./(2*sig.^2))./(sig*sqrt(2*pi));

end

function out = lorz_internal(Emat, center, fwhm)
% Calculates a Lorentzian function for disp2sqw.
    % fwhm should be scalar.
    out = abs(fwhm(1)/pi) ./ (bsxfun(@minus, center, Emat).^2 + fwhm(1)^2);
end

function out = voigt_internal(Emat, center, pars)
% Calculates a pseudo-Voigt function for disp2sqw.
    fwhm = pars(1);
    lorfrac = pars(2);
    sig = fwhm / sqrt(log(256));
    Ediff2 = bsxfun(@minus, center, Emat).^2;
    out = (abs(fwhm/pi) ./ (Ediff2 + fwhm^2)) .* lorfrac + ...
          (exp(-Ediff2 ./ (2*sig^2)) ./ (sig*sqrt(2*pi))) .* (1 - lorfrac);
end

function weight = sho_internal(en, om, sf, pars)
% Calculates weight from SHO model
    gam = 0.1;
    temp = 0;
    amp = 1;
    if numel(pars)>0; gam  = pars(1); end
    if numel(pars)>1; temp = pars(2); end
    if numel(pars)>2; amp  = pars(3); end
    
    nQ = numel(en);
    nM = size(om,2);
    Bose = ones(nQ,1);
    if abs(temp)>sqrt(eps) 
        Bose = en./ (1-exp(-11.602.*en./temp));
    end

    % Use damped SHO model to give intensity:
    enM2 = repmat(en.^2,[1,nM]);
    ok = om > 0 & isfinite(om);
    sho = zeros(nQ,nM);
    sho(ok) = 4*gam*om(ok).*sf(ok) ./ (pi*(enM2(ok)-om(ok).^2).^2 + 4*gam^2*enM2(ok));
    weight = sum(sho,2);
    weight = amp .* Bose .* weight;
end
