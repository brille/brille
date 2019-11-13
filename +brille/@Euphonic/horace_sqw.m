% Copyright 2019 Greg Tucker
%
% This file is part of brille.
%
% brille is free software: you can redistribute it and/or modify it under the
% terms of the GNU Affero General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% brille is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.
%
% See the GNU Affero General Public License for more details.
% You should have received a copy of the GNU Affero General Public License
% along with brille. If not, see <https://www.gnu.org/licenses/>.

function sqw = horace_sqw(obj,qh,qk,ql,en,varargin)
% Input:
% ------
%   qh,qk,ql,en Arrays containing points at which to evaluate sqw from the
%               broadened dispersion
%
%   pars        Arguments needed by the function.
%               - Should be a vector of parameters
%                 args = [scale_factor resolution_pars]
%               - The remaining parameters relate the convolution function
%                 defined by the 'resfun' option
%                 (either one or two parameters depending on function type)
%               - The last parameter is a scale factor for the intensity
%
%   kwds       A series of 'keywords' and parameters.
%
%              Some keywords control aspects of this function:
%              'coordtrans' - a matrix to transform the input coordinates
%                             (qh,qk,ql,en) before being sent to the
%                             py.brille.euphonic.FibEu object's method.
%                             [default: eye(4) % identity]
%
%              Any additional keyword parameters will be passed to SymSim
%              as a py.dict for processing.
%              Notably, one should provide at least
%              'resfun' - which determines the broadening function to
%                         convert S_i(Q(w_i)) into S(Q,E).
%                         This should be a string with one of the following
%                         values:
%                 'sho'   - Simple Harmonic Oscillator linewidth, with its
%                           own (optional) inclusion of the Bose thermal
%                           population factor. resolution_pars are
%                           [FWHM, sample temperature]
%                 'gauss' - gaussian with single fixed (fittable) FWHM
%                 'lor'   - lorentzian with single fixed (fittable) FWHM
%                 'voi'   - Faddeeva-based estimate to a true voigt
%                           function, requiring the FWHM of both the
%                           gaussian and lorentzian functions to be
%                           convoluted. resolution_pars are
%                           [gaussian_FWHM, lorentzian_FWHM]
%                         [default: 'gauss']
%
%              Other keyword-value pairs which are expected by SimPhony are
%              'scale'     - A scale factor, in addition to the one already
%                            required as part of args -- probably should
%                            not be used.
%                            [default: 1]
%              'dw_seed'   - Used in an optional Deby-Waller calculation
%                            [default: ?]
%              'dw_grid'   - Used in an optional Deby-Waller calculation
%                            [default: ?]
%              'calc_bose' - Whether the phonon structure factors should be
%                            corrected for the Bose factor. Note: if the
%                            simple harmonic oscillator broadening function
%                            is selected the *broadened lineshapes* will be
%                            corrected for the Bose factor.
%                            [default: true]
%              'T'         - The sample temperature used in the Deby-Waller
%                            and Bose factor corrections. For the simple
%                            harmonic oscillator function, the temperature
%                            provided with its parameters will override
%                            this keyword.
%                            [default: 5 K]
%
% Output:
% -------
%   sqw        Array with spectral weight at the Q,E points
%              [ size(sqw) == size(qh) ]
kdef = struct('coordtrans',eye(4));
[args,kwds] = brille.parse_arguments(varargin,kdef);
assert( numel(args)>0 && ismatrix(args{1}) && isnumeric(args{1}) ...
       ,'Parameter input is required')
input_pars = args{1};


nQ = numel(qh);
inshape = size(qh);
if size(qh,1) ~= nQ
    qh = qh(:);
    qk = qk(:);
    ql = ql(:);
    en = en(:);
end
% Transforms input coordinates if needed
if sum(sum(abs(kwds.coordtrans - eye(4)))) > 0
    qc = [qh qk ql en];
    qh = sum(bsxfun(@times, kwds.coordtrans(1,:), qc),2);
    qk = sum(bsxfun(@times, kwds.coordtrans(2,:), qc),2);
    ql = sum(bsxfun(@times, kwds.coordtrans(3,:), qc),2);
    en = sum(bsxfun(@times, kwds.coordtrans(4,:), qc),2);
    clear qc;
end
Q = brille.m2p(cat(2,qh,qk,ql));

if numel(args)>1
    dict = struct(args{2:end});
else
    dict = struct();
end
if isfield(dict,'scale')
    scale = dict.scale;
else
    scale = 1;
end
if isfield(dict,'T')
    temp = dict.T;
else
    temp = 5;
end
pars = [];
if isfield(dict,'resfun') && ischar(dict.resfun)
    switch lower(dict.resfun)
        case {'g','gauss','gaussian','l','lor','lorentz','lorentzian'}
            pars = input_pars(1);
            if numel(input_pars) > 1; scale = input_pars(2); end
        case {'v','voi','voigt'}
            pars = input_pars(1:2);
            if numel(input_pars) > 2; scale = input_pars(3); end
        case {'s','sho','simple harmonic oscillator'}
            pars = input_pars(1);
            temp = input_pars(2);
            if numel(args) > 2; scale = input_pars(3); end
    end
end
dict.scale = scale;
dict.T = temp;
if ~isempty(pars)
    dict.param = pars;
end
sqw = brille.p2m( obj.pyobj.s_qw(Q, brille.m2p(en), py.dict(dict)) );

if size(sqw) ~= inshape
  sqw = reshape(sqw, inshape);
end
end % horace_sqw
