function sqw = horace_sqw(obj,qh,qk,ql,en,pars,varargin) % split varagin into fill varargin, and interpreter varargin?
% Input:
% ------
%   qh,qk,ql,en Arrays containing points at which to evaluate sqw from the
%               broadened dispersion
%
%   pars        Arguments needed by the function.
%               - pars = [model_pars scale_factor resolution_pars]
%               - Should be a vector of parameters
%               - The first N parameters relate to the spin wave dispersion
%                 and correspond to spinW matrices in the order defined by
%                 the 'mat' option [N=numel(mat)]
%               - The next M parameters relate to the convolution parameters
%                 corresponding to the convolution function defined by the
%                 'resfun' option (either one or two parameters depending
%                 on function type.
%               - The last parameter is a scale factor for the intensity
%                 If this is omitted, a scale factor of 1 is used;
%
%   kwpars      - A series of 'keywords' and parameters. Specific to this
%                 function is:
%
%               - 'resfun' - determines the convolution / resolution
%                    function to get S(q,w). It can be either a string:
%                      'gauss' - gaussian with single fixed (fittable) FWHM
%                      'lor' - lorentzian with single fixed (fittable) FWHM
%                      'voigt' - pseudo-voigt with single fixed (fittable) FWHM
%                      @fun - a function handle satisfying the requirements of
%                             the 'fwhm' parameter of disp2sqw.
%                    NB. For 'gauss' and 'lor' only one fwhm parameter may be
%                    specified. For 'voigt', fwhm = [width lorz_frac]
%                    contains two parameters - the fwhm and lorentzian fraction
%                    [default: 'gauss']
%               - 'partrans' - a function to transform the fit parameters
%                    This transformation will be applied before each iteration
%                    and the transformed input parameter vector passed to
%                    spinW and the convolution function.
%                    [default: @(y)y  % identity operation]
%               - 'coordtrans' - a matrix to transform the input coordinates
%                    (qh,qk,ql,en) before being sent to SpinW.
%                    [default: eye(4) % identity]
%
%               In addition, the following parameters are used by this function
%                    and will also be passed on to spinw.matparser which will
%                    do the actual modification of spinW model parameters:
%
%               - 'mat' - A cell array of labels of spinW named 'matrix' or
%                    matrix elements. E.g. {'J1', 'J2', 'D(3,3)'}. These will
%                    be the model parameters to be varied in a fit, their
%                    order in this cell array will be the same as in the
%                    fit parameters vector.
%                    [default: [] % empty matrix - no model parameters]
%
%                 All other parameters will be passed to spinW. See the help
%                    for spinw/spinwave, spinw/matparser and spinw/sw_neutron
%                    for more information.
%
%   swobj       The spinwave object which defines the magnetic system to be
%               calculated.
%
% Output:
% -------
%   weight      Array with spectral weight at the q,e points
%               If q and en given:  weight is an nq x ne array, where nq
%                                   is the number of q points, and ne the
%                                   number of energy points
%               If qw given together: weight has the same size and dimensions
%                                     as q{1} i.e. qh
%

inpForm.fname  = {'partrans' 'coordtrans' 'mat'};
inpForm.defval = {@(y)y      eye(4)       []};
inpForm.size   = {[1 1]      [4 4]        [1 -1]};
inpForm.soft   = {false      false        false};

warnState = warning('off','sw_readparam:UnreadInput');
kwds = sw_readparam(inpForm, varargin{:});

% Prune varargin for passing to disp2sqw or horace:
passon = varargin;
for key = inpForm.fname
    keys = find(cellfun(@ischar, passon)); % possible keys in varargin
    if ~isempty(keys)
        iskey = cellfun(@(x) contains(key,x), passon(keys));
        keyloc = keys(iskey);
        if ~isempty(keyloc)
            passon = [passon(1:keyloc-1) passon(keyloc+2:end)];
        end
    end
end


% Sets the number of spinW model parameters. All others are interpreter pars.
n_horace_pars = numel(kwds.mat);
if isempty(n_horace_pars)
    n_horace_pars = 0;
end
fillerinpt = [{pars(1:n_horace_pars)} passon];
interpinpt = [{'pars'} {pars((1+n_horace_pars):end)} passon];

% The object holds one of:
%   one py.fibril._fibril.BZGridQ[,complex]
%   one py.fibril._fibril.BZGridQE[,complex]
%   a cell array of identical (but memory distinct) BZGridQs
%   a cell array of identical (but memory distinct) BZGridQEs
% plus a hash of the last parameters used, and functions to fill the
% grid(s) and convert interpolated values from the function outputs to
% what horace expects [a single intensity per (qh(i),qk(i),ql(i),en(i))
% tuple].
% This function needs to compare the inputs in varargin to the
% parameter hash, fill the grid(s) if the hash is different,
% interpolate the gridded Rank-N tensors ( 0<=N ) for the input Q [or
% (Q,E)] points, and then call the interpreter function(s) to convert
% the interpolated values into intensities for use by Horace/Tobyfit.

% First, refill the grid(s) if the hashes differ
obj.fill(fillerinpt{:});

% We have one or more filled BZGrids. We want to interpolate their
% stored values at the passed (qh,qk,ql) or (qh,qk,ql,en) points
assert( numel(qh)==numel(qk) && numel(qk)==numel(en)...
     && numel(qh)==numel(en),    'Expected matching numel arrays');
assert( all(size(qh)==size(qk)) && all(size(ql)==size(en)) ...
     && all(size(qh)==size(en)), 'Expected matching shaped arrays');
reshaped = false;
inshaped = size(qh);
if numel(qh) ~= size(qh,1)
    reshaped = true;
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

intres=obj.interpolate(qh,qk,ql,en);

% and then use the interpreter function(s) to convert this to S(Q,E)
for i=1:obj.nInt
    newintres = cell(1,obj.nRet(i));
    [newintres{:}] = obj.interpreter{i}(qh,qk,ql,en,intres{:},interpinpt{:});
    intres = newintres;
end
% If no one has changed the interpreter cells then intres is cell(1,1).
% Check to be sure
assert( numel(intres) == 1 );
sqw = intres{1};

if reshaped
    sqw = reshape(sqw,inshaped);
end

end
