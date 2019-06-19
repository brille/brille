function wq = w_q(obj,qh,qk,ql,varargin)
% Input:
% ------
%   qh,qk,ql   Arrays containing points at which to evaluate omega(q)
%
%   kwds       A series of 'keywords' and parameters.
%
%              Some keywords control aspects of this function:
%              'coordtrans' - a matrix to transform the input coordinates
%                             (qh,qk,ql,en) before being sent to the 
%                             py.symbz.simphony.SymSim object's method. 
%                             [default: eye(4) % identity]
%
%              Any additional keyword parameters will be passed to SymSim
%              as a py.dict for processing. 
%
% Output:
% -------
%   w(q)       Array with eigen energies at the Q points
%              [ size(w_q) == size(qh) ]
kdef = struct('coordtrans',eye(4));
[args,kwds] = symbz.parse_arguments(varargin,kdef);


nQ = numel(qh);
inshape = size(qh);
if size(qh,1) ~= nQ
    qh = qh(:);
    qk = qk(:);
    ql = ql(:);
end
% Transforms input coordinates if needed
if sum(sum(abs(kwds.coordtrans - eye(4)))) > 0
    qc = [qh qk ql 0*qh];
    qh = sum(bsxfun(@times, kwds.coordtrans(1,:), qc),2);
    qk = sum(bsxfun(@times, kwds.coordtrans(2,:), qc),2);
    ql = sum(bsxfun(@times, kwds.coordtrans(3,:), qc),2);
    clear qc;
end
Q = symbz.m2p(cat(2,qh,qk,ql));

if numel(args)>1
    dict = struct(args{:});
else
    dict = struct();
end
wq = symbz.p2m( obj.pyobj.w_q(Q, py.dict(dict)).magnitude );

wq_size = size(wq);
if wq_size(1:end-1) ~= inshape
    wq = reshape(wq, [inshape, wq_size(end)]);
end
end % horace_sqw