function [omega,S] = neutron_spinwave_intensity(obj,qh,qk,ql,en,omega,Sab,varargin)
assert( numel(qh)==numel(qk) && numel(qk)==numel(en)...
     && numel(qh)==numel(en),    'Expected matching numel arrays');
assert( all(size(qh)==size(qk)) && all(size(ql)==size(en)) ...
     && all(size(qh)==size(en)), 'Expected matching shaped arrays');

if numel(qh) ~= size(qh,1)
    qh = qh(:);
    qk = qk(:);
    ql = ql(:);
%     en = en(:);
end
nQ = numel(qh);
Q = cat(2,qh,qk,ql);
if obj.rluNeeded % qh,qk,ql are (or should be) in rlu, but we need absolute
    Bmatrix = double( obj.BZGrid.brillouinzone.lattice.get_B_matrix() );
    Q = permute(mtimesx_mex(Bmatrix,permute(Q,[2,1])),[2,1]); % Q is (nQ,3) but mtimesx needs (3,nQ)
end

assert(size(omega,1)==nQ, 'The number of q points does not match the omega. Is one a py.numpy.array?');
nMode = size(omega,2);
assert( all(size(Sab)==[nQ,3,3,nMode]) );
% Calculate the neutron scattering cross section (for unpolarized neutrons)
% given Q points and their corresponding spin-spin correlation tensors:

% The following method is lifted from SpinW/sw_neutron with very minor modifications

% Make sure we're only looking at the symmetric part of Sab
p = [1,3,2,4];
SabS = (Sab + permute(Sab,p))/2;

% Normalized scattering wavevector in xyz coordinate system.
Q2 = sqrt(sum(Q.^2,2));
Q2(Q2==0) = 1; % avoid dividing by zero for zero-length Q.
Qnorm = bsxfun(@rdivide,Q,Q2);

Ql = repmat(permute(Qnorm,[1 2 3]),[1 1 3]);
Qr = repmat(permute(Qnorm,[1 3 2]),[1 3 1]);

% Perpendicular part of the scattering wavevector. delta(a,b) - hatQa*hatQb
Qfact = repmat(permute(eye(3),[3 1 2])  ,[nQ,1,1])- Ql.*Qr;
Qfact = repmat(permute(Qfact ,[1 2 3 4]),[1 1 1 nMode]);

% Dynamical structure factor for neutron scattering
% S: nQ x nMode.
S = permute(sum(sum(Qfact.*SabS,2),3),[1 4 2 3]);
% end