function [dlat,rlat,trnm] = spinw2lat(sw,varargin)
% function [dlat,rlat, positions, types] = sw2lat(sw,varargin)
kdef = struct('k',NaN*[0;0;0],'nExt',NaN*[1;1;1]);
[~,kwds]=symbz.parse_arguments(varargin,kdef);

if all(isnan(kwds.k))
    k = sw.mag_str.k;
else
    k = kwds.k;
end
if all(isnan(kwds.nExt))
    nExt = double(sw.mag_str.nExt);
else
    nExt = kwds.nExt;
end
assert(numel(nExt)==3,'the number of unit cell extensions, nExt, must be (1,3) or (3,1)')
nExt = nExt(:);

if numel(k)==3
    % nExt and k should be compatible, with nExt a direct lattice "vector" and
    % k a reciprocal lattice vector
    [kn,kd]=rat(k(:),1e-5);
    % the rationalized denominator of k should (normally) be nExt
    if sum(abs(kd - nExt)) > sum(abs(kd + nExt))*eps()
        warning('k=(%d/%d,%d/%d,%d/%d) and nExt=(%d,%d,%d) are not compatible',...
            kn(1),kd(1),kn(2),kd(2),kn(3),kd(3),...
            nExt(1),nExt(2),nExt(3))
        nExt = max( kd, nExt);
        warning('Returning supercell lattice for nExt=(%d,%d,%d)',nExt(1),nExt(2),nExt(3))
    end
end

% the transformation matrix from units of Q in input lattice to those in
% the returned lattice (which SpinW expects)
trnm = diag( cat(1,nExt,1) ); 

lens = sw.lattice.lat_const(:) .* nExt;
angs = sw.lattice.angle(:); % SpinW stores angles in radian

[dlat,rlat]=symbz.lattice(lens,angs,'radian','direct');

% positions = sw.atom.r;  %(3,nAtoms)
% types = sw.atom.idx(:); %(nAtoms,1)
% if any(nExt > 1)
%     ijk = zeros(3, 1, prod(nExt));
%     l=1;
%     for i=1:nExt(1)
%     for j=1:nExt(2)
%     for k=1:nExt(3)
%         ijk(:,l) = [i;j;k];
%         l = l + 1;
%     end
%     end
%     end
%     positions = reshape( bsxfun(@plus,positions,ijk-1), [3, size(positions,2)*prod(nExt)]);
%     positions = bsxfun(@rdivide,positions,nExt);
%     types = repmat(types,[prod(nExt),1]);
% end