% Copyright 2019 Greg Tucker
%
% This file is part of fibril.
%
% fibril is free software: you can redistribute it and/or modify it under the
% terms of the GNU Affero General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% fibril is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.
%
% See the GNU Affero General Public License for more details.
% You should have received a copy of the GNU Affero General Public License
% along with fibril. If not, see <https://www.gnu.org/licenses/>.

function latmat = latmat(varargin)
% use varargin to define angle units
kdef = struct('degree',true,'radian',false);
[args,kwds]=parse_arguments(varargin,kdef,{'degree','radian'});

if numel(args)>=2
    lens = args{1};
    angs = args{2};
else
    latmat = args{1};
    assert( size(latmat,1)==3 && size(latmat,2)==3 && ismatrix(latmat) );
    return
end

assert(numel(lens)>=3 && numel(angs)>=3)
if kwds.degree && ~kwds.radian
    angs = angs / 180 * pi;
end

% spglib takes a matrix with columns defining the lattice vectors of
% the direct lattice. The first vector, a, is along the x axis of this
% orthonormal coordinate system. The second, b, is in the x-y plane; away
% from a by their mutual angle gamma. And the direction of the third is
% defined by alpha and beta.

xhat = [ 1;0;0];
yhat = [cos(angs(3));sin(angs(3));0];

cstar_hat = cross(xhat,yhat);
ccc2 = (cos(angs(2))+cos(angs(1))*cos(angs(3)))^2;
cz2 = ccc2*( 1/cos(angs(2))^2 -1) - (sin(angs(3))*cos(angs(1)))^2;
cz = sqrt(cz2)* cstar_hat/norm(cstar_hat);
zhat = xhat * cos(angs(2)) + yhat * cos(angs(1)) + cz;
zhat = zhat / norm(zhat);

latmat = cat(2,lens(1)*xhat,lens(2)*yhat,lens(3)*zhat);

end
