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

function isok = test()
isok=false;

[d,r] = brille.lattice([2*pi,2*pi,2*pi],[90,90,120],'direct');

if     ~isapprox(r.a, 2*pi/d.a/sin(d.gamma))
    return;
elseif ~isapprox(r.b, 2*pi/d.b/sin(d.gamma))
    return;
elseif ~isapprox(r.c, 2*pi/d.c)
    return;
end

if ~throws_message('A single py.brille._brille.Reciprocal lattice is required as input',@brille.brillouinzone,d)
    return;
end
bz = brille.brillouinzone(r);
if bz.faces.shape{1}~=8 || bz.faces.shape{2}~=3
    return;
elseif bz.faces_per_vertex.shape{1}~=12 || bz.faces_per_vertex.shape{2}~=3
    return;
elseif bz.vertices.shape{1}~=12||bz.vertices.shape{2}~=3
    return;
end

if ~throws_message('A single py.brille._brille.BrillouinZone is required as input',@brille.BZGridQ,r)
    return;
end
bzg = brille.BZGridQ(bz,'N',[5,5,5]);
if bzg.map.shape{1} ~=10 || bzg.map.shape{2} ~=10 || bzg.map.shape{3} ~= 10
    return;
elseif bzg.grid_rlu.shape ~= bzg.grid_invA.shape
    return;
elseif bzg.rlu.shape ~= bzg.invA.shape
    return;
end
% more tests?
[~,r]=brille.lattice(2*pi*[1,1,1],90*[1,1,1]);
bz=brille.brillouinzone(r);
bzg = brille.BZGridQ(bz,'N',[10,10,10]);
sq = @(Q)( cat(2, Q(:,1),Q(:,2)+0.5*Q(:,3)) ); % replace with any function linear in the components of Q
bzg.fill( py.numpy.array( sq(double(bzg.rlu)) ) );
qrand = (rand(30,3)-0.5);
intsq = double( bzg.interpolate_at( py.numpy.array(qrand) ) );
if ~all(all(abs(sq(qrand) - intsq) < 8*eps()*abs(sq(qrand) + intsq))) % 8*eps() since the interpolation is the sum of 8 terms
    return;
end

isok=true;
end

function tf=isapprox(a,b)
tf = abs(a-b) < abs(a+b)*eps();
end
function tf=throws_message(m,f,varargin)
    tf = false;
    try
        f(varargin{:});
        return;
    catch prob
        if strncmp(m,prob.message,min(length(m),length(prob.message)))
            tf=true;
        end
    end
end
