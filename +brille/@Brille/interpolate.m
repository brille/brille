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

function intres = interpolate(obj,qh,qk,ql,en)
% The python module expects an (N,3)
iat = cat(2,qh,qk,ql);
% or (N,4), if isQE is true
if obj.isQE
   iat = cat(2,iat,en);
end
s2 = size(iat,2);
trn = obj.Qtrans(1:s2,1:s2);
if sum(sum(abs(trn - eye(s2))))>0
    for i = 1:size(iat,1)
        iat(i,:) = permute( trn* permute(iat(i,:),[2,1]), [2,1]);
%         iat(i,:) = iat(i,:)/trn;
    end
end

% numpy.array as input to the interpolator
iat = brille.m2p(iat);

num = numel(qh);
numres = num * sum(cellfun(@prod,obj.shape));

% Do the actual interpolation
pyallres = obj.pygrid.interpolate_at(iat,true,obj.parallel);
allres = brille.p2m( pyallres );
assert( numel(allres) == numres )
% and then split-up the interpolated results into the expected outputs
intres = cell(1,obj.nFill);
if ismatrix(allres)
    offsets = cumsum( cat(2, 0, cellfun(@prod,obj.shape)) );
    for i=1:obj.nFill
        intres{i} = reshape( allres(:, (offsets(i)+1):offsets(i+1) ), cat(2,num,obj.shape{i}) );
    end
elseif ndims(allres)==3
    offsets = cumsum( cat(2, 0, obj.span) );
    for i=1:obj.nFill
        intres{i} = reshape( allres(:, (offsets(i)+1):offsets(i+1), :), cat(2,num,obj.shape{i}) );
    end
end


end
