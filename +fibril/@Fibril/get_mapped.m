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

function cellofvecs = get_mapped(obj)
% The python fibril mapped functions return either (N,3) or (N,4) arrays of
% (qh,qk,ql) or (qh,qk,ql,en).
if obj.rluNeeded
    QorQE = fibril.p2m(obj.pygrid.rlu);
else
    QorQE = fibril.p2m(obj.pygrid.invA);
end
s2 = size(QorQE,2);
trn = obj.Qtrans(1:s2,1:s2);
if sum(sum(abs(trn - eye(s2))))>0
    for i = 1:size(QorQE,1)
        QorQE(i,:) = permute( trn\permute(QorQE(i,:),[2,1]), [2,1]);
%         QorQE(i,:) = QorQE(i,:)*trn;
    end
end
% arrayfun splits the (N,3 or 4) array into a cell(1,3 or 4) of vectors
cellofvecs = arrayfun(@(i){QorQE(:,i)},1:s2);
end
