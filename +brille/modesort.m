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

function [omega,Sab] = modesort(omega,Sab)
% nP = size(omega,1);
% nM = size(omega,2);
% sz = size(Sab);
% assert( ismatrix(omega) && sz(1)==nP && sz(end)==nM);
%
% Sp = 1:ndims(Sab);
% Sp1 = circshift(Sp,1);
% Sp1([1,2])=Sp1([2,1]);
% Sab1 = permute( Sab, Sp1 );
%
% for i=1:nP
%     [omega(i,:),p] = sort( real(omega(i,:)) );
%     Sab1(i,:,:) = Sab1(i,p,:);
% end
%
% Spm1 = circshift(Sp,-1);
% Spm1([1,end])=Spm1([end,1]);
% Sab = permute( Sab1, Spm1);
