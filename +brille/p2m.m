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

function m = p2m(p)
if isa(p,'py.complex') % a scalar
    m = py.numpy.array(p).tolist(); % MATLAB converts the one element list to a complex number automatically
elseif contains(class(p),'uint','IgnoreCase',true)
    m = uint64(p);
elseif contains(class(p),'int','IgnoreCase',true)
    m = int64(p);
else
    if verLessThan('matlab','9.4') % before 2018a(?) For sure by 2018b=9.5
        warning('brille:p2m','Fast conversion of numpy.ndarrays not supported by this version of MATLAB. Consider upgrading.');
        ndim = int64(p.ndim);
        nmel = int64(p.size);
        if ndim>1
            toshape = zeros(1, ndim);
            for i=1:ndim
                toshape(i) = int64(p.shape{i});
            end
            p = p.reshape(nmel);
        else
            toshape = [1,nmel];
        end
        m = zeros(toshape);
        if contains(string(p.dtype.name),'complex')
            m = complex(m,0);
        end
        plist = p.tolist();
        for i=1:numel(m)
            m(i) = plist{i};
        end
    else
        eltype = lower(string(p.dtype.name));
        if contains(eltype,'complex')
            rp = py.numpy.array(p.real);
            ip = py.numpy.array(p.imag);
            eltype = lower(string(rp.dtype.name));
            if contains(eltype,'uint')
                rp = uint64(rp);
                ip = uint64(ip);
            elseif contains(eltype,'int')
                rp = int64(rp);
                ip = int64(ip);
            else
                rp = double(rp);
                ip = double(ip);
            end
            m = complex( rp, ip );
        elseif contains(eltype,'uint')
            m = uint64( p );
        elseif contains(eltype,'int')
            m = int64(p);
        elseif contains(eltype,'str')
            ndim = int64(p.ndim);
            nmel = int64(p.size);
            if ndim>1
                toshape = zeros(1,ndim);
                for i=1:ndim
                    toshape(i) = int64(p.shape{i});
                end
                p = p.reshape(nmel);
            else
                toshape = [1,nmel];
            end
            m = cell(toshape);
            plist = p.tolist();
            for i=1:numel(m)
                m{i} = char(plist{i});
            end
        else
            m = double(p);
        end
    end
end
end
