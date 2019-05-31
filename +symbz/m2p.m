function p = m2p(m,dtype,neversqueeze)
if nargin<3 || isempty(neversqueeze) || ~islogical(neversqueeze)
    neversqueeze = false;
end
if nargin<2 || isempty(dtype) || ~ischar(dtype)
    dtype = 'double'; 
end

if strncmpi(dtype,'complex',min(length(dtype),7)) || any(imag(m(:))) % the only way of testing for complex numbers in MATLAB?!
    if isscalar(m)
        p = py.complex(m);
    else
%         % ! Carefull not to use m(:)' here as ' is the complex conjugate
%         % transpose (you could use m(:).', the matrix transpose)
%         p = py.numpy.array( arrayfun(@py.complex, transpose(m(:)), 'UniformOutput',false) );
%         if ndims(m)>1
%             p=p.reshape( int64(size(m)) );
%         end
        rp = symbz.m2p( real(m) );
        ip = symbz.m2p( imag(m) );
        p = rp + 1i*ip;
    end
else
    switch lower(dtype)
        case {'int','int64','cint'}
            p = int64(m);
        case {'uint','uint64','size_t'}
            p = uint64(m);
        otherwise
            p = m;
    end
    if ~isscalar(m) && ~ischar(m)
        if verLessThan('matlab','9.5')
            p = py.numpy.array(transpose(p(:)));
            if neversqueeze || any( size(m)~=1 & size(m)~=numel(m) )
                p=p.reshape( int64(size(m)) );
            end
        else
            p = py.numpy.array(p);
        end
    end
end     
    
end