function intres = interpolate(obj,qh,qk,ql,en)
% The python module expects an (N,3)
iat = cat(2,qh,qk,ql);
% or (N,4), if isQE is true
if obj.isQE 
   iat = cat(2,iat,en);
end
% numpy.array as input to the interpolator
iat = symbz.m2p(iat);

num = numel(qh);
numres = num * sum(cellfun(@prod,obj.shape));

% Do the actual interpolation
allres = symbz.p2m( obj.BZGrid.interpolate_at(iat) );
assert( numel(allres) == numres )
% and then split-up the interpolated results into the expected outputs
offsets = cumsum( cat(2, 0, cellfun(@prod,obj.shape)) );
intres = cell(1,obj.nFill);
for i=1:obj.nFill
    intres{i} = reshape( allres(:, (offsets(i)+1):offsets(i+1) ), cat(2,num,obj.shape{i}) );
end

end