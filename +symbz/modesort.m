function [omega,Sab] = modesort(omega,Sab)
nP = size(omega,1);
nM = size(omega,2);
sz = size(Sab);
assert( ismatrix(omega) && sz(1)==nP && sz(end)==nM);

Sp = 1:ndims(Sab);
Sp1 = circshift(Sp,1);
Sp1([1,2])=Sp1([2,1]);
Sab1 = permute( Sab, Sp1 );

for i=1:nP
    [omega(i,:),p] = sort( real(omega(i,:)) );
    Sab1(i,:,:) = Sab1(i,p,:);
end

Spm1 = circshift(Sp,-1);
Spm1([1,end])=Spm1([end,1]);
Sab = permute( Sab1, Spm1);