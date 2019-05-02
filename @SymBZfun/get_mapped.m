function cellofvecs = get_mapped(obj)
% The python symbz mapped functions return either (N,3) or (N,4) arrays of
% (qh,qk,ql) or (qh,qk,ql,en).
if obj.rluNeeded
    QorQE = symbz.p2m(obj.BZGrid.mapped_rlu);
else
    QorQE = symbz.p2m(obj.BZGrid.mapped_invA);
end
% arrayfun splits the (N,3 or 4) array into a cell(1,3 or 4) of vectors
cellofvecs = arrayfun(@(i){QorQE(:,i)},1:size(QorQE,2));
end