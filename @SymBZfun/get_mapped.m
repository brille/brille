function cellofvecs = get_mapped(obj)
% The python symbz mapped functions return either (N,3) or (N,4) arrays of
% (qh,qk,ql) or (qh,qk,ql,en).
if obj.rluNeeded
    QorQE = symbz.p2m(obj.pygrid.rlu);
else
    QorQE = symbz.p2m(obj.pygrid.invA);
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
