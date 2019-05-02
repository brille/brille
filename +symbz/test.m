function isok = test()
isok=false;

[d,r] = symbz.lattice([2*pi,2*pi,2*pi],[90,90,120],'direct');

if     ~isapprox(r.a, 2*pi/d.a/sin(d.gamma))
    return;
elseif ~isapprox(r.b, 2*pi/d.b/sin(d.gamma))
    return;
elseif ~isapprox(r.c, 2*pi/d.c)
    return;
end

if ~throws_message('A singlepy.symbz._symbz.Reciprocal lattice is required as input',@symbz.brillouinzone,d)
    return;
end
bz = symbz.brillouinzone(r);
if bz.faces.shape{1}~=8 || bz.faces.shape{2}~=3
    return;
elseif bz.faces_per_vertex.shape{1}~=12 || bz.faces_per_vertex.shape{2}~=3
    return;
elseif bz.vertices.shape{1}~=12||bz.vertices.shape{2}~=3
    return;
end

if ~throws_message('A singlepy.symbz._symbz.BrillouinZone is required as input',@symbz.BZGridQ,r)
    return;
end
bzg = symbz.BZGridQ(bz,'N',[5,5,5]);
if bzg.map.shape{1} ~=10 || bzg.map.shape{2} ~=10 || bzg.map.shape{3} ~= 10
    return;
elseif bzg.rlu.shape ~= bzg.invA.shape
    return;
elseif bzg.mapped_rlu.shape ~= bzg.mapped_invA.shape
    return;
end
% more tests?
[~,r]=symbz.lattice(2*pi*[1,1,1],90*[1,1,1]);
bz=symbz.brillouinzone(r);
bzg = symbz.BZGridQ(bz,'N',[10,10,10]);
sq = @(Q)( cat(2, Q(:,1),Q(:,2)+0.5*Q(:,3)) ); % replace with any function linear in the components of Q
bzg.fill( py.numpy.array( sq(double(bzg.mapped_rlu)) ) );
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