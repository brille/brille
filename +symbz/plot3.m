function plot3(BrillouinZone_or_BZGrid,varargin)
    defs=struct('facecolor','none','facealpha',1,'edgecolor','k','edgealpha',1,'showgrid',true,'fullgrid',false);
    [~,kwds]=symbz.parse_arguments(varargin,defs,{'showgrid','fullgrid'});
    ph = ishold();
    hold on;
    
    bztype = 'py.symbz._symbz.BrillouinZone';
    bzgtype = 'py.symbz._symbz.BZGridQ';
    if isa(BrillouinZone_or_BZGrid,bztype)
        isgrid = false;
        bz = BrillouinZone_or_BZGrid;
    elseif isa(BrillouinZone_or_BZGrid,bzgtype)
        isgrid = true;
        bz = BrillouinZone_or_BZGrid.brillouinzone;
    end
    assert( exist('bz','var')==1, 'The first input must be either a BrillouinZone or BZGridQ object');
    
    % pull out MATLAB versions of python numpy arrays
    faces = double(bz.faces_invA)/2;
    verts = double(bz.vertices_invA);
    f_p_v = int32(bz.faces_per_vertex) +1 ; % +1 to convert from C to MATLAB indexing
    
    for i=1:size(faces,1)
        this_f_v = any( f_p_v == i, 2);
        corners = unique( verts(this_f_v,:), 'rows' );
        % to draw a patch, we need to sort the corners
        % we want dot( cross( corner(i)-dot(corner(i),faces(i)), corner(i+1)-dot(corner(i+1),faces(i)) ) , faces(i)) > 0 
        % for all corners
        n = faces(i,:) ./ norm(faces(i,:));
        % find the vectors from the face-centre to the corners
        perpc = corners - faces(i,:); 
        % if the cross product between each perpc is always positive
        % then we've fulfilled our goal
        nc = size(perpc,1);
        perm = 1:nc;
        for j=1:nc-1
            for k=j+1:nc
                if dot( cross( perpc(perm(j),:), perpc(perm(k),:)), n ) < 0
                    tmp = perm(j);
                    perm(j) = perm(k);
                    perm(k) = tmp;
                end
            end
        end
        patch('faces',perm,'Vertices',corners,...
              'facecolor',kwds.facecolor,'facealpha',kwds.facealpha,...
              'edgecolor',kwds.edgecolor,'edgealpha',kwds.edgealpha );
    end
    
    if isgrid && kwds.showgrid
        if kwds.fullgrid
            grid_points = double(BrillouinZone_or_BZGrid.invA);
        else
            grid_points = double(BrillouinZone_or_BZGrid.mapped_invA);
        end
        plotpoints3(grid_points,[],[],{'Q_x','Q_y','Q_z'});
    end
    
    view(3)
    
    if ~ph
        hold off;
    end
end
