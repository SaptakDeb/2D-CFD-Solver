function m = mesh_assemble(X, Y, id, isWall, periodic)
% MESH_ASSEMBLE  Turn a structured node array into the face-based mesh
% data that the solver uses.
%   X, Y      node coordinates, (NI nodes) x (NJ nodes)
%   id        node numbering; nodes that are physically the same point
%             share an id (this is how a C-grid wake cut is closed)
%   isWall    true for nodes lying on the body
%   periodic  true if the last i-cell wraps around to the first node (O-grid)
%
% Faces are found from the cell->node table, so the solver does not care
% whether the mesh came from an O-grid, a C-grid, or something else.

[NIn,NJn] = size(X);
NIc = NIn - ~periodic;
NJc = NJn - 1;

[uid,~,ic] = unique(id(:));
idc = reshape(ic,size(id));
Nn  = numel(uid);
xn  = zeros(Nn,1);  yn = zeros(Nn,1);  nodeWall = false(Nn,1);
xn(idc(:)) = X(:);  yn(idc(:)) = Y(:);  nodeWall(idc(:)) = isWall(:);

% Cells: quads (i,j) -> nodes (i,j) (i+1,j) (i+1,j+1) (i,j+1)
[I,J] = ndgrid(1:NIc,1:NJc);
I2 = I + 1;
if periodic, I2(I2 > NIn) = 1; end
lin = @(a,b) sub2ind([NIn NJn],a,b);
cells = [idc(lin(I,J))(:), idc(lin(I2,J))(:), idc(lin(I2,J+1))(:), idc(lin(I,J+1))(:)];
Nc = size(cells,1);

% Cell geometry
xq = xn(cells);  yq = yn(cells);
xr = circshift(xq,-1,2);  yr = circshift(yq,-1,2);
cr = xq.*yr - xr.*yq;
A  = 0.5*sum(cr,2);
if any(A > 0) && any(A < 0)
    error('mesh_assemble: folded cells (mixed orientation). Check clustering / domain size.');
end
vol = abs(A);
xc  = sum((xq + xr).*cr,2)./(6*A);
yc  = sum((yq + yr).*cr,2)./(6*A);

% Edges and face lists
e1 = [cells(:,1); cells(:,2); cells(:,3); cells(:,4)];
e2 = [cells(:,2); cells(:,3); cells(:,4); cells(:,1)];
ec = repmat((1:Nc)',4,1);
ea = min(e1,e2);  eb = max(e1,e2);
[ks,p] = sort(ea*(Nn+1) + eb);
dup    = ks(1:end-1) == ks(2:end);
k1     = find(dup);
isPair = false(size(ks));  isPair(k1) = true;  isPair(k1+1) = true;
k0     = find(~isPair);

% Interior faces
fA = p(k1);  fB = p(k1+1);
iL = ec(fA);  iR = ec(fB);
xa = xn(ea(fA));  ya = yn(ea(fA));  xb = xn(eb(fA));  yb = yn(eb(fA));
ilen = hypot(xb-xa, yb-ya);
inx  = (yb-ya)./ilen;   iny = -(xb-xa)./ilen;
flip = inx.*(xc(iR)-xc(iL)) + iny.*(yc(iR)-yc(iL)) < 0;
inx(flip) = -inx(flip);  iny(flip) = -iny(flip);
ixm = 0.5*(xa+xb);  iym = 0.5*(ya+yb);

% Neighbours: slot s of a cell is its edge s (linear index = edge index)
nbr = zeros(Nc,4);
nbr(fA) = iR;
nbr(fB) = iL;
nb0 = nbr;
self = repmat((1:Nc)',1,4);
nb0(nbr == 0) = self(nbr == 0);

% Boundary faces
fC = p(k0);
bL = ec(fC);
xa = xn(ea(fC));  ya = yn(ea(fC));  xb = xn(eb(fC));  yb = yn(eb(fC));
blen = hypot(xb-xa, yb-ya);
bnx  = (yb-ya)./blen;   bny = -(xb-xa)./blen;
bxm  = 0.5*(xa+xb);     bym = 0.5*(ya+yb);
flip = bnx.*(bxm-xc(bL)) + bny.*(bym-yc(bL)) < 0;
bnx(flip) = -bnx(flip);  bny(flip) = -bny(flip);
btype = 2*ones(size(bL));                              % 2 = far field
btype(nodeWall(ea(fC)) & nodeWall(eb(fC))) = 1;        % 1 = wall
wl = btype == 1;

% Least-squares gradient coefficients: grad q = sum_k c_k (q_nb(k) - q_c)
cx = zeros(Nc,4);  cy = zeros(Nc,4);
dxk = zeros(Nc,4); dyk = zeros(Nc,4); w2 = zeros(Nc,4);
for k = 1:4
    has = nbr(:,k) > 0;
    dxk(:,k) = (xc(nb0(:,k)) - xc).*has;
    dyk(:,k) = (yc(nb0(:,k)) - yc).*has;
    w2(:,k)  = has./max(dxk(:,k).^2 + dyk(:,k).^2, eps);      % weight 1/d^2
end
Axx = sum(w2.*dxk.^2,2);  Axy = sum(w2.*dxk.*dyk,2);  Ayy = sum(w2.*dyk.^2,2);
det = Axx.*Ayy - Axy.^2;
for k = 1:4
    cx(:,k) = w2(:,k).*( Ayy.*dxk(:,k) - Axy.*dyk(:,k))./det;
    cy(:,k) = w2(:,k).*(-Axy.*dxk(:,k) + Axx.*dyk(:,k))./det;
end

m = struct('xn',xn,'yn',yn,'cells',cells,'xc',xc,'yc',yc,'vol',vol, ...
           'iL',iL,'iR',iR,'inx',inx,'iny',iny,'ilen',ilen,'ixm',ixm,'iym',iym, ...
           'bL',bL,'bnx',bnx,'bny',bny,'blen',blen,'bxm',bxm,'bym',bym,'btype',btype, ...
           'nbr',nbr,'nb0',nb0,'cx',cx,'cy',cy);
m.wx = [xn(ea(fC(wl))), xn(eb(fC(wl)))].';      % wall segments, for drawing
m.wy = [yn(ea(fC(wl))), yn(eb(fC(wl)))].';
end
