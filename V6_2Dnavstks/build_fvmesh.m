function F = build_fvmesh(M,P)
%BUILD_FVMESH  Edge-based finite-volume data for a mixed triangle / quad mesh.
%
%  Every edge has a left cell L, a right cell R and a unit normal pointing from L to R.
%  Boundary edges get a ghost cell as their right cell (index nC + k, k = 1..nB).
%  Boundary edges are classified as wall (etype 1) or far field (etype 2) from the
%  physical-group names / tags of the mesh, from P.wallTags / P.farTags, or, if the
%  mesh carries no tags, from the geometry (edges on the bounding box = far field).

nodes = M.nodes;
nt = size(M.tris,1);
nq = size(M.quads,1);
nC = nt + nq;
if nC == 0
    error('build_fvmesh: the mesh has no cells');
end

C = [M.tris, M.tris(:,3); M.quads];
ncv = [3*ones(nt,1); 4*ones(nq,1)];
isT = ncv == 3;

%% orientation (counter-clockwise) and cell geometry
[a1,a2] = split_areas(nodes,C);
neg = (a1 + a2) < 0;
fT = neg & isT;
fQ = neg & ~isT;
C(fT,[2 3]) = C(fT,[3 2]);
C(fQ,:) = C(fQ,[1 4 3 2]);
C(isT,4) = C(isT,3);
[a1,a2] = split_areas(nodes,C);
Om = a1 + a2;
if any(Om <= 0)
    error('build_fvmesh: %d cells have zero or negative area (run inspect_mesh.py on the mesh)',sum(Om <= 0));
end
x = reshape(nodes(C,1),nC,4);
y = reshape(nodes(C,2),nC,4);
c1x = (x(:,1) + x(:,2) + x(:,3))/3;   c1y = (y(:,1) + y(:,2) + y(:,3))/3;
c2x = (x(:,1) + x(:,3) + x(:,4))/3;   c2y = (y(:,1) + y(:,3) + y(:,4))/3;
xc = [(a1.*c1x + a2.*c2x)./Om, (a1.*c1y + a2.*c2y)./Om];

%% edges
t = find(isT);
q = find(~isT);
E = [C(t,[1 2]); C(t,[2 3]); C(t,[3 1]); ...
     C(q,[1 2]); C(q,[2 3]); C(q,[3 4]); C(q,[4 1])];
cid = [t; t; t; q; q; q; q];

[uk,~,ic] = unique(sort(E,2),'rows');
nE = size(uk,1);

fwd = E(:,1) < E(:,2);
if any(accumarray(ic(fwd),1) > 1) || any(accumarray(ic(~fwd),1) > 1)
    error('build_fvmesh: inconsistent cell orientation or an edge shared by more than two cells');
end
Lrow = zeros(nE,1);
Rrow = zeros(nE,1);
Lrow(ic(fwd))  = find(fwd);
Rrow(ic(~fwd)) = find(~fwd);
nof = Lrow == 0;
Lrow(nof) = Rrow(nof);
Rrow(nof) = 0;

L = cid(Lrow);
R = zeros(nE,1);
R(Rrow > 0) = cid(Rrow(Rrow > 0));

ea = E(Lrow,1);
eb = E(Lrow,2);
tx = nodes(eb,1) - nodes(ea,1);
ty = nodes(eb,2) - nodes(ea,2);
A  = hypot(tx,ty);
nx =  ty./A;                  % outward normal of the left cell (counter-clockwise ordering)
ny = -tx./A;
xm = 0.5*(nodes(ea,:) + nodes(eb,:));

bE = find(R == 0);
nB = numel(bE);
isB = (R == 0);
R(bE) = nC + (1:nB)';

%% boundary classification
tag = zeros(nB,1);
if ~isempty(M.blines)
    key = sort(M.blines(:,1:2),2);
    [tf,loc] = ismember(uk(bE,:),key,'rows');
    tag(tf) = M.blines(loc(tf),3);
end

wt = P.wallTags;
ft = P.farTags;
wallPat = 'wall|body|airfoil|foil|surface|solid|cyl|profile|naca|object';
farPat  = 'far|inlet|outlet|outer|ext|free|inf|domain';
autoW = [];
autoF = [];
for k = 1:numel(M.physTag)
    nm = lower(M.physName{k});
    if ~isempty(regexp(nm,wallPat,'once'))
        autoW(end+1) = M.physTag(k);     %#ok<AGROW>
    elseif ~isempty(regexp(nm,farPat,'once'))
        autoF(end+1) = M.physTag(k);     %#ok<AGROW>
    end
end
if isempty(wt), wt = autoW; end
if isempty(ft), ft = autoF; end

bWall = ismember(tag,wt) & tag ~= 0;
bFar  = ismember(tag,ft) & tag ~= 0;
un = ~(bWall | bFar);
if any(un)
    x0 = min(nodes(:,1));  x1 = max(nodes(:,1));
    y0 = min(nodes(:,2));  y1 = max(nodes(:,2));
    tol = 1e-6*max(x1 - x0,y1 - y0);
    mb = xm(bE,:);
    onb = abs(mb(:,1)-x0) < tol | abs(mb(:,1)-x1) < tol | abs(mb(:,2)-y0) < tol | abs(mb(:,2)-y1) < tol;
    bFar(un & onb)   = true;
    bWall(un & ~onb) = true;
    if any(tag(un) ~= 0)
        fprintf('  note: boundary tags without a recognised name were classified by position\n');
    end
end
etype = zeros(nE,1);
etype(bE(bWall)) = 1;
etype(bE(bFar))  = 2;
if ~any(bWall)
    error('build_fvmesh: no wall boundary found (set P.wallTags to the physical tag of the body)');
end
if ~any(bFar)
    error('build_fvmesh: no far-field boundary found (set P.farTags)');
end

%% ghost cells
Lb = L(bE);
dn = (xm(bE,1) - xc(Lb,1)).*nx(bE) + (xm(bE,2) - xc(Lb,2)).*ny(bE);
xcg = [xc(Lb,1) + 2*dn.*nx(bE), xc(Lb,2) + 2*dn.*ny(bE)];
xcE = [xc; xcg];

%% sparse incidence matrices
ML = sparse(L,(1:nE)',1,nC,nE);
ie = find(~isB);
MR = sparse(R(ie),ie,1,nC,nE);

%% least-squares gradient weights
dx = xcE(R,1) - xcE(L,1);
dy = xcE(R,2) - xcE(L,2);
ds = hypot(dx,dy);
w  = 1./ds.^2;
Sxx = ML*(w.*dx.^2) + MR*(w.*dx.^2);
Sxy = ML*(w.*dx.*dy) + MR*(w.*dx.*dy);
Syy = ML*(w.*dy.^2) + MR*(w.*dy.^2);
det = Sxx.*Syy - Sxy.^2;

%% cell -> face / neighbour tables (padded to 4 entries)
FE = zeros(nC,4);
nT = numel(t);
nQ = numel(q);
if nT > 0
    FE(t,:) = [ic(1:nT), ic(nT+1:2*nT), ic(2*nT+1:3*nT), ic(1:nT)];
end
o = 3*nT;
if nQ > 0
    FE(q,:) = [ic(o+1:o+nQ), ic(o+nQ+1:o+2*nQ), ic(o+2*nQ+1:o+3*nQ), ic(o+3*nQ+1:o+4*nQ)];
end
NB = L(FE) + R(FE) - repmat((1:nC)',1,4);
dxf = reshape(xm(FE(:),1),nC,4) - repmat(xc(:,1),1,4);
dyf = reshape(xm(FE(:),2),nC,4) - repmat(xc(:,2),1,4);

%% pack
F.nC = nC;   F.nE = nE;   F.nB = nB;
F.nodes = nodes;   F.cells = C;   F.ncv = ncv;   F.isTri = isT;
F.xc = xc;   F.Om = Om;   F.h = sqrt(Om);
F.L = L;   F.R = R;   F.nx = nx;   F.ny = ny;   F.A = A;   F.xm = xm;
F.isB = isB;   F.etype = etype;   F.bE = bE;   F.bL = Lb;
F.bWall = find(bWall);   F.bFar = find(bFar);
F.wallE = bE(bWall);
F.ML = ML;   F.MR = MR;
F.wdx = w.*dx;   F.wdy = w.*dy;
F.Mi11 =  Syy./det;   F.Mi22 = Sxx./det;   F.Mi12 = -Sxy./det;
F.ex = dx./ds;   F.ey = dy./ds;   F.ds = ds;
F.dLx = xm(:,1) - xc(L,1);       F.dLy = xm(:,2) - xc(L,2);
F.dRx = xm(:,1) - xcE(R,1);      F.dRy = xm(:,2) - xcE(R,2);
F.FE = FE;   F.NB = NB;   F.dxf = dxf;   F.dyf = dyf;

% reference length and moment reference point from the wall nodes
wn = unique(uk(F.wallE,:));
wx = nodes(wn,1);
wy = nodes(wn,2);
if ~isempty(M.Lref)
    Lref = M.Lref;
else
    Lref = max(wx) - min(wx);
end
if ~isempty(P.Lref)
    Lref = P.Lref;
end
[xle,ile] = min(wx);
F.Lref = Lref;
F.xLE = [xle, wy(ile)];
if isempty(P.momentRef)
    F.momentRef = [xle + 0.25*Lref, wy(ile)];
else
    F.momentRef = P.momentRef;
end
end


function [a1,a2] = split_areas(nodes,C)
nC = size(C,1);
x = reshape(nodes(C,1),nC,4);
y = reshape(nodes(C,2),nC,4);
a1 = 0.5*((x(:,2)-x(:,1)).*(y(:,3)-y(:,1)) - (x(:,3)-x(:,1)).*(y(:,2)-y(:,1)));
a2 = 0.5*((x(:,3)-x(:,1)).*(y(:,4)-y(:,1)) - (x(:,4)-x(:,1)).*(y(:,3)-y(:,1)));
end
