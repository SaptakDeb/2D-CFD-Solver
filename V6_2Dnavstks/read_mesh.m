function M = read_mesh(P)
%READ_MESH  Read a 2D mesh and return it in one common format.
%
%  P.meshSource = 'file'     read P.meshFile
%                   - '*.msh'                 Gmsh ASCII, format 2.2 or 4.x
%                   - 'name' or 'name_nodes.txt'   text files written by mesh2d.py
%                                             (name_nodes.txt, _tri.txt, _quad.txt, _bnd.txt)
%  P.meshSource = 'builtin'  structured O-grid from mesh_<geometry>.m
%
%  Output struct M
%    nodes   N x 2
%    tris    nt x 3     node indices (rows of nodes)
%    quads   nq x 4
%    blines  nb x 3     [node1 node2 physical-tag] boundary segments (may be empty)
%    physTag, physName  names of the physical groups (cell array of strings), if known
%    Lref    reference length (empty if not known)

M = struct('nodes',zeros(0,2),'tris',zeros(0,3),'quads',zeros(0,4), ...
           'blines',zeros(0,3),'physTag',[],'physName',{{}},'Lref',[]);

if strcmpi(P.meshSource,'file')

    f = P.meshFile;
    [~,~,ext] = fileparts(f);
    if strcmpi(ext,'.msh')
        M = read_msh(f);
    elseif ~isempty(regexp(f,'_nodes\.txt$','once'))
        M = read_txt(regexprep(f,'_nodes\.txt$',''));
    elseif exist([f '_nodes.txt'],'file') == 2
        M = read_txt(f);
    else
        error('read_mesh: cannot interpret meshFile "%s" (use a .msh file or the base name of the mesher text files)',f);
    end

else

    meshFun = str2func(['mesh_' lower(P.geometry)]);
    [X,Y,Lref] = meshFun(P);
    M = grid_to_unstructured(X,Y);
    M.Lref = Lref;

end

end


%% ------------------------------------------------------------------------
function M = grid_to_unstructured(X,Y)
% structured O-grid (i around, periodic; j wall -> far field) -> quads
Ni = size(X,1) - 1;
Nj = size(X,2) - 1;

nodes = [reshape(X(1:Ni,:),[],1), reshape(Y(1:Ni,:),[],1)];
id = @(i,j) (j-1)*Ni + i;

[I,J] = ndgrid(1:Ni,1:Nj);
I = I(:);  J = J(:);
Ip = mod(I,Ni) + 1;
quads = [id(I,J), id(I,J+1), id(Ip,J+1), id(Ip,J)];      % counter-clockwise

i  = (1:Ni)';
ip = mod(i,Ni) + 1;
blines = [id(i,1),    id(ip,1),    ones(Ni,1); ...
          id(i,Nj+1), id(ip,Nj+1), 2*ones(Ni,1)];

M = struct('nodes',nodes,'tris',zeros(0,3),'quads',quads,'blines',blines, ...
           'physTag',[1;2],'physName',{{'wall','farfield'}},'Lref',[]);
end


%% ------------------------------------------------------------------------
function M = read_txt(base)
nodes = dlmread([base '_nodes.txt']);
tris  = read_opt([base '_tri.txt'],3);
quads = read_opt([base '_quad.txt'],4);
bl    = read_opt([base '_bnd.txt'],3);
M = struct('nodes',nodes(:,1:2),'tris',tris,'quads',quads,'blines',bl, ...
           'physTag',[1;2],'physName',{{'body','farfield'}},'Lref',[]);
end


function A = read_opt(fname,ncol)
A = zeros(0,ncol);
if exist(fname,'file') == 2
    d = dir(fname);
    if d.bytes > 0
        A = dlmread(fname);
        A = reshape(A,[],ncol);
    end
end
end


%% ------------------------------------------------------------------------
function M = read_msh(fname)
txt = fileread(fname);
txt = strrep(txt,char(13),'');
L = strsplit(txt,char(10));
nL = numel(L);

isS = false(1,nL);
for k = 1:nL
    s = L{k};
    isS(k) = ~isempty(s) && s(1) == '$';
end
sIdx  = find(isS);
sName = L(sIdx);

[a,~] = sec(sName,sIdx,'MeshFormat');
if a == 0
    error('read_msh: %s is not a Gmsh file ($MeshFormat missing)',fname);
end
fmt = sscanf(L{a},'%f');
ver = fmt(1);
if numel(fmt) >= 2 && fmt(2) ~= 0
    error('read_msh: binary .msh files are not supported, save the mesh as ASCII');
end

% physical names
physTag  = [];
physName = {};
[a,b] = sec(sName,sIdx,'PhysicalNames');
if a > 0
    for k = a+1:b
        tk = regexp(L{k},'^\s*(\d+)\s+(\d+)\s+"(.*)"','tokens','once');
        if ~isempty(tk)
            physTag(end+1,1) = str2double(tk{2});   %#ok<AGROW>
            physName{end+1}  = tk{3};               %#ok<AGROW>
        end
    end
end

if ver < 3
    % ---------------- format 2.x -----------------
    [a,~] = sec(sName,sIdx,'Nodes');
    n = sscanf(L{a},'%d');
    D = sscanf(strjoin(L(a+1:a+n),' '),'%f',[4 n]).';
    ids = D(:,1);
    xy  = D(:,2:3);
    map = zeros(max(ids),1);
    map(ids) = 1:n;

    [a,~] = sec(sName,sIdx,'Elements');
    ne = sscanf(L{a},'%d');
    tris = zeros(ne,3);   nt = 0;
    quads = zeros(ne,4);  nq = 0;
    bl = zeros(ne,3);     nb = 0;
    for k = 1:ne
        v = sscanf(L{a+k},'%d').';
        et = v(2);
        ntag = v(3);
        nd = v(4+ntag:end);
        switch et
            case 1
                nb = nb + 1;
                ph = 0;
                if ntag >= 1
                    ph = v(4);
                end
                bl(nb,:) = [map(nd(1)) map(nd(2)) ph];
            case 2
                nt = nt + 1;
                tris(nt,:) = map(nd(1:3)).';
            case 3
                nq = nq + 1;
                quads(nq,:) = map(nd(1:4)).';
        end
    end
    M = struct('nodes',xy,'tris',tris(1:nt,:),'quads',quads(1:nq,:),'blines',bl(1:nb,:), ...
               'physTag',physTag,'physName',{physName},'Lref',[]);

else
    % ---------------- format 4.x -----------------
    % curve entity -> first physical tag
    curvePhys = containers.Map('KeyType','double','ValueType','double');
    [a,~] = sec(sName,sIdx,'Entities');
    if a > 0
        h = sscanf(L{a},'%d').';
        np = h(1);  nc = h(2);
        ptr = a + 1 + np;
        for k = 1:nc
            v = sscanf(L{ptr},'%f').';
            if v(8) > 0
                curvePhys(v(1)) = v(9);
            end
            ptr = ptr + 1;
        end
    end

    [a,~] = sec(sName,sIdx,'Nodes');
    h = sscanf(L{a},'%d').';
    nb = h(1);  nn = h(2);
    ptr = a + 1;
    xy  = zeros(nn,2);
    ids = zeros(nn,1);
    cnt = 0;
    for b = 1:nb
        bh = sscanf(L{ptr},'%d').';
        ptr = ptr + 1;
        nin = bh(4);
        if nin > 0
            tg = sscanf(strjoin(L(ptr:ptr+nin-1),' '),'%d');
            ptr = ptr + nin;
            co = sscanf(strjoin(L(ptr:ptr+nin-1),' '),'%f',[3 nin]).';
            ptr = ptr + nin;
            ids(cnt+1:cnt+nin)  = tg;
            xy(cnt+1:cnt+nin,:) = co(:,1:2);
            cnt = cnt + nin;
        end
    end
    map = zeros(max(ids),1);
    map(ids) = 1:nn;

    [a,~] = sec(sName,sIdx,'Elements');
    h = sscanf(L{a},'%d').';
    nb = h(1);  ne = h(2);
    ptr = a + 1;
    tris = zeros(ne,3);   nt = 0;
    quads = zeros(ne,4);  nq = 0;
    bl = zeros(ne,3);     nbl = 0;
    npe = containers.Map('KeyType','double','ValueType','double');
    npe(1) = 2;  npe(2) = 3;  npe(3) = 4;
    for b = 1:nb
        bh = sscanf(L{ptr},'%d').';
        ptr = ptr + 1;
        etag = bh(2);  et = bh(3);  nin = bh(4);
        if nin > 0 && isKey(npe,et)
            k = npe(et);
            D = sscanf(strjoin(L(ptr:ptr+nin-1),' '),'%d',[k+1 nin]).';
            nd = map(D(:,2:end));
            switch et
                case 1
                    ph = -etag;
                    if isKey(curvePhys,etag)
                        ph = curvePhys(etag);
                    end
                    bl(nbl+1:nbl+nin,:) = [nd(:,1) nd(:,2) ph*ones(nin,1)];
                    nbl = nbl + nin;
                case 2
                    tris(nt+1:nt+nin,:) = nd;
                    nt = nt + nin;
                case 3
                    quads(nq+1:nq+nin,:) = nd;
                    nq = nq + nin;
            end
        end
        ptr = ptr + nin;
    end
    M = struct('nodes',xy,'tris',tris(1:nt,:),'quads',quads(1:nq,:),'blines',bl(1:nbl,:), ...
               'physTag',physTag,'physName',{physName},'Lref',[]);
end

end


%% ------------------------------------------------------------------------
function [a,b] = sec(sName,sIdx,name)
% first and last data line of the section $name ... $Endname
k1 = find(strcmp(strtrim(sName),['$' name]),1);
k2 = find(strcmp(strtrim(sName),['$End' name]),1);
if isempty(k1) || isempty(k2)
    a = 0;
    b = -1;
else
    a = sIdx(k1) + 1;
    b = sIdx(k2) - 1;
end
end
