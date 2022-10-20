function CH = homo3DY(lx,ly,lz,E,nu,Emin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lx       = Unit cell length in x-direction.
% ly       = Unit cell length in y-direction.
% lz       = Unit cell length in z-direction.
% lambda   = Lame's first parameter for solid materials.
% mu       = Lame's second parameter for solid materials.
% voxel    = Material indicator matrix. Used to determine nelx/nely/nelz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE
[nelx, nely, nelz] = size(E); %size of voxel model along x,y and z axis
% Stiffness matrix
dx = lx/nelx; dy = ly/nely; dz = lz/nelz;
nel = nelx*nely*nelz;
[ke, fe] = hexahedronY(dx/2,dy/2,dz/2,nu);
% Node numbers and element degrees of freedom for full (not periodic) mesh
nodenrs = reshape(1:(1+nelx)*(1+nely)*(1+nelz),1+nelx,1+nely,1+nelz);
edofVec = reshape(3*nodenrs(1:end-1,1:end-1,1:end-1)+1,nel,1);
addx = [0 1 2 3*nelx+[3 4 5 0 1 2] -3 -2 -1];
addxy = 3*(nely+1)*(nelx+1)+addx;
edof = repmat(edofVec,1,24) + repmat([addx addxy],nel,1);
%% IMPOOSE PERIODIC BOUNDARY CONDITIONS
% Use original edofMat to index into list with the periodic dofs
nn = (nelx+1)*(nely+1)*(nelz+1); % Total number of nodes
nnP = (nelx)*(nely)*(nelz);    % Total number of unique nodes
nnPArray = reshape(1:nnP, nelx, nely, nelz);
% Extend with a mirror of the back border
nnPArray(end+1,:,:) = nnPArray(1,:,:);
% Extend with a mirror of the left border
nnPArray(:, end+1, :) = nnPArray(:,1,:);
% Extend with a mirror of the top border
nnPArray(:, :, end+1) = nnPArray(:,:,1);
% Make a vector into which we can index using edofMat:
dofVector = zeros(3*nn, 1);
dofVector(1:3:end) = 3*nnPArray(:)-2;
dofVector(2:3:end) = 3*nnPArray(:)-1;
dofVector(3:3:end) = 3*nnPArray(:);
edof = dofVector(edof);
ndof = 3*nnP;
%% ASSEMBLE GLOBAL STIFFNESS MATRIX AND LOAD VECTORS
% Indexing vectors
iK = kron(edof,ones(24,1))';
jK = kron(edof,ones(1,24))';
% Material properties assigned to voxels with materials
% E = E*(voxel==1);
% The corresponding stiffness matrix entries
sK = ke(:)*E(:).';
K = sparse(iK(:), jK(:), sK(:), ndof, ndof);
K = 1/2*(K+K');
% Assembly three load cases corresponding to the three strain cases
iF = repmat(edof',6,1);
jF = [ones(24,nel); 2*ones(24,nel); 3*ones(24,nel);...
    4*ones(24,nel); 5*ones(24,nel); 6*ones(24,nel);];
sF = fe(:)*E(:).';
F  = sparse(iF(:), jF(:), sF(:), ndof, 6);
%% SOLUTION
% solve by PCG method, remember to constrain one node
activedofs = edof(E>Emin,:); activedofs = sort(unique(activedofs(:)));
X = zeros(ndof,6);
L = ichol(K(activedofs(4:end),activedofs(4:end)));
for i = 1:6
    X(activedofs(4:end),i) = pcg(K(activedofs(4:end),...
        activedofs(4:end)),F(activedofs(4:end),i),1e-10,300,L,L');
end
% X(activedofs(4:end),:) = K(activedofs(4:end),activedofs(4:end))...
%     \F(activedofs(4:end),:);    % Solving by direct method
%% HOMOGENIZATION
% The displacement vectors corresponding to the unit strain cases
X0 = zeros(nel, 24, 6);
% The element displacements for the six unit strains
X0_e = zeros(24, 6);
%fix degrees of nodes [1 2 3 5 6 12];
% ke = keMu + keLambda; % Here the exact ratio does not matter, because
% fe = feMu + feLambda; % it is reflected in the load vector
X0_e([4 7:11 13:24],:) = ke([4 7:11 13:24],[4 7:11 13:24])...
                           \fe([4 7:11 13:24],:);
X0(:,:,1) = kron(X0_e(:,1)', ones(nel,1)); % epsilon0_11 = (1,0,0,0,0,0)
X0(:,:,2) = kron(X0_e(:,2)', ones(nel,1)); % epsilon0_22 = (0,1,0,0,0,0)
X0(:,:,3) = kron(X0_e(:,3)', ones(nel,1)); % epsilon0_33 = (0,0,1,0,0,0)
X0(:,:,4) = kron(X0_e(:,4)', ones(nel,1)); % epsilon0_12 = (0,0,0,1,0,0)
X0(:,:,5) = kron(X0_e(:,5)', ones(nel,1)); % epsilon0_23 = (0,0,0,0,1,0)
X0(:,:,6) = kron(X0_e(:,6)', ones(nel,1)); % epsilon0_13 = (0,0,0,0,0,1)
CH = zeros(6);
DCH = cell(nely,nelx,nelz);
DCH(:,:,:) = {zeros(6)};
volume = lx*ly*lz;
for i = 1:6
    for j = 1:6
        sum_Y = ((X0(:,:,i) - X(edof+(i-1)*ndof))*ke).*...
            (X0(:,:,j) - X(edof+(j-1)*ndof));
        sum_Y = reshape(sum(sum_Y,2), nelx, nely, nelz);
        % Homogenized elasticity tensor
        CH(i,j) = 1/volume*sum(sum(sum(E.*sum_Y)));
    end
end
end