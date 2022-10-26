function [CH,DCH] = homogenize(lx, ly, E, nu, dE, phi)
    %% INITIALIZE
    % Deduce discretization
    [nely, nelx] = size(E);
    dx = lx/nelx; dy = ly/nely;
    nel = nely*nelx;
    [ke, fe] = elementMatVec(dx/2, dy/2, phi, nu);
    % Node number and element degrees of freedom for full(not periodic) mesh
    nodenrs = reshape(1:(nelx+1)*(nely+1),1+nely,1+nelx);
    edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nel,1);
    edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nel,1);
    %% IMPOSE PERIODIC BOUNDARY CONDITIONS
    % Use original edofMat to index into list with the periodic dofs
    nn = (nelx+1)*(nely+1); % Total number of nodes
    nnP = (nelx)*(nely); % total number of unique nodes
    nnPArray = reshape(1:nnP,nely ,nelx);
    % Extend with a mirror of the top border
    nnPArray(end+1,:) = nnPArray(1,:);
    % Extend with a mirror of the left border
    nnPArray(:,end+1) = nnPArray(:,1);
    % Make a vector into which we can index using edofMat:
    dofVector = zeros(2*nn,1);
    dofVector(1:2:end) = 2*nnPArray(:)-1;
    dofVector(2:2:end) = 2*nnPArray(:);
    edofMat = dofVector(edofMat);
    ndof = 2*nnP; % Number of dofs
    %% ASSEMBLE STIFFNESS MATRIX
    % Indexing vectors
    iK = kron(edofMat,ones(8,1))';
    jK = kron(edofMat,ones(1,8))';
    % The corresponding stiffness matrix entries
    sK = ke(:)*E(:).';
    K = sparse(iK(:), jK(:), sK(:), ndof, ndof);
    %% LOAD VECTORS AND SOLUTION
    % Assembly three load cases corresponding to the three strain cases
    sF = fe(:)*E(:).';
    iF = repmat(edofMat',3,1);
    jF = [ones(8,nel); 2*ones(8,nel); 3*ones(8,nel)];
    F = sparse(iF(:), jF(:), sF(:), ndof, 3);
    % Solve (remember to constrain one node)
    % K = 0.5*(K+K.');
    chi(3:ndof,:) = K(3:ndof,3:ndof)\F(3:ndof,:); % 1,2 dofs constrained
    %% HOMOGENIZATION
    % The displacement vector corresponding to the unit strain cases
    chi0 = zeros(nel, 8, 3);
    % The element displacements for the three unit strains
    chi0_e = zeros(8, 3);
    chi0_e([3 5:end],:) = ke([3 5:end],[3 5:end])\fe([3 5:end],:);
    % epsilon0_11 = (1, 0, 0)
    chi0(:,:,1) = kron(chi0_e(:,1)', ones(nel,1));
    % epsilon0_22 = (0, 1, 0)
    chi0(:,:,2) = kron(chi0_e(:,2)', ones(nel,1));
    % epsilon0_12 = (0, 0, 1)
    chi0(:,:,3) = kron(chi0_e(:,3)', ones(nel,1));
    CH = zeros(3);
    DCH = cell(nely,nelx);
    DCH(:,:) = {zeros(3)};
    cellVolume = lx*ly;
    for i =1:3
        for j = 1:3
            sumYoung = ((chi0(:,:,i) - chi(edofMat+(i-1)*ndof))*ke).*...
                (chi0(:,:,j) - chi(edofMat+(j-1)*ndof));
            sumYoung = reshape(sum(sumYoung,2), nely, nelx);
            % Homogenized elasticity tensor
            CH(i,j) = 1/cellVolume*sum(sum(E.*sumYoung));
            finalSum = dE.*sumYoung;
            for k=1:nely
                for l=1:nelx
                    DCH{k,l}(i,j) = 1/cellVolume*finalSum(k,l);
                end
            end
        end
    end
%     disp('--- Homogenized elasticity tensor ---'); disp(CH)
end