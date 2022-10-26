function UCOpt3D(lx,ly,lz,nelx,nely,nelz,nlx,nly,nlz,volfrac,penal,rmin,ft)
% USER-DEFINED LOOP PARAMETERS
maxloop = 200;    % Maximum number of iterations
tolx = 0.01;      % Termination criterion
displayflag = 1;  % Display structure flag
% USER-DEFINED MATERIAL PROPERTIES
E0 = 1;           % Young's modulus of solid material
Emin = 1e-9;      % Young's modulus of void-like material
nu = 0.3;         % Poisson's ratio
% USER-DEFINED LOAD DOFs
il = nelx; jl = 0; kl = 0:nelz;                         % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
% USER-DEFINED SUPPORT FIXED DOFs
[jf,kf] = meshgrid(1:nely+1,1:nelz+1);                  % Coordinates
fixednid = (kf-1)*(nely+1)*(nelx+1)+jf;                 % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-1,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
% KE = lk_H8(nu);
h1 = hexahedron(lx/nelx/2,ly/nely/2,lz/nelz/2);
h2 = hexahedron(0.000001*lx/nlx/2,0.000001*ly/nly/2,0.000001*lz/nlz/2);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = kron(edofMat,ones(24,1))';
jK = kron(edofMat,ones(1,24))';
% PREPARE FILTER
iH = ones(nlx*nly*nlz*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nlz
    for i1 = 1:nlx
        for j1 = 1:nly
            e1 = (k1-1)*nlx*nly + (i1-1)*nly+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nlz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nlx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nly)
                        e2 = (k2-1)*nlx*nly + (i2-1)*nly+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
% INITIALIZE ITERATION
x = repmat(ones,[nly,nlx,nlz]);
for i=1:nly
    for j=1:nlx
        for k=1:nlz
            if sqrt((i-nlx/2-0.5)^2+(j-nly/2-0.5)^2+(k-nlz/2-0.5)^2) < min(min(nlx,nly),nlz)/3
                x(i,j,k) = 0;
            end
        end
    end
end
xPhys = x; 
loop = 0; 
change = 1;
% START ITERATION
while change > tolx && loop < maxloop
    loop = loop+1;
    % INTERPOLATION
    [E,dE] = interpolate(x,E0,Emin,penal);
    % Homogenization
    [CH,DCH] = homo3DY(0.000001*lx,0.000001*ly,0.000001*lz,E,dE,nu,Emin,h2);
    % FE-ANALYSIS
    KE = h1.elemStiffness(CH);
    sK = reshape(KE(:)*ones(1,nely*nelx*nelz),576*nelx*nely*nelz,1);
    K = sparse(iK(:),jK(:),sK(:)); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    dKe = cellfun(@h1.elemStiffness,DCH,'UniformOutput',false);
    f = @(dKE) (reshape(sum((U(edofMat)*dKE).*U(edofMat),2),[nely,nelx,nelz]));
    ce = cellfun(f, dKe, 'UniformOutput',false);
    c = sum(sum(sum(cell2mat(cellfun(@sum, ce, 'UniformOutput',false)))));
    dc = - cell2mat(cellfun(@sum,cellfun(@sum,cellfun(@sum, ce, 'UniformOutput',false),'UniformOutput',false),'UniformOutput',false));
    dv = ones(nly,nlx,nlz);
    % FILTERING AND MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
    % OPTIMALITY CRITERIA UPDATE
    l1 = 0; l2 = 1e9; move = 0.2;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
        if ft == 1
            xPhys = xnew;
        elseif ft == 2
            xPhys(:) = (H*xnew(:))./Hs;
        end
        if sum(xPhys(:)) > volfrac*nlx*nly*nlz, l1 = lmid; else l2 = lmid; end
    end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    % PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
    % PLOT DENSITIES
    if displayflag, clf; display_3D(xPhys); end
end
clf; display_3D(xPhys);
end
% ===================== AUXILIARY FUNCTIONS ===============================
function [E,dE] = interpolate(x, E0, Emin, penal)
    E = Emin + x.^penal * (E0-Emin);
    dE = penal * x.^(penal-1) * (E0-Emin);
end
% DISPLAY 3D TOPOLOGY (ISO-VIEW)
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.5)  % User-defined display density threshold
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was modified by G. Kazakis and N.D. Lagaros and was 
% based on the top3d.m code written by K. Liu and A. Tovar, 
% in the paper "An efficient 3D topology optimization code written in 
% Matlab", E. K. Liu and A. Tovar, Struct Multidisc Optim, 2014
%
% The code is intended for educational purposes, extensions can be found in
% the paper "Topology optimization based material design for 3D domains
% using MATLAB"
%
% Disclaimer:                                                              
% The authors reserves all rights for the program.    
% The code may be distributed and used for educational purposes.       
% The authors do not guarantee that the code is free from errors, and
% they shall not be liable in any event caused by the use of the program.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%