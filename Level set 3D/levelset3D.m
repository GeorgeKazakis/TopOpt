% A 169 LINE 3D TOPOLOGY OPITMIZATION CODE BY LIU AND TOVAR (JUL 2013)
function levelset3D(nelx,nely,nelz,volfrac,stepLength,numReinit,topWeight)
% USER-DEFINED LOOP PARAMETERS
maxloop = 200;    % Maximum number of iterations
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
% KE
KE = lk_H8(nu);
% KTr
dx = 1; dy = 1; dz = 1;
[KTr] = D3_LB(dx, dy, dz, nu);
% lambda & mu
lambda = E0*nu/((1+nu)*(1-2*nu)); 
mu = E0/(2*(1+nu));
%%%%%%%%%%%%%%%%%%%
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
% INITIALIZE ITERATION
x = ones(nely,nelx,nelz); 
[lsf] = reinit(x);
loop = 0; 
% START ITERATION
while loop < maxloop
    loop = loop+1;
    % FE-ANALYSIS
    sK = KE(:)*(Emin+x(:)'*(E0-Emin));
    K = sparse(iK(:),jK(:),sK(:)); K = (K+K')/2;
    tolit=1e-8;
    maxit=8000;
    M=diag(diag(K(freedofs,freedofs)));
    U(freedofs,:) = pcg(K(freedofs,freedofs),F(freedofs,:),tolit,maxit,M);
    %U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    % Current Volume
    volCurr = sum(x(:))/(nelx*nely*nelz);
    % Set augmented Lagrangian parameters
    if loop == 1
        la = -0.01; La = 1000; alpha = 0.9;
    else
        la = la - 1/La * (volCurr - volfrac); La= alpha * La;
    end
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    ge = reshape(sum((U(edofMat)*KTr).*U(edofMat),2),[nely,nelx,nelz]);
    c(loop) = sum(sum(sum((Emin+x*(E0-Emin)).*ce)));
    dc = -((E0-Emin)*x+Emin).*ce - la + 1/La*(volCurr-volfrac);
    dg = x.*(pi*(lambda+2*mu)/mu/(9*lambda+14*mu)*(20*mu*ce+(3*lambda-2*mu)*ge))...
        + pi*(la - 1/La*(volCurr-volfrac));
    % Design update
    [x,lsf] = updateStep(lsf,dc,dg,stepLength,topWeight);
    % Reinitialize level-set function
    if ~mod(loop,numReinit)
        [lsf] = reinit(x);
    end
    % PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f\n',loop,c(loop),mean(x(:)));
    % PLOT DENSITIES
    if displayflag, clf; display_3D(x); end
end
clf; display_3D(x);
end
% ===================== AUXILIARY FUNCTIONS ===============================
% GENERATE ELEMENT STIFFNESS MATRIX
function [KE] = lk_H8(nu)
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/72*A'*[1; nu];
% GENERATE SIX SUB-MATRICES AND THEN GET KE MATRIX
K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
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

%%---- REINITIALIZATION OF LEVEL-SETFUNCTION ----
function [lsf] = reinit(struc)
strucFull = zeros(size(struc)+2);
strucFull(2:end-1,2:end-1,2:end-1) = struc;
% Use "bwdist" (Image Processing Toolbox)
lsf = (~strucFull).*(bwdist(strucFull)-0.5) - strucFull.*(bwdist(strucFull-1)-0.5);
end
%
%%----- DESIGN UPDATE ----
function [struc,lsf] = updateStep(lsf,shapeSens,topSens,stepLength,topWeight)
% Smooth the sensitivities
[shapeSens] = convn(padarray(shapeSens,[1,1,1],'replicate'),ConvolutionN(),'valid');
[topSens]   = convn(padarray(topSens,[1,1,1],'replicate'),ConvolutionN(),'valid');
% Load bearing pixels must remain solid -Bridge:
% edw nomizw einai karfota gia to paradeigma pou exei.
shapeSens(1:end,1,1:end) = 0;
shapeSens(end,end,1:end) = 0;
topSens(1:end,1,1:end) = 0;
topSens(end,end,1:end) = 0;
% Design update via evolution
[struc,lsf] = evolve(-shapeSens,topSens.*(lsf(2:end-1,2:end-1,2:end-1)<0),lsf,stepLength,topWeight);
end
%
%%---- EVOLUTION OF LEVEL-SET FUNCTION----
function [struc,lsf] = evolve(v,g,lsf,stepLength,w)
% Extend sensitivites using a zero border
vFull = zeros(size(v)+2); vFull(2:end-1,2:end-1,2:end-1) = v;
gFull = zeros(size(g)+2); gFull(2:end-1,2:end-1,2:end-1) = g;
% Choose time step for evolution based onCFL value
dt = 0.1/max(abs(v(:)));
% Evolve for total time stepLength * CFLvalue:
for i = 1:(10*stepLength)
    % Calculate derivatives on the grid
    dpx = circshift(lsf,[0,-1,0])-lsf;
    dmx = lsf - circshift(lsf,[0,1,0]);
    dpy = circshift(lsf,[-1,0,0]) - lsf;
    dmy = lsf - circshift(lsf,[1,0,0]);
    dpz = circshift(lsf,[0,0,-1]) - lsf;
    dmz = lsf - circshift(lsf,[0,0,1]);
    % Update level set function using anupwind scheme
    lsf = lsf - dt * min(vFull,0).* ...
    sqrt( min(dmx,0).^2+max(dpx,0).^2+min(dmy,0).^2+max(dpy,0).^2+min(dmz,0).^2+max(dpz,0).^2 ) ...
    - dt * max(vFull,0) .*...
    sqrt( max(dmx,0).^2+min(dpx,0).^2+max(dmy,0).^2+min(dpy,0).^2+max(dmz,0).^2+min(dpz,0).^2 )...
    - w*dt*gFull;
end
% New structure obtained from lsf
strucFull = (lsf<0); struc = strucFull(2:end-1,2:end-1,2:end-1);
end
%
function matric = ConvolutionN()
matric = zeros(3,3,3);
matric(:,:,1) = 1/6*[0 1 0; 1 2 1;0 1 0];
matric(:,:,2) = 1/6*[0 1 0; 1 2 1;0 1 0];
matric(:,:,3) = 1/6*[0 1 0; 1 2 1;0 1 0];
end