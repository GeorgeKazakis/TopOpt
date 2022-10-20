function levelset88(nelx,nely,volfrac,stepLength,numReinit,topWeight)
%% MATERIAL PROPERTIES
E0 = 1;
nu = 0.3;
Emin = 0.0001;
lambda = E0*nu/((1+nu)*(1-nu)); 
mu = E0/(2*(1+nu));
%% PREPARE FINITE ELEMENT ANALYSIS
% KE
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
% KTr
A11 = [4   3 -4   3;  3  4 -3  2; -4 -3  4 -3;  3  2 -3  4];
A12 = [-2 -3  2  -3; -3 -2  3 -4;  2  3 -2  3; -3 -4  3 -2];
KTr = 1/(1-nu)/12*([A11 A12;A12' A11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (Bridge)
F = sparse(2*(round(nelx/2)+1)*(nely+1),1,1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = [2*(nely+1)-1:2*(nely+1),2*(nelx+1)*(nely+1)-1:2*(nelx+1)*(nely+1)];
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% INITIALIZE ITERATION
x = ones(nely,nelx);
[lsf] = reinit(x); % Initialize levelset function
loop = 0;
%% START ITERATION
while loop < 200
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+x(:)'*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %%
  % Current Volume
  volCurr = sum(x(:))/(nelx*nely);
  % Set augmented Lagrangian parameters
  if loop == 1
      la = -0.01; La = 1000; alpha = 0.9;
  else
      la = la - 1/La * (volCurr - volfrac); La= alpha * La;
  end
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  ge = reshape(sum((U(edofMat)*KTr).*U(edofMat),2),nely,nelx);
  c(loop) = sum(sum((Emin+x*(E0-Emin)).*ce));
  dc = -((E0-Emin)*x+Emin).*ce - la + 1/La*(volCurr-volfrac);
  dg = x.*(pi/2*(lambda+2*mu)/mu/(lambda+mu)*(4*mu*ce+(lambda-mu)*ge))...
      + pi*(la - 1/La*(volCurr-volfrac));
  % Design update
  [x,lsf] = updateStep(lsf,dc,dg,stepLength,topWeight);
  % Reinitialize level-set function
  if ~mod(loop,numReinit)
      [lsf] = reinit(x);
  end
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f\n',loop,c(loop), ...
    mean(x(:)));
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-x); caxis([0 1]); axis equal; axis off; drawnow;
end
end
%
%%---- REINITIALIZATION OF LEVEL-SETFUNCTION ----
function [lsf] = reinit(struc)
strucFull = zeros(size(struc)+2);
strucFull(2:end-1,2:end-1) = struc;
% Use "bwdist" (Image Processing Toolbox)
lsf = (~strucFull).*(bwdist(strucFull)-0.5) - strucFull.*(bwdist(strucFull-1)-0.5);
end
%
%%----- DESIGN UPDATE ----
function [struc,lsf] = updateStep(lsf,shapeSens,topSens,stepLength,topWeight)
% Smooth the sensitivities
[shapeSens] = conv2(padarray(shapeSens,[1,1],'replicate'),1/6*[0 1 0; 1 2 1;0 1 0],'valid');
[topSens] = conv2(padarray(topSens,[1,1],'replicate'),1/6*[0 1 0; 1 2 1; 0 1 0],'valid');
% Load bearing pixels must remain solid -Bridge:
% edw nomizw einai karfota gia to paradeigma pou exei.
shapeSens(end,[1,round(end/2):round(end/2+1),end]) = 0;
topSens(end,[1,round(end/2):round(end/2+1),end]) = 0;
% Design update via evolution
[struc,lsf] = evolve(-shapeSens,topSens.*(lsf(2:end-1,2:end-1)<0),lsf,stepLength,topWeight);
end
%
%%---- EVOLUTION OF LEVEL-SET FUNCTION----
function [struc,lsf] = evolve(v,g,lsf,stepLength,w)
% Extend sensitivites using a zero border
vFull = zeros(size(v)+2); vFull(2:end-1,2:end-1) = v;
gFull = zeros(size(g)+2); gFull(2:end-1,2:end-1) = g;
% Choose time step for evolution based onCFL value
dt = 0.1/max(abs(v(:)));
% Evolve for total time stepLength * CFLvalue:
for i = 1:(10*stepLength)
    % Calculate derivatives on the grid
    dpx = circshift(lsf,[0,-1])-lsf;
    dmx = lsf - circshift(lsf,[0,1]);
    dpy = circshift(lsf,[-1,0]) - lsf;
    dmy = lsf - circshift(lsf,[1,0]);
    % Update level set function using anupwind scheme
    lsf = lsf - dt * min(vFull,0).* ...
    sqrt( min(dmx,0).^2+max(dpx,0).^2+min(dmy,0).^2+max(dpy,0).^2 ) ...
    - dt * max(vFull,0) .*...
    sqrt( max(dmx,0).^2+min(dpx,0).^2+max(dmy,0).^2+min(dpy,0).^2 )...
    - w*dt*gFull;
end
% New structure obtained from lsf
strucFull = (lsf<0); struc = strucFull(2:end-1,2:end-1);
end
%