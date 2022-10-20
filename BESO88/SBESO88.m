%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function SBESO88(lx,ly,nelx,nely,volfrac,penal,rmin,er)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
dx = lx/nelx;
dy = ly/nely;
r = dx/dy;
t = r^(-1);
%% PREPARE FINITE ELEMENT ANALYSIS
% A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
% A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
% B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
% B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
A11 = [8*t+4*r  3 -8*t+2*r -3; 3 8*r+4*t 3 4*r-4*t; -8*t+2*r 3 8*t+4*r -3; -3 4*r-4*t -3 8*r+4*t];
A12 = [-4*t-2*r -3  4*t-4*r  3; -3 -4*r-2*t -3 -8*r+2*t; 4*t-4*r -3 -4*t-2*r 3; 3 -8*r+2*t 3 -4*r-2*t];
B11 = [-4*r  3 -2*r  9;  3 -4*t -9  4*t; -2*r -9 -4*r -3;  9  4*t -3 -4*t];
B12 = [ 2*r -3  4*r -9; -3  2*t  9 -2*t;  4*r  9  2*r  3; -9 -2*t  3  2*t];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2*(nelx+1)*(nely+1),1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = (1:1:2*(nely+1));
alldofs = (1:2*(nely+1)*(nelx+1));
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = ones(nely,nelx);
vol = 1;
% xPhys = x;
loop = 0;
change = 1;
%% START ITERATION
while change > 0.001
  loop = loop + 1;
  vol = max(vol*(1-er),volfrac);
  if loop > 1; olddc = dc; end
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+x(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(0.5*(sum((U(edofMat)*KE).*U(edofMat),2)),nely,nelx);
  c(loop) = sum(sum((Emin+x.^penal*(E0-Emin)).*ce));
  dc = penal*(E0-Emin)*x.^(penal-1).*ce;
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  dc(:) = H*(dc(:))./Hs;
  if loop > 1; dc = (dc+olddc)/2.; end
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = min(min(dc)); l2 = max(max(dc)); 
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    x = max(0.001,sign(dc-lmid));
%     xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    if sum(x(:)) > vol*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  	if loop > 10
		change  = abs(sum(c(loop-9:loop-5))-sum(c(loop-4:loop)))/sum(c(loop-4:loop));
    end
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c(loop),mean(x(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-x); caxis([0 1]); axis equal; axis off; drawnow;
end
%
