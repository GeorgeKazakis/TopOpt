function c = UCOpt(lx,ly,nelx,nely,nlx,nly,volfrac,penal,rmin,ft)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2*(nelx+1)*(nely+1),1,-1,2*(nely+1)*(nelx+1),1);
fixeddofs = 1:2*(nely+1);
% F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
% fixeddofs = union(1:2:2*(nely+1),2*(nelx+1)*(nely+1));
U = zeros(2*(nely+1)*(nelx+1),1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nlx*nly*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nlx
  for j1 = 1:nly
    e1 = (i1-1)*nly+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nlx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nly)
        e2 = (i2-1)*nly+j2;
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
x = ones(nly,nlx);
for i = 1:nlx
    for j = 1:nly
        if sqrt((i-nlx/2-0.5)^2+(j-nly/2-0.5)^2) < min(nlx,nly)/3
            x(j,i) = 0;
        end
    end
end
xPhys = x;
loop = 0;
change = 1;
%% START ITERATION
while change > 0.01
  loop = loop + 1;
  %% INTERPOLATION
  [E, dE] = interpolate(x,E0,Emin,penal);
  %% HOMOGENIZATION
  [CH,DCH] = homogenize(0.000001*lx,0.000001*ly,E,nu,dE,90);
  %% FE-ANALYSIS
  KE = Q4elementStiffnessMatrix(lx/nelx,ly/nely,90,CH);
  sK = reshape(KE(:)*ones(1,nely*nelx),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  f = @(DCh) Q4elementStiffnessMatrix(lx/nelx,ly/nely,90,DCh);
  dKe = cellfun(f,DCH,'UniformOutput',false);
  f= @(dKE) reshape(sum((U(edofMat)*dKE).*U(edofMat),2),nely,nelx);
  ce = cellfun(f,dKe,'UniformOutput',false);
  c = sum(sum(cell2mat(cellfun(@sum, ce, 'UniformOutput',false))));
  dc = - cell2mat(cellfun(@sum,cellfun(@sum,cellfun(@sum, ce, 'UniformOutput',false),'UniformOutput',false),'UniformOutput',false));
  dv = ones(nly,nlx);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
    end
    if sum(xPhys(:)) > volfrac*nlx*nly, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
end

function [E,dE] = interpolate(x, E0, Emin, penal)
    E = Emin + x.^penal * (E0-Emin);
    dE = penal * x.^(penal-1) * (E0-Emin);
end
