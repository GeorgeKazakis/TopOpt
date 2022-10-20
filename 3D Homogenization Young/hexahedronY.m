%% COMPUTE ELEMENT STIFFNESS MATRIX AND LOAD VECTOR
function [ke, fe] = hexahedronY(a, b, c, nu)
% Constitutive matrix contributions
E = 1;
C = E/((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
    0 0 0 0 0 (1-2*nu)/2];
% Three Gauss points in both directions
% xx = [-1/sqrt(3), 1/sqrt(3)];
xx = [-sqrt(3/5), 0, sqrt(3/5)]; 
yy = xx; zz = xx;
ww = [5/9, 8/9, 5/9];
% Initialize
ke = zeros(24,24);
fe = zeros(24,6); 
for ii = 1:length(xx)
    for jj = 1:length(yy)
        for kk = 1:length(zz)
            %integration point
            x = xx(ii); y = yy(jj); z = zz(kk);
            %stress strain displacement matrix
            qx = [ -((y-1)*(z-1))/8, ((y-1)*(z-1))/8, -((y+1)*(z-1))/8,...
                ((y+1)*(z-1))/8, ((y-1)*(z+1))/8, -((y-1)*(z+1))/8,...
                ((y+1)*(z+1))/8, -((y+1)*(z+1))/8];
            qy = [ -((x-1)*(z-1))/8, ((x+1)*(z-1))/8, -((x+1)*(z-1))/8,...
                ((x-1)*(z-1))/8, ((x-1)*(z+1))/8, -((x+1)*(z+1))/8,...
                ((x+1)*(z+1))/8, -((x-1)*(z+1))/8];
            qz = [ -((x-1)*(y-1))/8, ((x+1)*(y-1))/8, -((x+1)*(y+1))/8,...
                ((x-1)*(y+1))/8, ((x-1)*(y-1))/8, -((x+1)*(y-1))/8,...
                ((x+1)*(y+1))/8, -((x-1)*(y+1))/8];
            % Jacobian
            J = [qx; qy; qz]*[-a a a -a -a a a -a; -b -b b b -b -b b b;...
                -c -c -c -c c c c c]';
            qxyz = J\[qx;qy;qz];
            B_e = zeros(6,3,8);
            for i_B = 1:8
                B_e(:,:,i_B) = [qxyz(1,i_B)   0             0;
                                0             qxyz(2,i_B)   0;
                                0             0             qxyz(3,i_B);
                                qxyz(2,i_B)   qxyz(1,i_B)   0;
                                0             qxyz(3,i_B)   qxyz(2,i_B);
                                qxyz(3,i_B)   0             qxyz(1,i_B)];
            end
            B = [B_e(:,:,1) B_e(:,:,2) B_e(:,:,3) B_e(:,:,4) B_e(:,:,5)...
                B_e(:,:,6) B_e(:,:,7) B_e(:,:,8)];
            % Weight factor at this point
            weight = det(J)*ww(ii) * ww(jj) * ww(kk);
            % Element matrices
            ke = ke + weight * B' * C * B;
            % Element loads
            fe = fe + weight * B' * C;       
        end
    end
end
end