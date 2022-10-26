function [KE,FE] = Q4elementStiffnessMatrix(a, b, phi, C)
%elementMatVec COMPUTE ELEMENT STIFFNESS MATRIX AND FORCE VECTOR (NUMERICALLY)
%   Compute element stiffness matrix numerically and break it into to parts
%   concerning the young modulus and poisson ratio
%   input
%   phi: elements angle (usually 90)
%   C: elasticity tensor
%   output
%   KE: element stiffness matrix

% Two Gauss points in both directions
xx = [-1/sqrt(3), 1/sqrt(3)]; yy = xx;
ww = [1,1];
% Initialize 
KE = zeros(8,8);
FE = zeros(8,3);
L = zeros(3,4); L(1,1) = 1; L(2,4) = 1; L(3,2:3) = 1;
    for ii=1:length(xx)
        for jj=1:length(yy)
            % Integration point 
            x = xx(ii); y = yy(jj);
            % Differentiated shape functions
            dNx = 1/4*[-(1-y) (1-y) (1+y) -(1+y)];
            dNy = 1/4*[-(1-x) -(1+x) (1+x) (1-x)];
            % Jacobian
            J = [dNx; dNy]*[-a a a+2*b/tan(phi*pi/180) 2*b/tan(phi*pi/180)-a;...
                -b -b b b]';
            detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1);
            ER = [J(2,2) -J(1,2); -J(2,1) J(1,1)];
            invJ = 1/detJ*ER;
            % Weight factor at this point 
            weight = ww(ii)*ww(jj)*detJ;
            % Strain-displacement matrix
            G = [invJ zeros(2); zeros(2) invJ];
            dN = zeros(4,8);
            dN(1,1:2:8) = dNx;
            dN(2,1:2:8) = dNy;
            dN(3,2:2:8) = dNx;
            dN(4,2:2:8) = dNy;
            B = L*G*dN;
            % Element matrices
            KE = KE + weight*(B' * C * B);
            % Element Loads
            FE = FE + weight*(B' * C * diag([1 1 1]));
        end
    end
end