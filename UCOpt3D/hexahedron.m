classdef hexahedron 
    %hexahedron hexahedron element stiffness matrix

    properties
        weights;
        Bv;
    end

    methods
        function obj = hexahedron(a,b,c)
            %hexahedron Construct an instance of this class
            %   compute weights and B for later use
            [weights,B] = obj.initialize(a,b,c);
            obj.weights = weights;
            obj.Bv = B;
        end

        function [ke,fe] = elemStiffness(obj,C)
            ke = zeros(24,24);
            fe = zeros(24,6); 
            for i = 1:size(obj.weights,1)
                ke = ke + obj.weights(i) * obj.Bv{i}' * C * obj.Bv{i};
                fe = fe + obj.weights(i) * obj.Bv{i}' * C;
            end
        end

        function [ke,fe] = elemStiffnessC(obj,E,nu)
            C = obj.elasticityTensor(E,nu);
            [ke,fe] = elemStiffness(obj,C);
        end
    end
    methods(Static)
        function [weights,Bv] = initialize(a,b,c)
            % Three Gauss points in both directions
            xx = [-sqrt(3/5), 0, sqrt(3/5)]; 
            yy = xx; zz = xx;
            ww = [5/9, 8/9, 5/9];
            % Initialize
            weights = zeros(length(xx) * length(yy) * length(zz),1);
            Bv = cell(length(xx) * length(yy) * length(zz),1);
            count = 0;
            for ii = 1:length(xx)
                for jj = 1:length(yy)
                    for kk = 1:length(zz)
                        count = count + 1;
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
                        % Save
                        weights(count) = weight;
                        Bv{count} = B;      
                    end
                end
            end
        end

        function [C] = elasticityTensor(E,nu)
            C = E/((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
                nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
                0 0 0 0 0 (1-2*nu)/2];
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was modified by G. Kazakis and N.D. Lagaros and was 
% based on the homo3D.m code written by G. Dong, Y. Tang, Y.F. Zhao, 
% in the paper "A 149 Line Homogenization Code for Three-Dimensional 
% Cellular Materials Written in matlab", G. Dong, Y. Tang, Y.F. Zhao,
% Journal of Engineering Materials and Technology, 2018
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