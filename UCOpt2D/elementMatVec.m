function [ke, fe] = elementMatVec(a, b, phi, nu)
        %elementMatVec COMPUTE ELEMENT STIFFNESS MATRIX AND FORCE VECTOR (NUMERICALLY)
        %   Compute element stiffness matrix numerically and break it into to parts
        %   concerning the young modulus and poisson ratio
        %   input
        %   a: x size of element / 2 (dx/2)
        %   b: y size of element / 2 (dy/2)
        %   phi: elements angle (usually 90)
        %   output
        %   ke:   ke connected with young modulus
        %   fe:   fe connected with young modulus
        A = [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
        C = 1/(1-nu^2)*A;
        [ke,fe] = Q4elementStiffnessMatrix(a, b, phi, C);
end