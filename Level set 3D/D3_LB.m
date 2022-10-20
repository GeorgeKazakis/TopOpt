function [K,B,E,V] = D3_LB (dx, dy, dz, miu) 

% x1 = -1; y1 = -1; z1 = -1;
% x2 =  1; y2 = -1; z2 = -1;
% x3 =  1; y3 =  1; z3 = -1;
% x4 = -1; y4 =  1; z4 = -1;
% x5 = -1; y5 = -1; z5 =  1;
% x6 =  1; y6 = -1; z6 =  1;
% x7 =  1; y7 =  1; z7 =  1;
% x8 = -1; y8 =  1; z8 =  1;
% Em = 1; miu = 0.3;
x1 = -dx/2; y1 = -dy/2; z1 = -dz/2;
x2 = dx/2; y2 = -dy/2; z2 = -dz/2;
x3 = dx/2; y3 = dy/2; z3 = -dz/2;
x4 = -dx/2; y4 = dy/2; z4 = -dz/2;
x5 = -dx/2; y5 = -dy/2; z5 = dz/2;
x6 = dx/2; y6 = -dy/2; z6 = dz/2;
x7 = dx/2; y7 = dy/2; z7 = dz/2;
x8 = -dx/2; y8 = dy/2; z8 = dz/2;
Em = 1;

% s=(1-miu);
E  = Em*(1-miu)/((1+miu)*(1-2*miu)).*...
         [   1   1  1   0   0  0;
             1   1  1   0   0  0;
             1   1  1   0   0  0;
             0   0  0   0   0  0;
             0   0  0   0   0  0;
             0   0  0   0   0  0];

K  = zeros(24,24); K2 = zeros(24,24); K3 = zeros(24,24);
B  = zeros(6,24);
V2 = 0; V3 = 0; V = 0;

int_point = [ -sqrt(3/5)  0  sqrt(3/5) ];


for i=1:3
   if i==1, xi1=int_point(1); end
   if i==2, xi1=int_point(2); end
   if i==3, xi1=int_point(3); end
   for j=1:3
      if j==1, xi2=int_point(1); end
      if j==2, xi2=int_point(2); end
      if j==3, xi2=int_point(3); end
      for k=1:3
         if k==1, xi3=int_point(1); end
         if k==2, xi3=int_point(2); end
         if k==3, xi3=int_point(3); end
         R = 1/8.*[  (-1+xi3+xi2-xi2*xi3) (1-xi3-xi2+xi2*xi3) (1-xi3+xi2-xi2*xi3) ...
               (-1+xi3-xi2+xi2*xi3) (-1-xi3+xi2+xi2*xi3) (1+xi3-xi2-xi2*xi3) ...
               (1+xi3+xi2+xi2*xi3) (-1-xi3-xi2-xi2*xi3);
            (-1+xi3+xi1-xi1*xi3) (-1+xi3-xi1+xi1*xi3) (1-xi3+xi1-xi1*xi3) ...
               (1-xi3-xi1+xi1*xi3) (-1-xi3+xi1+xi1*xi3) (-1-xi3-xi1-xi1*xi3) ...
               (1+xi3+xi1+xi1*xi3) (1+xi3-xi1-xi1*xi3);
            (-1+xi2+xi1-xi1*xi2) (-1+xi2-xi1+xi1*xi2) (-1-xi2-xi1-xi1*xi2) ...
               (-1-xi2+xi1+xi1*xi2) (1-xi2-xi1+xi1*xi2) (1-xi2+xi1-xi1*xi2) ...
               (1+xi2+xi1+xi1*xi2) (1+xi2-xi1-xi1*xi2)];

         J =  R * [ x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4; x5 y5 z5;
                    x6 y6 z6; x7 y7 z7; x8 y8 z8];
         if inv(J) <=0, disp('     error: negative Jacobian in element'); end
         % dN = inv(J)*R;
         dN = J\R;

         B1 = [dN(1,1) 0  0  dN(1,2) 0  0  dN(1,3) 0  0 dN(1,4) 0  0 dN(1,5) 0  0 ...
                 dN(1,6) 0  0 dN(1,7) 0  0 dN(1,8) 0  0;
               0 dN(2,1) 0  0 dN(2,2) 0  0 dN(2,3) 0 0 dN(2,4) 0 0 dN(2,5) 0 0 ...
                 dN(2,6) 0  0 dN(2,7) 0  0 dN(2,8) 0;
               0 0 dN(3,1) 0 0 dN(3,2) 0  0 dN(3,3) 0 0 dN(3,4) 0 0 dN(3,5) 0 0 ...
                 dN(3,6) 0  0 dN(3,7) 0  0 dN(3,8);
               dN(2,1) dN(1,1) 0 dN(2,2) dN(1,2) 0 dN(2,3) dN(1,3) ...
               0 dN(2,4) dN(1,4) 0 dN(2,5) dN(1,5) 0 dN(2,6) dN(1,6) 0 ...
                 dN(2,7) dN(1,7) 0 dN(2,8) dN(1,8) 0;
               0 dN(3,1) dN(2,1) 0 dN(3,2) dN(2,2) 0 dN(3,3) dN(2,3) ...
               0 dN(3,4) dN(2,4) 0 dN(3,5) dN(2,5) 0 dN(3,6) dN(2,6) 0 ...
                 dN(3,7) dN(2,7) 0 dN(3,8) dN(2,8);
               dN(3,1) 0 dN(1,1) dN(3,2) 0 dN(1,2) dN(3,3) 0 dN(1,3) ...
               dN(3,4) 0 dN(1,4) dN(3,5) 0 dN(1,5) dN(3,6) 0 dN(1,6) ...
                 dN(3,7) 0 dN(1,7) dN(3,8) 0 dN(1,8)];
         
         Vo = det(J);
         K1 = B1'*E*B1.*Vo ;

         if k==1, V2 = V2 + 5/9.*Vo; end
         if k==2, V2 = V2 + 8/9.*Vo; end
         if k==3, V2 = V2 + 5/9.*Vo; end

         if (xi1==0 && xi2==0 && xi3==0), B(1:6,1:24,1) = B1; end
         
         if k==1, K2 = K2 + 5/9.*K1; end
         if k==2, K2 = K2 + 8/9.*K1; end
         if k==3, K2 = K2 + 5/9.*K1; end
      end

      if j==1, V3 = V3 + 5/9.*V2; end 
      if j==2, V3 = V3 + 8/9.*V2; end
      if j==3, V3 = V3 + 5/9.*V2; end
      
      if j==1, K3 = K3 + 5/9.*K2; end
      if j==2, K3 = K3 + 8/9.*K2; end
      if j==3, K3 = K3 + 5/9.*K2; end
      K2 = K2.*0; V2 = 0;
   end
   
   if i==1, V = V + 5/9.*V3; end
   if i==2, V = V + 8/9.*V3; end
   if i==3, V = V + 5/9.*V3; end
   
   if i==1, K = K + 5/9.*K3; end
   if i==2, K = K + 8/9.*K3; end
   if i==3, K = K + 5/9.*K3; end
   K3 = K3.*0; V3 = 0;
end

nd_shape = [1 -1 -1 -1; 2 +1 -1 -1; 3 +1 +1 -1; 4 -1 +1 -1; 5 -1 -1 +1; 6 +1 -1 +1; 7 +1 +1 +1; 8 -1 +1 +1];

for i=1:8
    xi1  = nd_shape(i,2);
    xi2  = nd_shape(i,3);
    xi3  = nd_shape(i,4);

    R = 1/8.*[  (-1+xi3+xi2-xi2*xi3) (1-xi3-xi2+xi2*xi3) (1-xi3+xi2-xi2*xi3) ...
            (-1+xi3-xi2+xi2*xi3) (-1-xi3+xi2+xi2*xi3) (1+xi3-xi2-xi2*xi3) ...
            (1+xi3+xi2+xi2*xi3) (-1-xi3-xi2-xi2*xi3);
        (-1+xi3+xi1-xi1*xi3) (-1+xi3-xi1+xi1*xi3) (1-xi3+xi1-xi1*xi3) ...
            (1-xi3-xi1+xi1*xi3) (-1-xi3+xi1+xi1*xi3) (-1-xi3-xi1-xi1*xi3) ...
            (1+xi3+xi1+xi1*xi3) (1+xi3-xi1-xi1*xi3);
        (-1+xi2+xi1-xi1*xi2) (-1+xi2-xi1+xi1*xi2) (-1-xi2-xi1-xi1*xi2) ...
            (-1-xi2+xi1+xi1*xi2) (1-xi2-xi1+xi1*xi2) (1-xi2+xi1-xi1*xi2) ...
            (1+xi2+xi1+xi1*xi2) (1+xi2-xi1-xi1*xi2)];

    J =  R * [ x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4; x5 y5 z5;
        x6 y6 z6; x7 y7 z7; x8 y8 z8];
    if inv(J) <=0, disp('     error: negative Jacobian in element'); end
%     dN = inv(J)*R;
    dN = J\R;    

         Bi = [dN(1,1) 0  0  dN(1,2) 0  0  dN(1,3) 0  0 dN(1,4) 0  0 dN(1,5) 0  0 ...
                 dN(1,6) 0  0 dN(1,7) 0  0 dN(1,8) 0  0;
               0 dN(2,1) 0  0 dN(2,2) 0  0 dN(2,3) 0 0 dN(2,4) 0 0 dN(2,5) 0 0 ...
                 dN(2,6) 0  0 dN(2,7) 0  0 dN(2,8) 0;
               0 0 dN(3,1) 0 0 dN(3,2) 0  0 dN(3,3) 0 0 dN(3,4) 0 0 dN(3,5) 0 0 ...
                 dN(3,6) 0  0 dN(3,7) 0  0 dN(3,8);
               dN(2,1) dN(1,1) 0 dN(2,2) dN(1,2) 0 dN(2,3) dN(1,3) ...
               0 dN(2,4) dN(1,4) 0 dN(2,5) dN(1,5) 0 dN(2,6) dN(1,6) 0 ...
                 dN(2,7) dN(1,7) 0 dN(2,8) dN(1,8) 0;
               0 dN(3,1) dN(2,1) 0 dN(3,2) dN(2,2) 0 dN(3,3) dN(2,3) ...
               0 dN(3,4) dN(2,4) 0 dN(3,5) dN(2,5) 0 dN(3,6) dN(2,6) 0 ...
                 dN(3,7) dN(2,7) 0 dN(3,8) dN(2,8);
               dN(3,1) 0 dN(1,1) dN(3,2) 0 dN(1,2) dN(3,3) 0 dN(1,3) ...
               dN(3,4) 0 dN(1,4) dN(3,5) 0 dN(1,5) dN(3,6) 0 dN(1,6) ...
                 dN(3,7) 0 dN(1,7) dN(3,8) 0 dN(1,8)];
        
    B(1:6,1:24,i+1) = Bi;
end

