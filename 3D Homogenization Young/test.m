clc; clear; close all;

lx = 1; ly = lx; lz = lx;
E = 1;
nu = 0.2;

[voxel,Density] = GenerateVoxel(40,'topology/x.txt',0.1);

E = E*(voxel==1);
CH = homo3DY(lx,ly,lz,E,nu,10e-3);
visual(CH);