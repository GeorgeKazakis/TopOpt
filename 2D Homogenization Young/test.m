%% Test for homogenization
clc; clear;
x = randi([1 2],10);
E = 1 * (x==1) + 0.0001 * (x==2);
homogenizeY(1,1,E,0.2,90);