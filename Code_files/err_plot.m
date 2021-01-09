%%error calculator
%% clear workspace
clear;
clc;
%% input loading matrics
load('T_gs_10');
load('T_gs_50');
load('T_gs_100');
load('T_gs_200');
load('T_gs_500');
load('T_gs_1000');

A=T_gs_10;
B=T_gs_100;
C=T_gs_200;
D=T_gs_500;
E=T_gs_1000;


I=E(1:100:1001,1:100:1001);
H=D(1:50:501,1:50:501);
G=C(1:20:201,1:20:201);
J=B(1:10:101,1:10:101);


err1=abs(max(max(A-I)))/36;
err2=abs(max(max(A-H)))/18;
err3=abs(max(max(A-G)))/9;
err4=abs(max(max(A-J)))*3;

% plot error magnitude against the mesh size

ERR=[err4,err3,err2,err1];
SIZE=[0.1,0.01,0.005,0.002];
loglog(SIZE,ERR,'r','linewidth',1);
xlabel('Mesh Size');
ylabel('Error Max Magnitude');
title('Maximum error magnitude-mesh size');
grid on





