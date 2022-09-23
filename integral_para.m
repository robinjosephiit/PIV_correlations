clc;clear all;close all;
avg=dlmread('X125_Average.dat');
xy=dlmread('X125_xy.dat');
y=xy(1:100:10000,2);
U=avg(50:100:10000,1);


y=flip(y);U=flip(U);
ff=fit(y(1:4),U(1:4),'poly1');
y=y+ff.p2/ff.p1;

[disp_t,mom_t,H,Uinf,Delta]=inte(y,U);