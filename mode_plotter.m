clc;clear all;
temp=importdata('Energy_fraction_eigen_values.dat');
modes=temp.data();

figure(1)
hold on;
subplot(1,2,1)
semilogx(modes(:,3),modes(:,1),'-or');
xlabel('Mode number');

legend('Mode energy');

subplot(1,2,2)
semilogx(modes(:,3),modes(:,2),'-dk');
set(gca,'FontSize',12);
legend('Cumulative energy');
xlabel('Mode number');