clc;clear all;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%% Reading data from file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sprintf('IS THE INPUT FILE datamean PRESENT IN THE FOLDER')
sprintf('ARE I,J VALUES WHILE WRITING DATA CORRECTED ACCORDING TO MESH SIZE')

zones=input('\n Enter the number of zones: ');
rows=input('\n Enter the number of rows: ');
columns=input('\n Enter the number of columns: ');
data=load('UV_clean.dat');                                %Enter the filename containing Instataneous velocity%
a=data(:,1);                                              %The first column contains u velocity%
b=data(:,2);                                              %The second column contains v velocity%  

clear data;

Uinst = zeros(rows,columns,zones);
Vinst = zeros(rows,columns,zones);
u = zeros(rows,columns,zones);
v = zeros(rows,columns,zones);
m=1;
for k=1:zones                                             %i------->X position, j-------> Y position, k------->No of zones%
    for j=1:columns
        for i=1:rows
            Uinst(i,j,k)=a(m);                                %u=u%
            Vinst(i,j,k)=b(m);                                %v=v%
            m=m+1;
        end
    end
end
clear a; clear b; clear m;
for k=1:zones 
Uinst(:,:,k)=flip(transpose(Uinst(:,:,k)));
Vinst(:,:,k)=flip(transpose(Vinst(:,:,k)));
end



umean = zeros(rows,columns);
vmean = zeros(rows,columns);

for j=1:columns
    for i=1:rows
        umean(i,j)= mean(Uinst(i,j,:));                       %Mean at a given position over all zones in u%
        vmean(i,j)= mean(Vinst(i,j,:));                       %Mean at a given position over all zones in v%

    end
end


for k=1:zones                                             
    for j=1:columns
        for i=1:rows
            u(i,j,k)=(Uinst(i,j,k)-umean(i,j));
            v(i,j,k)=(Vinst(i,j,k)-vmean(i,j));                                
        end
    end
end
%% 

loc_x=50;
loc_y=[5 10 20 30 40 50 60 70 80 90];
loc_y=flip(loc_y);

for kk=1:length(loc_y)

tempu_1=reshape(u(loc_y(kk),loc_x,:),[length(u(loc_y(kk),loc_x,:)),1]);

for i=1:100
    for j=1:100
    tempu_2=reshape(u(i,j,:),[length(u(i,j,:)),1]);
    c_tempuu=corrcoef(tempu_1,tempu_2);
    c_uu(i,j)=c_tempuu(1,2);
    if(c_uu(i,j)<0.3)
        c_uu(i,j)=0;
    end
    end
end

tempu_1=reshape(u(loc_y(kk),loc_x,:),[length(u(loc_y(kk),loc_x,:)),1]);
for i=1:100
    for j=1:100
    tempu_2=reshape(v(i,j,:),[length(v(i,j,:)),1]);
    c_tempuv=corrcoef(tempu_1,tempu_2);
    c_uv(i,j)=c_tempuv(1,2);
%     if(c_uv(i,j)<0.3)
%         c_uv(i,j)=0;
%     end
    end
end




tempv_1=reshape(v(loc_y(kk),loc_x,:),[length(v(loc_y(kk),loc_x,:)),1]);
for i=1:100
    for j=1:100
    tempv_2=reshape(v(i,j,:),[length(v(i,j,:)),1]);
    c_tempvv=corrcoef(tempv_1,tempv_2);
    c_vv(i,j)=c_tempvv(1,2);
    if(c_vv(i,j)<0.3)
        c_vv(i,j)=0;
    end
    end
end

tempv_1=reshape(v(loc_y(kk),loc_x,:),[length(v(loc_y(kk),loc_x,:)),1]);
for i=1:100
    for j=1:100
    tempv_2=reshape(u(i,j,:),[length(u(i,j,:)),1]);
    c_tempvu=corrcoef(tempv_1,tempv_2);
    c_vu(i,j)=c_tempvu(1,2);
%     if(c_vu(i,j)<0.3)
%         c_vu(i,j)=0;
%     end
    end
end



xy=dlmread('xy_dt_5p9mm.dat');x=xy(1:100,1);y=xy(1:100:10000,2);y=y/5.9;
[X,Y]=meshgrid(flip(x),flip(y));
mymap=[254 255 177
    254 197 97
    252 79 42
    220 22 29
    177 0 38];
mymap=mymap/255;
xx=linspace(min(x),max(x),50);yy=ones(1,50);

figure(1)
subplot(length(loc_y)/2,2,kk); 
contourf(X,Y,c_uu,[0 0.3 0.6 0.9 1.0]);ylabel('y/\delta');set(gca,'FontSize',5);hold on;
plot(xx,yy,'--k');title(num2str(kk))
colormap(mymap)

figure(2)
subplot(length(loc_y)/2,2,kk)
contourf(X,Y,c_vv,[0 0.3 0.6 0.9 1.0]);ylabel('y/\delta');set(gca,'FontSize',5);hold on;
plot(xx,yy,'--k');title(num2str(kk))
set(gca,'Xticklabel',[]);
colormap(mymap)

figure(3)
subplot(length(loc_y)/2,2,kk); 
contourf(X,Y,c_uv);ylabel('y/\delta');set(gca,'FontSize',5);hold on;
plot(xx,yy,'--k');title(num2str(kk))
% colormap(mymap)

figure(4)
subplot(length(loc_y)/2,2,kk)
contourf(X,Y,c_vu);ylabel('y/\delta');set(gca,'FontSize',5);hold on;
plot(xx,yy,'--k');title(num2str(kk))
set(gca,'Xticklabel',[]);
% colormap(mymap)

end

