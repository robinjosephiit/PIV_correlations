clc;clear all;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%% Reading data from file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


zones=input('\n Enter the number of zones: ');
rows=100;
columns=100;
data=load('Recon3.dat');                                %Enter the filename containing Instataneous velocity%
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

for k=1:zones 
Uinst(:,:,k)=fliplr(Uinst(:,:,k));
Vinst(:,:,k)=fliplr(Vinst(:,:,k));
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
loc_y=[5 10 20 30 40 50 60 90];


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
%      if(abs(c_vu(i,j))<0.1)
%          c_vu(i,j)=0;
%      end
    end
end



xy=dlmread('X125_xy.dat');x=xy(1:100,1);y=xy(1:100:10000,2);
disp_t=1.78;delta=5.2;

y=y+0.0; % Wall correction
y1=flip(y);
x1=flip(x);
[X,Y]=meshgrid(flip(x),flip(y)/delta);
mymap=[254 255 177
    254 197 97
    252 79 42
    220 22 29
    177 0 38];
mymap=mymap/255;
mymap2=[126 5 2
    208 47 5
    251 125 33
    238 207 58
    165 252 59
    50 242 151
    40 189 234
    71 111 231
    48 15 59];
mymap2=mymap2/255;
xx=linspace(min(x),max(x),50);yy=ones(1,50);

figure(1)
subplot(length(loc_y)/2,2,kk); 
contourf(X,Y,c_uu,[0 0.3 0.6 0.9 1.0]);ylabel('y/\delta');set(gca,'FontSize',5);hold on;
plot(xx,yy,'.-k','linewidth',1);title(num2str(kk),'fontsize',5); 
plot(x1(50),y1(loc_y(kk))/delta,'-ow','linewidth',0.5,'MarkerFaceColor','w')
plot(x1(50),y1(loc_y(kk))/delta,'-+k','linewidth',1.5); hold on;
set(gca,'FontSize',12);
caxis([0 1]);
colormap(mymap);
ch=colorbar;set(ch,'FontSize',8);

figure(2)
subplot(length(loc_y)/2,2,kk)
contourf(X,Y,c_vv,[0 0.3 0.6 0.9 1.0]);ylabel('y/\delta');set(gca,'FontSize',5);hold on;
plot(xx,yy,'.-k','linewidth',1);title(num2str(kk),'fontsize',5); 
plot(x1(50),y1(loc_y(kk))/delta,'-ow','linewidth',0.5,'MarkerFaceColor','w')
plot(x1(50),y1(loc_y(kk))/delta,'-+k','linewidth',1.5); hold on;
set(gca,'FontSize',12);
caxis([0 1]);
colormap(mymap);
ch=colorbar;set(ch,'FontSize',8);


figure(3)
subplot(length(loc_y)/2,2,kk); 
contourf(X,Y,c_uv,[0.4 0.3 0.2 0.1 0 -0.1 -0.2 -0.3 -0.4], 'LineStyle','none');ylabel('y/\delta');set(gca,'FontSize',5);hold on;
plot(xx,yy,'.-k','linewidth',1);title(num2str(kk),'fontsize',5); 
plot(x1(50),y1(loc_y(kk))/delta,'-ow','linewidth',0.5,'MarkerFaceColor','w')
plot(x1(50),y1(loc_y(kk))/delta,'-+k','linewidth',1.5); hold on;
set(gca,'FontSize',12);
caxis([-0.4 0.4]);
colormap(flipud(mymap2));
ch=colorbar;set(ch,'FontSize',8);


figure(4)
subplot(length(loc_y)/2,2,kk)

contourf(X,Y,c_vu,[0.4 0.3 0.2 0.1 0 -0.1 -0.2 -0.3 -0.4,'ShowText','on'], 'LineStyle','none');ylabel('y/\delta');set(gca,'FontSize',5);hold on;
plot(xx,yy,'.-w','linewidth',1);title(num2str(kk),'fontsize',5); 
plot(x1(50),y1(loc_y(kk))/delta,'-ow','linewidth',0.5,'MarkerFaceColor','w')
plot(x1(50),y1(loc_y(kk))/delta,'-+k','linewidth',1.5); hold on;
set(gca,'FontSize',12);
caxis([-0.4 0.4]);
colormap(flipud(mymap2));
ch=colorbar;set(ch,'FontSize',8);

end

