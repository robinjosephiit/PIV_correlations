clc;clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%% Reading data from file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zones=input('\n Enter the number of zones: ');
rows=input('\n Enter the number of rows: ');
columns=input('\n Enter the number of columns: ');

data=load('X125_clean.dat');                                %Enter the filename containing Instataneous velocity%
a=data(:,1);                                              %The first column contains u velocity%
b=data(:,2);                                              %The second column contains v velocity%  

clear data;

u = zeros(rows,columns,zones);
v = zeros(rows,columns,zones);
m=1;
for k=1:zones                                             %i------->X position, j-------> Y position, k------->No of zones%
    for j=1:columns
        for i=1:rows
            u(i,j,k)=a(m);                                %u=u%
            v(i,j,k)=b(m);                                %v=v%
            m=m+1;
        end
    end
end
clear a; clear b; clear m;

sprintf('\n \t CALCULATIONS IN PROGRESS')

%%%%%%%%%%%%%%%%%%%%%%% Averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uv    = zeros(rows,columns,zones);
umean = zeros(rows,columns);
vmean = zeros(rows,columns);
urms  = zeros(rows,columns);
vrms  = zeros(rows,columns);
rs    = zeros(rows,columns);

for j=1:columns
    for i=1:rows
        umean(i,j)= mean(u(i,j,:));                       %Mean at a given position over all zones in u%
        vmean(i,j)= mean(v(i,j,:));                       %Mean at a given position over all zones in v%
        urms(i,j) = std(u(i,j,:),1);                      %Standard deviation in u%
        vrms(i,j) = std(v(i,j,:),1);                      %Standard deviation in v%
    end
end

for k=1:zones                                             
    for j=1:columns
        for i=1:rows
            uv(i,j,k)=(u(i,j,k)-umean(i,j))*(v(i,j,k)-vmean(i,j));                                
        end
    end
end

for j=1:columns
    for i=1:rows
        rs(i,j)=mean(uv(i,j,:));                            
    end
end

fid=fopen('Average.dat','w');
fprintf(fid,'TITLE     = "2C PIV DATASET" \n');
fprintf(fid,'VARIABLES = "Umean" \n');
fprintf(fid,'"Vmean" \n');
fprintf(fid,'"Urms" \n');
fprintf(fid,'"Vrms" \n');
fprintf(fid,'"UV" \n');
fprintf(fid,'ZONE T="Program" \n');
fprintf(fid,'I=100, J=100, K=1, ZONETYPE=Ordered \n');        %I and J has to modified according to mesh size%
fprintf(fid,'DATAPACKING=POINT \n');
fprintf(fid,'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE ) \n');
for j=1:columns
    for i=1:rows
        fprintf(fid,'%.9f \t %.9f \t %.9f \t %.9f \t %.9f',umean(i,j),vmean(i,j),urms(i,j),vrms(i,j),rs(i,j));    
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end
fclose(fid);
