%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PROGRAM FOR OBTAINING POD BASED ON MEYER (2007) FORMULATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;close all;

zones   = input('\n Enter the number of zones: ');
rows    = input('\n Enter the number of rows: ');
columns = input('\n Enter the number of columns: ');
data=load('UV_clean.dat'); 
                                   
a=data(:,1);                                              %instantaneous fluctuation, u' velocity%
b=data(:,2);                                              %instantaneous fluctuation, v' velocity%
clear data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Writing data in the format specified by Meyer (2007) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=1;
u = zeros(rows*columns,zones);  v = zeros(rows*columns,zones);
for k=1:zones                                             %i------->x position, j-------> y position, k------->No of zones FOR X-Y PLANE PIV DATA%                                    
    for i=1:rows*columns 
        u(i,k) = a(m);                                
        v(i,k) = b(m);      
        m=m+1;
    end
end
clear a; clear b; clear m;
Uall=[u;v];                                             % Creating matrix will all fluctuating velocity components for each snapshot in a column
clear u; clear v;
%Active variables are Uall, rows, columns, zones.

R=Uall'*Uall;                                           % Autocovariance matrix

[eV,D]=eig(R);                                          % Columns of eV are eigenvectors, D is a diagonal matrix whose diagonal elements as eigenvalues

[L,I]=sort(abs(diag(D)));                               % Sorting of the index of the eigen values in ascending order. 
for i=1:length(I)
    L(i)=D(I(i),I(i));                                  % L(i) contains the sign of eigen value also.
end

for i=1:length(D)
eValue(length(D)+1-i) = L(i);                           % Eigenvalues sorted in descending order of their magnitude
eVec(:,length(D)+1-i) = eV(:,I(i));                     % Eigenvectors sorted in the same order
end
eValue(length(eValue))=0;                               % last eigenvalue should be zero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Basis functions and their coeffieients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:zones
    tmp=Uall*eVec(:,i);                                  % find mode
     if (norm(tmp)<= 10^-5)                               % normalize mode
         phi(:,i)=0;                                         
     else
        phi(:,i)=tmp/norm(tmp);  
                            

     end
end

C=phi'*Uall;                                              % The columns corresponds to coefficient of modes of snapshots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Writing basis functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phiu=zeros(rows*columns,zones);                           % Separating the basis functions of u and v fluctuations
phiv=zeros(rows*columns,zones);
for k=1:zones
    phiu(:,k)=phi(1:(rows*columns),k);
    phiv(:,k)=phi((rows*columns)+1:length(phi),k);
end
sprintf('\t WRITING THE BASIS FUNCTION')
N_m=input('\n Enter the number of modes to write: ');
sprintf('\t \t Is I,J for writing data corrected according to rows and columns')
fid=fopen('Basis_function.dat','w');
for k=1:N_m
    if (k==1)
        fprintf(fid,'TITLE     = "2C PIV DATASET" \n');
        fprintf(fid,'VARIABLES = "phiu" \n');
        fprintf(fid,'"phiv" \n');
        fprintf(fid,'ZONE T="Basis function" \n');                    
        fprintf(fid,'I=100, J=100, K=1, ZONETYPE=Ordered \n');            %I,J has to modified according to number of rows and columns.
        fprintf(fid,'DATAPACKING=POINT \n');
        fprintf(fid,'DT=(SINGLE SINGLE) \n');
        for i=1:(rows*columns)
            fprintf(fid,'%.9f \t %.9f',phiu(i,k),phiv(i,k));    
            fprintf(fid,'\n');
        end
    else
        fprintf(fid,'ZONE T="Basis function" \n');                    
        fprintf(fid,'I=100, J=100, K=1, ZONETYPE=Ordered \n');            %I,J has to modified according to rows and columns.
        fprintf(fid,'DATAPACKING=POINT \n');
        fprintf(fid,'DT=(SINGLE SINGLE) \n');
        for i=1:(rows*columns)
            fprintf(fid,'%.9f \t %.9f',phiu(i,k),phiv(i,k));    
            fprintf(fid,'\n');
        end
    end   
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Writing Coefficient matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sprintf('\t WRITING THE COEFFICIENT MATRIX')
sprintf('\t \t Is I and J modified according to number of zones')
fid=fopen('Coefficients.dat','w');
fprintf(fid,'TITLE     = "2C PIV DATASET" \n');
fprintf(fid,'VARIABLES = "C" \n');
fprintf(fid,'ZONE T="Coefficients" \n');                    
fprintf(fid,'I=996, J=996, K=1, ZONETYPE=Ordered \n');            %I,J has to modified according to number of zones.
fprintf(fid,'DATAPACKING=POINT \n');
fprintf(fid,'DT=(SINGLE) \n');
for j=1:zones
    for i=1:zones
        fprintf(fid,'%.9f \t',C(i,j));    
    end
    fprintf(fid,'\n');
end
sprintf('\t \t The columns contains the coefficients of modes for a snapshot')
fclose(fid);


program = input('\n ENTER PROGRAM 1 TO OBTAIN ENERGY SPECTRUM CORRESPONDING TO EIGEN VALUES ONLY: ');
if (program == 1)
    KE_t = sum(abs(eValue));                                        % KE_t:   Total mean energy.
    KE_f = (abs(eValue)/KE_t)*100;                                  % KE_f:   % fraction of kinetic energy contributed by each mode.
    C_KE_t = cumsum(KE_f);                                          % C_KE_t: Cumulative sum of the % fraction of kinetic energy.

    sprintf('\t \t Is J corrected according to number of zones')
    fid=fopen('Energy_fraction_eigen_values.dat','w');
    fprintf(fid,'TITLE     = "2C PIV DATASET" \n');
    fprintf(fid,'VARIABLES = "KE_f" \n');
    fprintf(fid,'"C_KE_t" \n');
    fprintf(fid,'"Mode num" \n');
    fprintf(fid,'ZONE T="Kinetic energy fraction" \n');                    
    fprintf(fid,'I=1, J=996, K=1, ZONETYPE=Ordered \n');            %J has to modified according to number of zones.
    fprintf(fid,'DATAPACKING=POINT \n');
    fprintf(fid,'DT=(SINGLE SINGLE SINGLE) \n');
    for i=1:zones
        fprintf(fid,'%.9f \t %.9f \t %.9f',KE_f(i),C_KE_t(i),i); 
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    clear KE_t; clear KE_f; clear C_KE_t;
end

program = input('\n ENTER PROGRAM 2 TO OBTAIN ENERGY SPECTRUM CORRESPONDING TO COEFFICIENTS OF EIGEN VALUES: ');
if (program == 2)
    E    = (1/2)*(C'.^2);                                           % C': The rows corresponds to coefficient of modes of snapshots, E: Kinetic energy per unit mass of the modes in snapshots
    KE_m = sum(E);                                                  % KE_m: It is row matrix, which contains total energy of a mode in all snapshots.
    KE_t = sum(sum(KE_m));                                          % KE_t: Total kinetic energy of all modes.
    KE_f = (KE_m/KE_t)*100;                                         % KE_f: Fraction of kinetic energy contributed by each mode.
    C_KE_t = cumsum(KE_f);                                          % C_KE_t: Cumulative sum of the kinetic energy.

    sprintf('\t \t Is J corrected according to number of zones')
    fid=fopen('Energy_fraction_coefficients.dat','w');
    fprintf(fid,'TITLE     = "2C PIV DATASET" \n');
    fprintf(fid,'VARIABLES = "KE_f" \n');
    fprintf(fid,'"C_KE_t" \n');
    fprintf(fid,'"Mode num" \n');
    fprintf(fid,'ZONE T="Kinetic energy fraction" \n');                    
    fprintf(fid,'I=1, J=996, K=1, ZONETYPE=Ordered \n');            %J has to modified according to number of zones.
    fprintf(fid,'DATAPACKING=POINT \n');
    fprintf(fid,'DT=(SINGLE SINGLE SINGLE) \n');
    for i=1:zones
        fprintf(fid,'%.9f \t %.9f \t %.9f',KE_f(i),C_KE_t(i),i); 
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    clear E; clear KE_m; clear KE_t; clear KE_f; clear C_KE_t;
end

sprintf(' Are all the eigen values positive? Are the eigen vectors real or complex? ')
sprintf(' Is the energy obtained from the eigen values and coefficients same? ')
sprintf('Check whether spectrum changes for norm values < 10^-5')


%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program = input('\n ENTER PROGRAM 3 TO RECONSTRUCT A GIVEN ZONE: ');
error = 0; RI = 0;                                                   
while (program == 3)
    clear N_M; clear N_Z; clear p; clear q; clear Recon;  
    N_M1=input('\n Starting number of mode: ');
    N_M2=input('\n Ending number of mode ');
    N_Z=1;
    for N_Z=1:zones
    for k=N_M1:N_M2
        p(:,k)=phi(:,k);   
        q(k)  =C(k,N_Z);
    end
    
    Recon = p*q';
    Recon1(:,1)=Recon(1:10000);Recon1(:,2)=Recon(10001:20000);
    dlmwrite('Recon3.dat',Recon1,'-append','Delimiter',' ');  
    end
    program=4;
%     %%%%%%%%%%%%%%%% Computing error due to reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     clear r1; clear r2;
%     r1 = norm(Uall(:,N_Z)-Recon(:,1),'fro');                          %Absolute error in Recon
%     r2 = norm(Uall(:,N_Z),'fro');
%     error(N_M,N_Z) = r1/r2;                                              %Relative error in Recon
%     sprintf('error(%d,%d) = %f',N_M,N_Z,error(N_M,N_Z))     
%     
%     %%%%%%%%%%%%%%%% Finding similarity between structures (Chen,2012)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     clear r1; clear r2;
%     r1 = Uall(:,N_Z); 
%     r2 = Recon(:,1); 
%     RI(N_M,N_Z) = sum(r1.*r2)/(norm(r1,'fro')*norm(r2,'fro'));           %RI is Relevence Index; RI=1 if structures are similar
%     sprintf('RI(%d,%d) = %f',N_M,N_Z,RI(N_M,N_Z)) 
%     
%         %%%%%%%%%%%%%%%% Writing reconstructed velocity field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Ru=zeros(rows*columns,1);                           % Separating the Recon u and v fluctuations
%     Rv=zeros(rows*columns,1);
%     Ru(:,1)=Recon(1:(rows*columns),1);
%     Rv(:,1)=Recon((rows*columns)+1:length(Recon),1);
%     
%     
%     fid=fopen('Recon_zone383_60modes.dat','w');
%         fprintf(fid,'TITLE     = "2C PIV DATASET" \n');
%         fprintf(fid,'VARIABLES = "Ru" \n');
%         fprintf(fid,'"Rv" \n');
%         fprintf(fid,'ZONE T="Recon" \n');                    
%         fprintf(fid,'I=100, J=100, K=1, ZONETYPE=Ordered \n');            %I,J has to modified according to number of rows and columns.
%         fprintf(fid,'DATAPACKING=POINT \n');
%         fprintf(fid,'DT=(SINGLE SINGLE) \n');
%         for i=1:(rows*columns)
%             fprintf(fid,'%.9f \t %.9f',Ru(i,1),Rv(i,1));    
%             fprintf(fid,'\n');
%         end
%     fclose(fid);
%     sprintf('\n Rename the reconstructed file if you want to reconstruct an other one')
%     program = input('\n Enter 3 to reconstruct again: ');
end
