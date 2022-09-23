clear all;

a=input('\n No of lines to skip at entrance:');
b=input('\n No of lines to skip in between:');
c=input('\n No of zones:');
d=input('\n No of data lines:');

fid1=fopen('X125_971zones.dat','r');                                   
fid2=fopen('X125_clean.dat','w');



for i=1:a                                
    e=fgetl(fid1);                     
end

for i=1:c
    disp(i)
    for j=1:d
        e=fgetl(fid1);
        fwrite(fid2,e);
        fprintf(fid2,'\n');
    end
    for j=1:b
        e=fgetl(fid1);
    end
end 

fclose(fid1);
fclose(fid2);