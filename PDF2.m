% clear all
clc
format long g

%auto populate 3D sector limits based on division matrices
%Mach start,Mach end,Alt start,Alt end, Weight start, Weight end
Md= [0.74,  0.767,  0.793,  0.82];
Ad= [31000, 34333.33, 37666.67, 41000];
Wd= [27700, 35233.33, 42766.67, 50300];
sec=zeros(3,3,3,7);
for i=1:3
    for j=1:3
        for k=1:3            
            sec(i,j,k,:)=[Md(i), Md(i+1), Ad(j), Ad(j+1), Wd(k), Wd(k+1),0];
        end        
    end    
end

%--Mach PDF
pdm = makedist('Normal','mu',0.78,'sigma',0.015);
%--Alt PDF
pda = makedist('Normal','mu',36000,'sigma',1750);
%--Weight PDF
pdw = makedist('Normal','mu',39100,'sigma',4000);

%--wing area
Sref=92.5;
Mnom=0.78;
Cltars=0.4;

%for info arrays
Mmax=100;
Wmax=100;
Amax=100;
Mach=zeros(Mmax,1);
Alt=zeros(Amax,1);
rho=zeros(Amax,1);
a=zeros(Amax,1);
W=zeros(Wmax,1);

%array of normally distributed mach numbers
for i=1:Mmax
    Mach(i)=random(pdm);    
end

%array of normally distributed altitudes
for i=1:Amax    
    Alt(i)=random(pda);     
    [rho(i),a(i),~,~,~,~] = atmosphere(Alt(i));
end

%array of normally distributed weights
for i=1:Wmax
    W(i)=random(pdw);    
end

%Divide info array into sectors
for i=1:Mmax
    for j=1:Amax
        for k=1:Wmax
            
            %%%%%%%%
            for ii=1:3
                if(Mach(i)>sec(ii,1,1,1) && Mach(i)<=sec(ii,1,1,2))
                    isec=ii;
                end
            end
            for jj=1:3
                if(Alt(j)>sec(1,jj,1,3) && Alt(j)<=sec(1,jj,1,4))
                    jsec=jj;
                end
            end
            for kk=1:3
                if(W(k)>sec(1,1,kk,5) && W(k)<=sec(1,1,kk,6) )
                    ksec=kk;
                end
            end
            %%%%%%%%
            
            %increase count of this sector
            sec(isec,jsec,ksec,7)=sec(isec,jsec,ksec,7)+1;
            
        end
    end
    if(mod(i,100)==0)
        i
    end
end
%Normalize the values
sec(:,:,:,7)=sec(:,:,:,7)/(Mmax*Amax*Wmax);
figure
sec(:,:,:,7)
Mp= [0.74,  0.78,  0.82];
Ap= [31000, 36000, 41000];
Wp= [27700, 39000, 50300];
contour(Mp,Ap,sec(:,:,2,7),'ShowText','on')

%Plot sector information
fileID = fopen('oppts1.dat','wt');
for i=1:3
    %Set Up of Tecplot Block 
    fprintf(fileID,'TITLE = "Cruise Flight"\n');
    fprintf(fileID,'VARIABLES="M","A","W","Freq"\n');
    fprintf(fileID,'ZONE T="');
    fprintf(fileID,'%.f',3); 
    fprintf(fileID,'", I=   ');
    fprintf(fileID,'%.f',3);
    fprintf(fileID,', J=   ');
    fprintf(fileID,'%.f',3);
    fprintf(fileID,', F=BLOCK\n');
    
    %Print M Values
    for j=1:3
        for k=1:3
            fprintf(fileID,'%4f\n',0.5*(sec(i,j,k,1)+sec(i,j,k,2)));
        end
    end
    
    %Print A Values
    for j=1:3
        for k=1:3
            fprintf(fileID,'%4f\n',0.5*(sec(i,j,k,3)+sec(i,j,k,4)));
        end
    end
    
    %Print W Values
    for j=1:3
        for k=1:3
            fprintf(fileID,'%4f\n',0.5*(sec(i,j,k,5)+sec(i,j,k,6)));
        end
    end
        
    %Print Frequency Values
    for j=1:3
        for k=1:3
            fprintf(fileID,'%4f\n',sec(i,j,k,7));
        end
    end
end

% vel=zeros(1000,1000);
% CL=zeros(1000,1000,1000);
% for i=1:1000
%     for j=1:1000
%         vel(i,j)=Mach(i)*a(i);
%         for k=1:1000
%             CL(i,j,k)=2*9.81*W(k)/(rho(j)*Sref*vel(i,j)^2);            
%         end      
%     end
%     if(mod(i,100)==0)
%         i
%     end
% end