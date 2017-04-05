clear all
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
[W,Wmax]=fracW();
% Wmax=100;
Amax=3*Wmax;
Mach=zeros(Mmax,1);
Alt=zeros(Wmax,Mmax,3);
rho=zeros(Amax,1);
a=zeros(Amax,1);
%W=zeros(Wmax,1);

%array of normally distributed mach numbers
for i=1:Mmax
    Mach(i)=random(pdm);    
end

%array of normally distributed altitudes
% for i=1:Amax    
%     Alt(i)=random(pda);     
%     [rho(i),a(i),~,~,~,~] = atmosphere(Alt(i));
% end

%array of normally distributed weights

%%%%%%initial weight%%%%%%%
j=1;
for j=1:Mmax
    %target rho*v^2
    rhoV2=2*W(i)*9.81/(Sref*Cltars);
    %bisection method to find altitude
    Alo=0;
    Ahi=41000;
    Alti=(Alo+Ahi)/2;
    res=1;
    tol=1e-5;
    while(abs(res)>tol)
        [rhoi,ai,~,~,~,~] = atmosphere(Alti);
        rs=rhoi*(Mach(j)*ai)^2;
        res=rhoV2-rs;
        if (res<0)
            Alo=Alti;
        elseif (res>=0)
            Ahi=Alti;
        end
        Alti=(Alo+Ahi)/2;
        
        %--if reaching one of the limits, exit and keep value as limit
        if(abs(Alo-Ahi)<tol)
            res=0;
        end
        
    end
    %introduce altitude variation and
    %Assign rho and a with correct altitude
    for t=1:3
        if (t==2)
            inc=1000;
        elseif (t==3)
            inc=-1000;
        else
            inc=0;
        end
        Alt(1,j,t)=Alti+inc;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%all other weights in array
for i=2:Wmax
    %     W(i)=random(pdw);
    if (W(i)~=W(i-1))
        for j=1:Mmax
            %target rho*v^2
            rhoV2=2*W(i)*9.81/(Sref*Cltars);
            
            %bisection method to find altitude
            Alo=0;
            Ahi=41000;
            Alti=(Alo+Ahi)/2;
            res=1;
            tol=1e-5;
            while(abs(res)>tol)
                [rhoi,ai,~,~,~,~] = atmosphere(Alti);
                rs=rhoi*(Mach(j)*ai)^2;
                res=rhoV2-rs;
                if (res<0)
                    Alo=Alti;
                elseif (res>=0)
                    Ahi=Alti;
                end
                Alti=(Alo+Ahi)/2;
                
                %--if reaching one of the limits, exit and keep value as limit
                if(abs(Alo-Ahi)<tol)
                    res=0;
                end
                
            end
            
            %introduce altitude variation and
            %Assign rho and a with correct altitude
            for t=1:3
                if (t==2)
                    inc=1000;
                elseif (t==3)
                    inc=-1000;
                else
                    inc=0;
                end
                Alt(i,j,t)=Alti+inc;
            end
            
        end     
        
    else        
        for j=1:Mmax
            for t=1:3
                Alt(i,j,t)=Alt(i-1,j,t);
            end         
        end
        
    end
    i
end

%Divide info array into sectors
for i=1:Wmax
    for j=1:Mmax
        for t=1:3%ensure variable altitude stays with respective weight
                        
            %%%%%%%%
            for ii=1:3
                if(Mach(j)>sec(ii,1,1,1) && Mach(j)<=sec(ii,1,1,2))
                    isec=ii;
                end
            end
            for jj=1:3
                if(Alt(i,j,t)>sec(1,jj,1,3) && Alt(i,j,t)<=sec(1,jj,1,4))
                    jsec=jj;
                end
            end
            for kk=1:3
                if(W(i)>sec(1,1,kk,5) && W(i)<=sec(1,1,kk,6) )
                    ksec=kk;
                end
            end
            %%%%%%%%
            
            %increase count of this sector
            sec(isec,jsec,ksec,7)=sec(isec,jsec,ksec,7)+1;
            
        end
    end
    i
    if(mod(i,100)==0)
        i
    end
end
%Normalize the values
sec(:,:,:,7)=sec(:,:,:,7)/(Mmax*3*Wmax);

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
