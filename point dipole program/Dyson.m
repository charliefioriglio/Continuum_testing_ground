function [orbL,orbR,llim]=Dyson(CuO_NTOs,ps,step,X,Y,Z,tol,Dystol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%       Input from QChem Output, name of which is 'file'          %%%%% 
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%              Last Modified:    Feb 20 2023                      %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A,N,CZ,wL,wR,llim]=CuO_NTOs();
a0=0.5291172109;                                               %Bohr Radius
A=A/a0;            %convert atom coordinates from Angstroms to atomic units
ats=size(A,1);                                             %number of atoms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%         Translate the Input into Dyson Orbitals on grid         %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npts=ps+1;
for l=0:llim
f2(1,l+1)=dfac(2*l-1);
end
rmat(:,:,:,1)=(X.^2+Y.^2+(Z-A(1,3)).^2); % r^2 from each atom center
rmat(:,:,:,2)=(X.^2+Y.^2+(Z-A(2,3)).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%          Load Regular Solid Harmonics for each atom, r          %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars S
S(:,:,:,:,:,1)=Ylm(ps,llim,X,Y,Z,A(1,:));
S(:,:,:,:,:,2)=Ylm(ps,llim,X,Y,Z,A(2,:));
for a=1:ats
    o(a)=CZ{a}(1,1);                  % # of atomic orbitals for atom a 
    p(a)=CZ{a}(1,2);               % max # of Gaussian Primitives for an AO
end
nmax=max(o);
mmax=max(p);
l(nmax,ats)=zeros;                              % array for l values of AOs
for i=0:llim
    f2p(1,i+1)=dfac(2*i+1);
end
% % % v(4,2)=zeros;                          % array for number of AOs for each l
v(llim + 1,ats)=zeros;                           %modified by annie to accomodate l=4,modified by Pratichhya to accomodate # of atoms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                Set up l values and # of AOS per l               %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for a=1:ats
    row=1;
for i=1:CZ{a}(row,1)
    row=row+1;
    l(i,a)=CZ{a}(row,1);
    v(l(i,a)+1,a)=v(l(i,a)+1,a)+1;
for k=1:CZ{a}(row,2)
    row=row+1;
end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%      Read exponents and contraction coefficients. Normalize     %%%%%
%%%%%      each primitive and multiply by contraction coeffcient      %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sl=max(max(v));
c(mmax,sl,llim+1,ats)=zeros;      % placeholders for normalized contraction
ze(mmax,sl,llim+1,ats)=zeros;                        %and zeta coefficients
clearvars v
v(llim+1,ats)=zeros;
m(llim+1,sl,ats)=zeros;           % array for # of contraction coefficients
for a=1:ats
    row=1;
for i=1:CZ{a}(row,1)
    row=row+1;
    v(l(i,a)+1,a)=v(l(i,a)+1,a)+1;
    m(l(i,a)+1,v(l(i,a)+1,a),a)=CZ{a}(row,2);
    for k=1:CZ{a}(row,2)
        row=row+1;
    ze(k,v(l(i,a)+1,a),l(i,a)+1,a)=CZ{a}(row,1);
    zeta=CZ{a}(row,1);
    c(k,v(l(i,a)+1,a),l(i,a)+1,a)=CZ{a}(row,2)*sqrt((2*zeta/pi)^(1/2)*((2*zeta)^(l(i,a)+1))*2^((l(i,a)+2))/f2p(l(i,a)+1));
end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%               Normalize each basis function AO                  %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars nm
nm(llim+1,sl,ats)=zeros;          % modified by Pratichhya to accomodate # of atoms   
for a=1:ats
    for l=0:llim
        for j=1:v(l+1,a)
            i=j;
            for k=1:m(l+1,j,a)
                for kd=1:m(l+1,i,a)
nm(l+1,j,a)=nm(l+1,j,a)+c(k,j,l+1,a)*c(kd,i,l+1,a)*f2p(l+1)*sqrt(pi)/(2^(l+2))/(ze(k,j,l+1,a)+ze(kd,i,l+1,a))^(l+1.5);
                end
            end
            nm(l+1,j,a)=1/sqrt(nm(l+1,j,a));
        end
    end
end

clearvars AOSL AOSR
AOSL(npts,npts,npts)=zeros;                             %left Dyson orbital
AOSR(npts,npts,npts)=zeros;                            %right Dyson orbital
AOC=0;
for a=1:ats
    for l=0:llim %was l=0:3 modified by annie for g fn
        for j=1:v(l+1,a)
            clearvars addmat
            addmat(npts,npts,npts)=zeros;        
            for k=1:m(l+1,j,a)
                addmat(:,:,:)=addmat(:,:,:)+nm(l+1,j,a)*c(k,j,l+1,a)*exp(-ze(k,j,l+1,a).*rmat(:,:,:,a));
            end
            for ml=-l:l
                AOC=AOC+1;
                if abs(wR(AOC))>Dystol || abs(wL(AOC))>Dystol
                    AOSL=AOSL+wL(AOC)*addmat.*S(:,:,:,l+1,ml+llim+1,a);
                    AOSR=AOSR+wR(AOC)*addmat.*S(:,:,:,l+1,ml+llim+1,a);
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                    Centre the Dyson Orbitals                    %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
message1='\t normalizing and centering \n';
fprintf(message1);
orbL=AOSL;
orbR=AOSR; %%%%%%%%%%%%%%
HsqL=orbL.^2;
HsqR=orbR.^2; %%%%%%%%%%%%%%%%
clearvars normL normR
normL=1/sqrt(simp3D(HsqL,ps,step)); %normalization constant for L
normR=1/sqrt(simp3D(HsqR,ps,step)); %normalization constant for R
message1='Dyson Left norm %4f, \t Dyson Right norm %4f \n';
fprintf(message1,[normL,normR])
HZL=Z.*HsqL.*normL.^2; %integrand for <zL>=zcl
HYL=Y.*HsqL.*normL.^2; %integrand for <yL>=ycl
HXL=X.*HsqL.*normL.^2; %integrand for <xL>=xcl
HZR=Z.*HsqR.*normR.^2; %integrand for <zR>=zcr
HYR=Y.*HsqR.*normR.^2; %integrand for <yR>=ycr
HXR=X.*HsqR.*normR.^2; %integrand for <xR>=xcr
%calculate center of density offset from grid origin
zcl=round(simp3D(HZL,ps,step),tol);
xcl=round(simp3D(HXL,ps,step),tol);
ycl=round(simp3D(HYL,ps,step),tol);
zcr=round(simp3D(HZR,ps,step),tol);
xcr=round(simp3D(HXR,ps,step),tol);
ycr=round(simp3D(HYR,ps,step),tol);
message2='Initial Left Dyson centroid \t x = %.5f \t y = %.5f \t z = %.5f \t \n';
fprintf(message2,[xcl,ycl,zcl]);
message2='Initial Right Dyson centroid \t x = %.5f \t y = %.5f \t z = %.5f \t \n';
fprintf(message2,[xcr,ycr,zcr]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%            Shift to set centre of density at origin             %%%%%
%%%%%               for Left and Right Dyson orbitals                 %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Recentering The Left and Right Dyson Orbitals')    
D=griddedInterpolant(X,Y,Z,AOSL,'spline');
F=griddedInterpolant(X,Y,Z,AOSR,'spline');
orbL=normL*D(X+xcl,Y+ycl,Z+zcl);
orbR=normR*F(X+xcr,Y+ycr,Z+zcr);
end
