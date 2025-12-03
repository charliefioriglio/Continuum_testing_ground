function [output]=pointdpmatrixc(lmax,dp,lam,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%% 
%%%%% Determination of eigenvectors and eigenvalues of point dipole   %%%%%
%%%%% continuum electron wave functions (P.Rev.A. 23, 632, 1981).     %%%%% 
%%%%% dp is the dipole moment, lam is the projection of the angular   %%%%%
%%%%% momentum onto the dipole axis, size represents the number of    %%%%%
%%%%% (orbital) angular momentum basis functions to include in the    %%%%%
%%%%% expansion, number is the number of eigenvectors to calculate    %%%%%
%%%%% and w determines whether to output the eigenvalues (1) or       %%%%%
%%%%% eigenvectors (2). The eigenvalues associated with the angular   %%%%%
%%%%% part of the Schrodinger eq are N(N+1) and the output values N   %%%%%
%%%%% are determined by solving the quadratic formula. At present     %%%%%
%%%%% this limits the range of applicable dipole moments to the first %%%%%
%%%%% critial limit for dipole bound states.                          %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%              Last Modified:    Jan 29 2022                      %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dipop(lmax+1,lmax+1)=zeros;
es(lmax+1,1)=zeros;
for lc=1:lmax+1
    for lr=1:lmax+1
        dipop(lc,lr)=0;
        if (lc>abs(lam))
        if (lc==lr)
        dipop(lc,lr)=(lc-1)*(lc-1+1);
        else
        dipop(lc,lr)=-2*dp*sqrt(4*pi/3)*sqrt((2*(lc-1)+1)*3*(2*(lr-1)+1)/4/pi)*(-1^lam)*wigner3j((lc-1),1,(lr-1),0,0,0)*wigner3j((lc-1),1,(lr-1),-lam,0,lam);
        end
        end
        end
end

[v,e]=eigs(dipop,lmax+1,0);

for lc=1:lmax+1
    for lr=1:lmax+1
        if ((lc-1)<abs(lam))
        v(lc,lr)=0;
        end
    end
end

if (w==1)
    es=diag((-1+sqrt(1+4*e))/2);
    output=es;
end
if (w==2)
    output=v;
end