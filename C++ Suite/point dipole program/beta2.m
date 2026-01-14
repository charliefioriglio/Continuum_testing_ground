function b=beta2(CG,B,Bd,eigv,eign,SPV,lmax,kr,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%         Calculation of beta values using photoelectron          %%%%%
%%%%%                   matrix element integrals                      %%%%%
%%%%%  CG - Clebsch-Gordan Coefficients                               %%%%%
%%%%%  B - integrals for Left Dyson orbital                           %%%%%
%%%%%  Bd - integrals for Right Dyson orbitals                        %%%%%
%%%%%  eigv - point dipole eienvectors                                %%%%%
%%%%%  eign - point dipole eigenvectors                               %%%%%
%%%%%  SPV - parallel/perpendicular values for spherical harmonics    %%%%%
%%%%%  lmax - highest l in continuum                                  %%%%%
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
sumpar=0;
sumper=0;
for N=0:lmax
    for Nd=0:lmax
        for lam=-N:N
            for lamd=-Nd:Nd
                for m2=-1:1
                    for m2d=-1:1
In=B(N+1,lam+lmax+1,m2+2)*Bd(Nd+1,lamd+lmax+1,m2d+2);
if ne(In,0)
if (lam+m2)==(lamd+m2d)    
                        for l=abs(lam):lmax
                            for ld=abs(lamd):lmax
if ne(eigv(N+1,l+1,lam+lmax+1)*eigv(Nd+1,ld+1,lamd+lmax+1),0)
                                for m=-l:l
                                    for md=-ld:ld
if m==md
                        for L=abs(l-1):l+1                    
                            for Ld=abs(ld-1):ld+1
                                if L==Ld
cleb=CG(l+1,L+1,m2+2,lam+lmax+1,(m2+lam)+lmax+2)*CG(l+1,L+1,2,m+lmax+1,m+lmax+2)*CG(ld+1,Ld+1,m2d+2,lamd+lmax+1,(m2d+lamd)+lmax+2)*CG(ld+1,Ld+1,2,md+lmax+1,md+lmax+2);
                         if ne(cleb,0)
sumpar=sumpar+eigv(N+1,l+1,lam+lmax+1)*eigv(Nd+1,ld+1,lamd+lmax+1)*In*real(exp(pi/2*1i*(eign(N+1,lam+lmax+1)-eign(Nd+1,lamd+lmax+1))))*8/3*8*pi^2/(2*L+1)*cleb*2*(-1)^m*SPV(l+1,-m+lmax+1,1)*SPV(ld+1,md+lmax+1,1);
sumper=sumper+eigv(N+1,l+1,lam+lmax+1)*eigv(Nd+1,ld+1,lamd+lmax+1)*In*real(exp(pi/2*1i*(eign(N+1,lam+lmax+1)-eign(Nd+1,lamd+lmax+1))))*8/3*8*pi^2/(2*L+1)*cleb*2*(-1)^m*SPV(l+1,-m+lmax+1,2)*SPV(ld+1,md+lmax+1,2);                            
                         end
                         end                                        
                            end                                    
                        end    
end                            
                                    end
                                end
end
                            end
                        end
end
end
                    end                    
                end
            end
        end
    end
end
b=round(2*(sumpar-sumper)/(sumpar+2*sumper),4);
end
