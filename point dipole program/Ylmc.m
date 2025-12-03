function YLM=Ylmc(ps,llmax,X,Y,Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%            Program to Generate Complex Harmonics                %%%%%
%%%%%        Complex Harmonics are output on Cartesian grids          %%%%%
%%%%%                                                                 %%%%%
%%%%%  Complex Harmonics indexed as l,m where the l index is l + 1    %%%%%
%%%%%  and the m index is m+lmax+1.                                   %%%%%
%%%%%  Uses module thph                                               %%%%%
%%%%%  Spherical Harmonics are normalized over theta and phi          %%%%%
%%%%%  and include Condon-Shortley phase                              %%%%%
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
clearvars Ph Th YLM l m
Th=thph(X,Y,Z,1,ps);
Ph=thph(X,Y,Z,2,ps);
PHI(ps+1,ps+1,ps+1,2*llmax+1)=zeros;                        % m*Phi on grid
for m=-llmax:llmax
    PHI(:,:,:,llmax+m+1)=exp(1i*m*Ph);
end
YLM(ps+1,ps+1,ps+1,llmax+1,2*llmax+1)=zeros;% Array for spherical harmonics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%  Calculation starts from generating function for Associated     %%%%%
%%%%%  Legendre functions.                                            %%%%%
%%%%%                                                                 %%%%%
%%%%%   m        m  l    2 (m/2)   l            (k-m)( l) ((l+k-1)/2) %%%%%
%%%%%  P (x)=(-1)  2 (1-x )     SUM  k!/(k-m)! x(    ( k )(    l    ) %%%%%
%%%%%   l                          k=m                                %%%%%
%%%%%                                                                 %%%%%
%%%%%  (l)   and ((l+k-1)/2)                                          %%%%%  
%%%%%  (k )      (    l    )   are binomial coefficients.             %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l=0:llmax
    for m=0:l
        for k=m:l
           bc1=0;
           if le(k,l)
              bc1=nchoosek(l,k);    
           end
           if bc1==0
              continue
           end
            k2=(l+k-1)/2;
            if rem((l+k-1),2)==0
                bc2=0;
                if le(l,k2)
                    bc2=nchoosek(k2,l);
                end
            end
            if ne(rem(l+k-1,2),0)
               bc2=1;
               if l>0
               for a=0:l-1
               bc2=bc2*(k2-a);
               end
               bc2=bc2/factorial(l);
               end
            end
        if ne(bc1*bc2,0) 
YLM(:,:,:,l+1,m+llmax+1)=YLM(:,:,:,l+1,m+llmax+1)+bc1*bc2*factorial(k)/factorial(k-m).*(cos(Th)).^(k-m);
        end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%  Construct remaining part of spherical harmonics and include    %%%%%
%%%%%  negative m's.                                                  %%%%%
%%%%%                                                                 %%%%%
%%%%% m                                               m          imph %%%%%
%%%%%Y (th,ph)=sqrt((2*l+1)/(4pi))sqrt((l-m)!/(l+m)!)P (cos(th))e     %%%%%
%%%%% l                                               l               %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
if m==0
YLM(:,:,:,l+1,m+llmax+1)=(2^l)*sqrt((2*l+1)/4/pi)*sqrt(factorial(l-m)/factorial(l+m))*YLM(:,:,:,l+1,m+llmax+1);
end
if m>0
YLM(:,:,:,l+1,-m+llmax+1)=(2^l)*sqrt((2*l+1)/4/pi)*sqrt(factorial(l-m)/factorial(l+m))*YLM(:,:,:,l+1,m+llmax+1).*(sin(Th)).^m.*PHI(:,:,:,llmax-m+1);    
YLM(:,:,:,l+1,m+llmax+1)=((-1)^m)*(2^l)*sqrt((2*l+1)/4/pi)*sqrt(factorial(l-m)/factorial(l+m))*YLM(:,:,:,l+1,m+llmax+1).*(sin(Th)).^m.*PHI(:,:,:,llmax+m+1);
end
    end
end
end
