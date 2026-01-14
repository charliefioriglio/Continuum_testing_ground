function SPV=parper(lmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%      Spherical Harmonic parallel and perperndicular values      %%%%%
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
SPV(lmax+1,2*lmax+1,2)=zeros;
CS(lmax+1,2*lmax+1)=zeros;
NRM(lmax+1,2*lmax+1)=zeros;
SPV(1,lmax+1,1)=1;
SPV(2,lmax+1,1)=1;
SPV(1,lmax+1,2)=1;
SPV(2,lmax+1,2)=0;
NRM(1,lmax+1)=sqrt((1)/4/pi*factorial(0)/factorial(0));
NRM(2,lmax+1)=sqrt(3/4/pi*factorial(1)/factorial(1));
for l=0:lmax
    for m=-lmax:lmax
        CS(l+1,m+lmax+1)=1;
        if m<0
            CS(l+1,m+lmax+1)=(-1)^m;
        end
    end
end
for l = 0:lmax
    if l>1
    SPV(l+1,lmax+1,1)=((2*(l-1)+1)*SPV(l,lmax+1,1)-(l-1)*SPV(l-1,lmax+1,1))/l;
    NRM(l+1,lmax+1)=sqrt((2*l+1)/4/pi*factorial(l)/factorial(l));
    end
    for m = 0:lmax
        if (l==m)
        SPV(l+1,m+lmax+1,2)=((-1)^(l)*dfac(2*l-1));
        SPV(l+1,-m+lmax+1,2)=((-1)^(l)*dfac(2*l-1));
        NRM(l+1,m+lmax+1)=sqrt((2*l+1)/4/pi*factorial(l-abs(m))/factorial(l+abs(m)));
        NRM(l+1,-m+lmax+1)=sqrt((2*l+1)/4/pi*factorial(l-abs(m))/factorial(l+abs(m)));
        end
        if ((l-1)==m)
        SPV(l+1,m+lmax+1,2)=0;
        NRM(l+1,m+lmax+1)=sqrt((2*l+1)/4/pi*factorial(l-abs(m))/factorial(l+abs(m)));
        SPV(l+1,-m+lmax+1,2)=0;
        NRM(l+1,-m+lmax+1)=sqrt((2*l+1)/4/pi*factorial(l-abs(m))/factorial(l+abs(m)));

        end
        if ((l-1)>m)
        SPV(l+1,m+lmax+1,2)=-(l+m-1)*SPV(l-1,m+lmax+1,2)/(l-m);
        SPV(l+1,-m+lmax+1,2)=-(l+m-1)*SPV(l-1,m+lmax+1,2)/(l-m);
        NRM(l+1,m+lmax+1)=sqrt((2*l+1)/4/pi*factorial(l-abs(m))/factorial(l+abs(m)));
        NRM(l+1,-m+lmax+1)=sqrt((2*l+1)/4/pi*factorial(l-abs(m))/factorial(l+abs(m)));        
        end
    end
end
SPV(:,:,1)=SPV(:,:,1).*NRM(:,:).*CS(:,:);
SPV(:,:,2)=SPV(:,:,2).*NRM(:,:).*CS(:,:);
end