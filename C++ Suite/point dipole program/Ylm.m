function YLM=Ylm(ps,llim,X,Y,Z,A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%               Program to Generate Cubic Harmonics               %%%%%
%%%%%       Real, Cartesian representation Spherical Harmonics        %%%%%
%%%%%         Cubic Harmonics are output on Cartesian grids           %%%%%
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
clearvars Ph YLM
Ph=thph(X,Y,Z,2,ps);
fs={    0,    0,    0,    0, @C00,    0,    0,    0,    0;
        0,    0,    0, @C11, @S11, @C10,    0,    0,    0;
        0,    0, @S22, @S21, @C20, @C21, @C22,    0,    0;
        0, @S33, @S32, @S31, @C30, @C31, @C32, @C33,    0;
     @S44, @S43, @S42, @S41, @C40, @C41, @C42, @C43, @C44};
lsmax=size(fs,1);
lcen=5;
clearvars YLM
YLM(ps+1,ps+1,ps+1,llim+1,2*llim+1)=zeros;
cen=ps/2+1;                                                    %grid center
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%   Cubic Harmonics indexed as l,m where the l index is l + 1     %%%%%
%%%%%   and the m index is m+lmax+1                                   %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for l=0:llim
    for m=-l:l
YLM(:,:,:,l+1,m+llim+1)=fs{l+1,lcen+m}(X,Y,Z,A,Ph);
        if ne(l,0)
YLM(cen,cen,cen,l+1,m+llim+1)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c00=C00(X,Y,Z,A,~)
c00=1*prefactor(0,0);
end

function s11=S11(X,Y,Z,A,Ph)
    z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);          %Radial distance from atom position
    s11=sqrt(1-(z./r).^2).*sin(Ph).*r*prefactor(1,1);
end
function c10=C10(X,Y,Z,A,~)
    z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
c10=z./r.*r*prefactor(1,0);
end
function c11=C11(X,Y,Z,A,Ph)
    z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    c11=sqrt(1-(z./r).^2).*cos(Ph).*r*prefactor(1,1);
end

function s22=S22(X,Y,Z,A,Ph)
    z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    s22=3.*(r.^2).*(1-(z./r).^2).*sin(2*Ph)*prefactor(2,2);
end
function s21=S21(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);

s21=3*z./r.*sqrt(1-(z./r).^2).*sin(Ph).*r.^2*prefactor(2,1);
end
function c20=C20(X,Y,Z,A,~)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);

c20=0.5*(3*((z./r).^2)-1).*r.^2*prefactor(2,0);
end
function c21=C21(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
c21=3*z./r.*sqrt(1-(z./r).^2).*cos(Ph).*r.^2*prefactor(2,1);
end
function c22=C22(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
c22=3.*(r.^2).*(1-(z./r).^2).*cos(2*Ph)*prefactor(2,2);
end
function s33=S33(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    s33=15*(r.^3).*(((1-(z./r).^2)).^(3/2)).*sin(3*Ph)*prefactor(3,3);
end
function s32=S32(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    s32=15*z./r.*((1-(z./r).^2)).*(r.^3).*sin(2*Ph)*prefactor(3,2);
end
function s31=S31(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    s31=3/2*(5*(z./r).^2-1).*sqrt((1-(z./r).^2)).*(r.^3).*sin(Ph)*prefactor(3,1);
end
function c30=C30(X,Y,Z,A,~)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
 c30=0.5*(5*(z./r).^3-3*z./r).*r.^3*prefactor(3,0);   
end
function c31=C31(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    c31=3/2*(5*(z./r).^2-1).*sqrt((1-(z./r).^2)).*(r.^3).*cos(Ph)*prefactor(3,1);
end
function c32=C32(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    c32=15*z./r.*((1-(z./r).^2)).*(r.^3).*cos(2*Ph)*prefactor(3,2);
end
function c33=C33(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    c33=15*(r.^3).*((1-(z./r).^2)).^(3/2).*cos(3*Ph)*prefactor(3,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s44=S44(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    s44=sqrt(35/64)*r.^4.*((1-(z./r).^2)).^2.*sin(4*Ph);
end
function s43=S43(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    s43=-sqrt(35/8)*r.^4.*z./r.*(((1-(z./r).^2)).^(3/2)).*sin(3*Ph);
end
function s42=S42(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    s42=sqrt(5/16).*r.^4.*(7*(z./r).^2-1).*((1-(z./r).^2)).*sin(2*Ph);
end
function s41=S41(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    s41=-sqrt(5/8)*r.^4.*(7*(z./r).^3-3*z./r).*sqrt((1-(z./r).^2)).*sin(Ph);
end
function c40=C40(X,Y,Z,A,~)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    c40=1/8*(35*(z./r).^4-30*(z./r).^2+3).*r.^4;
end
function c41=C41(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    c41=sqrt(5/8)*r.^4.*(7*(z./r).^3-3*z./r).*sqrt((1-(z./r).^2)).*cos(Ph);
end
function c42=C42(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    c42=sqrt(5/16).*r.^4.*(7*(z./r).^2-1).*((1-(z./r).^2)).*cos(2*Ph);
end
function c43=C43(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    c43=sqrt(35/8)*r.^4.*z./r.*((1-(z./r).^2)).^(3/2).*cos(3*Ph);
end
function c44=C44(X,Y,Z,A,Ph)
        z=Z-A(1,3);
    r=sqrt(X.^2+Y.^2+z.^2);
    c44=sqrt(35/64)*r.^4.*((1-(z./r).^2)).^2.*cos(4*Ph);
end
end