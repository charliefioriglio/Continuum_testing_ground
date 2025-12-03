%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                Matlab Implementation of ezDyson                 %%%%%
%%%%%      approach for beta values from linearly polarized           %%%%%
%%%%%      single photon detachment of anions. Treatment is           %%%%%
%%%%%      extended to include a point dipole description of          %%%%%
%%%%%      the continuum.                                             %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%     Last Modified:    March 2023 to include point dipole        %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clf
clearvars -except CG Sc orbL orbR grid0 ps0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%          User Defined Parameters for Calculation                %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qcoutput=@()CuO_NTOs;     %name of dyson orbital file
lmax=6;         %was 5                        %highest l partial wave to consider
E=[0.1,1.0];                                %first and last eKE values (eV)
Estep=0.2;                                         %step size (eV) for beta
dstart=0.639;          %first dipole moment (atomic units) value in calculation
dend=0.639;    %end dipole moment to calculate - (do not exceed critical limit)
dstep=0.1;                  %interval between successive dipole moment values
gridmax=18.897260;                              %Grid dimension (+/-) in a0
ps=100;                              %number of points on the integral grid
Dystol=1E-8;       %Max. magnitude of Dyson orbital coefficients to include
isoval=0.0010;        %isovalue at which to display parent orbital (L Dyson)
tol=10;                              %tolerance for rounding integrals to 0
redo=[1,1,1,1];                           %Flags for What to Recalculate/Do
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%   redo(1)   recalculate the Clebsch Gordan Coefficients         %%%%%
%%%%%   redo(2)   reload the parent orbital                           %%%%%
%%%%%   redo(3)   plot the parent orbital                             %%%%%
%%%%%   redo(4)   repopulate spherical harmonic arrays                %%%%%  
%%%%%   New grid dimensions require (2) and (4) to be set to 1        %%%%%
%%%%%   If lmax is changed then (1) and (4) need to be set to 1       %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%            Constant, Controls and Derived Variables             %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a0=0.5291172109;                                               %Bohr radius
npts=ps+1;                            %number of gridpoints (per dimension)
stepsz=2*gridmax/ps;                              %stepsize for integration
cen=ps/2+1;                                                    %grid center
x1=-gridmax:stepsz:gridmax;                      %3D (cubic) grid dimension
[X,Y,Z]=ndgrid(x1);
r(npts,npts,npts)=zeros;
r(:,:,:)=sqrt(X.^2+Y.^2+Z.^2);                      %(r values on the grid)
nbeta=round((E(2)-E(1))/Estep,0)+1;     %number of beta values to calculate
Eval(nbeta)=zeros;                                                         
Eval(1:nbeta)=E(1)+([1:nbeta]-1)*Estep;                                    
qval(1:nbeta)=sqrt(2*Eval(1:nbeta)/27.211);                                
ndps=round((dend-dstart)/dstep,0)+1;
betas(nbeta,ndps)=zeros;                  %array for calculated beta values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%         Recalculate the Clebsch-Gordan coefficients             %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  redo(1)==1 || exist('CG','var')==0
    clearvars CG
    CG(lmax+1,lmax+1,3,2*lmax+1,2*lmax+3)=zeros;
    CG=clebarray1(lmax);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%      Spherical Harmonics and Parallel/Perpendicular Values      %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if redo(4)==1 ||  exist('Sc','var')==0 || ( ne(ps,ps0)) ...
   || ne(ps,ps0) || (exist('grid0','var')==0 && ne(grid0,gridmax))
message1='Loading (Complex) Spherical Harmonics (up to l = %d )\n';
fprintf(message1,lmax);
clearvars Sc
Sc(npts,npts,npts,lmax+1,2*lmax+1)=zeros;
Sc=Ylmc(ps,lmax,X,Y,Z);
end
message1='Generating Parallel and Perpendicular Values for Spherical Harmonics (up to l = %d )\n';
fprintf(message1,lmax);
SPV(lmax+1,2*lmax+1,2)=zeros;
SPV=parper(lmax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                    Reload The Dyson Orbital                     %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if redo(2)==1 || exist('orbL','var')==0 || exist('orbR','var')==0 ...
 || ( ne(ps,ps0)) ...
 || ne(ps,ps0) || (exist('grid0','var')==0 && ne(grid0,gridmax))
message1='Loading the L & R Dyson Orbitals \n \t including all Dyson orbital coefficients > %1.2e \n';
fprintf(message1,Dystol);
clearvars orbL orbR
orbL(npts,npts,npts)=zeros;
orbR(npts,npts,npts)=zeros;
[orbL,orbR]=Dyson(qcoutput,ps,stepsz,X,Y,Z,tol,Dystol);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%             Plot the Left and Right Dyson Orbitals              %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if redo(3)==1
figplot(X,Y,Z,orbL,gridmax/2,0.015,2,1,1);
xlabel('x','FontSize',12,'FontWeight','bold','Color','r');
ylabel('y');
axis on;
hold on
figplot(X,Y,Z,orbR,gridmax/2,0.015,2,2,1);
drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%   Save parameters to determine if next calculation has          %%%%%
%%%%%   same grid size or # of points                                 %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid0=gridmax;
ps0=ps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                Calculation for Each dipole moment               %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dind=0;
for dp=dstart:dstep:dend
    dind=dind+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%      Calculate Eigenvalues and Eigenvectors for Point Dipole      %%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
disp('Setting Up and Normalizing Eigenfunctions and Eigenvalues');
eign(lmax+1,2*lmax+1)=zeros;
eigv(lmax+1,lmax+1,2*lmax+1)=zeros;
for lam=-lmax:lmax
    eign(:,lam+lmax+1)=pointdpmatrixc(lmax,dp,lam,1);
    eigv(:,:,lam+lmax+1)=pointdpmatrixc(lmax,dp,lam,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%               Normalize the Continuum Functions                 %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for lam=-lmax:lmax
for N=abs(lam):lmax
    norm=0;
    parfor l=abs(lam):lmax
        norm=norm+eigv(N+1,l+1,lam+lmax+1).^2;
    end
    for l=abs(lam):lmax
        eigv(N+1,l+1,lam+lmax+1)=eigv(N+1,l+1,lam+lmax+1)/sqrt(norm);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                Load Continuum Functions onto grids              %%%%%
%%%%%                 (Required every time eKE changes)               %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nbeta
    clearvars qrgrid F
    F(npts,npts,npts,lmax+1,2*lmax+1)=zeros;
    Fd(npts,npts,npts,lmax+1,2*lmax+1)=zeros;
    message2='Calculating Continuum Functions for k = %.4f a.u., \t E = %.2f eV';
    fprintf(message2,[qval(i), (E(1)+((i-1)*Estep))])
for N=0:lmax
    for lam=-N:N
 F(:,:,:,N+1,lmax+lam+1)=continuum(N,lam,lmax,eign,eigv,qval(i),r,Sc,npts);
 Fd(:,:,:,N+1,lmax+lam+1)=((-1)^lam)*continuum(N,-lam,lmax,eign,eigv,qval(i),r,Sc,npts);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%                Set up and Evaluate Integrals on grid              %%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Integrals(lmax+1,2*lmax+1,3,nbeta+1)=zeros; 
Integralsd(lmax+1,2*lmax+1,3,nbeta+1)=zeros;
for lam=-lmax:1:lmax
    for N=abs(lam):lmax
        clearvars Integrand Integrandd
        Integrand(npts,npts,npts,3)=zeros;
        Integrandd(npts,npts,npts,3)=zeros;
        Integrand(:,:,:,1)=Sc(:,:,:,2,lmax+0).*orbL.*F(:,:,:,N+1,lmax+lam+1);
        Integrand(:,:,:,2)=Sc(:,:,:,2,lmax+1).*orbL.*F(:,:,:,N+1,lmax+lam+1);
        Integrand(:,:,:,3)=Sc(:,:,:,2,lmax+2).*orbL.*F(:,:,:,N+1,lmax+lam+1);
        Integrandd(:,:,:,1)=-Sc(:,:,:,2,lmax+2).*orbR.*Fd(:,:,:,N+1,lmax+lam+1);
        Integrandd(:,:,:,2)=+Sc(:,:,:,2,lmax+1).*orbR.*Fd(:,:,:,N+1,lmax+lam+1);
        Integrandd(:,:,:,3)=-Sc(:,:,:,2,lmax+0).*orbR.*Fd(:,:,:,N+1,lmax+lam+1);
        Integrals(N+1,lam+lmax+1,1,i+1) =round(simp3D(Integrand(:,:,:,1),npts-1,stepsz),tol);
        Integrals(N+1,lam+lmax+1,2,i+1) =round(simp3D(Integrand(:,:,:,2),npts-1,stepsz),tol);
        Integrals(N+1,lam+lmax+1,3,i+1) =round(simp3D(Integrand(:,:,:,3),npts-1,stepsz),tol);
        Integralsd(N+1,lam+lmax+1,1,i+1)=round(simp3D(Integrandd(:,:,:,1),npts-1,stepsz),tol);
        Integralsd(N+1,lam+lmax+1,2,i+1)=round(simp3D(Integrandd(:,:,:,2),npts-1,stepsz),tol);
        Integralsd(N+1,lam+lmax+1,3,i+1)=round(simp3D(Integrandd(:,:,:,3),npts-1,stepsz),tol);     
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                     Calculate beta value                        %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betas(i,1)=Eval(i);
betas(i,dind+1)=beta2(CG,Integrals(:,:,:,i+1),Integralsd(:,:,:,i+1),eigv,eign,SPV,lmax);
fprintf(repmat('\b',1,length(message2)+1))
message1='Dipole moment %.2f a.u., \t eKE %.2f eV, \t beta = %.4f \n';
fprintf(message1,[dp,betas(i,1),betas(i,dind+1)]);
end
end
