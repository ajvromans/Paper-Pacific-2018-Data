function timeend  = NewGypsum4P(phii10,phii30,del,epsilon,KK,phirct,JJ,phires,AA,gamma,simulation)
%--------------------------------------------------------------------------
%This code is designed to execute the Rothe method for the PDE model as 
%written in Vromans, van de Ven & Muntean (2017) titled "EXISTENCE OF WEAK 
%SOLUTIONS FOR A PSEUDO-PARABOLIC SYSTEM COUPLING CHEMICAL REACTIONS, 
%DIFFUSION AND MOMENTUM EQUATIONS", but applied to the Gypsum model as 
%introduced by Fons van de Ven.
%--------------------------------------------------------------------------
%Version 1.4 (BVP5C implementation + Trapezium integrator)
%Upgrade vs Version 1.0: using sol structure + deval in BVP4C ODEFUN.
%Upgrade vs Version 1.1: put v eqn out of BVP4C + corrected errors
%Upgrade vs Version 1.2: solve sequentially. First w eqns, then phi.
%Upgrade vs Version 1.3: Replace BVP4C by BVP5C. Evaluate derivatives.
%--------------------------------------------------------------------------
%METHOD OF NUMERICAL CALCULATION:
%First we test whether or not the initial conditions are physical and 
%whether or not the parameters yield a hyperbolic system. Then we use the 
%BVP5C algorithm to calculate the initial function of v(z). This function 
%is then used to test whether or not V^2/Delta t is large enough.
%Then we fill the necessary recursive variables with the new values to be 
%able to start the recursion.
%
%The main part of the programme is an iterative application of the BVP4C
%algorithm to obtain (almost) all 0th, 1st  and 3rd order derivatives in 
%space. Then we add the new functions to the plotting algorithm.
%We apply this method separately on the w en phi equations to prevent any
%unwanted problems like singular matrices.
%First we integrate the w eqns with BVP5C, then the phi eqns with the w
%solution filled in at the correct spots.
%
%Lastly we redo the physical checks and the V^2/Delta t check to determine
%the optimal stopping time. Additionally, we check if the height function 
%is decreasing.Furthermore, we fill the recursive variables with their new 
%values if necessary.
%--------------------------------------------------------------------------
%We use globals to transfer fixed values at a timestep between different
%functions. Usually this are parameters, but sometimes it must be applied
%to vectors and matrices as well in order to keep the structure of the in
%build functions of MATLAB.

%phii10  = 0.3301;
%phii30  = 0.33;
%del    = 0;
%epsilon    = 0;
%KK     = 0;
%phirct = 0;
%JJ     = 0;
%phires = 0;
%AA     = 0;
%gamma  = 0;

global delta1 delta2 DT eps chi1 chi2 gamma1 gamma2 Ehat1 Ehat2 phiS phiT
global kappa1 kappa3 J2 J3 phi2res phi3res Veloc TimeEND tol N yinit FFF
global xx phi10 phi30 phi20 A1 A2 FF CC Wstart Wprev Vtotal solw winit

%Intrinsic parameters of the system
delta1 = 1.00*10^(del);
delta2 = 1.00*10^(del);
eps = 0.0014*10^(0.5*epsilon);
kappa1 = 23.00*10^(KK);
kappa3 = 13.50*10^(KK);
gamma1 = 0.50*10^(gamma);
gamma2 = 0.50*10^(gamma);
chi1 = 01.00; %By definition either chi1 or chi2 is equal to 1
chi2 = 01.00; %and the other smaller or equal.
Ehat1 = 0.380; %By definition either Ehat1 or Ehat2 is equal to 1 
Ehat2 = 1.000; %and the other smaller or equal.
phiS = 1.00*10^(-0.5*phirct);
phiT = 0.00;
J2 = 0.4*10^(JJ);
J3 = 2.0*10^(JJ);
phi2res = 1.0*10^(-0.5*phires);
phi3res = 1.0*10^(-0.5*phires);
A1 = 0.50*10^(AA);
A2 = 0.50*10^(AA);
Wstart = 0.0;

%Analytical control parameters
Veloc = 1000000.00; %The arbitrary upper bound to L^2H^1 of v
tol = 0.00001; %The value of phimin. Basically the resolution of 0.

%Numerical parameters
DT = 0.001; %The timestep
TimeEND = 0.500; %Arbitrary maximum time to ensure programme ending
N = 300; %N+1 = number of points in position space z. DZ = 1/N.

%Initial condition parameters
phi10 = 1/3;
phi30 = 1/3;
phi20 = 1-phi10-phi30;
Wprev = Wstart;

%The total amount of integration steps must be an integer
M = round(TimeEND/DT);

%The output surface plot initialisation
 Z1 = zeros(M+1,N+1);    %phi1     (volume fraction of gypsum)
 Z2 = zeros(M+1,N+1);    %phi2     (volume fraction of kalk)
 Z3 = zeros(M+1,N+1);    %phi3     (volume fraction of sulphuric acid)
 Z4 = zeros(M+1,N+1);    %v        (displacement velocity of acid)
 Z5 = zeros(M+1,N+1);    %w1       (displacement of gypsum)
 Z6 = zeros(M+1,N+1);    %w2       (displacement of kalk)
 Z7 = zeros(M+1,N+1);    %dzphi1   (1st derivative)
 Z8 = zeros(M+1,N+1);    %dzphi3   (1st derivative)
 Z9 = zeros(M+1,N+1);    %dzv      (1st derivative)
Z10 = zeros(M+1,N+1);    %dzw1     (1st derivative)
Z11 = zeros(M+1,N+1);    %dzw2     (1st derivative)
Z12 = zeros(M+1,N+1);    %dzdzphi1 (3rd derivative)
Z13 = zeros(M+1,N+1);    %dzdzphi3 (3rd derivative)
Z14 = zeros(M+1,N+1);    %dzdzw1   (3rd derivative)
Z15 = zeros(M+1,N+1);    %dzdzw2   (3rd derivative)
Z16 = zeros(M+1,N+1);    %Freact   (Chemical reaction speed)
Z17 = zeros(M+1,1);      %W        (modified height function)
Z18 = zeros(M+1,1);      %Vtotal   (integrated L^2H^1 velocity)

%The domain of the system is a one-dimensional linear vector space
xx = linspace(0,1,N+1);

%Check whether the initial conditions are physical.
%if (((phi10 < tol)||(phi30 < tol))||((1-phi10-phi30) < tol))
if (min([phi10,phi30,1-phi10-phi30]) < tol)
    disp('Initial concentration are within tolerance level of vanishing or having unphysical values. Program stopped.');
    timeend = 0;
    return
else
    %Check whether the initial conditions are in Steinmetz Solid and thus
    %that all future volume fractions are physical.
    tolerance = (1-tol)*(1-tol);
    if (max([phi10*phi10+phi30*phi30,phi10*phi10+phi20*phi20,phi20*phi20+phi30*phi30])>tolerance/4)
        disp('Initial concentration outside Steinmetz Solid. Program continues.');
        disp(phi10*phi10+phi30*phi30>tolerance/4);
        disp(phi10*phi10+phi20*phi20>tolerance/4);
        disp(phi20*phi20+phi30*phi30>tolerance/4);
    end
end

%Check whether the system parameters give a hyperbolic system.
if ((gamma1/chi1+tol<=0)||(gamma2/chi2+tol<=0))||((Ehat1*DT>gamma1+tol)||(Ehat2*DT>gamma2+tol))
    disp('Either gamma close to wrong sign or system is close to hyperbolic for this timestep. Program stopped.');
    timeend = 0;
    return
end

%-------------------------------------------------------------------------
%Determining the correct constants for the integration.
FFP = (phiS-phi10)*(phiS>phi10)*(phi30-phiT)*(phiT<phi30);
FF = 1/phi30*FFP;
FFF = FFP*ones(1,N+1);
CC = J2/phi30*(phi2res-phi20)*(phi2res>phi20);
%Defining the differential system and BC for the BVP4c integration.
function dydx = odev(x,y)
    dydx = [[1,0,0,0];
            [0,1,0,0];
            [0,0,gamma1*(1-phi10)/chi1,-gamma2*phi10/chi1];
            [0,0,-gamma1*phi20/chi2,gamma2*(1-phi20)/chi2]]\[y(3);
        y(4);
        (1-phi20)/phi30*y(1)+phi20/phi30*y(2)-FF*x-CC;
        chi2*(1-phi10)/phi30*y(2)+chi2*phi10/phi30*y(1)-FF*x-CC];
end 
function res = bcv(ya,yb)
    res = [ya(1);
        ya(2)-phi30/phi20*CC;
        yb(3)+A1*(FF+CC-(1-phi20)/phi30*yb(1)-phi20/phi30*yb(2)+J3*(phi3res-phi30)...
            *(phi3res>phi30)/phi30);
        yb(4)+A2*(FF+CC-(1-phi10)/phi30*yb(2)-phi10/phi30*yb(1)+J3*(phi3res-phi30)...
            *(phi3res>phi30)/phi30)];
end
%Determination of the auxilliary velocities with the help BVP4c.
solinit = bvpinit(xx,[1;1;1;1]);
options = bvpset();
solv = bvp5c(@odev,@bcv,solinit,options);
vv = deval(solv,xx);
%Determination of the initial condition of v with the auxilliary velocities
v0 = FF*xx+CC-phi10/phi30*vv(1,:)-phi20/phi30*vv(2,:);
v0 = v0-v0(1)*ones(1,N+1);
v1 = FF-phi10/phi30*vv(3,:)-phi20/phi30*vv(4,:);
%Determination of L2H1 norm of v with trapezium integration.
Vcorr = v0(1)*v0(1)+v0(N+1)*v0(N+1)+v1(1)*v1(1)+v1(N+1)*v1(N+1);
Vtotal = DT*((v0*v0')+(v1*v1')-Vcorr/2)/N;
%Check whether the chosen upper bound for the velocity is large enough or
%that the integration time step is too small for the velocity upper bound.
if (Vtotal>Veloc^2)
    disp(sprintf('V or 1/DT is too small. V^2 must be at least ||dzv||^2DT = %f.',Vtotal));
    timeend = 0;
    return
end
%-------------------------------------------------------------------------
%Creating the initial conditions by inserting them in the vector containing
%the previous time step values.
yprev = zeros(15,N+1);
yprev( 1,:) = phi10*ones(1,N+1);    %phi1
yprev( 2,:) = phi30*ones(1,N+1);    %phi3
yprev( 3,:) = v0;                   %v
yprev( 4,:) = zeros(1,N+1);         %w1
yprev( 5,:) = zeros(1,N+1);         %w2
yprev( 6,:) = zeros(1,N+1);         %dzphi1
yprev( 7,:) = zeros(1,N+1);         %dzphi3
yprev( 8,:) = v1;                   %dzv
yprev( 9,:) = zeros(1,N+1);         %dzw1
yprev(10,:) = zeros(1,N+1);         %dzw2
yprev(11,:) = zeros(1,N+1);         %dzdzphi1
yprev(12,:) = zeros(1,N+1);         %dzdzphi3
yprev(13,:) = zeros(1,N+1);         %dzdzw1
yprev(14,:) = zeros(1,N+1);         %dzdzw2
yprev(15,:) = (phiS-phi10)*(phiS>phi10)*(phi30-phiT)*(phiT<phi30)*ones(1,N+1);  %Freact
winit = [yprev(4,1);
         yprev(5,1);
         yprev(9,1);
         yprev(10,1)];   %create w initial guess vector for bvp4c
yinit = [yprev(1,1);
         yprev(2,1);
         yprev(6,1);
         yprev(7,1);];   %create phi initial guess vector for bvp4c
%-------------------------------------------------------------------------     
%Insert the initial conditions in the sol structure called start. In this
%way we can use the deval function that allows for arbitrary interpolation
%evaluations. This type of interpolation is done automatically by bvp4c,
%but not for our previous time step terms. Hence we need to make the format
%such that it can be done again automatically by bvp4c.
start.solver = 'bvp4c'; %must be bvp4c, because deval fails when bvp5c
%Attach domain vector.
start.x = xx;
%Attach previous function
start.y = [yprev( 1,:);
            yprev( 2,:);
            yprev( 3,:);
            yprev( 4,:);
            yprev( 5,:);
            yprev( 6,:);
            yprev( 7,:);
            yprev( 9,:);
            yprev(10,:)];
%Attach previous function derivative
start.yp = [yprev( 6,:);
           yprev( 7,:);
           yprev( 8,:);
           yprev( 9,:);
           yprev(10,:);
           yprev(11,:);
           yprev(12,:);
           yprev(13,:);
           yprev(14,:)];
%Attach the statistics of the structure as BVP4c does in output.
start.stats = struct('nmeshpoints',1+N,'maxres',0.0001,...
                   'nODEevals',1,'nBCevals',1);

%Insert the initial conditions in the data matrices for plotting at the
%end.
 Z1(1,:)=yprev( 1,:);            %phi1     (volume fraction of gypsum)
 Z2(1,:)=1-yprev(1,:)-yprev(2,:);%phi2     (volume fraction of kalk)
 Z3(1,:)=yprev( 2,:);            %phi3     (volume fraction of sulphuric acid)
 Z4(1,:)=yprev( 3,:);            %v        (displacement velocity of acid)
 Z5(1,:)=yprev( 4,:);            %w1       (displacement of gypsum)
 Z6(1,:)=yprev( 5,:);            %w2       (displacement of kalk)
 Z7(1,:)=yprev( 6,:);            %dzphi1   (1st derivative)
 Z8(1,:)=yprev( 7,:);            %dzphi3   (1st derivative)
 Z9(1,:)=yprev( 8,:);            %dzv      (1st derivative)
Z10(1,:)=yprev( 9,:);            %dzw1     (1st derivative)
Z11(1,:)=yprev( 10,:);           %dzw2     (1st derivative)
Z12(1,:)=yprev( 11,:);           %dzdzphi1 (3rd derivative)
Z13(1,:)=yprev( 12,:);           %dzdzphi3 (3rd derivative)
Z14(1,:)=yprev( 13,:);           %dzdzw1   (3rd derivative)
Z15(1,:)=yprev( 14,:);           %dzdzw2   (3rd derivative)
Z16(1,:)=yprev( 15,:);           %Freact   (Chemical reaction speed)
Z17(1,1)=Wprev;                  %W        (modified height function)
Z18(1,1)=Vtotal;                 %Vtotal   (integrated L^2H^1 velocity)
%-------------------------------------------------------------------------
%We introduce the functions for the system equations and the BC.
%First the functions of the w eqns.
function dydx = odew(x,y)
    %Obtain the automatic interpolation evaluation of previous time step.
    [ystar, ypstar] = deval(start,x);
    %Specify the ODE of the w equations.
    dydx = [[1,0,0,0];
            [0,1,0,0];
            [0,0,gamma1*(1-ystar(1,1))+Ehat1*DT,-gamma2*ystar(1,1)];
            [0,0,-gamma1*(1-ystar(1,1)-ystar(2,1)),gamma2*(ystar(1,1)+ystar(2,1))+Ehat2*DT]]\[y(3);
         y(4);
         gamma1*(1-ystar(1,1))*ypstar(8,1)-gamma2*ystar(1,1)*ypstar(9,1)+...
            chi1*(y(1)-ystar(4,1)-DT*ystar(3,1))+Ehat1*DT*(ystar(6,1)*ystar(8,1)+ystar(1,1)*ypstar(8,1))+...
            Ehat2*DT*(ystar(6,1)*ystar(9,1)+ystar(1,1)*ypstar(9,1));
         -gamma1*(1-ystar(1,1)-ystar(2,1))*ypstar(8,1)+gamma2*(ystar(1,1)+ystar(2,1))*ypstar(9,1)+...
            chi2*(y(2)-ystar(5,1)-DT*ystar(3,1))+Ehat1*DT*((1-ystar(1,1)-ystar(2,1))*ypstar(8,1)-(ystar(6,1)+ystar(7,1))*ystar(8,1))+...
            Ehat2*DT*((1-ystar(1,1)-ystar(2,1))*ypstar(9,1)-(ystar(6,1)+ystar(7,1))*ystar(9,1))];
end 
function res = bcw(ya,yb)
    %Specify the BC of the w equations
    res = [ya(1)-yprev(4,1);
        ya(2)-yprev(5,1)-J2*DT*(phi2res-1+yprev(1,1)+yprev(2,1))*(phi2res+...
            yprev(1,1)+yprev(2,1)>1)/(1-yprev(1,1)-yprev(2,1));
        yb(3)-A1*(yb(1)-Wprev-yprev(3,1+N)*DT-J3*DT*(phi3res-...
            yprev(2,1+N))*(phi3res>yprev(2,1+N))/yprev(2,1+N));
        yb(4)-A2*(yb(2)-Wprev-yprev(3,1+N)*DT-J3*DT*(phi3res-...
            yprev(2,1+N))*(phi3res>yprev(2,1+N))/yprev(2,1+N));];
end
%Then the functions of the phi eqns.
function dydx = odey(x,y)
    %Create values of special constants
    k1  = eps*kappa1/delta1;
    k3  = eps*kappa3/delta2;
    %Obtain the automatic interpolation evaluation of previous time step.
    [ystar, ypstar] = deval(start,x);
    %Obtain the automatic interpolation evaluation of w of current time step.
    [ww, wwx] = deval(solw,x);
    %Determine special function (the chemical reaction speed) 
    FF1 = (phiS-ystar(1,1))*(phiS>ystar(1,1))*(ystar(2,1)-phiT)*...
        (ystar(2,1)>phiT);
    %Specify the ODE of the phi equations.
    dydx = [y(3);
         y(4);
         (y(1)-ystar(1,1)+eps*(ystar(1,1)*(wwx(1)-ypstar(4,1))+ystar(6,1)*(ww(1)-ystar(4,1))))/(DT*delta1)-k1*FF1;
         (y(2)-ystar(2,1))/(DT*delta2)+eps*(ystar(2,1)*ypstar(3,1)+ystar(7,1)*ystar(3,1))/delta2+k3*FF1;];
end 
function res = bcy(ya,yb)
    %Specify the BC of the phi equations
    res = [ya(3);ya(4);yb(3);yb(4)];
end
%----------------------------------------------------------------------
%Starting the system integration
%Execute for number of time steps that need integration
for t = 1:M
    disp(t);
    %----------------------------------------------------------------------
    %The actual integration routine for w eqns
    solwinit = bvpinit(xx,winit);
    options = bvpset();
    solw = bvp5c(@odew,@bcw,solwinit,options);
    [yw, ywx] = deval(solw,xx);%evaluate the solution at our preferred grid
    %----------------------------------------------------------------------
    %The actual integration routine for phi eqns
    solyinit = bvpinit(xx,yinit);
    options = bvpset();
    soly = bvp5c(@odey,@bcy,solyinit,options);
    [yy, yyx] = deval(soly,xx);%evaluate the solution at our preferred grid
    %----------------------------------------------------------------------
    %Create special constants and matrices
    C = J2*(phi2res+yprev(1,1)+yprev(2,1)-1).*(phi2res+yprev(1,1)+yprev(2,1)>1)./(1-yprev(1,1)-yprev(2,1));
    MM = toeplitz([1/2,ones(1,N)],[1/2,zeros(1,N)]);%Trapezium intergrator
    MM(:,1) = MM(:,1)-1/2;                          %matrix initialisation
    %Determining reaction function F from phi functions.
    FFF = (phiS-yy(1,:)).*(phiS>yy(1,:)).*(yy(2,:)-phiT).*(yy(2,:)>phiT);
    %Determining v from the w function.
    vy = ((MM*FFF')'/N-yprev(1,:).*(yw(1,:)-yprev(4,:))/DT...
        -(1-yprev(1,:)-yprev(2,:)).*(yw(2,:)-yprev(5,:))/DT)./yprev(2,:);
    vy = vy-vy(1)*ones(1,N+1);
    vyx = (FFF-vy.*yprev(7,:)...
        -yprev(1,:).*(ywx(1,:)-yprev(9,:))/DT...
        -yprev(6,:).*(yw(1,:)-yprev(4,:))/DT...
        -(1-yprev(1,:)-yprev(2,:)).*(ywx(2,:)-yprev(10,:))/DT...
        +(yprev(6,:)+yprev(7,:)).*(yw(2,:)-yprev(5,:))/DT)./yprev(2,:);
    %----------------------------------------------------------------------
    %Updating the vector for previous timestep with new values
    yprev( 1,:) = yy(1,:);    %phi1
    yprev( 2,:) = yy(2,:);    %phi3
    yprev( 3,:) = vy;         %v
    yprev( 4,:) = yw(1,:);    %w1
    yprev( 5,:) = yw(2,:);    %w2
    yprev( 6,:) = yyx(1,:);   %dzphi1
    yprev( 7,:) = yyx(2,:);   %dzphi3
    yprev( 8,:) = vyx;        %dzv
    yprev( 9,:) = yw(3,:);    %dzw1
    yprev(10,:) = yw(4,:);    %dzw2
    yprev(11,:) = yyx(3,:);   %dzdzphi1
    yprev(12,:) = yyx(4,:);   %dzdzphi3
    yprev(13,:) = ywx(3,:);   %dzdzw1
    yprev(14,:) = ywx(4,:);   %dzdzw2
    yprev(15,:) = FFF;        %Freact
    winit = [yprev(4,1);
             yprev(5,1);
             yprev(9,1);
             yprev(10,1)];   %update initial w guess for bvp4c
    yinit = [yprev(1,1);
             yprev(2,1);
             yprev(6,1);
             yprev(7,1);];   %update initial phi guess for bvp4c
    Wtemp = Wprev;           %update the W value
    Wprev = Wtemp+yprev(3,1+N)*DT+J3*DT*(phi3res-yprev(2,1+N))...
        *(phi3res>yprev(2,1+N))/yprev(2,1+N);
    Vtemp = Vtotal;          %update the L2H1 norm of v (trapezium integration)
    Vcorr = vy(1)*vy(1)+vyx(1)*vyx(1)+vy(1+N)*vy(1+N)...
        +vyx(1+N)*vyx(1+N);
    Vtotal = Vtemp + DT*(vy*vy'+vyx*vyx'-Vcorr/2)/N;
    %-------------------------------------------------------------------------     
    %Update the start structure with the new time step values
    start.solver = 'bvp4c'; %must be bvp4c, because deval fails when bvp5c
    %Attach domain vector.
    start.x = xx; 
    %Attach previous functions
    start.y = [yprev( 1,:);
            yprev( 2,:);
            yprev( 3,:);
            yprev( 4,:);
            yprev( 5,:);
            yprev( 6,:);
            yprev( 7,:);
            yprev( 9,:);
            yprev(10,:)];
    %Attach previous function derivative
    start.yp = [yprev( 6,:);
           yprev( 7,:);
           yprev( 8,:);
           yprev( 9,:);
           yprev(10,:);
           yprev(11,:);
           yprev(12,:);
           yprev(13,:);
           yprev(14,:)];
    %Attach the statistics of the structure as BVP4c does in output.
    start.stats = struct('nmeshpoints',1+N,'maxres',0.0001,...
                   'nODEevals',1,'nBCevals',1);
    %----------------------------------------------------------------------
    %Insert the new values in the data matrices for plotting at the
    %end.
     Z1(1+t,:)=yprev( 1,:);            %phi1     (volume fraction of gypsum)
     Z2(1+t,:)=1-yprev(1,:)-yprev(2,:);%phi2     (volume fraction of kalk)
     Z3(1+t,:)=yprev( 2,:);            %phi3     (volume fraction of sulphuric acid)
     Z4(1+t,:)=yprev( 3,:);            %v        (displacement velocity of acid)
     Z5(1+t,:)=yprev( 4,:);            %w1       (displacement of gypsum)
     Z6(1+t,:)=yprev( 5,:);            %w2       (displacement of kalk)
     Z7(1+t,:)=yprev( 6,:);            %dzphi1   (1st derivative)
     Z8(1+t,:)=yprev( 7,:);            %dzphi3   (1st derivative)
     Z9(1+t,:)=yprev( 8,:);            %dzv      (1st derivative)
    Z10(1+t,:)=yprev( 9,:);            %dzw1     (1st derivative)
    Z11(1+t,:)=yprev( 10,:);           %dzw2     (1st derivative)
    Z12(1+t,:)=yprev( 11,:);           %dzdzphi1 (3rd derivative)
    Z13(1+t,:)=yprev( 12,:);           %dzdzphi3 (3rd derivative)
    Z14(1+t,:)=yprev( 13,:);           %dzdzw1   (3rd derivative)
    Z15(1+t,:)=yprev( 14,:);           %dzdzw2   (3rd derivative)
    Z16(1+t,:)=yprev( 15,:);           %Freact   (Chemical reaction speed)
    Z17(1+t,1)=Wprev;                  %W        (modified height function)
    Z18(1+t,1)=Vtotal;                 %Vtotal   (integrated L^2H^1 velocity)
    %---------------------------------------------------------------------
    %Test whether we can go another time step.
    %First we test volume fractions
    if (min([min(Z1(1+t,:)),min(Z2(1+t,:)),min(Z3(1+t,:))]) < tol)
        disp('Concentrations are within tolerance level of vanishing or having unphysical values.');
        disp(sprintf('Program stopped at time step %u out of %u planned.',t,M));
        %------------------------------------------------------------------
        %Create the current surface plots.
        %subplot(3,3,1);surf(Z1,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('phi1');
        %subplot(3,3,2);surf(Z2,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('phi2');
        %subplot(3,3,3);surf(Z3,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('phi3');
        %subplot(3,3,4);surf(Z4,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('V3');
        %subplot(3,3,5);surf(Z5,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('w1');
        %subplot(3,3,6);surf(Z6,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('w2');
        %subplot(3,3,7);surf(Z16,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('F');
        %subplot(3,3,8);plot(Z18);title('Vtotal');
        %subplot(3,3,9);plot(Z17);title('W-height');
        timeend = t;
        
        %Create output files in specific directory
        filedirectory = 'D:\Simulations\3rd-set\Pressure\'; %specific directory
        names = cellstr(['phi1.txt    ';'phi2.txt    ';'phi3.txt    ';'V3.txt      ';'w1.txt      ';'w2.txt      ';'F.txt       ';'Vtotal.txt  ';'W-height.txt']); %specific filenames
        variables = cellstr(['Z1 ';'Z2 ';'Z3 ';'Z4 ';'Z5 ';'Z6 ';'Z16';'Z17';'Z18']); %variable names associated with the specific filenames
        ZZ = {Z1;Z2;Z3;Z4;Z5;Z6;Z16;Z17;Z18}; %variables associated with the variable names
        Zxx = {xx;xx;xx;xx;xx;xx;xx;'';''}; %domains associated with the variables
        for i = 1:9
            mkdir(strcat(filedirectory,simulation,'\'));
            filepath = strcat(filedirectory,simulation,'\',names(i)); %make specific location for files
            fileName = fopen(char(filepath),'w'); %create txt file at specific location
            fmt1 = ['%8s\t|\t',repmat('% 7f\t',1,size(ZZ{i}',1)),'\r\n\v\r\n'];
            fmt3 = ['%6.6f\t|\t',repmat('% 7f\t',1,size(ZZ{i}',1)),'\r\n'];
            fprintf(fileName,fmt1,char(variables(i)),Zxx{i}); %fill the txt file with the correct variable
            fprintf(fileName,fmt3,[DT*[0:M];ZZ{i}']); %fill the txt file with the correct variable
        end
        fclose('all');
        return 
    else
        %Then we test if the L2H1 norm of v has become too large
        if (Vtotal>Veloc^2)
            disp('V or DT is too small. V^2/DT must be at least ||dzv||^2.');
            disp(sprintf('Program stopped at time step %u out of %u planned.',t,M));
            %--------------------------------------------------------------
            %Create the current surface plots.
            %subplot(3,3,1);surf(Z1,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('phi1');
            %subplot(3,3,2);surf(Z2,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('phi2');
            %subplot(3,3,3);surf(Z3,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('phi3');
            %subplot(3,3,4);surf(Z4,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('V3');
            %subplot(3,3,5);surf(Z5,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('w1');
            %subplot(3,3,6);surf(Z6,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('w2');
            %subplot(3,3,7);surf(Z16,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('F');
            %subplot(3,3,8);plot(Z18);title('Vtotal');
            %subplot(3,3,9);plot(Z17);title('W-height');
            timeend = t;
            
            %Create output files in specific directory
            filedirectory = 'D:\Simulations\3rd-set\Pressure\'; %specific directory
            names = cellstr(['phi1.txt    ';'phi2.txt    ';'phi3.txt    ';'V3.txt      ';'w1.txt      ';'w2.txt      ';'F.txt       ';'Vtotal.txt  ';'W-height.txt']); %specific filenames
            variables = cellstr(['Z1 ';'Z2 ';'Z3 ';'Z4 ';'Z5 ';'Z6 ';'Z16';'Z17';'Z18']); %variable names associated with the specific filenames
            ZZ = {Z1;Z2;Z3;Z4;Z5;Z6;Z16;Z17;Z18}; %variables associated with the variable names
            Zxx = {xx;xx;xx;xx;xx;xx;xx;'';''}; %domains associated with the variables
            for i = 1:9
                mkdir(strcat(filedirectory,simulation,'\'));
                filepath = strcat(filedirectory,simulation,'\',names(i)); %make specific location for files
                fileName = fopen(char(filepath),'w'); %create txt file at specific location
                fmt1 = ['%8s\t|\t',repmat('% 7f\t',1,size(ZZ{i}',1)),'\r\n\v\r\n'];
                fmt3 = ['%6.6f\t|\t',repmat('% 7f\t',1,size(ZZ{i}',1)),'\r\n'];
                fprintf(fileName,fmt1,char(variables(i)),Zxx{i}); %fill the txt file with the correct variable
                fprintf(fileName,fmt3,[DT*[0:M];ZZ{i}']); %fill the txt file with the correct variable
            end
            fclose('all');
            return
        else
             %Then we test for non-decreasing in height function
            if (Z17(1+t,1)<Z17(t,1))
                disp('Decreasing height!');
                disp(sprintf('Program stopped at time step %u out of %u planned.',t,M));
                %--------------------------------------------------------------
                %Create the current surface plots.
                %subplot(3,3,1);surf(Z1,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('phi1');
                %subplot(3,3,2);surf(Z2,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('phi2');
                %subplot(3,3,3);surf(Z3,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('phi3');
                %subplot(3,3,4);surf(Z4,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('V3');
                %subplot(3,3,5);surf(Z5,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('w1');
                %subplot(3,3,6);surf(Z6,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('w2');
                %subplot(3,3,7);surf(Z16,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('F');
                %subplot(3,3,8);plot(Z18);title('Vtotal');
                %subplot(3,3,9);plot(Z17);title('W-height');
                timeend = t;
                
                %Create output files in specific directory
                filedirectory = 'D:\Simulations\3rd-set\Pressure\'; %specific directory
                names = cellstr(['phi1.txt    ';'phi2.txt    ';'phi3.txt    ';'V3.txt      ';'w1.txt      ';'w2.txt      ';'F.txt       ';'Vtotal.txt  ';'W-height.txt']); %specific filenames
                variables = cellstr(['Z1 ';'Z2 ';'Z3 ';'Z4 ';'Z5 ';'Z6 ';'Z16';'Z17';'Z18']); %variable names associated with the specific filenames
                ZZ = {Z1;Z2;Z3;Z4;Z5;Z6;Z16;Z17;Z18}; %variables associated with the variable names
                Zxx = {xx;xx;xx;xx;xx;xx;xx;'';''}; %domains associated with the variables
                for i = 1:9
                    mkdir(strcat(filedirectory,simulation,'\'));
                    filepath = strcat(filedirectory,simulation,'\',names(i)); %make specific location for files
                    fileName = fopen(char(filepath),'w'); %create txt file at specific location
                    fmt1 = ['%8s\t|\t',repmat('% 7f\t',1,size(ZZ{i}',1)),'\r\n\v\r\n'];
                    fmt3 = ['%6.6f\t|\t',repmat('% 7f\t',1,size(ZZ{i}',1)),'\r\n'];
                    fprintf(fileName,fmt1,char(variables(i)),Zxx{i}); %fill the txt file with the correct variable
                    fprintf(fileName,fmt3,[DT*[0:M];ZZ{i}']); %fill the txt file with the correct variable
                end
                fclose('all');
                return
            end
        end
    end
end
timeend = t;
%----------------------------------------------------------------------
%Create the final surface plots.
%subplot(3,3,1);surf(Z1,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('phi1');
%subplot(3,3,2);surf(Z2,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('phi2');
%subplot(3,3,3);surf(Z3,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('phi3');
%subplot(3,3,4);surf(Z4,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('V3');
%subplot(3,3,5);surf(Z5,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('w1');
%subplot(3,3,6);surf(Z6,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('w2');
%subplot(3,3,7);surf(Z16,'EdgeColor','none','LineStyle','none','FaceLighting','phong');title('F');
%subplot(3,3,8);plot(Z18);title('Vtotal');
%subplot(3,3,9);plot(Z17);title('W-height');

%Create output files in specific directory
filedirectory = 'D:\Simulations\3rd-set\Pressure\'; %specific directory
names = cellstr(['phi1.txt    ';'phi2.txt    ';'phi3.txt    ';'V3.txt      ';'w1.txt      ';'w2.txt      ';'F.txt       ';'W-height.txt';'Vtotal.txt  ']); %specific filenames
variables = cellstr(['Z1 ';'Z2 ';'Z3 ';'Z4 ';'Z5 ';'Z6 ';'Z16';'Z17';'Z18']); %variable names associated with the specific filenames
ZZ = {Z1;Z2;Z3;Z4;Z5;Z6;Z16;Z17;Z18}; %variables associated with the variable names
Zxx = {xx;xx;xx;xx;xx;xx;xx;'';''}; %domains associated with the variables
for i = 1:9
    mkdir(strcat(filedirectory,simulation,'\'));
    filepath = strcat(filedirectory,simulation,'\',names(i)); %make specific location for files
    fileName = fopen(char(filepath),'w'); %create txt file at specific location
    fmt1 = ['%8s\t|\t',repmat('% 7f\t',1,size(ZZ{i}',1)),'\r\n\v\r\n'];
    fmt3 = ['%6.6f\t|\t',repmat('% 7f\t',1,size(ZZ{i}',1)),'\r\n'];
    fprintf(fileName,fmt1,char(variables(i)),Zxx{i}); %fill the txt file with the correct variable
    fprintf(fileName,fmt3,[DT*[0:M];ZZ{i}']); %fill the txt file with the correct variable
end
fclose('all');
return
end