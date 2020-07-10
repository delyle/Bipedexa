function output = Bipedexa(auxdata)

%-------------------------------------------------------------------%
%-------------------- Data Required by Problem ---------------------%
%-------------------------------------------------------------------%
g= auxdata.g; %gravity
lmax=auxdata.lmax; %max leg length
d=auxdata.d; %step length
D=auxdata.D; %sride length
Fmax=auxdata.Fmax;
Taumax=auxdata.Taumax;
I=auxdata.I;
r=auxdata.r;
T=auxdata.T;
% specify auxdata if not already done

%-------------------------------------------------------------------%
%------------------------- Variable Bounds -------------------------%
%-------------------------------------------------------------------%
% ----- PHASE 1 ----- %
i = 1;
bounds.phase(i).initialtime.lower = 0;              % scalar
bounds.phase(i).initialtime.upper = 0;              % scalar
bounds.phase(i).finaltime.lower = T ;                % scalar
bounds.phase(i).finaltime.upper = T ;                % scalar
%States
%6 kinematic states
%3 Force states
%3 Torque states
%1 Integrated Force
%1 Integrated Torque
 xlow = 0;
 xupp = D;
 ylow = 0;
 yupp = Inf;
 Flow= zeros(1,3);
 Fupp= [1 1 1]*Fmax;
 Taulow= [1 1 1]*(-Taumax);
 Tauupp= [1 1 1]*Taumax;
bounds.phase(i).initialstate.lower = [xlow,ylow,-Inf,-Inf,-pi,-Inf,Flow,Taulow,0,0];           % row vector, length = numstates
bounds.phase(i).initialstate.upper = [xlow,yupp,Inf,Inf,Inf,pi,Fupp,Tauupp,0,Inf];           % row vector, length = numstates
bounds.phase(i).state.lower = [xlow,ylow,-Inf,-Inf,-pi,-Inf,Flow,Taulow,0,0];             % row vector, length = numstates
bounds.phase(i).state.upper = [xupp,yupp,Inf,Inf,pi,Inf,Fupp,Tauupp,Inf,Inf];                 % row vector, length = numstates
bounds.phase(i).finalstate.lower = [xupp,ylow,-Inf,-Inf,-pi,-Inf,Flow,Taulow,0,0];             % row vector, length = numstates
bounds.phase(i).finalstate.upper = [xupp,yupp,Inf,Inf,pi,Inf,Fupp,Tauupp,Fmax*T,Inf];             % row vector, length = numstates
% 3 Time derivative of force controls
% 3 time derivative of torque controls
neg=[1 1 1]*(-Inf);
pos= [1 1 1]*(Inf);
bounds.phase(i).control.lower = [neg,neg];                % row vector, length = numstates
bounds.phase(i).control.upper = [pos,pos];                % row vector, length = numstates
% ???
bounds.phase(i).integral.lower = 0;                 % row vector, length = numintegrals
bounds.phase(i).integral.upper = Inf;                 % row vector, length = numintegrals
% no parameters introduced
%bounds.parameter.lower = ;                      % row vector, length = numintegrals
%bounds.parameter.upper = ;                      % row vector, length = numintegrals

% Endpoint constraints (if required)

bounds.eventgroup.lower = [0,0,D,0,0,0,0,0]; % row vector
bounds.eventgroup.upper = [0,0,D,0,0,0,0,0]; % row vector

% Path constraints (if required)
% 6 complimentary limb length constaints 
% 2 complimentary exclusion constaints 

% ----- PHASE 1 ----- %
i = 1;
bounds.phase(i).path.lower = zeros(1,8); % row vector, length = number of path constraints in phase


bounds.phase(i).path.upper =[inf inf inf inf inf inf 0 0]; % row vector, length = number of path constraints in phase
%-------------------------------------------------------------------------%
%---------------------------- Provide Guess ------------------------------%
%-------------------------------------------------------------------------%
% ----- PHASE 1 ----- %
i = 1;
guess.phase(i).time    = [0;T];                % column vector, min length = 2
guess.phase(i).state   = [0,lmax+r,D/T,0,pi/2,0,Fmax*[1 1 1],[0 0 0],0,0;...
                          D,lmax+r,D/T,0,pi/2,0,Fmax*[1 1 1],[0 0 0],Fmax*T,0];% array, min numrows = 2, numcols = numstates
guess.phase(i).control = zeros(2,6);               % array, min numrows = 2, numcols = numcontrols
guess.phase(i).integral = 1;               % scalar

%guess.parameter = [];                    % row vector, numrows = numparams


%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
setup.mesh.maxiterations=2;

% not required

%-------------------------------------------------------------------%
%--------------------------- Problem Setup -------------------------%
%-------------------------------------------------------------------%
setup.name                        = 'bipedalprob';
setup.functions.continuous        = @Continuous;
setup.functions.endpoint          = @Endpoint;
setup.auxdata                     = auxdata; % not necessary
setup.bounds                      = bounds;
setup.nlp.solver= 'snopt';
setup.guess                       = guess;

setup.derivatives.derivativelevel = 'first';


%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);
end


function phaseout = Continuous(input)

% extract data
t = input.phase(1).time;
X = input.phase(1).state;
U = input.phase(1).control;

auxdata = input.auxdata;
g= auxdata.g; %gravity
lmax=auxdata.lmax; %max leg length
d=auxdata.d; %step length
D=auxdata.D; %sride length
Fmax=auxdata.Fmax;
Taumax=auxdata.Taumax;
I=auxdata.I;
r=auxdata.r;
T=auxdata.T;
c1=auxdata.c1;
c2=auxdata.c2;

%P = input.phase(1).parameter;

x=X(:,1);
y=X(:,2);
theta=X(:,5);
xdot = X(:,3); % provide derivative
ydot = X(:,4);
thetadot=X(:,6);
F=X(:,7:9); %Collecting Forces 
Tau=X(:,10:12); %Collecting Torques
P= X(:,13);
Q= X(:,14);

Ftr=F(:,1);
Flead=F(:,2);
Fref=F(:,3);
Tautr=Tau(:,1);
Taulead=Tau(:,2);
Tauref=Tau(:,3);
Tautrsqr=Tautr.^2;
Tauleadsqr=Taulead.^2;
Taurefsqr=Tauref.^2;

Fdot=U(:,1:3); 
Taudot=U(:,4:6);

Pdot= Flead;
Qdot= Taulead.^2;

%ntime=length(x);
zs=zeros(size(x));
os=ones(size(x));
temp=[os,zs,zs];
tempcol=os;
dvec=temp*d;
Dvec=temp*D;
dveccol=tempcol*d;
Dveccol=tempcol*D;

rvec=-r.*[cos(theta),sin(theta),zs];
xc=[x,y,zs];
ltr = xc +rvec;
llead=(xc+rvec)-dvec;
lref=(xc+rvec)-Dvec;

magnitudeltr=sqrt(dot(ltr,ltr,2));
magnitudellead=sqrt(dot(llead,llead,2));
magnitudelref= sqrt(dot(lref,lref,2));
ultr=ltr./magnitudeltr;
ullead=llead./magnitudellead;
ulref =lref./magnitudelref;

Ftrvec= Ftr.*ultr;
Fleadvec = Flead.*ullead;
Frefvec = Fref.*ulref;

crossFtr=cross(rvec,Ftrvec);
crossFtrz=crossFtr(:,3);  %Extracting z column 
crossFlead=cross(rvec,Fleadvec);
crossFleadz=crossFlead(:,3); %Extracting z column 
crossFref=cross(rvec,Frefvec);
crossFrefz=crossFref(:,3); %Extracting z column

xddot= Ftr.*(x./magnitudeltr)+Fref.*((x-dveccol)./magnitudelref)+Flead.*((x-Dveccol)./magnitudellead);
yddot= Ftr.*(y./magnitudeltr)+Fref.*(y./magnitudelref)+Flead.*(y./magnitudellead)-g;
thetaddot=(Tautr+Taulead+Tauref+crossFtrz+crossFleadz+crossFrefz)./I;

phaseout.dynamics = [xdot,ydot,xddot,yddot,thetadot,thetaddot,Fdot,Taudot,Pdot,Qdot];


%what is c1 and c2 respectively?
phaseout.integrand = c1*(Ftr.^2+Fref.^2+Flead.^2)+c2*(Tautr.^2+Tauref.^2+Taulead.^2);
%Path constraint
Ftrllc= Ftr.*(lmax-magnitudeltr);
%magFtrllc= sqrt(dot(Ftrllc,Ftrllc,2));
Fleadllc= Flead.*(lmax-magnitudellead);
%magFleadllc= sqrt(dot(Fleadllc,Fleadllc,2));
Frefllc= Fref.*(lmax-magnitudelref);
%magFrefllc=sqrt(dot(Frefllc,Frefllc,2));
Tautrllc= Tautrsqr.*(lmax-magnitudeltr);
%magTautrllc=sqrt(dot(Tautrllc,Tautrllc,2));
Tauleadllc=Tauleadsqr.*(lmax-magnitudellead);
%magTauleadllc=sqrt(dot(Tauleadllc,Tauleadllc,2));
Taurefllc= Taurefsqr.*(lmax-magnitudelref);
%magTaurefllc=sqrt(dot(Taurefllc,Taurefllc,2));
Fxc=P.*Ftr;
Tauxc=Q.*Tautr;
phaseout.path = [Ftrllc,Fleadllc,Frefllc,Tautrllc,Tauleadllc,Taurefllc,Fxc,Tauxc]; % path constraints, matrix of size num collocation points X num path constraints
end

function output = Endpoint(input)

Finalstates =input.phase(1).finalstate;
Initialstates= input.phase(1).initialstate;

Ftr=Initialstates(7); 
Flead=Finalstates(8);

Ttr=Initialstates(10);
Tlead=Finalstates(11);

xbeg=Initialstates(1); %diff of d
xend=Finalstates(1);

ybeg=Initialstates(2); %equal
yend=Finalstates(2);

xdotbeg=Initialstates(3);
xdotend=Finalstates(3);

ydotbeg=Initialstates(4);
ydotend=Finalstates(4);

thetabeg=Initialstates(5);
thetaend=Finalstates(5);

omegabeg=Initialstates(6);
omegaend=Finalstates(6);

output.eventgroup.event = [(Ftr-Flead) (Ttr-Tlead) (xend-xbeg) (ybeg-yend) (xdotbeg-xdotend) (ydotbeg-ydotend) (thetabeg-thetaend) (omegabeg-omegaend)];% event constraints (row vector)

J = input.phase(1).integral(1);
output.objective = J; % objective function (scalar)

end