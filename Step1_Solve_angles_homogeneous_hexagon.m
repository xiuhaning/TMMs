clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Blue Triangle 
%Define dimensions of blue triangle
global a_b; global b_b; global c_b;
a_b=0.5;b_b=0.7;c_b=1;
% a_b=0.8;b_b=1;c_b=0.4;

global psi_ab; global psi_bb; global psi_cb;
psi_ab=acos((a_b^2-b_b^2-c_b^2)/(-2*b_b*c_b));
psi_bb=acos((b_b^2-a_b^2-c_b^2)/(-2*a_b*c_b));
psi_cb=acos((c_b^2-b_b^2-a_b^2)/(-2*b_b*a_b));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Red Triangle
%Define dimensions of red triangle
global a_r; global b_r; global c_r;
a_r=0.4;b_r=0.8;c_r=1;
% a_r=1;b_r=0.5;c_r=0.7;

global psi_ar; global psi_br; global psi_cr;
psi_ar=acos((a_r^2-b_r^2-c_r^2)/(-2*b_r*c_r));
psi_br=acos((b_r^2-a_r^2-c_r^2)/(-2*a_r*c_r));
psi_cr=acos((c_r^2-b_r^2-a_r^2)/(-2*b_r*a_r));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha0=100/180*pi; %(must be greater than 18.3 degree, less than 2pi-psi_cb-psi_ar)
theta0=170/180*pi;gamma0=120/180*pi; % initial guess for theta and gamma
% theta must be less than 2pi-psi_cr-psi_bb, and gamma must be below 2pi-psi_ab-psi_br

%First Shot
R0=[theta0,gamma0]';   
F_coorD = solve_coordinate_D(alpha0,R0(1),R0(2));
f_D=[F_coorD(1),F_coorD(2)]';
J_f_D=[F_coorD(3) F_coorD(4); F_coorD(5) F_coorD(6)];
R1=R0-J_f_D\f_D;
% i=0;
%Netow's method
while norm(R1-R0)>10^-9  %&&  rcond(J_f)>10^-10
    i=i+1
    R0=R1;
    F_coorD = solve_coordinate_D(alpha0,R0(1),R0(2));
    f_D=[F_coorD(1),F_coorD(2)]';
    J_f_D=[F_coorD(3) F_coorD(4); F_coorD(5) F_coorD(6)];
    R1=R0-J_f_D\f_D;
end

theta0=R1(1); gamma0=R1(2); %solve for theta and gamma 
l_s=sqrt((a_b/2)^2+(b_r/2)^2-1/2*a_b*b_r*cos(alpha0+psi_ar));
rest_length_srping0=l_s/((a_b+b_r)/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%
%Find alpha, theta, gamma when 0<alpha<120/180*pi(alpha0)
alpha=alpha0-0.005:-0.005:0;   
aa=0;
theta=theta0;gamma=gamma0;

while aa<length(alpha) 
    aa=aa+1;
    %First Shot
    R0=[theta,gamma]';   
    F_coorD = solve_coordinate_D(alpha(aa),R0(1),R0(2));
    f_D=[F_coorD(1),F_coorD(2)]';
    J_f_D=[F_coorD(3) F_coorD(4); F_coorD(5) F_coorD(6)];
    R1=R0-J_f_D\f_D;
    i=0;
    %Netow's method
    while norm(R1-R0)>10^-9  %&&  rcond(J_f)>10^-10
        i=i+1
        R0=R1;
        F_coorD = solve_coordinate_D(alpha(aa),R0(1),R0(2));
        f_D=[F_coorD(1),F_coorD(2)]';
        J_f_D=[F_coorD(3) F_coorD(4); F_coorD(5) F_coorD(6)];
        R1=R0-J_f_D\f_D;
    end
    
    if R1(1)>=0 && R1(1)<=2*pi-psi_cr-psi_bb && R1(2)>=0 && R1(2)<=2*pi-psi_ab-psi_br  % constraints for angles
        alpha1(aa)=alpha(aa); theta1(aa)=R1(1); gamma1(aa)=R1(2); %solve for theta(aa) and gamma(aa) for alpha1(aa) 
        theta=R1(1); gamma=R1(2); %initial guess for theta(aa+1) and gamma(aa+1）
        l_s=sqrt((a_b/2)^2+(b_r/2)^2-1/2*a_b*b_r*cos(alpha1(aa)+psi_ar));
        rest_length_srping1(aa)=l_s/((a_b+b_r)/2);
    else
        aa=length(alpha);
    end
end

%%
%Find alpha, theta, gamma when 120/180*pi(alpha0)<alpha<2pi-psi_cb-psi_ar
alpha=alpha0+0.005:0.005:2*pi;   
aa=0;
theta=theta0;gamma=gamma0;

while aa<length(alpha) 
    aa=aa+1;
    %First Shot
    R0=[theta,gamma]';   
    F_coorD = solve_coordinate_D(alpha(aa),R0(1),R0(2));
    f_D=[F_coorD(1),F_coorD(2)]';
    J_f_D=[F_coorD(3) F_coorD(4); F_coorD(5) F_coorD(6)];
    R1=R0-J_f_D\f_D;
    i=0;
    %Netow's method
    while norm(R1-R0)>10^-9  %&&  rcond(J_f)>10^-10
        i=i+1
        R0=R1;
        F_coorD = solve_coordinate_D(alpha(aa),R0(1),R0(2));
        f_D=[F_coorD(1),F_coorD(2)]';
        J_f_D=[F_coorD(3) F_coorD(4); F_coorD(5) F_coorD(6)];
        R1=R0-J_f_D\f_D;
    end
    
    if R1(1)>=0 && R1(1)<=2*pi-psi_cr-psi_bb && R1(2)>=0 && R1(2)<=2*pi-psi_ab-psi_br  % constraints for angles
        alpha2(aa)=alpha(aa); theta2(aa)=R1(1); gamma2(aa)=R1(2); %solve for theta(aa) and gamma(aa) for alpha1(aa) 
        theta=R1(1); gamma=R1(2); %initial guess for theta(aa+1) and gamma(aa+1）
        l_s=sqrt((a_b/2)^2+(b_r/2)^2-1/2*a_b*b_r*cos(alpha2(aa)+psi_ar));
        rest_length_srping2(aa)=l_s/((a_b+b_r)/2);
    else
        aa=length(alpha);
    end
end

%%
Alpha=[flip(alpha1),alpha0,alpha2]; %combine all angles found by Newton's method
Theta=[flip(theta1),theta0,theta2];
Gamma=[flip(gamma1),gamma0,gamma2];
rest_length_srping=[flip(rest_length_srping1),rest_length_srping0,rest_length_srping2];

figure;plot(Alpha,Theta,'k--','linewidth', 2)
hold on; plot(Alpha,Gamma,'r-','linewidth', 2)
ylabel('\theta and \gamma')
xlabel('\alpha')
legend('\theta','\gamma','Location','east')

% save('Homogeneous_lattice_angles_new.mat','Alpha','Theta','Gamma') 
save('Homogeneous_lattice_angles_bistable_11.mat','Alpha','Theta','Gamma','rest_length_srping') 

figure;plot(Alpha,rest_length_srping,'k--','linewidth', 2)


    
    