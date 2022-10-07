clear all
clc
load Homogeneous_lattice_angles.mat
load Variables_names_m_by_n_lattice_30.mat
load Variables_names_unit_cell_30.mat
% load Variables_names_m_by_n_lattice.mat
% load Variables_names_unit_cell.mat

%Present n (row) by m (column) unit cells, as well as n-1 by m-1 hexagons
n=30;m=30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Blue Triangle
global a_b; global b_b; global c_b;
a_b=0.5;b_b=0.7;c_b=1;

global psi_ab; global psi_bb; global psi_cb;
psi_ab=acos((a_b^2-b_b^2-c_b^2)/(-2*b_b*c_b));
psi_bb=acos((b_b^2-a_b^2-c_b^2)/(-2*a_b*c_b));
psi_cb=acos((c_b^2-b_b^2-a_b^2)/(-2*b_b*a_b));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Red Triangle
global a_r; global b_r; global c_r;
a_r=0.4;b_r=0.8;c_r=1;

global psi_ar; global psi_br; global psi_cr;
psi_ar=acos((a_r^2-b_r^2-c_r^2)/(-2*b_r*c_r));
psi_br=acos((b_r^2-a_r^2-c_r^2)/(-2*a_r*c_r));
psi_cr=acos((c_r^2-b_r^2-a_r^2)/(-2*b_r*a_r));

global l_s; l_s=(a_b+b_r)/2*0.95;

Edge_solver=[2]; %1 or 2    
%if Edge_solver==1, Left, bottom and right boundaries of the lattice are fixed.
%u_bottom, u_left, and u_right ar equal to zero. 
%Calculate edge stiffness of the top boundary.
%if Edge_solver==2, Left, top and right boundaries of the lattice are fixed.
%u_top, u_left, and u_right ar equal to zero.
%Calculate edge stiffness of the bottom boundary

global k_bond; global k_spring; 
k_bond=1; k_spring=0.01;

u00=0.01;
% Point load is applied at the mid-point of the test edge and make the
% point have an initial movement of u0.


%%

%Choose the angle of the homogeneous lattice (1<=i_alpha<=176), (0.3344<=alpha<=3.8344)
% for i_alpha=1:length(Alpha)
i_alpha=89;
alpha=Alpha(i_alpha);
gamma=Gamma(i_alpha);
theta=Theta(i_alpha);
kappa_horizontal=psi_ab+gamma-pi;
shift_horizontal=[c_b+a_r*cos(kappa_horizontal),a_r*sin(kappa_horizontal)];

vertical1(1)=1/sqrt(1+shift_horizontal(1)^2/shift_horizontal(2)^2);
vertical1(2)=-shift_horizontal(1)/shift_horizontal(2)*vertical1(1);

vertical2(1)=-1/sqrt(1+shift_horizontal(1)^2/shift_horizontal(2)^2);
vertical2(2)=-shift_horizontal(1)/shift_horizontal(2)*vertical2(1);
Vector=vertical1-shift_horizontal;
if atan2d(Vector(2),Vector(1))>0
    vertical=vertical1';
else
    vertical=vertical2';
end

%%
%Calculate K, K_b, K_r matrices of any unit cell
K_b=solve_Kb_unit_cell();
K_r=solve_Kr_unit_cell(theta);
[K,D]=solve_KD_bond_unit_cell(alpha,theta);
% rank(K)
% eig(K)
% eig(K_r)
% eig(K_b)

K_entire=zeros(length(U_entire_name)); %creates the entire K square matrix with length(U_entire) rows and columns
D_entire=zeros(length(U_entire_name),1); %creates the entire D 1-D matrix with length(U_entire) columns

%%
%Incorporate elments of K_b for the first row, i=1, into the entire K_entire matrix
for j=1:m
    %find U_entire(p_entire,q_entire) which is equal to U_b(p,q), and
    %incorporate K_b(p,q) into the corresponding K_entire(p_entire,q_entire).
    for p=1:length(U_b_name(:,j))
        p_entire=find(abs(U_entire_name-U_b_name(p,j))<1e-5); %find U_entire(p_entire)==U_b(p)
        for q=1:length(U_b_name(:,j))
            q_entire=find(abs(U_entire_name-U_b_name(q,j))<1e-5);%find U_entire(q_entire)==U_b(q)
            K_entire(p_entire,q_entire)=K_entire(p_entire,q_entire)+K_b(p,q);
        end
    end
end


%%
for j=1:m
    %find U_entire(p_entire,q_entire) which is equal to U_r(p,q), and
    %incorporate K_r(p,q) into the corresponding K_entire(p_entire,q_entire).
    for p=1:length(U_r_name(:,j))
        p_entire=find(abs(U_entire_name-U_r_name(p,j))<1e-5); %find U_entire(p_entire)==U_r(p)
        for q=1:length(U_r_name(:,j))
            q_entire=find(abs(U_entire_name-U_r_name(q,j))<1e-5); %find U_entire(q_entire)==U_r(q)
            K_entire(p_entire,q_entire)=K_entire(p_entire,q_entire)+K_r(p,q);
        end
    end
end


%%
%Incorporate elments of K for 2nd to nth row, i=2:n, into the entire K_entire matrix
for i=2:n
    for j=1:m
        %find U_entire(p_entire,q_entire) which is equal to U(p,q), and
        %incorporate K(p,q) and D(p) into the corresponding K_entire(p_entire,q_entire).
        for p=1:length(U_name(i,j,:))
            p_entire=find(abs(U_entire_name-U_name(i,j,p))<1e-5); %find U_entire(p_entire)==U(p)
            for q=1:length(U_name(i,j,:))
                q_entire=find(abs(U_entire_name-U_name(i,j,q))<1e-5); %find U_entire(q_entire)==U(q)
                K_entire(p_entire,q_entire)=K_entire(p_entire,q_entire)+K(p,q);
            end
            D_entire(p_entire)=D_entire(p_entire)+D(p);
        end
    end
end
% M_entire=eye(length(U_entire_name));
% [vector,omega]=eig(K_entire,M_entire);
% omega=diag(omega);
omega=eig(K_entire);

%%
for i_edge=1:length(Edge_solver)
%Remove boundary corresponding rows and columns of the dynamic matrix for solve other variables.
%Remove variable names of corresponding boundaries.
u0=u00*vertical;
U_entire_movable_name=U_entire_name;
K_entire_movable=K_entire;
D_entire_movable=D_entire;
for p=1:length(U_left_name)%Left boundary is fixed
    p_entire_movable1(p)=find(U_entire_name==U_left_name(p)); %find U_entire(p_entire)==U_left_name(p)
end
p_entire_movable=p_entire_movable1;
for p=1:length(U_right_name)%Right boundary is fixed
    p_entire_movable1(p)=find(U_entire_name==U_right_name(p)); %find U_entire(p_entire)==U_right_name(p)
end
p_entire_movable=[p_entire_movable,p_entire_movable1];
if Edge_solver(i_edge)==1 %Bottom boundary are fixed
    for p=1:length(U_bottom_name)
        p_entire_movable1(p)=find(U_entire_name==U_bottom_name(p)); %find U_entire(p_entire)==U_bottom_name(p)
    end
elseif Edge_solver(i_edge)==2 %Left, top and right boundaries are fixed
    for p=1:length(U_top_name)
        p_entire_movable1(p)=find(U_entire_name==U_top_name(p)); %find U_entire(p_entire)==U_top_name(p)
    end
end

p_entire_movable=[p_entire_movable,p_entire_movable1];
K_entire_movable(p_entire_movable,:)=[];
K_entire_movable(:,p_entire_movable)=[];
U_entire_movable_name(p_entire_movable)=[];
D_entire_movable(p_entire_movable)=[];

% eig(K_entire_movable)
% eig(K_entire)

U0_entire_movable = pcg(K_entire_movable, D_entire_movable, 1e-6,length(D_entire_movable)*10,[],[],D_entire_movable);


%%
%Give a point load at the middle of the test edge
if Edge_solver(i_edge)==1 %point load at mid-top
    i=n+1;j=round(m/2)+1;
    u0=-u0;
elseif Edge_solver(i_edge)==2 %point load at mid-bottom
    i=1; j=round(m/2);
    u0=u0;
end
  
result = strcat(num2str(i),num2str(j),num2str(1),num2str(3));  u0_name(1) = str2num(result);
result = strcat(num2str(i),num2str(j),num2str(2),num2str(3));  u0_name(2) = str2num(result);

if j<10
    u0_name=u0_name/1e3;
elseif j<100
    u0_name=u0_name/1e4;
else
    u0_name=u0_name/1e5;
end    

K_entire_movable2=K_entire_movable;
%Find and remove rows and columns of the K_entire_movable corresponding 
%to variables named form u0_name
for i=1:length(u0_name)
p_entire_movable2(i)=find(U_entire_movable_name==u0_name(i));
end

K_entire_movable2(p_entire_movable2,:)=[];
K_entire_movable2(:,p_entire_movable2)=[];

%%
%Use conjugate gradient method to solve dynamic variables

%f(u)=u'*K*u, df/du=2K*u. When df/du=0, f(u) is minimum, and K*u=0. u0 are two known variables from
%u, the other unknown variables can be solve based on u0 using conjugate gradient method .

% U_entire_movable=1e-3*ones(length(K_entire_movable2),1);
b=zeros(length(K_entire_movable),1);
u0=u0+U0_entire_movable(p_entire_movable2);
for p=1:length(u0_name)
    b= b-K_entire_movable(:,p_entire_movable2(p))*u0(p);
end
b=b-D_entire_movable;
b(p_entire_movable2)=[];

%Precondition conjugate gradient method
U_entire_movable = pcg(K_entire_movable2, b, 1e-6,length(b)*10,[],[],b);
% [U_entire_movable, Rsnew, iteration] = conjugate_gradient(K_entire_movable2, b, U_entire_movable);
% iteration

for p=1:length(p_entire_movable2)
    U_entire_movable=[U_entire_movable(1:p_entire_movable2(p)-1);...
        u0(p);U_entire_movable(p_entire_movable2(p):end)];
end
U_entire_movable=U_entire_movable-U0_entire_movable;
% K_entire_movable*U_entire_movable

%Store all dynamic variables including boundaries
U_entire=zeros(length(U_entire_name),1);
for p=1:length(U_entire_movable_name)%Left boundary is fixed
    p_entire_movable(p)=find(U_entire_name==U_entire_movable_name(p)); %find U_entire(p_entire)==U_left_name(p)
end
U_entire(p_entire_movable)=U_entire_movable;

Force_entire=K_entire*U_entire;  

for i=1:length(u0_name)
p_entire_movable3(i)=find(U_entire_name==u0_name(i));
end
Force=norm(Force_entire(p_entire_movable3));
Edge_stiffness(i_alpha,i_edge)=Force/u00;
end
% end

% p_edge=find(Edge_stiffness(:,1)>0); 

Edge_stiffness=abs(Edge_stiffness);
% % save('Dynamic_variables','U_entire')
figure;
semilogy(Alpha,Edge_stiffness,'k:','linewidth',2)
xlabel('\alpha')
ylabel('Edge stiffness')
title(['Edge stiffness against \alpha of a homogeneous Lattice '])

% THETA=2*pi-psi_cb-psi_ar-Alpha;
% figure;
% semilogy(THETA,Edge_stiffness,'k:','linewidth',2)

