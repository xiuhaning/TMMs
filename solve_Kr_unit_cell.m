function [K_r]=solve_Kr_unit_cell(theta)
%%
K_r=zeros(6,6);
K=zeros(10,10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global a_b; global b_b; global c_b;
global psi_ab; global psi_bb; global psi_cb;
global a_r; global b_r; global c_r;
global psi_ar; global psi_br; global psi_cr;
global k_bond; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=k_bond; 
t_bb=theta+psi_bb;
t_bbcr=theta+psi_bb+psi_cr;


%%
%Square terms k_ii
Kr11=k1*cos(t_bb)^2+k1*cos(t_bbcr)^2;
K(2,2)=k1*cos(t_bb)^2+k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))^2/...
    ((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
Kr33=k1*cos(t_bbcr)^2+k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))^2/((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
Kr66=k1*sin(t_bb)^2+k1*sin(t_bbcr)^2;
K(7,7)=k1*sin(t_bb)^2+k1*(a_r*sin(t_bb)-b_r*sin(t_bbcr))^2/((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
Kr88=k1*sin(t_bbcr)^2+k1*(a_r*sin(t_bb)-b_r*sin(t_bbcr))^2/((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));



%%
%Quadratic terms k_ij, i!=j
%coefficients for u_x1*R, i=1
K(1,2)=-k1*cos(t_bb)^2;
K(1,3)=-k1*cos(t_bbcr)^2;

Kr16=-k1*cos(t_bb)*sin(t_bb)-k1*cos(t_bbcr)*sin(t_bbcr);
K(1,7)=k1*cos(t_bb)*sin(t_bb);
K(1,8)=k1*cos(t_bbcr)*sin(t_bbcr);

%coefficients for u_y1*R, i=6
K(2,6)=k1*cos(t_bb)*sin(t_bb);
K(3,6)=k1*cos(t_bbcr)*sin(t_bbcr);
K(6,7)=-k1*sin(t_bb)^2;
K(6,8)=-k1*sin(t_bbcr)^2;

%coefficients for u_x2*R, i=2
K(2,3)=-k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))^2/((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
K(2,7)=-k1*cos(t_bb)*sin(t_bb)+...
    k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))*(-a_r*sin(t_bb)+b_r*sin(t_bbcr))/...
    ((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
K(2,8)=k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))*(a_r*sin(t_bb)-b_r*sin(t_bbcr))/...
    ((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));

%coefficients for u_y2*R, i=7
K(3,7)=k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))*(a_r*sin(t_bb)-b_r*sin(t_bbcr))/...
    ((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
K(7,8)=-k1*(a_r*sin(t_bb)-b_r*sin(t_bbcr))^2/((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));

%coefficients for u_x3*R, i=3
Kr38=-k1*cos(t_bbcr)*sin(t_bbcr)+k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))*(-a_r*sin(t_bb)+b_r*sin(t_bbcr))/...
   ((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));



K_r(1,1)=Kr11;K_r(1,2)=K(1,2);K_r(1,3)=K(1,3);K_r(1,4)=Kr16;K_r(1,5)=K(1,7);K_r(1,6)=K(1,8);
K_r(2,2)=K(2,2);K_r(2,3)=K(2,3);K_r(2,4)=K(2,6);K_r(2,5)=K(2,7);K_r(2,6)=K(2,8);
K_r(3,3)=Kr33;K_r(3,4)=K(3,6);K_r(3,5)=K(3,7);K_r(3,6)=Kr38;
K_r(4,4)=Kr66;K_r(4,5)=K(6,7);K_r(4,6)=K(6,8);
K_r(5,5)=K(7,7);K_r(5,6)=K(7,8);
K_r(6,6)=Kr88;
K_r = triu(K_r,0) + tril(K_r',-1);







