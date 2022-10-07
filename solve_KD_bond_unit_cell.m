function [K,D]=solve_KD_bond_unit_cell(alpha,theta)
%%
K=zeros(10,10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global a_b; global b_b; global c_b;
global psi_ab; global psi_bb; global psi_cb;
global a_r; global b_r; global c_r;
global psi_ar; global psi_br; global psi_cr;
global k_bond; global k_spring; global l_s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=k_bond; k2=k_spring;
t_bb=theta+psi_bb;
t_bbcr=theta+psi_bb+psi_cr;
at_arbbcr=alpha+theta+psi_ar+psi_bb+psi_cr;
at_arbbcrcb=alpha+theta+psi_ar+psi_bb+psi_cr+psi_cb;
a_ar=alpha+psi_ar;

%%
%Square terms k_ii
K(1,1)=k1*cos(t_bb)^2+k1*cos(t_bbcr)^2+k2*(b_r*cos(t_bbcr)-a_b*cos(at_arbbcr))^2/...
    (4*(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)));
K(2,2)=k1*cos(t_bb)^2+k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))^2/...
    ((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
K(3,3)=k1*cos(t_bbcr)^2+k1*cos(at_arbbcr)^2+k1*cos(at_arbbcrcb)^2+...
    k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))^2/((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
K(4,4)=k1*cos(at_arbbcr)^2+k2*(b_r*cos(t_bbcr)-a_b*cos(at_arbbcr))^2/...
    (4*(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)))+...
    k1*(a_b*cos(at_arbbcr)-b_b*cos(at_arbbcrcb))^2/((a_b^2+b_b^2-2*a_b*b_b*cos(psi_cb)));
K(5,5)=k1*cos(at_arbbcrcb)^2+k1*(a_b*cos(at_arbbcr)-b_b*cos(at_arbbcrcb))^2/...
    ((a_b^2+b_b^2-2*a_b*b_b*cos(psi_cb)));
K(6,6)=k1*sin(t_bb)^2+k1*sin(t_bbcr)^2+k2*(b_r*sin(t_bbcr)-a_b*sin(at_arbbcr))^2/...
    (4*(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)));

K(7,7)=k1*sin(t_bb)^2+k1*(a_r*sin(t_bb)-b_r*sin(t_bbcr))^2/((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
K(8,8)=k1*sin(t_bbcr)^2+k1*sin(at_arbbcr)^2+k1*sin(at_arbbcrcb)^2+...
    k1*(a_r*sin(t_bb)-b_r*sin(t_bbcr))^2/((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
K(9,9)=k1*sin(at_arbbcr)^2+k2*(b_r*sin(t_bbcr)-a_b*sin(at_arbbcr))^2/...
    (4*(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)))+...
    k1*(a_b*sin(at_arbbcr)-b_b*sin(at_arbbcrcb))^2/...
    ((a_b^2+b_b^2-2*a_b*b_b*cos(psi_cb)));
K(10,10)=k1*sin(at_arbbcrcb)^2+k1*(a_b*sin(at_arbbcr)-b_b*sin(at_arbbcrcb))^2/...
    ((a_b^2+b_b^2-2*a_b*b_b*cos(psi_cb)));

%%
%Quadratic terms k_ij, i!=j
%coefficients for u_x1*R, i=1
K(1,2)=-k1*cos(t_bb)^2;
K(1,3)=-k1*cos(t_bbcr)^2;
K(1,4)=-k2*(b_r*cos(t_bbcr)-a_b*cos(at_arbbcr))^2/(4*(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)));
K(1,6)=-k1*cos(t_bb)*sin(t_bb)-k1*cos(t_bbcr)*sin(t_bbcr)-...
    k2*(b_r*cos(t_bbcr)-a_b*cos(at_arbbcr))*(b_r*sin(t_bbcr)-a_b*sin(at_arbbcr))/...
    (4*(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)));
K(1,7)=k1*cos(t_bb)*sin(t_bb);
K(1,8)=k1*cos(t_bbcr)*sin(t_bbcr);
K(1,9)=k2*(b_r*cos(t_bbcr)-a_b*cos(at_arbbcr))*(b_r*sin(t_bbcr)-a_b*sin(at_arbbcr))/...
    (4*(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)));

%coefficients for u_y1*R, i=6
K(2,6)=k1*cos(t_bb)*sin(t_bb);
K(3,6)=k1*cos(t_bbcr)*sin(t_bbcr);
K(4,6)=k2*(b_r*cos(t_bbcr)-a_b*cos(at_arbbcr))*(b_r*sin(t_bbcr)-a_b*sin(at_arbbcr))/...
    (4*(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)));
K(6,7)=-k1*sin(t_bb)^2;
K(6,8)=-k1*sin(t_bbcr)^2;
K(6,9)=-k2*(b_r*sin(t_bbcr)-a_b*sin(at_arbbcr))^2/(4*(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)));

%coefficients for u_x2*R, i=2
K(2,3)=-k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))^2/((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
K(2,7)=-k1*cos(t_bb)*sin(t_bb)+k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))*(-a_r*sin(t_bb)+b_r*sin(t_bbcr))/...
    ((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
K(2,8)=k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))*(a_r*sin(t_bb)-b_r*sin(t_bbcr))/...
    ((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));

%coefficients for u_y2*R, i=7
K(3,7)=k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))*(a_r*sin(t_bb)-b_r*sin(t_bbcr))/...
    ((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
K(7,8)=-k1*(a_r*sin(t_bb)-b_r*sin(t_bbcr))^2/((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));

%coefficients for u_x3*R, i=3
K(3,8)=-k1*cos(t_bbcr)*sin(t_bbcr)-k1*cos(at_arbbcr)*sin(at_arbbcr)-k1*cos(at_arbbcrcb)*sin(at_arbbcrcb)+...
    k1*(a_r*cos(t_bb)-b_r*cos(t_bbcr))*(-a_r*sin(t_bb)+b_r*sin(t_bbcr))/((a_r^2+b_r^2-2*a_r*b_r*cos(psi_cr)));
K(3,4)=-k1*cos(at_arbbcr)^2;
K(3,5)=-k1*cos(at_arbbcrcb)^2;
K(3,9)=k1*cos(at_arbbcr)*sin(at_arbbcr);
K(3,10)=k1*cos(at_arbbcrcb)*sin(at_arbbcrcb);

%coefficients for u_y3*R, i=8
K(4,8)=k1*cos(at_arbbcr)*sin(at_arbbcr);
K(5,8)=k1*cos(at_arbbcrcb)*sin(at_arbbcrcb);
K(8,9)=-k1*sin(at_arbbcr)^2;
K(8,10)=-k1*sin(at_arbbcrcb)^2;

%coefficients for u_x'1*R, i=4
K(4,5)=-k1*(a_b*cos(at_arbbcr)-b_b*cos(at_arbbcrcb))^2/((a_b^2+b_b^2-2*a_b*b_b*cos(psi_cb)));
K(4,9)=-k1*cos(at_arbbcr)*sin(at_arbbcr)-...
    k2*(b_r*cos(t_bbcr)-a_b*cos(at_arbbcr))*(b_r*sin(t_bbcr)-a_b*sin(at_arbbcr))/...
    (4*(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)))-...
    k1*(a_b*cos(at_arbbcr)-b_b*cos(at_arbbcrcb))*(a_b*sin(at_arbbcr)-b_b*sin(at_arbbcrcb))/...
    ((a_b^2+b_b^2-2*a_b*b_b*cos(psi_cb)));
K(4,10)=k1*(a_b*cos(at_arbbcr)-b_b*cos(at_arbbcrcb))*...
    (a_b*sin(at_arbbcr)-b_b*sin(at_arbbcrcb))/...
    ((a_b^2+b_b^2-2*a_b*b_b*cos(psi_cb)));

%coefficients for u_y'1*R, i=9
K(5,9)=-k1*(a_b*cos(at_arbbcr)-b_b*cos(at_arbbcrcb))*...
    (a_b*sin(at_arbbcr)-b_b*sin(at_arbbcrcb))/...
    ((a_b^2+b_b^2-2*a_b*b_b*cos(psi_cb)));
K(9,10)=-k1*(a_b*sin(at_arbbcr)-b_b*sin(at_arbbcrcb))^2/...
    ((a_b^2+b_b^2-2*a_b*b_b*cos(psi_cb)));

%coefficients for u_x''1*R, i=5
K(5,10)=-k1*cos(at_arbbcrcb)*sin(at_arbbcrcb)-...
    k1*(a_b*cos(at_arbbcr)-b_b*cos(at_arbbcrcb))*...
    (a_b*sin(at_arbbcr)-b_b*sin(at_arbbcrcb))/...
    ((a_b^2+b_b^2-2*a_b*b_b*cos(psi_cb)));

K = triu(K,0) + tril(K',-1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D=zeros(10,1);
D(1)=k2*(-2*l_s+sqrt(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)))*...
    (-b_r*cos(a_ar)+a_b*cos(at_arbbcr))/(2*sqrt(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)));
D(4)=k2*(-2*l_s+sqrt(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)))*...
    (b_r*cos(a_ar)-a_b*cos(at_arbbcr))/(2*sqrt(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)));
D(6)=k2*(-2*l_s+sqrt(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)))*...
    (b_r*sin(a_ar)-a_b*sin(at_arbbcr))/(2*sqrt(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)));
D(9)=k2*(-2*l_s+sqrt(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)))*...
    (-b_r*sin(a_ar)+a_b*sin(at_arbbcr))/(2*sqrt(a_b^2+b_r^2-2*a_b*b_r*cos(a_ar)));





