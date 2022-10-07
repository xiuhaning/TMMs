function [K_b]=solve_Kb_unit_cell()
%%
K_b=zeros(6,6);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global a_b; global c_b;
global psi_bb; 
global k_bond; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=k_bond; 


%%
%Square terms k_ii
K_b(1,1)=k1+k1*cos(psi_bb)^2;
K_b(2,2)=k1+k1*(c_b-a_b*cos(psi_bb))^2/(a_b^2+c_b^2-2*a_b*c_b*cos(psi_bb));
K_b(3,3)=k1*cos(psi_bb)^2+k1*(c_b-a_b*cos(psi_bb))^2/(a_b^2+c_b^2-2*a_b*c_b*cos(psi_bb));
K_b(4,4)=k1*sin(psi_bb)^2;
K_b(5,5)=k1*(a_b*sin(psi_bb))^2/(a_b^2+c_b^2-2*a_b*c_b*cos(psi_bb));
K_b(6,6)=k1*sin(psi_bb)^2+k1*(a_b*sin(psi_bb))^2/(a_b^2+c_b^2-2*a_b*c_b*cos(psi_bb));

%%
K_b(1,2)=-k1;
K_b(1,3)=-k1*cos(psi_bb)^2;
K_b(1,4)=-k1*cos(psi_bb)*sin(psi_bb);
K_b(1,6)=k1*cos(psi_bb)*sin(psi_bb);

K_b(2,3)=-k1*(c_b-a_b*cos(psi_bb))^2/(a_b^2+c_b^2-2*a_b*c_b*cos(psi_bb));
K_b(2,5)=k1*(c_b-a_b*cos(psi_bb))*a_b*sin(psi_bb)/(a_b^2+c_b^2-2*a_b*c_b*cos(psi_bb));
K_b(2,6)=-k1*(c_b-a_b*cos(psi_bb))*a_b*sin(psi_bb)/(a_b^2+c_b^2-2*a_b*c_b*cos(psi_bb));

K_b(3,4)=k1*cos(psi_bb)*sin(psi_bb);
K_b(3,5)=-k1*(c_b-a_b*cos(psi_bb))*a_b*sin(psi_bb)/(a_b^2+c_b^2-2*a_b*c_b*cos(psi_bb));
K_b(3,6)=-k1*cos(psi_bb)*sin(psi_bb)+...
    k1*(c_b-a_b*cos(psi_bb))*a_b*sin(psi_bb)/(a_b^2+c_b^2-2*a_b*c_b*cos(psi_bb));

K_b(4,6)=-k1*sin(psi_bb)^2;
K_b(5,6)=-k1*(a_b*sin(psi_bb))^2/(a_b^2+c_b^2-2*a_b*c_b*cos(psi_bb));

K_b = triu(K_b,0) + tril(K_b',-1);









