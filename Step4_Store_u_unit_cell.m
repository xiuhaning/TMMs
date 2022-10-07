clear all
clc

%Present n (row) by m (column) unit cells, as well as n-1 by m-1 hexagons
n=2;m=1;

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


%%
%Incorporate elments of K_b for the first row, i=1, into the entire K_entire matrix
U_b_name=zeros(6,m);
i=1;
for j=1:m
    %Store names of variables of a unit cell only including blue, on the bottom of thr lattice
    %u_ijx3, u_i+1jx1, u_i+1j+1x2, u_ijy3, u_i+1jy1, u_i+1j+1y2  
    result = strcat(num2str(i+1),num2str(j),num2str(1),num2str(1));U_b_name(1,j) = str2num(result);
    result = strcat(num2str(i+1),num2str(j+1),num2str(1),num2str(2));U_b_name(2,j) = str2num(result);
    result = strcat(num2str(i),num2str(j),num2str(1),num2str(3));U_b_name(3,j) = str2num(result);
    result = strcat(num2str(i+1),num2str(j),num2str(2),num2str(1));U_b_name(4,j) = str2num(result);
    result = strcat(num2str(i+1),num2str(j+1),num2str(2),num2str(2));U_b_name(5,j) = str2num(result);
    result = strcat(num2str(i),num2str(j),num2str(2),num2str(3));U_b_name(6,j) = str2num(result);
    if j<10
        U_b_name(:,j)=U_b_name(:,j)/1e3;
        if j==9
            U_b_name(2,j)=U_b_name(2,j)/10;U_b_name(5,j)=U_b_name(5,j)/10;
        end
    elseif j<100
        U_b_name(:,j)=U_b_name(:,j)/1e4;
    else
        U_b_name(:,j)=U_b_name(:,j)/1e5;
    end
end


%%
%Incorporate elments of K_r for the last row, i=n+1, into the entire K_entire matrix
U_r_name=zeros(6,m);
i=n+1;
for j=1:m
    %Store names of variables of a unit cell only including red, on the top of thr lattice
    %u_ijx1, u_ijx2, u_ijx3, u_ijy1, u_ijy2, u_ijy3
    result = strcat(num2str(i),num2str(j),num2str(1),num2str(1));U_r_name(1,j) = str2num(result);
    result = strcat(num2str(i),num2str(j),num2str(1),num2str(2));U_r_name(2,j) = str2num(result);
    result = strcat(num2str(i),num2str(j),num2str(1),num2str(3));U_r_name(3,j) = str2num(result);
    result = strcat(num2str(i),num2str(j),num2str(2),num2str(1));U_r_name(4,j) = str2num(result);
    result = strcat(num2str(i),num2str(j),num2str(2),num2str(2));U_r_name(5,j) = str2num(result);
    result = strcat(num2str(i),num2str(j),num2str(2),num2str(3));U_r_name(6,j) = str2num(result);
    if j<10
        U_r_name(:,j)=U_r_name(:,j)/1e3;
    elseif j<100
        U_r_name(:,j)=U_r_name(:,j)/1e4;
    else
        U_r_name(:,j)=U_r_name(:,j)/1e5;
    end
end


%%
%Incorporate elments of K for 2nd to nth row, i=2:n, into the entire K_entire matrix
U_name=zeros(n,m,10);
for i=2:n
    for j=1:m
        %Store names of variables of an entire unit cell: x and y replaced by 1 and 2,respectively, in the vectors
        %u_ijx1, u_ijx2, u_ijx3, u_i+1jx1, u_i+1j+1x2; 
        %u_ijy1, u_ijy2, u_ijy3, u_i+1jy1, u_i+1j+1y2
        result = strcat(num2str(i),num2str(j),num2str(1),num2str(1));U_name(i,j,1) = str2num(result);
        result = strcat(num2str(i),num2str(j),num2str(1),num2str(2));U_name(i,j,2) = str2num(result);
        result = strcat(num2str(i),num2str(j),num2str(1),num2str(3));U_name(i,j,3) = str2num(result);
        result = strcat(num2str(i+1),num2str(j),num2str(1),num2str(1));U_name(i,j,4) = str2num(result);
        result = strcat(num2str(i+1),num2str(j+1),num2str(1),num2str(2));U_name(i,j,5) = str2num(result);
        result = strcat(num2str(i),num2str(j),num2str(2),num2str(1));U_name(i,j,6) = str2num(result);
        result = strcat(num2str(i),num2str(j),num2str(2),num2str(2));U_name(i,j,7) = str2num(result);
        result = strcat(num2str(i),num2str(j),num2str(2),num2str(3));U_name(i,j,8) = str2num(result);
        result = strcat(num2str(i+1),num2str(j),num2str(2),num2str(1));U_name(i,j,9) = str2num(result);
        result = strcat(num2str(i+1),num2str(j+1),num2str(2),num2str(2));U_name(i,j,10) = str2num(result);
        if j<10
            U_name(i,j,:)=U_name(i,j,:)/1e3;
            if j==9
                U_name(i,j,5)=U_name(i,j,5)/10;U_name(i,j,10)=U_name(i,j,10)/10;
            end
        elseif j<100
            U_name(i,j,:)=U_name(i,j,:)/1e4;
        else
            U_name(i,j,:)=U_name(i,j,:)/1e5;
        end
    end
end


save('Variables_names_unit_cell_2','U_b_name','U_r_name','U_name')




