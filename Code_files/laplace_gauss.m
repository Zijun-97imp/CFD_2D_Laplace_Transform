% Gauss Seidel to solve Laplace Eqn
%     T_xx + T_yy = 0
% with Dirichlet BCs.

%% clear the workspace
clc
clear

% Parameters
nx = 51;
ny = nx;                       % number of space steps
x=linspace(0,1,nx);            % x range
y=linspace(0,1,ny);            % y range
dx = 0.02;
dy = dx;
B=(dx/dy)^2;
T_gs=zeros(nx,ny);             % initialize T for Point Jacobi
tol=1e-5;                      % error tolerance
err=1;                         % error
k=1;                           % iteration index
omega=1.95;                    % (optimal) relaxation parameter
%%
% Boundary Conditions
T_gs(1,:) = 1;%left
T_gs(nx,:) = cos(6*1.5*pi*y)+1;%right
T_gs(:,1) = 1+x;%bottom
T_gs(:,end) = 1;%top
Texact=T_gs;     %analytical 

%%
% Time loop
i=1;j=1;
rms=0;
while err>tol
    T_gsold=T_gs;
    for i=2:nx-1
        for j=2:ny-1
            %Gauss-Seidel
            T_gs(i,j)= (1-omega)*T_gsold(i,j) +(omega/(2*(1+B)))*(T_gs(i-1,j)+T_gsold(i+1,j)+B*(T_gs(i,j-1)+T_gsold(i,j+1)));
        end
    end
    %boundary conditions
    T_gs(1,:)= 1;%left
    T_gs(nx,:) = cos(6*1.5*pi*y)+1;%right
    T_gs(:,1) = 1+x;%bottom
    T_gs(:,end) = 1;%top
    T_gs(25,25)=1.5;
    T_gs(10,10)=0.5;
    
    err= max(max(abs(T_gs-T_gsold)));
    k=k+1;
end
%%
% analytical solution
for iii=1:nx
    for jjj=1:ny  
        A=0;
        for n=1:101
            if mod(n,2)==1
                A = A +((n*pi)^-2 * csch(2*n*pi) * sinh(n*pi*x(iii)) * cos(n*pi*y(jjj)));
            end 
        end
        Texact(iii,jjj)=(x(iii)/4)-(4*A);
    end
end


% plot
figure(1)
figure,surf(x,y,T_gs');
xlabel('x');ylabel('y');zlabel('Temperature');
axis square;
title('T(x,y) Gauss-Seidel');
legend('Gauss Seidel','Analytical','location','sw');
colorbar
pts={[0.5 0.25] [0.5 0.75] [1.5 0.25] [1.5 0.75]};

figure(2)
contourf(x,y,T_gs');
xlabel('x');ylabel('y');
axis square;
title('T(x,y) Gauss-Seidel');
legend('Gauss-Seidel','Analytical','location','best');
colorbar



