% Point Jacobi Method to solve Laplace Eqn
%     T_xx + T_yy = 0
% with Dirichlet BCs.
%% clear the workspace
clc
clea
%%
clearvars
% Parameters
nx = 51; 
ny = nx;                       % number of space steps
x=linspace(0,1,nx);            % x range
y=linspace(0,1,ny);            % y range

dx = 0.02;
dy = dx;                       % each unit length
B=(dx/dy)^2;
T_pj=zeros(nx,ny);             % initialize T 
tol=1e-5;                      % error tolerance for convergence
err=1;                         % to check convergence
k=1;                           % iteration index

%%
% Boundary Conditions
T_pj(1,:) = 1;%left
T_pj(nx,:) = cos(6*1.5*pi*y)+1;%right
T_pj(:,1) = 1+x;%bottom
T_pj(:,end) = 1;%top
Texact=T_pj;     % initialize analytical solution


%%
while err>tol
    T_pjold=T_pj;
    
    for i=2:nx-1
        for j=2:ny-1
            %Point Jacobi
            T_pj(i,j)=(1/(2*((1+B))))*(T_pjold(i-1,j)+T_pjold(i+1,j)+B*(T_pjold(i,j-1)+T_pjold(i,j+1)));
        end
    end
    %boundary conditions
    T_pj(1,:) = 1;%left
    T_pj(nx,:) = cos(6*1.5*pi*y)+1;%right
    T_pj(:,1) = 1+x;%bottom
    T_pj(:,end) = 1;%top
    T_pj(25,25)=1.5;
    T_pj(10,10)=0.5;
   
    err= norm(T_pj(:)-T_pjold(:),Inf); 
    k=k+1;  
end
rmspj=rms(T_pj(:)-T_pjold(:));
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
figure,surf(x,y,T_pj');
xlabel('x');ylabel('y');zlabel('Temperature');
axis square;
title('T(x,y) Point Jacobi');
legend('Point Jacobi','Analytical','location','sw');
colorbar
pts={[0.5 0.25] [0.5 0.75] [1.5 0.25] [1.5 0.75]};

figure(2)
contourf(x,y,T_pj');
hold on
xlabel('x');ylabel('y');
axis square;
title('T(x,y) Jacobi');
legend('Point-Jacobi','Analytical','location','best');
colorbar
