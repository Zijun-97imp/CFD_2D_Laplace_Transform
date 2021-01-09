%truncation error plotting
%Error through the step number
%step number from: 0 to 10^3

%% clear the workspace
clc
clear





%% 10

% Parameters
nx = 11;
ny = nx;                       % number of space steps
x=linspace(0,1,nx);            % x range
y=linspace(0,1,ny);            % y range
dx = 0.1;
dy = dx;
B=(dx/dy)^2;
T_gs_10=zeros(nx,ny);             % initialize T for Point Jacobi
tol=1e-5;                      % error tolerance
err=1;                         % error
k=1;                           % iteration index
omega=1.95;                    % (optimal) relaxation parameter
%%
% Boundary Conditions
T_gs_10(1,:) = 1;%left
T_gs_10(nx,:) = cos(6*1.5*pi*y)+1;%right
T_gs_10(:,1) = 1+x;%bottom
T_gs_10(:,end) = 1;%top
Texact_10=T_gs_10;     %analytical 

%%
% Time loop
i=1;j=1;
rms=0;
while err>tol
    T_gsold_10=T_gs_10;
    for i=2:nx-1
        for j=2:ny-1
            %Gauss-Seidel
            T_gs_10(i,j)= (1-omega)*T_gsold_10(i,j) +(omega/(2*(1+B)))*(T_gs_10(i-1,j)+T_gsold_10(i+1,j)+B*(T_gs_10(i,j-1)+T_gsold_10(i,j+1)));
        end
    end
    %boundary conditions
    T_gs_10(1,:)= 1;%left
    T_gs_10(nx,:) = cos(6*1.5*pi*y)+1;%right
    T_gs_10(:,1) = 1+x;%bottom
    T_gs_10(:,end) = 1;%top
    T_gs_10(5,5)=1.5;
    T_gs_10(2,2)=0.5;
    
    err_10= max(max(abs(T_gs_10-T_gsold_10)));
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
figure,surf(x,y,T_gs_10');
xlabel('x');ylabel('y');zlabel('Temperature');
axis square;
title('T(x,y) Gauss Seidel');
legend('Gauss Seidel','Analytical','location','sw');
colorbar
pts={[0.5 0.25] [0.5 0.75] [1.5 0.25] [1.5 0.75]};









%% 100
% Parameters
nx = 101;
ny = nx;                       % number of space steps
x=linspace(0,1,nx);            % x range
y=linspace(0,1,ny);            % y range
dx = 0.01;
dy = dx;
B=(dx/dy)^2;
T_gs_100=zeros(nx,ny);             % initialize T for Point Jacobi
tol=1e-5;                      % error tolerance
err=1;                         % error
k=1;                           % iteration index
omega=1.95;                    % (optimal) relaxation parameter
%%
% Boundary Conditions
T_gs_100(1,:) = 1;%left
T_gs_100(nx,:) = cos(6*1.5*pi*y)+1;%right
T_gs_100(:,1) = 1+x;%bottom
T_gs_100(:,end) = 1;%top
Texact_100=T_gs_100;     %analytical 

%%
% Time loop
i=1;j=1;
rms=0;
while err>tol
    T_gsold_100=T_gs_100;
    for i=2:nx-1
        for j=2:ny-1
            %Gauss-Seidel
            T_gs_100(i,j)= (1-omega)*T_gsold_100(i,j) +(omega/(2*(1+B)))*(T_gs_100(i-1,j)+T_gsold_100(i+1,j)+B*(T_gs(i,j-1)+T_gsold_100(i,j+1)));
        end
    end
    %boundary conditions
    T_gs_100(1,:)= 1;%left
    T_gs_100(nx,:) = cos(6*1.5*pi*y)+1;%right
    T_gs_100(:,1) = 1+x;%bottom
    T_gs_100(:,end) = 1;%top
    T_gs_100(50,50)=1.5;
    T_gs_100(20,20)=0.5;
    
    err_100= max(max(abs(T_gs_100-T_gsold_100)));
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
figure,surf(x,y,T_gs');
xlabel('x');ylabel('y');zlabel('Temperature');
axis square;
title('T(x,y) Gauss Seidel');
legend('Gauss Seidel','Analytical','location','sw');
colorbar
pts={[0.5 0.25] [0.5 0.75] [1.5 0.25] [1.5 0.75]};










%% 500

% Parameters
nx = 501;
ny = nx;                       % number of space steps
x=linspace(0,1,nx);            % x range
y=linspace(0,1,ny);            % y range
dx = 0.002;
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
    T_gs(250,250)=1.5;
    T_gs(100,100)=0.5;
    
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
figure,surf(x,y,T_gs');
xlabel('x');ylabel('y');zlabel('Temperature');
axis square;
title('T(x,y) Gauss Seidel');
legend('Gauss Seidel','Analytical','location','sw');
colorbar
pts={[0.5 0.25] [0.5 0.75] [1.5 0.25] [1.5 0.75]};

