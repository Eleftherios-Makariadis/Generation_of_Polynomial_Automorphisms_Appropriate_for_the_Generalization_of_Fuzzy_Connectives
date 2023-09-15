%Creating the points that will be used in the Lagrange
%Interpolation as X and Y vectors
X = zeros([1,3]);
Y = zeros([1,3]);
X(1)=0;
X(3)=1;
flag=true;
while (flag)
X(2)=input("Please provide an abscissa value between 0 and 1 to be used in the Lagrange Interpolation (0<x<1):  ");
flag= (X(2)<=0)||(X(2)>=1);
end
for i=1:3
    Y(i)=X(i)*sqrt(X(i));
end
%Lagrange Interpolation algorithm
L=zeros(3, 3);
for k = 1: 3
         V =1;
        for j=1: 3
            if k~=j
            V=conv (V,poly (X(j))) /(X(k) -X(j));
            end
        end
        L(k,:)=V;
end
% Calculation of the coefficients of the Lagrange interpolating polynomial
C = Y*L;
% Display of messages to user
disp("The following vector contains the list of abscissas used in the Lagrange Interpolation:");
disp (X);
disp("The following vector contains a list of ordinates used in the Lagrange Interpolation:");
disp (Y);
disp("The following vector contains the coefficients of the Lagrange interpolatory polynomial:");
disp (C);
% Creation of Lagrange interpolating polynomial
syms f(x);
f(x)= poly2sym([C]);
disp("The Lagrange interpolatory polynomial is:");
% Display Lagrange interpolating polynomial
disp (f(x));
% Display of automorphism properties to user
disp("The automorphism properties are:")
fprintf("f(0)=");
disp(f(0));
fprintf("f(1)=");
disp(f(1));
syms df(x);
df(x)= diff(f(x));
figure('Name','Graph of the derivative of the interpolatory polynomial','NumberTitle','off');
fplot (df(x), [0 1], "b");
% Creation of Inverse Lagrange interpolating polynomial
syms g(x);
syms y
g(x) = finverse(f);
% Display of the inverse Lagrange interpolating polynomial
disp("The inverse Lagrange interpolatory polynomial is:");
disp (g(x));
% Ploting of the polynomial automorphism and it's inverse
figure('Name',' Graph of the interpolatory and inverse interpolatory polynomials','NumberTitle','off');
fplot (f(x), [0 1], "b");
axis([0 1 0 1]);
hold on
fplot (g(x), [0 1], "g");
hold off
grid on
% Display of the menu to the user
disp("Select one of the following fuzzy connectives for generalization:");
disp("1. Negation");
disp("2. T-Norm");
disp("3. S-Norm");
disp("4. I-Implication");
disp("5. Exit");
flag= true;
while (flag)
ch = input("Provide your selection: ");
flag=~((ch==1)||(ch==2)||(ch==3)||(ch==4)||(ch==5));
end
switch ch
    case 1
        % Creation of natural negation
        syms h(x);
        h(x)= 1-x ;
        % Creation of generalized strong fuzzy negation
        syms N(x);
        N = compose( g, compose(h,f));
        % Ploting of strong fuzzy negation
        figure('Name','Graph of the strong fuzzy negation','NumberTitle','off');
        fplot (N, [0 1], "b");
        % Creation of involutive property
        syms N1(x);
        N1= compose( N , N );
        % Ploting of involutive property
        figure('Name','Graph of the involutive property','NumberTitle','off');
        fplot (N1, [0 1], "r");
    case 2
        % Creation of Gödel t-Norm
        syms t(x,y);
        t(x,y)=min(x,y);
        % Creation of generalized T-Norm
        syms T(x,y)
        T(x,y)= compose(g,compose(t(x,y),f));
        % Ploting of generalized T-Norm
        [x, y] = meshgrid(0:0.01:1, 0:0.01:1);
        z = zeros(size(x));
        for i = 1:numel(x)
             z(i) = T(x(i),y(i));
        end
        figure('Name','Graph of T-Norm','NumberTitle','off');
        surf(x, y, z);
    case 3
        % Creation of Maximum s-Norm
        syms s(x,y);
        s(x,y)=max(x,y);
        % Creation of generalized S-Norm
        syms S(x,y)
        S(x,y)= compose(g,compose(s(x,y),f));
        % Ploting of generalized S-Norm
        [x, y] = meshgrid(0:0.01:1, 0:0.01:1);
        z = zeros(size(x));
        for i = 1:numel(x)
             z(i) = S(x(i),y(i));
        end
        figure('Name','Graph of S-Norm','NumberTitle','off');
        surf(x, y, z);
    case 4
         % Creation of Łukasiewicz I-Implication
        syms L(x,y);
        L(x,y)=min(1,1-x+y);
        % Creation of generalized I-Implication
        syms I(x,y)
        I(x,y)= compose(g,compose(L(x,y),f));
        % Ploting of generalized I-Implication
        [x, y] = meshgrid(0:0.01:1, 0:0.01:1);
        z = zeros(size(x));
        for i = 1:numel(x)
             z(i) = I(x(i),y(i));
        end
        figure('Name','Łukasiewicz I-Implication','NumberTitle','off');
        surf(x, y, z);
end