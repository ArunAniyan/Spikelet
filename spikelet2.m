%% Matching Wavelet - Spikelet Algorithm - Generate Filter Coefficients
% Algorithm - Prof. Dr. Rodrigo Capobianco Guido, University of Sao Paulo,Brazil
% Matlab Program - Arun Kumar A , Santhom Computing Facility
% aka.bhagya@gmail.com
% 31-05-10
% Filter support - 6

%% Initializations

clc;
clear all               % Clear the Workspace


%% Load Data File

% Data file should contain a data vector of any format csv or text

file_name=input('Enter filename with extension : ','s');
data=load(file_name);                      % Data loaded into variable "data"

% Choice for saving coefficients and wavelet

choice=input('Do you want to save the generated coefficients and wavelets (Y/N) : ' ,'s');

if (choice=='Y')
    hname=input('Enter filename to store coefficients (with extension) : ','s');
    wtname=input('Enter filename to store wavelet (with extension) : ','s');
end


%% Declarations and Initial Variables


filt_sup=6;            % Filter Support

l=length(data);        % Length of data vector  



%% Set filter coefficients for Daubechies ( Now set for testing  ) - Stage A1-A3

% h0
  a1 = 1;
  a2 = 0;
  a3 = 0;
  a4 = 1;
  
% h1
  b1 = -1;
  b2 = -1;
  b3 = -1;
  b4 = 1;
  
% h2
  c1 = 1;
  c2 = 2;
  c3 = 4;
  c4 = 1;
  
% h3
  d1 = -1;
  d2 = -3;
  d3 = -9;
  d4 = 1;
  
% h4
  e1 = 1;
  e2 = 4;
  e3 = 16;
  e4 = 1;

% h5
  f1 = -1;
  f2 = -5;
  f3 = -25;
  f4 = 1;
  
  
  
cof1=[f1,e1,d1,c1,b1,a1;f2,e2,d2,c2,b1,a2;f3,e3,d3,c3,b3,a3;f4,e4,d4,c4,b4,a4];

%% Stage B1

norm_fac=abs(min(data));     % Normalisation factor

norm_vec= data/norm_fac;     % Normalisation

f=norm_vec;


%% Stages B2 & B3

cof2=zeros(1,filt_sup);

for j=0:filt_sup-1
    
    for k=0:filt_sup-1
        
        for i=0:round(l/2)-1
            
            in1=rem((i*2+j),(l-1));
            in2=rem((i*2+k),(l-1));
                      
            r(i+1) = f(in1+1)*f(in2+1);
          
         end
       
      cof2(k+1,j+1)=sum(r);
    end

end



  
%% Stage C1

for i=0:filt_sup-1
    cof2(:,i+1)=((-1)^i*cof2(:,i+1));
end



B=vertcat(cof2,cof1);


C=[0;0;0;0;0;0;0;0;0;2];

%% Stage  C2

A=horzcat(cof2,cof1');


%% Stage C3

% Solve for Coefficients

h=flipud(B\C);

g=flipud(h);

for j=0:filt_sup-1
    g(j+1)=((-1)^j)*g(j+1);
    
end

disp('Low Pass Filter Coefficients')
disp('        ')
for i=1:filt_sup
    disp(sprintf(' h%i = %f',i-1,h(i)))
end

disp('-------------------')
disp('High Pas Filter Coefficients')
disp('           ')

for i=1:filt_sup
    disp(sprintf(' g%i = %f',i-1,g(i)))
end

%% Generate Scaling Function


phi(1)=0;
phi(7)=0;

A1=[h(2)-1 h(1);h(4) h(3)-1; 1 1];
B2=[0;0;1];

C=A1\B2;

phi(3)=C(1);
phi(5)=C(2);

phi(2)=h(1)*phi(3);
phi(4)=(h(2)*phi(5))+(h(3)*phi(3));
phi(6)=h(4)*phi(5);
%----------------------------------------------------------|
% In matlab, indexing different that from other languages. |
% In the actual paper the index values differ from what is | 
% written in the program. The actual indexes are as below. |
% The first phi is the notation in the code and the 2nd of |
% the paper.                                               |  
% phi(1) => phi(0), phi(2) => phi(1/2), phi(3) => phi(1)   | 
% phi(4) => phi(3/2), phi(5) => phi(2), phi(6) => phi(5/2) |
% phi(7) => phi(3) .                                       |
% ---------------------------------------------------------|


%% Wavelet Function

psi(1)=g(1)*phi(1);
psi(3)=(g(1)*phi(5))+(g(2)*phi(3))+(g(3)*phi(1));
psi(5)=(g(2)*phi(7))+(g(3)*phi(5))+(g(4)*phi(3));
psi(7)=g(4)*phi(7);


psi(2)=g(1)*phi(3);
psi(4)=(g(2)*phi(5))+(g(3)*phi(3));
psi(6)=g(4)*phi(5);


%% Plot all important stuff

subplot(2,2,1),plot(data),title('Original data');
subplot(2,2,2),plot(f),title('Normalised data');
subplot(2,2,3),plot(phi),title('Scaling Function - phi');
subplot(2,2,4),plot(psi),title('Wavelet - psi');


%% Write Coefficients and Wavelet to file
if ( choice=='Y')
        
    dlmwrite(hname,h,'delimiter',',');
    disp(sprintf('Created file %s',hname));
    dlmwrite(wtname,psi,'delimiter',',');
    disp(sprintf('Created file %s',wtname));
    disp('All task Completed Successfully !');
    
else
    disp('All task Completed Successfully !');
    
end

    
    
    
%% End

