clear all
data=open('scaled.dat');
data=data.scaled;
wn=data(:,1);
G=data(:,2:3);

beta=round(pi/wn(1));

nt=500;
t=linspace(eps,beta-eps,nt)';


%numerial part
Gnt=zeros(nt,1);
for j=1:nt
   Gnt(j)=sum(G(:,1).*cos(wn(:)*t(j))+G(:,2).*sin(wn(:)*t(j))); 
end
Gnt=2/beta*Gnt;

%asymptotic part: G(i*w_n) \approx  -i/w_n. 
%Add this contribution from Matsubara frequencies higher than have in input file
%G(t)= 2/beta*Sum_{n=n0}^{\inf} Re[G_n]*cos(wn*t)+Im[G_n]*sin(wn*t) \approx
% =  2/beta*Sum_{n=n0}^{\inf} -sin(wn*t)/wn =
% = -2/beta*Sum_{n=n0}^{\inf}  sin(wn*t)/wn =
% = -2/beta*(Sum_{n=1}^{\inf} sin(wn*t)/wn - Sum_{n=1}^{n0-1} sin(wn*t)/wn )
% = -2/beta*(1I*beta/(2*pi)*(Arctanh(exp(-1I*pi*t/beta))-Arctanh(exp(1I*pi*t/beta))) - Sum_{n=1}^{n0-1} sin(wn*t)/wn )                         )
Gasinf=zeros(nt,1);
Gasn0=zeros(nt,1);
n0=size(wn,1)+1; 
for j=1:nt
    Gasinf(j)=1i*beta/(2*pi)*( atanh(exp(-1i*pi*t(j)/beta))-atanh(exp(1i*pi*t(j)/beta)) );
    for k=1:n0-1
        Gasn0(j)=Gasn0(j)+sin((2*k-1)*pi/beta*t(j))/((2*k-1)*pi/beta);
    end
end
Gasinf=-2/beta*Gasinf;
Gasn0=-2/beta*Gasn0;
Gas=Gasinf-Gasn0;
Gas=real(Gas); %to remove round off errors

%Sum of domain contributions: 1:n0-1 and n0:inf
Gt=Gnt+Gas;

close all
figure(1)
plot(t,Gnt,'-o')
hold on;
plot(t,Gas,'-o')
plot(t,Gt,'-o');
legend('n=1:n0-1','n=n0:inf, approx','n=1:inf, approx')


%save t and G(t)
header1 = '#   tau';
header2 = '                G(tau)';
fn='scaledImagTime.dat';
fid=fopen(fn,'w');
fprintf(fid, [ header1 ' ' header2 '\n']);
%fprintf(fid, '%d %d \n', [t Gt]');
fclose(fid);
x=[t,Gt];
save(fn,'x','-ascii','-double','-append')


