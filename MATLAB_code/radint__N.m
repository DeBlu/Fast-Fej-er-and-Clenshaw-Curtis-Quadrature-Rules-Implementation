function radint_spd2_vsN_a_NEW
format long e
hold off
myfile = fopen('radint_sp.dat', 'wt');
exact=(1.d0+1.d0/sqrt(10.d0)+1.d0/sqrt(100.d0))*sqrt(pi)/4.d0;
% nucrxspeed2 is versus the scale factor for 3 values of N - fixed b
yerr=[];
int=[];
str1='-ok';
str2='-^k';
str3='-sk';
str4='-dk';
symbol=[str1 str2 str3 str4];
%sc=[.5 .6 .7 .8];
%sc=[.6 .7 .8 1];
sc=[1 1.25 1.5 1.75];
%sc=[.3 .4 .45 .5];
% Coefficients
A = [1 1 1 1 10 10 10 10 10 10 10 100 100 100];
B = [0.5 0.5 1 50 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 50 50];
C = [0.1 10 1 100 0.1 0.1 1 1 10 10 0.5 0.1 1 1];
p = [2 1 2 2 1 2 1 2 1 2 1 1 1 2];
q = [4 4 4 4 4 4 4 4 4 4 4 4 4 4];
r = [6 7 6 6 7 6 7 6 7 6 7 7 7 6];
kmax=0;
%This is the small loop through 4 scales
yns=[];
for k=1:1:4
    kmin=kmax+1;
    kmax=kmin+2;
    nx=[];
    int=[];
    yerr=[];
    scale=sc(k);
% This is the larger loop through N
%    for m=4:1:16
    %for m=4:1:16
    for m=4:1:30
        nx=[nx m];
%    fprintf(myfile,'%8.2f\n',scale);
    ptwt=abspeed2(m);
x=ptwt(:,1);
w=ptwt(:,2);
n1=length(x);
s=0.;
bigw=w./(x.^2.*exp(-x.*x));
ws=scale*bigw;
xs=scale*x;
%fx=xs.^2.*(exp(-xs.^2)+10*exp(-10*xs.^2)+100*exp(-100*xs.^2));
%fx=xs.^2.*(exp(-xs.^2)+10*exp(-10*xs.^2));
% Choose scheme
sch = 13;   
for ii = sch:sch
    fx=xs.^2.*(A(ii)*exp(-(xs-p(ii)).^2)+B(ii)*exp(-(xs-q(ii)).^2)+C(ii)*exp(-(xs-r(ii)).^2));
%fx=xs.^2.*(A(1)*exp(-(xs-p(1)).^2)+B(1)*exp(-(xs-q(1)).^2)+C(1)*exp(-(xs-r(1)).^2));
%fx=xs.^2.*(A(1)*exp(-(xs-p(1)).^2)+B(1)*exp(-(xs-q(1)).^2)+C(1)*exp(-(xs-r(1)).^2));
%fx=xs.^2.*(A(2)*exp(-(xs-p(2)).^2)+B(2)*exp(-(xs-q(2)).^2)+C(2)*exp(-(xs-r(2)).^2));
%fx=xs.^2.*(A(3)*exp(-(xs-p(3)).^2)+B(3)*exp(-(xs-q(3)).^2)+C(3)*exp(-(xs-r(3)).^2));
%fx=xs.^2.*(A(4)*exp(-(xs-p(4)).^2)+B(4)*exp(-(xs-q(4)).^2)+C(4)*exp(-(xs-r(4)).^2));
%fx=xs.^2.*(A(5)*exp(-(xs-p(5)).^2)+B(5)*exp(-(xs-q(5)).^2)+C(5)*exp(-(xs-r(5)).^2));
%fx=xs.^2.*(A(6)*exp(-(xs-p(6)).^2)+B(6)*exp(-(xs-q(6)).^2)+C(6)*exp(-(xs-r(6)).^2));
%fx=xs.^2.*(A(7)*exp(-(xs-p(7)).^2)+B(7)*exp(-(xs-q(7)).^2)+C(7)*exp(-(xs-r(7)).^2));
%fx=xs.^2.*(A(8)*exp(-(xs-p(8)).^2)+B(8)*exp(-(xs-q(8)).^2)+C(8)*exp(-(xs-r(8)).^2));
%fx=xs.^2.*(A(9)*exp(-(xs-p(9)).^2)+B(9)*exp(-(xs-q(9)).^2)+C(9)*exp(-(xs-r(9)).^2));
%fx=xs.^2.*(A(10)*exp(-(xs-p(10)).^2)+B(10)*exp(-(xs-q(10)).^2)+C(10)*exp(-(xs-r(10)).^2));
%fx=xs.^2.*(A(11)*exp(-(xs-p(11)).^2)+B(11)*exp(-(xs-q(11)).^2)+C(11)*exp(-(xs-r(11)).^2));
%fx=xs.^2.*(A(12)*exp(-(xs-p(12)).^2)+B(12)*exp(-(xs-q(12)).^2)+C(12)*exp(-(xs-r(12)).^2));
%fx=xs.^2.*(A(13)*exp(-(xs-p(13)).^2)+B(13)*exp(-(xs-q(13)).^2)+C(13)*exp(-(xs-r(13)).^2));
%fx=xs.^2.*(A(14)*exp(-(xs-p(14)).^2)+B(14)*exp(-(xs-q(14)).^2)+C(14)*exp(-(xs-r(14)).^2));
s=sum(ws.*fx);
int=[int s];
%exact=(1.d0+1.d0/sqrt(10.d0)+1.d0/sqrt(100.d0))*sqrt(pi)/4.d0;
%exact=(1.d0+1.d0/sqrt(10.d0))*sqrt(pi)/4.d0;
exact = [29.06790389 894.6209168 101.9157572 7939.706686 49.73155391 100.8492314 128.6943730 159.0743404 918.3225636 741.3254305 84.82614016 286.7480217 1813.362524 2324.539299];
%err=log10(abs(s-exact(7))/exact(7));
%fprintf('%i %16.8f\n',m,s-exact(7))

err=log10(abs(s-exact(ii))/exact(ii));
fprintf('%i %16.8f\n',m,s-exact(ii))
end
%pause
yerr=[yerr err];
    end
    %ns=length(s)
%line=exact*ones(1,ns);
pcoeff=polyfit(nx,yerr,1)
for kkk=1:2
fprintf(myfile,'%13.5f %13.5\n',pcoeff(kkk)); end
%pause
yn=polyval(pcoeff,nx);
yns=[yns,yn'];
plot(nx,yerr,symbol(kmin:kmax),'markersize',10,'markerfacecolor','k')
%plot(nx,yerr,'-k')
hold on
xlabel('${\rm N}$','Interpreter','latex','fontsize',32)
ylabel('$\log_{10}[{\rm Relative}\;\; {\rm Error}]$','Interpreter','latex','fontsize',32)
%ylabel('Integral','Interpreter','latex','fontsize',20)
%axis([4 16 -9 0])
axis([4 30 -12 1])
set(gca,'FontSize',36)
%set(gca,'Ytick',[-3:2:1],'linewidth',1.6)
set(gca,'Ytick',[-12:2:1],'linewidth',1.6)
%set(gca,'Xtick',[4:2:24],'linewidth',1.6)
set(gca,'Xtick',[4:2:30],'linewidth',1.6)
end
plot(nx,line,'--k','linewidth',1.6)
%str={'N = 20'}; 
%text(0.6,.55,str,'Interpreter','latex','fontsize',24)
%title('fx=x^2*(exp(-x^2)+10*exp(-10*x^2)+100*exp(-100*x^2)','fontsize',14)
%legend('s=0.6','s=0.7','s=0.8','s=1','Location','SouthWest')
legend(strcat('s=',num2str(sc(1))),...
    strcat('s=',num2str(sc(2))),...
    strcat('s=',num2str(sc(3))),...
    strcat('s=',num2str(sc(4))),...
    'Location',...
    'SouthWest')
for ii = sch:sch
str={strvcat(strcat('A=',num2str(A(ii))), strcat('B=',num2str(B(ii))),...
    strcat('C=',num2str(C(ii))),...
strcat('a=',num2str(p(ii))),...
strcat('b=',num2str(q(ii))),...
strcat('c=',num2str(r(ii))))};
end
%str={'(B)'};
text(24,-2.9,str,'Interpreter','Latex','Fontsize',36) 

plot(nx,yns(:,1),'--k','linewidth',1.6)
plot(nx,yns(:,2),'--k','linewidth',1.6)
plot(nx,yns(:,3),'--k','linewidth',1.6)
plot(nx,yns(:,4),'--k','linewidth',1.6)
set(gcf, 'Units','centimeters','Papersize',[36,36])
set(gcf, 'Units','centimeters','Position',[3 3 24 20])
%
% Speed points and weights for p = 1
%
function ptwt = abspeed2(m)
load abspd.dat; n=90;
a=abspd(1:m,1);
b=abspd(1:m,2);
rtb=sqrt(b);
rtb(m)=[];
t=diag(rtb,-1)+diag(a)+diag(rtb,1);
[f,lambda]=eig(t);
pt=diag(lambda);
for i=1:m
wt(i)=sqrt(pi)*f(1,i)^2/4.d0; end
ptwt=[pt,wt'];