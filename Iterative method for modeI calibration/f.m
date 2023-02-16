
function [S_1,S_2,S_3,S_4]=f(K_O,K_S,K_UAP,K_NO,k_STO,n_NO,b_HO,Y_STOO,i_NSS,k_BAPO)
%%%%%%%%%%%%输入动力学参数的典型值%%%%%%%%%%%%%
% K_O范围0.05-0.4%%%;
% K_O=0.2;
% K_S范围0.5-5%%%;
% K_S=2;
% K_UAP范围50-300%%%;
% K_UAP=100;
K_BAP=85;
K_STO=1;
u_H=2;
K_NH=0.01;
K_ALK=0.1;
% K_NO范围0.1-1%%%;
% K_NO=0.5;
k_H=2;
K_X=1;
% k_STO范围0.5-10%%%;
% k_STO=5;
k_USTO=1.53;
k_BSTO=0.085;
% n_NO范围0.2-1%%%;
% n_NO=0.6;
% b_HO范围0.1-1%%%;
% b_HO=0.2;
b_HNO=0.05;
b_STOO=0.1;
b_STONO=0.05;
u_A=0.35;
K_ANH=1;
K_AO=0.5;
K_ANO=0.5;
K_AALK=0.5;
b_AO=0.05;
b_ANO=0.02;
k_a=0.08;

%%%%%%%%%%%%输入化学计量学系数的典型值%%%%%%%%%%%%%
%Y_STOO=0.85;%
%Y_STOO=0.935;%
% Y_STOO范围0.5-5%%%;
% Y_STOO=0.85;
Y_STONO=0.80;
%Y_STONO=0.88;%
Y_HO=0.63;
%Y_HO=0.693;%
Y_HNO=0.54;
%Y_HNO=0.594;%
%Y_A=0.24;%
%Y_A=0.264;%
Y_A=0.24;
% i_NSS范围0.01-1%%%;
% i_NSS=0.03;
%i_NSS=0.033;%
i_NXI=0.02;
%i_NXI=0.022;%
i_NBM=0.07;
%i_NBM=0.077;%
i_TSBM=0.9;
%i_TSBM=0.99;%
i_TSSTO=0.60;
%i_TSSTO=0.66;%
i_NUAP=0.03;
%i_NUAP=0.033;%
i_NBAP=0.02;
%i_NBAP=0.022;%
k_UAPO=0.12;
%k_UAPO=0.132;%
k_UAPNO=0.12;
%k_UAPNO=0.132;%
% k_BAPO范围0.01-2%%%;
% k_BAPO=0.09;
%k_BAPO=0.099;%
k_BAPNO=0.09;
%k_BAPNO=0.099;%
k_UAPAO=0.1;
%k_UAPAO=0.11;%
k_BAPAO=0.09;
%k_BAPAO=0.099;%
k_BAPANO=0.09;
%k_BAPANO=0.099;%
%%%%%%%%%%%%计算未知化学计量学系数的典型值%%%%%%%%%%%%%
x_2=Y_STOO-1;
x_21=Y_STOO-1;
x_22=Y_STOO-1;
x_3=(Y_STONO-1)/2.86;
x_31=(Y_STONO-1)/2.86;
x_32=(Y_STONO-1)/2.86;
x_4=1-1/Y_HO+k_UAPO;
x_5=(k_UAPNO+1-1/Y_HNO)/2.86;
x_6=k_BAPO-1;
x_7=(k_BAPNO-1)/2.86;
x_8=-1;
x_9=-1/2.86;
x_10=k_UAPAO+1-4.57/Y_A;
x_11=k_BAPAO-1;
x_12=(k_BAPANO-1)/2.86;

y_2=i_NSS;
y_21=i_NUAP;
y_22=i_NBAP;
y_3=i_NSS;
y_31=i_NUAP;
y_32=i_NBAP;
y_4=-i_NBM-i_NUAP*k_UAPO;
y_5=-i_NBM-i_NUAP*k_UAPNO;
y_6=i_NBM-i_NBAP*k_BAPO;
y_7=i_NBM-i_NBAP*k_BAPNO;
y_8=0;
y_9=0;
y_10=-1/Y_A-i_NBM-i_NUAP*k_UAPAO;
y_11=i_NBM-i_NBAP*k_BAPAO;
y_12=i_NBM-i_NBAP*k_BAPANO;

z_2=y_2/14;
z_21=y_21/14;
z_22=y_22/14;
z_3=(y_3-x_3)/14;
z_31=(y_31-x_31)/14;
z_32=(y_32-x_32)/14;
z_4=y_4/14;
z_5=(y_5-x_5)/14;
z_6=y_6/14;
z_7=(y_7-x_7)/14;
z_8=0;
z_9=-x_9/14;
z_10=(y_10-1/Y_A)/14;
z_11=y_11/14;
z_12=(y_12-x_12)/14;

%%%%%%%%%%%%%输入反应器运行初始值%%%%%%%%%%%%
%%%%溶解氧(mg/L)%%%
S_O=[6.6 7 6.8 6.7 7.5];
%%基质浓度(mg/L)%%
S_S=[9.69 7.93 7.33 6.92 6.69];
%%%进水UAP浓度%%%%%
S_UAP=[0.98 1.65 2.30 3.06 4.11];
%%进水BAP浓度(mg/L)%%
S_BAP=[0.89 0.91 0.93 0.95 0.95];
%%%氨氮浓度%%%
S_NH=[8.41 8.14 7.40 6.24 5.02];
%%%亚硝氮硝氮浓度%%%
S_NO=[0.81 1.18 2.14 3.45 4.86];
%%%碱度%%%
S_ALK=[0.01 0.03 0.05 0.05 0.1];
%%异养菌(gX/m3)%%
X_H=[5938.5832 4375.3694 1345.9254 266.19069 56.494469];
%%%自养菌((mgX/L))%%
X_A=[4937.1518 5695.37 4402.0802 3992.0263 2792.156];
%%%胞内贮存物浓度%%%%
X_STO=[7.82 5.36 4.11 2.90 1.62];
%X_STO=[5.5 5.36 4.11 3 2.5];%
%%%溶解性有机氮浓度%%%
S_ND=[0.91 1.00 1.23 1.32 1.36];
%%%水力停留时间(d)%%%
T=0.00212;

%%%%各组分的动力学模型%%%%%
p_2=k_STO.*S_O.*S_S.*X_H./(K_O+S_O)./(K_S+S_S);
p_21=k_USTO.*S_O.*S_UAP.*X_H./(K_O+S_O)./(K_UAP+S_UAP);
p_22=k_BSTO.*S_O.*S_BAP.*X_H./(K_O+S_O)./(K_BAP+S_BAP);
p_3=k_STO*n_NO.*K_O.*S_NO.*S_S.*X_H./(K_O+S_O)./(K_NO+S_NO)./(K_S+S_S);
p_31=k_USTO.*n_NO.*K_O.*S_NO.*S_UAP.*X_H./(K_O+S_O)./(K_NO+S_NO)./(K_UAP+S_UAP);
p_32=k_BSTO*n_NO.*K_O.*S_NO.*S_BAP.*X_H./(K_O+S_O)./(K_NO+S_NO)./(K_BAP+S_BAP);
p_4=u_H.*S_O.*S_NH.*S_ALK.*X_STO./(K_O+S_O)./(K_NH+S_NH)./(K_ALK+S_ALK)./(K_STO+X_STO./X_H);
p_5=u_H*n_NO.*S_O.*S_NO.*S_NH.*S_ALK.*X_STO./(K_O+S_O)./(K_NO+S_NO)./(K_NH+S_NH)./(K_ALK+S_ALK)./(K_STO+X_STO./X_H);
p_6=b_HO.*S_O.*X_H./(K_O+S_O);
p_7=b_HNO.*K_O.*S_NO.*X_H./(K_O+S_O)./(K_NO+S_NO);
p_8=b_STOO.*S_O.*X_STO./(K_O+S_O);
p_9=b_STONO.*K_O.*S_NO.*X_STO./(K_O+S_O)./(K_NO+S_NO);
p_10=u_A.*S_O.*S_NH.*S_ALK.*X_A./(K_AO+S_O)./(K_ANH+S_NH)./(K_AALK+S_ALK);
p_11=b_AO.*S_O.*X_A./(K_AO+S_O);
p_12=b_ANO.*K_O.*S_NO.*X_A./(K_AO+S_O)./(K_ANO+S_NO);
p_13=k_a.*S_ND;


%%%%%%%%%%%%各组分的生成速率1-COD,2-NH,3-NO,4-ND%%%%%%%%%%%
%%%Ui=k_UAP*q_m.*S_BOM.*X./ai-q_UAP.*S_UAP.*X_H./bi;%%%%
%%%%Bi=k_BAP.*X-q_BAP.*S_BAP.*X_H./ci;%%%%
r_1=(Y_STOO-1).*p_2+(Y_STOO-1).*p_21+(Y_STOO-1).*p_22+(Y_STONO-1).*p_3+(Y_STONO-1).*p_31+(Y_STONO-1).*p_32+(k_UAPO+1-1/Y_HO).*p_4+(k_UAPNO+1-1/Y_HNO).*p_5+(k_BAPO-1).*p_6+(k_BAPNO-1).*p_7-p_8-p_9+(k_UAPAO+1).*p_10+(k_BAPAO-1).*p_11+(k_BAPANO-1).*p_12;
r_2=y_2.*p_2+y_21.*p_21+y_22.*p_22+y_3.*p_3+y_31.*p_31+y_32.*p_32+y_4.*p_4+y_5.*p_5+y_6.*p_6+y_7.*p_7+y_10.*p_10+y_11.*p_11+y_12.*p_12+p_13;
r_3=x_3.*p_3+x_31.*p_31+x_32.*p_32+x_5.*p_5+x_7.*p_7+x_9.*p_9+p_10./Y_A+x_12.*p_12;
r_4=(-p_21-p_31+k_UAPO.*p_4+k_UAPNO.*p_5+k_UAPAO.*p_10-p_22-p_32+k_BAPO.*p_6+k_BAPNO.*p_7+k_BAPAO.*p_11+k_BAPANO.*p_12).*0.7-p_13;

%%%%反应器模型%%%%
R_1=r_1.*T;
R_2=r_2.*T;
R_3=r_3.*T;
R_4=r_4.*T;
S_0C=11.13;
S_0NH=7.67;
S_0NO=0.68;
S_0ND=0.67;

S_1C=S_0C+R_1(1);
S_2C=S_1C+R_1(2);
S_3C=S_2C+R_1(3);
S_4C=S_3C+R_1(4);
S_5C=S_4C+R_1(5);
S_1=[S_1C S_2C S_3C S_4C S_5C];

S_1NH=S_0NH+R_2(1);
S_2NH=S_1NH+R_2(2);
S_3NH=S_2NH+R_2(3);
S_4NH=S_3NH+R_2(4);
S_5NH=S_4NH+R_2(5);
S_2=[S_1NH S_2NH S_3NH S_4NH S_5NH];

S_1NO=S_0NO+R_3(1);
S_2NO=S_1NO+R_3(2);
S_3NO=S_2NO+R_3(3);
S_4NO=S_3NO+R_3(4);
S_5NO=S_4NO+R_3(5);
S_3=[S_1NO S_2NO S_3NO S_4NO S_5NO];

S_1ND=S_0ND+R_4(1);
S_2ND=S_1ND+R_4(2);
S_3ND=S_2ND+R_4(3);
S_4ND=S_3ND+R_4(4);
S_5ND=S_4ND+R_4(5);
S_4=[S_1ND S_2ND S_3ND S_4ND S_5ND];

end


