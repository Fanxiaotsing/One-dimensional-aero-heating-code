function [ rho,valueC,valueK,value_dC_dT,value_dK_dT] = Multilinear_rhoCK( temp,name )
%defining material
   value_dC_dT=0;
   value_dK_dT=0;

if ( strcmp(name,'G1'))
    %Material='Al2O3';       
   valueC=600;
 %%
    k_T=[
        20,	0.021;
       200,	0.025;
       400,	0.031;
       600,	0.041;
       800,	0.052;
      1000,	0.066;
      1200,	0.082;
      1400,	0.099;
      1600,	0.12
     ];
    [mk,nk]=size(k_T); 
    temp=temp-273.15;

    if( temp<=k_T(1,1) ) 
        valueK=k_T(1,2);
        
    end
    if ( temp>=k_T(mk,1) )
        valueK=k_T(mk,2);
    end
    
    for i=1:mk-1
        if(temp>k_T(i,1)&&temp<k_T(i+1,1) )
         value_dK_dT=( k_T(i+1,2)-k_T(i,2) )/( k_T(i+1,1)-k_T(i,1) );
         valueK     =k_T(i,2)+( temp-k_T(i,1) )*value_dK_dT;          
        end
        
    end
 %% 
    rho=350;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'G2'))
    %Material='SiO2-A';
valueC=750; 
 %%
    k_T=[
        20,	0.015;
       300,	0.022;
       400,	0.025;
       500,	0.029;
       600, 0.034;
       700,	0.04;
       800,	0.046;
       900,	0.053;
      1000,	0.062
     ];
    [mk,nk]=size(k_T); 
   temp=temp-273.15;

    if( temp<=k_T(1,1) ) 
        valueK=k_T(1,2);
        
    end
    if ( temp>=k_T(mk,1) )
        valueK=k_T(mk,2);
    end
    
    for i=1:mk-1
        if(temp>k_T(i,1)&&temp<k_T(i+1,1) )
         value_dK_dT=( k_T(i+1,2)-k_T(i,2) )/( k_T(i+1,1)-k_T(i,1) );
         valueK     =k_T(i,2)+( temp-k_T(i,1) )*value_dK_dT;          
        end
        
    end
 %% 
    rho=250;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'G3'))
    %Material='SiO2-B';
valueC=800; 
 %%
     k_T=[
     200,0.038;
     400,0.042;
     600,0.062;
     800,0.09;
    1000,0.12;
    1200,0.15;
    1400,0.18;
    1600,0.2;
    1800,0.22;
    2000,0.24;
    2200,0.27
     ];
    [mk,nk]=size(k_T); 
    temp=temp-273.15;

    if( temp<=k_T(1,1) ) 
        valueK=k_T(1,2);
        
    end
    if ( temp>=k_T(mk,1) )
        valueK=k_T(mk,2);
    end
    
    for i=1:mk-1
        if(temp>k_T(i,1)&&temp<k_T(i+1,1) )
         value_dK_dT=( k_T(i+1,2)-k_T(i,2) )/( k_T(i+1,1)-k_T(i,1) );
         valueK     =k_T(i,2)+( temp-k_T(i,1) )*value_dK_dT;          
        end
        
    end   
 %% 
    rho=420; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'G4'))
    %Material='石英酚醛-251厂';
valueC=1130;
valueK=0.453;
 %% 
    rho=1500; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'G5'))
    %Material='高硅氧酚醛-王思民论文';
    c_T=[
     20 ,1761;
    100,1746;
    200,1727;
    380,1692;
    400,1687;
    500,1590;
    600,1409;
    727,1143;
    800,1129.8;
    2000,1129.8
     ];
    [mc,nc]=size(c_T);
    %linear interpolation 
    temp=temp-273.15;

    if( temp<=c_T(1,1) ) 
        valueC=c_T(1,2);
        
    end
    if ( temp>=k_T(mc,1) )
        valueC=c_T(mc,2);
    end
    
    for i=1:mc-1
        if(temp>c_T(i,1)&&temp<c_T(i+1,1) )
         value_dC_dT=( c_T(i+1,2)-c_T(i,2) )/( c_T(i+1,1)-c_T(i,1) );
         valueC     =c_T(i,2)+( temp-c_T(i,1) )*value_dC_dT;          
        end
        
    end   
 %%
 valueK=0.453;
 %% 
    rho=1500; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'G6'))
    %Material='高硅氧酚醛-时圣波博士论文';
    c_T=[
     20 ,1761;
    100,1746;
    200,1727;
    380,1692;
    400,1687;
    500,1590;
    600,1409;
    727,1143;
    800,1129.8;
    2000,1129.8
     ];
    [mc,nc]=size(c_T);
    %linear interpolation 
    temp=temp-273.15;

    if( temp<=c_T(1,1) ) 
        valueC=c_T(1,2);
        
    end
    if ( temp>=k_T(mc,1) )
        valueC=c_T(mc,2);
    end
    
    for i=1:mc-1
        if(temp>c_T(i,1)&&temp<c_T(i+1,1) )
         value_dC_dT=( c_T(i+1,2)-c_T(i,2) )/( c_T(i+1,1)-c_T(i,1) );
         valueC     =c_T(i,2)+( temp-c_T(i,1) )*value_dC_dT;          
        end
        
    end 
 %%
valueK=0.453;
 %% 
    rho=1500; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'G7'))
    %Material='高硅氧酚醛-703所数据';
%%
valueC=1227;
valueK=0.55;
%% 
rho=1710; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'G8'))
    %Material='石英玻璃钢-703所数据';
%%设定炭化温度/摄氏度
T_cr=700;
temp=temp-273.15;
rhock_T=[
       100,1710,1100,0.50;
       150,1710,1130,0.52;
       260,1710,1200,0.56;
      2000,1710,1560,1.0
     ];
    [mck,nck]=size(rhock_T);
    %linear interpolation  

    if( temp<=rhock_T(1,1) )
        rho   =rhock_T(1,2);
        valueC=rhock_T(1,3);
        valueK=rhock_T(1,4);         
    end
    if( temp>=rhock_T(mck,1) )
        rho   =rhock_T(mck,2);   
        valueC=rhock_T(mck,3);
        valueK=rhock_T(mck,4);        
    end
    
    for i=1:mck-1
        if(temp>rhock_T(i,1)&&temp<rhock_T(i+1,1) )
         rho   =rhock_T(i,2);   
         value_dC_dT=( rhock_T(i+1,3)-rhock_T(i,3) )/( rhock_T(i+1,1)-rhock_T(i,1) );   
         valueC=rhock_T(i,3)+( temp-rhock_T(i,1) )*value_dC_dT;
         value_dK_dT=( rhock_T(i+1,4)-rhock_T(i,4) )/( rhock_T(i+1,1)-rhock_T(i,1) );
         valueK=rhock_T(i,4)+( temp-rhock_T(i,1) )*value_dK_dT; 
         
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'G9'))
    %Material='C/C-ZrC-SiC';
    ck_T=[
        25,670,12.98;
       200,1016,16.61;
       400,1275,18.36;
       600,1420,18.65;
       800,1450,18.17;
      1000,1450,21.523;
      1400,1450,23.166;
      1800,1450,25.822        
     ];
    [mck,nck]=size(ck_T);
    %linear interpolation 
  temp=temp-273.15;

    if( temp<=ck_T(1,1) ) 
        valueC=ck_T(1,2);
        valueK=ck_T(1,3);
        
    end
    if ( temp>=ck_T(mck,1) )
        valueC=ck_T(mck,2);
        valueK=ck_T(mck,3);
    end
    
    for i=1:mck-1
        if(temp>ck_T(i,1)&&temp<ck_T(i+1,1) )
         value_dC_dT=( ck_T(i+1,2)-ck_T(i,2) )/( ck_T(i+1,1)-ck_T(i,1) );   
         valueC=ck_T(i,2)+( temp-ck_T(i,1) )*value_dC_dT;
         value_dK_dT=( ck_T(i+1,3)-ck_T(i,3) )/( ck_T(i+1,1)-ck_T(i,1) );
         valueK=ck_T(i,3)+( temp-ck_T(i,1) )*value_dK_dT; 
         
        end
        
    end
 %% 
    rho=1824; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'G10'))
    %Material='改性C/SiC-一室';
valueC=1100;
valueK=9;
 %% 
 rho=2150; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'G11'))
    %Material='气凝胶-一室';
valueC=800;
valueK=0.35;
 %% 
 rho=650; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'G12'))
    %Material='炭气凝胶-306';
valueC=550; 
 %%
     k_T=[
     25,0.051;    
    1000,0.085;   
    2000,0.1101 
     ];
    [mk,nk]=size(k_T); 
    temp=temp-273.15;

    if( temp<=k_T(1,1) ) 
        valueK=k_T(1,2);
        
    end
    if ( temp>=k_T(mk,1) )
        valueK=k_T(mk,2);
    end
    
    for i=1:mk-1
        if(temp>k_T(i,1)&&temp<k_T(i+1,1) )
         value_dK_dT=( k_T(i+1,2)-k_T(i,2) )/( k_T(i+1,1)-k_T(i,1) );
         valueK     =k_T(i,2)+( temp-k_T(i,1) )*value_dK_dT;          
        end
        
    end   
 %% 
    rho=500; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'G13'))
    %Material='SiO2气凝胶-306';
valueC=800; 
 %%
     k_T=[
     25,0.017;
    300,0.038;
    400,0.039;
    500,0.0404;
    600,0.0424;   
    700,0.045 
     ];
    [mk,nk]=size(k_T); 
    temp=temp-273.15;

    if( temp<=k_T(1,1) ) 
        valueK=k_T(1,2);
        
    end
    if ( temp>=k_T(mk,1) )
        valueK=k_T(mk,2);
    end
    
    for i=1:mk-1
        if(temp>k_T(i,1)&&temp<k_T(i+1,1) )
         value_dK_dT=( k_T(i+1,2)-k_T(i,2) )/( k_T(i+1,1)-k_T(i,1) );
         valueK     =k_T(i,2)+( temp-k_T(i,1) )*value_dK_dT;          
        end
        
    end   
 %% 
    rho=250; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'T1'))
    %Material='ZrO2-陈剑';
   valueC=4500;
   valueK=0.429;
 %% 
    rho=4000; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'T2'))
    %Material='ZrO2-过程所';
    ck_T=[
        25,250.7,0.5659;
       400,309.1,0.5442;
       600,337.3,0.5179;
       800,378.6,0.5714;
      1200,418.8,0.639;
      1400,446.3,0.7981;
      1600,535.1,1.4624        
     ];
    [mck,nck]=size(ck_T);
    %linear interpolation 
  temp=temp-273.15;

    if( temp<=ck_T(1,1) ) 
        valueC=ck_T(1,2);
        valueK=ck_T(1,3);
        
    end
    if ( temp>=ck_T(mck,1) )
        valueC=ck_T(mck,2);
        valueK=ck_T(mck,3);
    end
    
    for i=1:mck-1
        if(temp>ck_T(i,1)&&temp<ck_T(i+1,1) )
         value_dC_dT=( ck_T(i+1,2)-ck_T(i,2) )/( ck_T(i+1,1)-ck_T(i,1) );   
         valueC=ck_T(i,2)+( temp-ck_T(i,1) )*value_dC_dT;
         value_dK_dT=( ck_T(i+1,3)-ck_T(i,3) )/( ck_T(i+1,1)-ck_T(i,1) );
         valueK=ck_T(i,3)+( temp-ck_T(i,1) )*value_dK_dT; 
         
        end
        
    end
 %% 
    rho=5520; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'T3'))
    %Material='ZrO2-金属所';
    ck_T=[
        20,350,0.75;
       200,447,0.868;      
       600,568,0.866;      
      1000,548,0.705;
      1400,645,0.949     
     ];
    [mck,nck]=size(ck_T);
    temp=temp-273.15;
    if( temp<=ck_T(1,1) ) 
        valueC=ck_T(1,2);
        valueK=ck_T(1,3);
    end
    if ( temp>=ck_T(mck,1) )
        valueC=ck_T(mck,2);
        valueK=ck_T(mck,3);
    end
    
    for i=1:mck-1
        if(temp>ck_T(i,1)&&temp<ck_T(i+1,1) )
         value_dC_dT=( ck_T(i+1,2)-ck_T(i,2) )/( ck_T(i+1,1)-ck_T(i,1) );    
         valueC=ck_T(i,2)+( temp-ck_T(i,1) )*value_dC_dT;
         value_dK_dT=( ck_T(i+1,3)-ck_T(i,3) )/( ck_T(i+1,1)-ck_T(i,1) );
         valueK=ck_T(i,3)+( temp-ck_T(i,1) )*value_dK_dT;
        end
    end
        

 %% 
    rho=5400; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(name,'J9'))
    %Material='GH3128';
ck_T=[100,515.1,11.3;
      200,512.8,12.56;
      300,533.4,14.24;
      400,518.7,15.49;
      500,515.2,16.75;
      600,538.9,18.42;
      700,537  ,19.68;
      800,618.2,21.35;
      900,628.1,23.02;
      950,651  ,23.86];
 [mck,nck]=size(ck_T);
    %linear interpolation 
        temp=temp-273.15;
    if( temp<=ck_T(1,1) ) 
        valueC=ck_T(1,2);
        valueK=ck_T(1,3);
    end
    if ( temp>=ck_T(mck,1) )
        valueC=ck_T(mck,2);
        valueK=ck_T(mck,3);
    end
    
    for i=1:mck-1
        if(temp>ck_T(i,1)&&temp<ck_T(i+1,1) )    
         value_dC_dT=( ck_T(i+1,2)-ck_T(i,2) )/( ck_T(i+1,1)-ck_T(i,1) );    
         valueC=ck_T(i,2)+( temp-ck_T(i,1) )*value_dC_dT;
         value_dK_dT=( ck_T(i+1,3)-ck_T(i,3) )/( ck_T(i+1,1)-ck_T(i,1) );
         valueK=ck_T(i,3)+( temp-ck_T(i,1) )*value_dK_dT;
        end
    end
 
 %% 
 rho=8810;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( strcmp(name,'J6'))
    %Material='TA15';
ck_T=[100,545,8.8;
      200,587,10.2;
      300,628,10.9;
      400,670,12.2;
      500,712,13.8;
      600,755,15.1;
      700,838,16.8;
      800,880,18.0;
      900,922,19.7 ];
 [mck,nck]=size(ck_T);
    %linear interpolation 
        temp=temp-273.15;
    if( temp<=ck_T(1,1) ) 
        valueC=ck_T(1,2);
        valueK=ck_T(1,3);
    end
    if ( temp>=ck_T(mck,1) )
        valueC=ck_T(mck,2);
        valueK=ck_T(mck,3);
    end
    
    for i=1:mck-1
        if(temp>ck_T(i,1)&&temp<ck_T(i+1,1) )    
         value_dC_dT=( ck_T(i+1,2)-ck_T(i,2) )/( ck_T(i+1,1)-ck_T(i,1) );    
         valueC=ck_T(i,2)+( temp-ck_T(i,1) )*value_dC_dT;
         value_dK_dT=( ck_T(i+1,3)-ck_T(i,3) )/( ck_T(i+1,1)-ck_T(i,1) );
         valueK=ck_T(i,3)+( temp-ck_T(i,1) )*value_dK_dT;
        end
    end
 
 %% 
 rho=4450;
end

if ( strcmp(name,'J2'))
    %Material='不锈钢';
ck_T=[100,502,16.3;
      200,502,17.6;
      300,502,18.8;
      400,502,20.5;
      500,502,21.8;
      600,502,23.5;
      700,502,24.7;
      800,502,28.5 ];
 [mck,nck]=size(ck_T);
    %linear interpolation 
        temp=temp-273.15;
    if( temp<=ck_T(1,1) ) 
        valueC=ck_T(1,2);
        valueK=ck_T(1,3);
    end
    if ( temp>=ck_T(mck,1) )
        valueC=ck_T(mck,2);
        valueK=ck_T(mck,3);
    end
    
    for i=1:mck-1
        if(temp>ck_T(i,1)&&temp<ck_T(i+1,1) )    
         value_dC_dT=( ck_T(i+1,2)-ck_T(i,2) )/( ck_T(i+1,1)-ck_T(i,1) );    
         valueC=ck_T(i,2)+( temp-ck_T(i,1) )*value_dC_dT;
         value_dK_dT=( ck_T(i+1,3)-ck_T(i,3) )/( ck_T(i+1,1)-ck_T(i,1) );
         valueK=ck_T(i,3)+( temp-ck_T(i,1) )*value_dK_dT;
        end
    end
 
 %% 
 rho=7850;
end




end

