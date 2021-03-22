function [h,dh_dT]=Enthalpy(T,material_flag)

if ( strcmp(material_flag,'dry_air'))
a=[297.84357;
  998.37261; 
   27.75219; 
 -189.63421;
  264.98218; 
 -192.38503; 
   86.80328 ;
  -25.05508 ;
    4.50058 ;
   -0.45827 ;
    0.0202];
b=[-301.304;
   1043.76;
   -206.916;
    388.667;
   -220.448;
     29.94557;
     26.75992;
    -15.8513;
      3.83286;
     -0.45165;
      0.021144];
  
  sum=0;sum1=0;
for i=2:11
    sum=sum+b(i)*(T/1000)^(i-1);sum1=sum1+(i-1)*b(i)*(T/1000)^(i-2);
end
sum=sum+b(1);
h=sum*1000;
dh_dT=sum1;
end

if ( strcmp(material_flag,'full_quan') )
h=T;
end



end

    








    





