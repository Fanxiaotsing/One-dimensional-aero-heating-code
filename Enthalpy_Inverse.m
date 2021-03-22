function [T]=Enthalpy_Inverse(h,material_flag)
eps=100;
if ( strcmp(material_flag,'dry_air'))

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
  
syms x
y=0;
for i=1:11;
    y=y+b(i)*x^(i-1);
end
diff_y_x=diff(y,'x');
%二分法求得满足h的T
%Tmax 4000K
xmax=4000/1000;
%Tmin 100K
xmin=100/1000;
hmax=1000*subs(y,'x',xmax);
hmin=1000*subs(y,'x',xmin);

x_ave=0.5*xmax+0.5*xmin;
h_ave=1000*subs(y,'x',x_ave);
k=1000*subs(diff_y_x,'x',x_ave);
x1=x_ave+(h-h_ave)/k;
h1=h_ave;
while(abs(h-h1)>eps)

h1=1000*subs(y,'x',x1);
%tangent
k1=1000*subs(diff_y_x,'x',x1);
x1=x1+(h-h1)/k1;
end

T=1000*x1;  


    


end

end

    








    





