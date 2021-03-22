function [St,dSt_T_dT1_temp]=laminar(ps,R,Ue,L,Pr,T_ref,model_flag)

  syms T  
  %gas 粘性 萨瑟兰公式 见安德森书
  v=1.7894*10^-5*(T/288)^1.5*398/(T+110);
  %Re
  rho_ref=ps/R/T;
  Re   =rho_ref*Ue*L/v; 
  
if ( strcmp(model_flag,'model_la1'))  
  lg_Re=log10(Re);
  cf=0.37/ ((lg_Re)^2.584 ); 
  St=cf/2/( Pr^(2/3) );
end

if ( strcmp(model_flag,'model_la2'))
  cf=0.332/ ((Re)^0.5 ); 
  St=cf/( Pr^(2/3) );
end

  dSt_T_dT=diff(St/T,'T');
  dSt_T_dT1_temp=0.5*subs(dSt_T_dT,'T',T_ref);
  St=subs(St,'T',T_ref); 

















end

  