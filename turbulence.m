function  [St,dSt_T_dT1_temp]=turbulence(ps,R,Ue,L,Pr,T_ref,model_flag)


 syms T
  %gas 粘性 萨瑟兰公式 见安德森书
  v=1.7894*10^-5*(T/288)^1.5*398/(T+110);
  %Re
  rho_ref=ps/R/T;
  Re   =rho_ref*Ue*L/v; 
  
if ( strcmp(model_flag,'model_tu1'))  
  cf=0.0287/ ((Re)^0.2 ); 
  St=cf/( Pr^(2/5) ); 
end

if ( strcmp(model_flag,'model_tu2'))

end

  dSt_T_dT=diff(St/T,'T');
  dSt_T_dT1_temp=0.5*subs(dSt_T_dT,'T',T_ref);
  St=subs(St,'T',T_ref);















end

