%multilayer_gere program
clc 
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relax_number=0.8;
 % gas  input 
    %fluent_position(mm),Ma,total temperture(K),total pressure(KPa),workingtime(s)     
 LMTPt=importdata('gas.txt');
 [mm,nn]=size(LMTPt);
 %reading material and thickness from input
material_and_thickness=importdata('input.txt');
 %m=numbers of material
[m0,n0]=size(material_and_thickness.data);

 s=0;  
 % prepare to write 
 fid=fopen('T-t.txt','wt'); 

sum_num =0;
for j=1:m0
    n(j)       =material_and_thickness.data(j,2);
    sum_num    =sum_num+n(j);
end
 %initiate temperature
 Ta=300;%(degree_K)
T=Ta*ones(sum_num+1,1)      ;
  %outlayer exchange temperature coeffients
 alpha=0.0;
  %outlayer  temperature
 T_out=Ta;%(degree_K)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 for ii=1:mm
  L=LMTPt(ii,1)/1000;
 Ma=LMTPt(ii,2); 
 T0=LMTPt(ii,3);  
 pt=LMTPt(ii,4)*1000;
 t1=LMTPt(ii,5);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %gas output
 gamma=1.4;R=287;
 sigema=5.67*10^-8;
  %static temperature
 %Ma= (ps/pt)^(-(gamma-1)/gamma)-1;
 %Ma=2/(gamma-1)*Ma;
 %Ma=sqrt(Ma);
 Te=T0/( 1+(gamma-1)/2*Ma*Ma );
 %T0=Te*( 1+(gamma-1)/2*Ma*Ma );
 %pt=ps*( 1+(gamma-1)/2*Ma*Ma )^(gamma/(gamma-1));
 ps=pt*( 1+(gamma-1)/2*Ma*Ma )^(-gamma/(gamma-1));
 % gas velosity 
  Ue=Ma*sqrt(gamma*R*Te); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %recover temperature coeffients 
  [mu,~,~,Pr]=hotgasprop(Te);
  L=1.5;
  Re=ps/R/Te*Ue*L/mu;
  if (Re<=3E3)
  r=Pr^(1/2);
  end
  if (Re>3E3)
  r=Pr^(1/3);
  end
  
  [he,~]=Enthalpy(Te,'dry_air');  
  [h0,~]=Enthalpy(T0,'dry_air');  
    
  % ¾øÈÈ Enthalpy
  haw=he+r*(h0-he); 
  Taw=Enthalpy_Inverse(haw,'dry_air');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% start temperature loop
 %time step and total time and number of time step  
deta_t=0.1;t=t1;n_k=t/deta_t;  
 %wilson seta
 seta=3/3;

 %time marching
  for k=1:n_k
 % initiate global matrix 
C=zeros(sum_num+1,sum_num+1);
K=zeros(sum_num+1,sum_num+1);
P=zeros(sum_num+1,1)        ;

Ce=zeros(sum_num,2,2);
Ke=zeros(sum_num,2,2);

sum_num=0;
%loop  material
for j=1:m0
   
   material(j)=material_and_thickness.textdata(j);  
   H(j)       =material_and_thickness.data(j,1)/1000; 
   %n(j)       =material_and_thickness.data(j,2);
   sum_num    =sum_num+n(j);
   eth(j)     =H(j)/n(j);
  
    
   %loop elements
  for i=sum_num-n(j)+1:sum_num 
    %elemental Ce Ke     
    [~ , tempC1 tempK1,~,~]=Multilinear_rhoCK(T(i),  material(j));
    [rho tempC2 tempK2,~,~]=Multilinear_rhoCK(T(i+1),material(j));
    
    Ce(i,1,1)=1/4 *tempC1+1/12*tempC2;
    Ce(i,1,2)=1/12*tempC1+1/12*tempC2;
    Ce(i,2,1)=1/12*tempC1+1/12*tempC2;
    Ce(i,2,2)=1/12*tempC1+1/4 *tempC2;
    
    Ce(i,:,:)=eth(j)*rho* Ce(i,:,:);    
    
    Ke(i,1,1)= 1/2*tempK1+1/2*tempK2;
    Ke(i,1,2)=-1/2*tempK1-1/2*tempK2;
    Ke(i,2,1)=-1/2*tempK1-1/2*tempK2;
    Ke(i,2,2)= 1/2*tempK1+1/2*tempK2;
    
    Ke(i,:,:)=1/eth(j)*   Ke(i,:,:);
    
    %global matrix C K¡¡P
    for i1=i:i+1
      for j1=i:i+1
        C(i1,j1)=Ce(i,i1-i+1,j1-i+1)+C(i1,j1);
        K(i1,j1)=Ke(i,i1-i+1,j1-i+1)+K(i1,j1);
      end
    end   
    
  end   

    
end
 %if there exist outerlayer convection just plus it    
    K(sum_num+1,sum_num+1)=K(sum_num+1,sum_num+1)+alpha;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
  %wall Enthalpy
  [h_wall,~]=Enthalpy(T(1),'dry_air');
  
 %ref Enthalpy
  h_ref=0.5* (h_wall+he)+0.22* (haw-he);  
  T_ref=Enthalpy_Inverse(h_ref,'dry_air');
  
   %gas densetiy p/R/T_ref
  rho_ref=ps/R/T_ref;
  %%choosing different model of St
  if (Re>3E3)  
  [St,~]=turbulence(ps,R,Ue,L,Pr,T_ref,'model_tu1');
  end  
  if (Re<=3E3)
  [St,~]=laminar   (ps,R,Ue,L,Pr,T_ref,'model_la1');
  end      
  %gas convection flux 
  qc=rho_ref*Ue*St*(haw-h_wall);
  
  %gas radiation  flux  
   %gas emission coeffients
  e_em_g=0.0;
   %solid emission coeffients
  e_em_s=0.0;
  hr_min=4*sigema*e_em_g*T(1)^3;
  hr_max=4*sigema*e_em_g*Te^3;
  qr=sigema*(e_em_g*Te^4-e_em_s*T(1)^4);    
  % total gas flux
  q=qr+qc;  
  hc=q/(Taw-T(1));
  %force vector
     P(1)        =q;
     P(sum_num+1)=alpha*T_out;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  %N-L iteration  
   %T_temp   
     eps=1.0*10^-2;
     deta_T_temp=Ta*ones(sum_num+1,1);
     %initial temperature of k+1
     T_temp=T;         
     
     while( norm(deta_T_temp,2)>eps )    
     % initial global matrix 
   sum_num1=0; 
     C_temp=zeros(sum_num+1,sum_num+1);
     K_temp=zeros(sum_num+1,sum_num+1);
     P_temp=zeros(sum_num+1,1); 
     
     dC_dT_temp=zeros(sum_num+1,sum_num+1);
     dK_dT_temp=zeros(sum_num+1,sum_num+1);
     
     Ce_temp=zeros(sum_num,2,2);
     Ke_temp=zeros(sum_num,2,2);
     
     dCe_dT_temp=zeros(sum_num,2,2);
     dKe_dT_temp=zeros(sum_num,2,2);
     
     
for j=1:m0  
   %material(j)=material_and_thickness.textdata(j);  
   %H(j)       =material_and_thickness.data(j,1)/1000; 
   %n(j)       =material_and_thickness.data(j,2);
   sum_num1    =sum_num1+n(j);
   %eth(j)     =H(j)/n(j);
    
   %loop elements
  for i=sum_num1-n(j)+1:sum_num1
      
    [~ , tempC1 tempK1 dC1_dT1_temp dK1_dT1_temp]=Multilinear_rhoCK(T_temp(i),  material(j));
    [rho tempC2 tempK2 dC2_dT2_temp dK2_dT2_temp]=Multilinear_rhoCK(T_temp(i+1),material(j));
    
    Ce_temp(i,1,1)=1/4 *tempC1+1/12*tempC2;
    Ce_temp(i,1,2)=1/12*tempC1+1/12*tempC2;
    Ce_temp(i,2,1)=1/12*tempC1+1/12*tempC2;
    Ce_temp(i,2,2)=1/12*tempC1+1/4 *tempC2;
    
    Ce_temp(i,:,:)=eth(j)*rho* Ce_temp(i,:,:); 
    
    dCe_dT_temp(i,1,1)=1/4 *dC1_dT1_temp*T_temp(i)+1/12*dC2_dT2_temp*T_temp(i+1);
    dCe_dT_temp(i,1,2)=1/12*dC1_dT1_temp*T_temp(i)+1/12*dC2_dT2_temp*T_temp(i+1);
    dCe_dT_temp(i,2,1)=1/12*dC1_dT1_temp*T_temp(i)+1/12*dC2_dT2_temp*T_temp(i+1);
    dCe_dT_temp(i,2,2)=1/12*dC1_dT1_temp*T_temp(i)+1/4 *dC2_dT2_temp*T_temp(i+1);

    Ke_temp(i,1,1)= 1/2*tempK1+1/2*tempK2;
    Ke_temp(i,1,2)=-1/2*tempK1-1/2*tempK2;
    Ke_temp(i,2,1)=-1/2*tempK1-1/2*tempK2;
    Ke_temp(i,2,2)= 1/2*tempK1+1/2*tempK2;
    
    Ke_temp(i,:,:)=1/eth(j)*   Ke_temp(i,:,:);
    
    dKe_dT_temp(i,1,1)=1/4 *dK1_dT1_temp*T_temp(i)+1/12*dK2_dT2_temp*T_temp(i+1);
    dKe_dT_temp(i,1,2)=1/12*dK1_dT1_temp*T_temp(i)+1/12*dK2_dT2_temp*T_temp(i+1);
    dKe_dT_temp(i,2,1)=1/12*dK1_dT1_temp*T_temp(i)+1/12*dK2_dT2_temp*T_temp(i+1);
    dKe_dT_temp(i,2,2)=1/12*dK1_dT1_temp*T_temp(i)+1/4 *dK2_dT2_temp*T_temp(i+1);    
    
    
     %global matrix C K¡¡P
     for i1=i:i+1
       for j1=i:i+1
        C_temp(i1,j1)=Ce_temp(i,i1-i+1,j1-i+1)+C_temp(i1,j1);
        K_temp(i1,j1)=Ke_temp(i,i1-i+1,j1-i+1)+K_temp(i1,j1);
        dC_dT_temp(i1,j1)=dCe_dT_temp(i,i1-i+1,j1-i+1)+dC_dT_temp(i1,j1);
        dK_dT_temp(i1,j1)=dKe_dT_temp(i,i1-i+1,j1-i+1)+dK_dT_temp(i1,j1);
        
       end
     end 
    
  end  

    
    
end
 %if there exist outerlayer convection just plus it    
    K_temp(sum_num1+1,sum_num1+1)=K_temp(sum_num1+1,sum_num1+1)+alpha; 
    
  [h_wall_temp,dh_dT1_temp]=Enthalpy(T_temp(1),'dry_air');    
     %ref Enthalpy
  h_ref_temp=0.5* (h_wall_temp+he)+0.22* (haw-he);
  T_ref_temp=Enthalpy_Inverse(h_ref_temp,'dry_air');
  
  %gas densetiy p/R/T_ref
  rho_ref_temp=ps/R/T_ref_temp;
    
  if (Re>3E3)  
  [St_temp,dSt_T_dT1_temp]=turbulence(ps,R,Ue,L,Pr,T_ref_temp,'model_tu1');
  end  
  if (Re<=3E3)
  [St_temp,dSt_T_dT1_temp]=laminar   (ps,R,Ue,L,Pr,T_ref_temp,'model_la1');
  end   
  
  %gas convection flux
  qc_temp=rho_ref_temp*Ue*St_temp*(haw-h_wall_temp);
  %gas radiation  flux 
  qr_temp=sigema*(e_em_g*Te^4-e_em_s*T_temp(1)^4);      
  % total gas flux
  q_temp=qr_temp+qc_temp;  
  %force vector
     P_temp(1)  =q_temp;
     P_temp(sum_num1+1)=alpha*T_out;     
   
   %tan Matrix 
    %dP/dT
    dPdT_temp=zeros(sum_num1+1,sum_num1+1);      

     %dq/dT1
    dPdT_temp(1,1)=ps/R*Ue*( haw-h_wall_temp )*dSt_T_dT1_temp;
    dPdT_temp(1,1)=dPdT_temp(1,1)-ps/R/T_ref_temp*Ue*St_temp* dh_dT1_temp;
    dPdT_temp(1,1)=dPdT_temp(1,1)-4*e_em_s*sigema*T_temp(1)^3;
   %dF
    dF_temp=1/deta_t*dC_dT_temp+seta*dK_dT_temp+1/deta_t*C_temp+seta*K_temp-seta*dPdT_temp;
 
    F=(1.0/deta_t*C_temp+seta*K_temp)*T_temp-( 1/deta_t*C-(1.0-seta)*K )*T-(1-seta)*P-seta*P_temp;
    deta_T_temp=-dF_temp\F;
    T_temp=T_temp+relax_number*deta_T_temp;    
   
     end     
     T=T_temp;  
  
    %time out
    fprintf('Time ');
    fprintf('%f  ',s+(k)*deta_t);
    fprintf('s \n');
    
      fprintf(fid,'%f   ', s+(k)*deta_t); 
      for out=1:sum_num+1
      fprintf(fid,'%f   ', T(out)-273.15); 
      end
    
      fprintf(fid,'%f   \n', q/10^6);    
  end
   s=s+t; 
 end  
 fclose(fid);

 
 
 
 
 

 
 
 

 
  
  
  
    
    
    
   
     

 
 
 





