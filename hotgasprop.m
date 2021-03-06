function [mu,cp,k,pr]=hotgasprop(T)
    mu=(1.0E-6)*1.458*(T^1.5)/(T+110.4);
    cp=10*(0.021144)*((0.001*T)^9)+9*(-0.45165)*((0.001*T)^8)...
         +8*(3.83286)*((0.001*T)^7)+7*(-15.8513)*((0.001*T)^6)...
         +6*(26.75992)*((0.001*T)^5)+5*(29.9455575)*((0.001*T)^4)...
         +4*(-220.448)*((0.001*T)^3)+3*(388.667)*((0.001*T)^2)...
         +2*(-206.916)*((0.001*T)^1)+1043.76;   
    k=(1.0E-3)*2.648151*(T^1.5)/( T+( 245.4*(10^(-12/T)) ) );
   pr=mu*cp/k;