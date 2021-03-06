function [fc]= stress_unconfined_conc(ec, e_pc, ecu,  f_pc, ne, lambda_conc, rc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Inputs%%%%%%%%
%ec: concrete strain 
%e_pc: strain at the unconfined compressive strength 
%ecu: crushing strain of unconfined concrete
%ne: nonlinearity coefficent of pre-peak stress-strain curve
%lambda_conc: regulariation coefficient eq. 10b
%rc: Post-Peak non-linear function coefficient

%%%%%%%Outputs%%%%%%%%
%fc: unconfined cocnrete stress at inputed strain

if  e_pc < ec
    fc = -f_pc*(1-(1-ec/e_pc)^ne);
elseif e_pc > ec && ec >ecu
    fc = -(1+1/lambda_conc*(ec/e_pc-1))/(rc-1+(1+1/lambda_conc*(ec/e_pc-1))^rc)*rc*f_pc;
else
    fc=0;
end