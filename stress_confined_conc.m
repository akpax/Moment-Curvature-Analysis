function [fcc] = stress_confined_conc(ecc, epcc, fpcc, lambda_conc, eccu, rcc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Inputs%%%%%%%%%%%%%%%%%%%%%%
%ecc: strain of confined cocnrete
%epcc: confined concrete strain
%fpcc: compressive strength of confined concrete
%lambda_conc: regulariation coefficient
%ecu: regularied crushing strain of unconfined concrete
%rcc: Non-linearity function parameter
%%%%%%%%%%%Output %%%%%%%%%%%%%%%%%%%%%%%%
%unconfined concrete stress at ecc

if ecc>=epcc
    fcc=-(ecc/epcc)/(rcc-1+(ecc/epcc)^rcc)*rcc*fpcc;
elseif ecc > eccu && ecc< epcc
    fcc=-(1+1/lambda_conc*(ecc/epcc-1))/(rcc-1+(1+1/lambda_conc*(ecc/epcc-1))^rcc)*rcc*fpcc;
else
    fcc = 0;
end
