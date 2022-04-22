function[ftc]=tens_conc(Ec,fcr,ec,ecr)
% Obtains tensile stress of concrete model
%% Input
% Ec is the Young's modulus of concrete
% fcr is the cracking strength of concrete
% ec is the strain of concrete
% ecr is the cracking strain of concrete
%% Output
% ftc is the stress at tensile strain ec
%%
if ec>ecr
    ftc=0.7*fcr/(1+sqrt(500*ec));
else, ftc=ec*Ec;
end
end