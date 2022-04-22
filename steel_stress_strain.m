function[fs]=steel_stress_strain(Es,fy,P,esh,esu,es,fsu)
% Obtains steel stress using Mander's reinforcing steel model.
%% Input
% Es is the Young's modulus of steel
% fy is the yield stress of steel
% P is the power of the function for the strain hardening region
% esh is the strain at onset of strain hardening
% esu is the uniform strain
% es is the strain at the fiber
% fsu is the ultimate yield strength
%% Output
% fs is the stress at strain es
%%
if es>0
    if es<=fy/Es
        fs=es*Es;
    else
        if es<=esh
            fs=fy;
        else
            if es<esu
                fs=fsu-(fsu-fy)*(((esu-es)/(esu-esh)))^P;
            else, fs=0;
            end
        end
    end
else
    if es > -fy/Es
        fs=es*Es;
    else, fs=-fy;
    end
end
end