clear all
format long

onskad_x = 7.5;  
onskad_y = 15;     
steg = 0.5;
tolerans = 0.01;   

bast_F = 0;
bast_Fi = 0;
bast_fel = 100;

while true
    F_range = 0.4:steg:1.5;   
    Fi_range = 1.3:steg:1.5;  
    
    for F = F_range
        for Fi = Fi_range
            [x, max_hojd] = simulera_raket2(F, Fi);
            fel = sqrt((x - onskad_x)^2 + (max_hojd - onskad_y)^2);
            if fel < bast_fel
                bast_fel = fel;
                bast_F = F;
                bast_Fi = Fi;
            end
            if bast_fel < tolerans
                break;
            end
        if bast_fel < tolerans
            break;
        end
        end
    if bast_fel < tolerans
        break;
    end
    end
    if bast_fel < tolerans
        break;
    end
    steg = steg/2;
end

disp(['Kraften ska vara: ', num2str(bast_F)]);
disp(['Vinkeln ska vara: ', num2str(rad2deg(bast_Fi))]);
disp(['X-position: ', num2str(x)]);
disp(['Y-position: ', num2str(max_hojd)]);

 