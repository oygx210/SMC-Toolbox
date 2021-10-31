function SMCreplace(Shs, int, Fil)

switch Shs

    case 'on'
        
        if strcmp(get_param([gcb '/s(t)'], 'BlockType'), 'Terminator')
        replace([gcb '/s(t)'], 'built-in/Outport');
        set_param([gcb '/s(t)'],'Port',num2str(2));
        end
        
    case 'off'
        
        if strcmp(get_param([gcb '/s(t)'], 'BlockType'),'Outport')
        replace([gcb '/s(t)'], 'built-in/Terminator');
        end
    
end

switch int
    
    case 'on'
        
        switch Fil
            
            case 'on'
                
                replace([gcb '/r(t)'], 'built-in/Ground');
                
                if strcmp(get_param([gcb '/R(t)'], 'BlockType'), 'Ground')
                replace([gcb '/R(t)'], 'built-in/Inport');
                set_param([gcb '/R(t)'],'Port',num2str(2));
                end
                
            case 'off'
                
                replace([gcb '/R(t)'], 'built-in/Ground');
             
                
                if strcmp(get_param([gcb '/r(t)'], 'BlockType'), 'Ground')
                replace([gcb '/r(t)'], 'built-in/Inport');
                set_param([gcb '/r(t)'],'Port',num2str(2));
                end
        end
        
    case 'off'
        
        replace([gcb '/r(t)'], 'built-in/Ground');
        replace([gcb '/R(t)'], 'built-in/Ground');
end

if strcmp(Fil, 'on') && strcmp(FRT, 'on')
     if strcmp(get_param([gcb '/FR(t)'], 'BlockType'), 'Terminator')
        replace([gcb '/FR(t)'], 'built-in/Outport');
         set_param([gcb '/FR(t)'],'Port',num2str(2));
     end
else
    replace([gcb '/FR(t)'], 'built-in/Terminator')
end

end