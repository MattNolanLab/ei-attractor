function [sheet_size alpha] = parseOptions(options)
    % Parse options written into matlab file by Brian
    % We don't expect any errors in the input string
    
    sheet_size = extractOptionNum('sheet_size', options);
    alpha = extractOptionNum('alpha', options);
   
end