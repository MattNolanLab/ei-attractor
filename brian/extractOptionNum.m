function [optVal] = extractOptionNum(optName, str)
    [s e] = regexp(str, ['''' optName ''': \d+\.*\d*'], 'start', 'end');
    if numel(s) ~= 0
        optStr = str(s(1):e(1));
   
        [s e] = regexp(optStr, '\d+\.*\d*', 'start', 'end');
        optStr = optStr(s(1):e(1));
        
        optVal = str2num(optStr);
    else
        optVal = nan;
    end
end
