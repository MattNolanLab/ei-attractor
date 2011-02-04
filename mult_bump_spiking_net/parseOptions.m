function [optOut] = parseOptions(options)
    % Parse options written into matlab file by Brian
    % We don't expect any errors in the input string
    
    optOut.sheet_size = extractOptionNum('sheet_size', options);
    optOut.alpha = extractOptionNum('alpha', options);
    optOut.a = extractOptionNum('a', options);
    optOut.sim_dt = extractOptionNum('sim_dt', options);
    optOut.lambda_net = extractOptionNum('lambda_net', options);
    optOut.job_num = extractOptionNum('job_num', options);
    optOut.l = extractOptionNum('l', options);
    optOut.connMult = extractOptionNum('connMult', options);
    optOut.time = extractOptionNum('time', options);
    optOut.update_interval = extractOptionNum('update_interval', options);
    optOut.input = extractOptionNum('input', options);
    optOut.taum =  extractOptionNum('taum', options);
    optOut.taui = extractOptionNum('taui', options);
    optOut.threshold = extractOptionNum('threshold', options);
end
