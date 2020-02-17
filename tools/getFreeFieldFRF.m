function G = getFreeFieldFRF(ES_loc,Mic_loc,k)

M = size(Mic_loc,1);
L = size(ES_loc,1);
G = zeros(M,L);

% Loop over reconstruction virtual microphones
for m = 1 : M    
    
    % location of the current reconstructio mic
    rm = Mic_loc(m,:);    
    
    % Loop over each equivalent source    
     for l = 1 : L
        
        % Location of the current equivalent source
        r0 = ES_loc(l,:);
        R = norm(r0 - rm);
        G(m,l) = exp(-1i * k * R)/(4*pi*R);         
        
    end
end
