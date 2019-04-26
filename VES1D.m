function con_a = VES1D(con, thick, ab2)
%   Calculate forward VES response with half the current electrode 
%   spacing for a 1D layered earth.
%
%   Parameters
%   ----------
%   con : Electrical conductivity of n layers (S/m), 
%         last layer n is assumed to be infinite
%
%   thick : Thickness of n-1 layers (m), last layer n 
%   is assumed to be infinite and does not require a thickness
%
%   ab2 : Half the current (AB/2) electrode spacing (m)
%
%
%   Returns
%   -------
%   app_con : Apparent half-space electrical conductivity (S/m)
%
%
%   References
%   ----------
%   Ekinci, Y. L., Demirci, A., 2008. A Damped Least-Squares Inversion Program 
%   for the Interpretation of Schlumberger Sounding Curves, Journal of Applied Sciences, 8, 4070-4078.
%
%   Koefoed, O., 1970. A fast method for determining the layer distribution from the raised
%   kernel function in geoelectrical soundings, Geophysical Prospection, 18, 564-570.
%
%   Nyman, D. C., Landisman, M., 1977. VES Dipole-dipole filter coefficients,
%   Geophysics, 42(5), 1037-1044.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Profile characteristics
    res = 1./con;
    lays = length(res);
    
    % Constants
    LOG = log(10);
    COUNTER = 1 + (2 * 13 - 2);
    UP = exp(0.5 * LOG / 4.438);
    
    % Filter integral variable
    up = ab2 * exp(-10 * LOG / 4.438);

    % Loop for digital filter
    for ii = 1:COUNTER
        
        % Resistivity of last layer equals T
        ti1 = res(lays);
        
            % Loop for recursive formula of Gosh
            lay = lays;
            while lay > 1
                lay = lay - 1;
                tan_h = tanh(thick(lay) / up);
                ti1 = (ti1 + res(lay) * tan_h) / ...
                    (1 + ti1 * tan_h / res(lay));
            end
        
        % Set layer to previous
        ti(ii) = ti1;
        
        % Update digital filter
        up = up * UP;
    end

    % Set to top layer
    ii = 1;
    
    % Apply digital filter
    res_a = (105 * ti(ii) - 262 * ti(ii + 2) + ...
        416 * ti(ii + 4) - 746 * ti(ii + 6) + 1605 * ti(ii + 8) ...
        - 4390 * ti(ii + 10) + 13396 * ti(ii + 12) - 27841 * ...
        ti(ii + 14) + 16448 * ti(ii + 16) + 8183 * ti(ii + 18) + ...
        2525 * ti(ii + 20) + 336 * ti(ii + 22) + 225 * ti(ii + 24)) / 1e4;
    
    % Res to con
    con_a = 1/res_a;
    
return