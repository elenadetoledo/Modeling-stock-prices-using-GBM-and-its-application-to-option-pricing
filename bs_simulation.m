function [V, Vex, ster, S0, actual_ST, error, vol, cpt, end_date] = OPT_bs_simulation_v2(start_date_test, ndaystraining, t0, K, T_days, r, M)
    
    V = [];
    Vex = [];
    ster = [];
    S0 = [];
    actual_ST = [];
    error = [];
    vol = [];
    cpt = [];
    end_date = [];

    try

        % Start timing the execution
        ndayspredicting = T_days;
        T_days = T_days / 365;
        warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    
        filename_5years = ''; %Include here the path to the file containing daily data
        dataTable_total = readtable(filename_5years);
        closeData_total = dataTable_total.Close;
        datesData_total = dataTable_total.Date;
    
        test_date = datetime(start_date_test, 'InputFormat', 'dd-MMM-yyyy', 'Locale', 'en_US');
        idx_startTest = find(datesData_total == test_date);
        traindates = datesData_total(idx_startTest - ndaystraining : idx_startTest - 1);
        St = closeData_total(idx_startTest - ndaystraining : idx_startTest - 1);
        ratio = St(2:end) ./ St(1:end-1); %because we can fit distribution with Δt = 1
        logreturns = log(ratio);
        actual_data = closeData_total(idx_startTest : idx_startTest + ndayspredicting - 1); 
        preddates = datesData_total(idx_startTest : idx_startTest + ndayspredicting - 1);
        end_date = datestr(preddates(end), 'dd-mmm-yyyy');
        dates_logreturns = traindates(2:end);
        stock_return = (St(2:end) - St(1:end-1)) ./ St(1:end-1);
        vol = sqrt(1 / (length(stock_return)-1) * sum((stock_return - mean(stock_return)).^2));
    
        S0 = actual_data(1);
    
        tic;
        % M - number of trajectories
        % n - number of time steps
        n = T_days- t0;
    
        % Initialize variables
        S_T = zeros(M, 1); % Matrix to store stock prices
    
        % Vectorized simulation of the value of the stock at maturity (S_t)
        S_T = S0 * exp((r - 0.5 * vol^2) * n + vol * sqrt(n) * randn(M, 1));
    
        S_T(S_T <= 0) = 0; 
    
        % Compute the payoff for each simulation at maturity (T)
        payoff = max(S_T - K, 0); % Payoff of a call option
    
        % Calculate the option value (V) as the present value of the average payoff
        V = exp(-r * (T_days - t0)) * mean(payoff);
    
        % Calculate the standard error (ster) of the simulation
        ster = std(payoff) * exp(-r * (T_days - t0)) / sqrt(M);
    
        % % Compare with exact solution
        [Vex, ~] = blsprice(S0, K, r, T_days, vol);
        %fprintf('blsprice %f\n', Vex);
        error = V - Vex; % Difference between Monte Carlo and exact solution
        actual_ST = actual_data(end);
        cpt = toc; % End timing the execution
    
        % Display results
        % fprintf('European call with K=%g, T=%g days, vol=%g, r=%g, N=%d, computation time=%g seconds\n', K, T_days, vol, r, M, cpt);
        % fprintf('Spot price: S0=%f\n', S0);
        % fprintf('V(t=%g,S=%g)= %g (exact) %g +/- %g (Monte Carlo). Error= %g\n', t0, S0, Vex, V, ster, abs(error));
        % fprintf('Actual price of the stock at maturity: %g\n\n', actual_ST);
    catch ME
        %warning('An error occurred: %s', ME.message);
        % Return empty values
        V = [];
        Vex = [];
        ster = [];
        S0 = [];
        actual_ST = [];
        error = [];
        vol = [];
        cpt = [];
        end_date = [];
        
    end
end

