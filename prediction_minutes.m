function [MAPE_onesimulation, MAPE_manysimulations, MAPE_mean, MAPE_median, MAPE_mode] = prediction_minutes(filename, start_min, nminstraining, nminspredicting, M, callibration_option,  Z, z)
    MAPE_onesimulation = [];
    MAPE_manysimulations = [];
    MAPE_mean = [];
    MAPE_median = [];
    MAPE_mode = [];

    %Prediction strategy: initial step S(i-1). (t-T0)=1, muo & sigmao are constant
    % For the callibration option, choose 1 for using fitdist, choose 2 to use
    %the methods explained in Farida (sample mean, sample standard deviation)
    %Show images flag. 
    images = 0 ; % Change to 1 to display images
    final_image=1; % Change to 1 if you want to display the final image
    display = 0; % Change to 1 if you want printfs

    if Z == 0
        Z = randn(M,1);
    end
    if z == 0
        z = randn;
    end
    if display == 1
        fprintf('\nMAKING PREDICTIONS USING S(i-1) AS BASIS AND CONSTANT MUO AND SIGMA\n');
    end

    try 
        
        dataTable_total = readtable(filename);
        closeData_total = dataTable_total.Cierre;
        datesData_total = dataTable_total.Date;
        

        if display == 1
        fprintf('You entered the following startdate: %s\n', start_min);
        end
        
        test_date = datetime(start_min, 'InputFormat', 'dd-MMM-yyyy HH:mm', 'Locale', 'en_US');
        idx_startTest = find(datesData_total == test_date);

        if display == 1
            fprintf('Position in the data: %d\n', idx_startTest);
        end
       
        traindates = datesData_total(idx_startTest - nminstraining : idx_startTest - 1);
        St = closeData_total(idx_startTest - nminstraining : idx_startTest - 1);
        ratio=St(2:end)./St(1:end-1); % Computing minutely logarithmic returns
        logreturns = log(ratio);
        
        actual_data = closeData_total(idx_startTest : idx_startTest + nminspredicting - 1); 
        preddates = datesData_total(idx_startTest : idx_startTest + nminspredicting - 1);

        if display == 1
            fprintf('Train data: %s - %s (%d data points) , Test data: %s - %s (%d data points)\n', datestr(traindates(1), 'dd-MMM-yyyy HH:mm'), datestr(traindates(end), 'dd-MMM-yyyy HH:mm'), length(traindates), datestr(preddates(1), 'dd-MMM-yyyy HH:mm'), datestr(preddates(end), 'dd-MMM-yyyy HH:mm'), length(preddates));
            fprintf('TRAINING: Hypothesis testing\n');
            fprintf('%s - %s: Sample size Stocks: %d, Max: %f, Min: %f, Mean:%f, Std dev:%f\n', datestr(traindates(1), 'dd-MMM-yyyy HH:mm'), datestr(traindates(end), 'dd-MMM-yyyy HH:mm'),length(St), max(St), min(St), mean(St), std(St));
        end

        dates_logreturns = traindates(2:end);
        if images == 1
            figure();
            plot(traindates, St);
            xlim([traindates(1), traindates(end)]);
            title(sprintf('AAPL Price %s - %s', datestr(traindates(1), 'dd-MMM-yyyy HH:mm'), datestr(traindates(end), 'dd-MMM-yyyy HH:mm')));
            ylabel('Closing Price ($)');
            xlabel('Date');
        end
        
        % INDEPENDENCE
        %1 - Visually assessing independence of log returns
        if images == 1
            figure();
            plot(dates_logreturns, logreturns);
            xlim([dates_logreturns(1), dates_logreturns(end)]);
            title(sprintf('Log Returns %s - %s', datestr(traindates(1), 'dd-MMM-yyyy HH:mm'), datestr(traindates(end), 'dd-MMM-yyyy HH:mm')));
            ylabel('Log return');
            xlabel('Date');
        
            figure;
            autocorr(logreturns, ceil(length(logreturns)/4));
            title(sprintf('ACF AAPL Log Returns %s - %s', datestr(traindates(1), 'dd-MMM-yyyy HH:mm'), datestr(traindates(end), 'dd-MMM-yyyy HH:mm')));
            xlim([1, ceil(10*log10(length(logreturns)))]);
            ylabel('ACF');
        end
        
        %Analytical tests
        [h_jb,p_jb] = jbtest(logreturns);
        [h_lbq, ~, ~, ~] = lbqtest(logreturns);
        if display == 1
            fprintf('%s - %s: Sample size logreturns: %d, Max: %f, Min: %f, Mean:%f, Std dev:%f, Skewness:%f, Jarque-Bera Test:%d, p-value:%f, Ljung-Box Q-test:%d\n', datestr(traindates(1), 'dd-MMM-yyyy HH:mm'), datestr(traindates(end), 'dd-MMM-yyyy HH:mm'), length(logreturns), max(logreturns), min(logreturns), mean(logreturns), std(logreturns), skewness(logreturns), h_jb, p_jb, h_lbq);
        end
        if display == 1
            if h_jb == 0
                fprintf('** LOGRETURNS ARE NORMAL\n');
            end
            
            if h_lbq == 0
                fprintf('** LOGRETURNS ARE UNCORRELATED\n');
            end
        end
        
        %Plotting the histogram
        if images == 1
            ceil(sqrt(length(logreturns)));
            figure;
            histogram(logreturns, ceil(sqrt(length(logreturns))));
            title(sprintf('Histogram of Log Returns %s - %s', datestr(traindates(1), 'dd-MMM-yyyy HH:mm'), datestr(traindates(end), 'dd-MMM-yyyy HH:mm')));
            xlabel('Value');
            ylabel('Frequency');
        end
        
        %Q-Q plot
        if images == 1
            figure;
            qqplot(logreturns);
            title(sprintf('Q-Q Plot of Data vs. Standard Normal %s - %s', datestr(traindates(1), 'dd-MMM-yyyy HH:mm'), datestr(traindates(end), 'dd-MMM-yyyy HH:mm')));
            xlabel('Theoretical Quantiles');
            ylabel('Sample Quantiles');
        end
        
        %Visualization of the log returns in a histogram
        if images == 1
            figure;
            histfit(logreturns);
        end
        
        if callibration_option == 1
            pd=fitdist(logreturns,'Normal');
            sigmao = pd.sigma; % Volatility of the stock
            muo = pd.mu + 1/2*(sigmao)^2; % Drift of the stock
            if display == 1
                fprintf('TRAINING: Model calibration using fitdist (MLE)\n');
            end
        else %callibration_option == 2
            stock_return = (St(2:end) - St(1:end-1))./St(1:end-1);
            muo = 1/length(stock_return) * sum(stock_return);
            sigmao = sqrt(1/(length(stock_return)-1)*sum((stock_return - mean(stock_return)).^2));
            if display == 1
                fprintf('TRAINING: Model calibration using sample mean and sample std dev of arithmetic returns\n');
            end
        end
        if display == 1
            fprintf('Calibrated parameters: mu=%f sigma=%f\n', muo, sigmao);
        end
        
        % PARAMETERS OF THE SIMULATION
        T = length(actual_data)-1;             % (*) Time at which we want to assess the price (in minutes)
        N = length(actual_data)-1;           %  (*) Number of time steps
        if display == 1
            fprintf('TESTING: Model creation\n');
            fprintf('Parameters of testing: T=%d, N=%d, M=%d, S0=%f$\n\n', T, N, M, actual_data(1));
        end
        % Time discretization
        dt = 1;
        t = 0:dt:T;
        
        if display == 1
            fprintf('Using one simulation\n');
        end
        S = zeros(1, N+1);
        S(:,1) = actual_data(1);
        for i = 1:N
             S(:,i+1) = S(:,i) .* exp((muo - 0.5 *sigmao^2) * dt + sigmao .* sqrt(dt) .* z);
        end
        if display == 1
            fprintf('The MAPE for one simulation is %f %%\n\n', mape(S',actual_data));
        end
        S_onesimulation = S;
        MAPE_onesimulation = mape(S_onesimulation', actual_data);
        if display == 1
            fprintf('Using various simulations and applying mean, median and mode\n')
        end
        S = zeros(1, N+1);
        S(1) = actual_data(1);
        for i = 1:N
            S_day = S(i).* exp((muo - 0.5 *sigmao^2) * dt + sigmao .* sqrt(dt) .* Z);
            S(i+1) = mean(S_day);
        end
        S_manysimulations = S;
        if display == 1
            fprintf('The MAPE for the mean of many simulations is %f %% \n', mape(S_manysimulations', actual_data));
        end
        MAPE_manysimulations = mape(S_manysimulations', actual_data);
        
        if display == 1
            fprintf('Using closed form solution of mean, median and mode\n');
        end
        S = zeros(1, N+1);
        S(:,1) = actual_data(1);
        for i = 1:N
             S(:, i+1) = S(:,i) .*exp((muo)*dt);
        end
        if display == 1
            fprintf('The MAPE for the closed form of the mean (expectation) is %f %% \n', mape(S', actual_data));
        end
        S_mean = S;
        MAPE_mean = mape(mean(S_mean)', actual_data);
        S = zeros(1, N+1);
        S(:,1) = actual_data(1);
        for i=1:N 
            S(:, i+1) = S(:,i) .*exp((muo - 1/2*sigmao^2)*dt);
        end
        if display == 1
            fprintf('The MAPE for the closed form of the median is %f %% \n', mape(S', actual_data));
        end
        S_median = S;
        MAPE_median = mape(mean(S_median)', actual_data);
        S = zeros(1, N+1);
        S(:,1) = actual_data(1);
        for i=1:N 
           S(:, i+1) = S(:,i) .*exp((muo - 3/2*sigmao^2)*dt);
        end
        if display == 1
            fprintf('The MAPE for the closed form of the mode is %f %% \n', mape(S', actual_data));
        end
        S_mode = S;
        MAPE_mode = mape(mean(S_mode)', actual_data);
        if images == 1 || final_image==1
            
            figure;
            p1=plot(preddates, actual_data', 'LineWidth', 1,'Marker', 's' );
            hold on;
            p2=plot(preddates, S_onesimulation', 'LineWidth', 1,'Marker', 'x');
            p3=plot(preddates, mean(S_manysimulations)', 'LineWidth', 1, 'Marker', 'x');
            p4=plot(preddates, median(S_manysimulations)', 'LineWidth', 1,'Marker', 'x');
            p6=plot(preddates, S_mean', 'LineWidth', 1,'Marker', 'x');
            p7=plot(preddates, S_median', 'LineWidth', 1, 'Marker', 'x');
            p8=plot(preddates, S_mode', 'LineWidth', 1, 'Marker', 'x');
            entry1 = 'Actual data';
            entry2 = ['One simulation, MAPE = ', num2str(mape(S_onesimulation', actual_data)),'%'];
            entry3 = ['Average of M simulations, MAPE = ', num2str(mape(mean(S_manysimulations)', actual_data)),'%'];
            entry4 = ['Median of many simulations, MAPE = ', num2str(mape(median(S_manysimulations)', actual_data)),'%'];
            entry6 = ['Closed form mean, MAPE = ', num2str(mape(S_mean', actual_data)),'%'];
            entry7 = ['Closed form median, MAPE = ', num2str(mape(S_median', actual_data)),'%'];
            entry8 = ['Closed form mode, MAPE = ', num2str(mape(S_mode', actual_data)),'%'];
            legend([p1, p2, p3, p4, p6, p7, p8], entry1, entry2, entry3, entry4, entry6, entry7, entry8);
            xlabel('Time (minutes)');
            ylabel('Closing price ($)');
            xlim([preddates(1), preddates(end)]);
            title(sprintf('Prediction vs Real AAPL stock data. Training data: %s - %s Test data: %s - %s \n ', datestr(traindates(1), 'dd-MMM-yyyy HH:mm'), datestr(traindates(end), 'dd-MMM-yyyy HH:mm'), datestr(preddates(1), 'dd-MMM-yyyy HH:mm'), datestr(preddates(end), 'dd-MMM-yyyy HH:mm')));
        end
    catch ME
        warning('An error occurred: %s', ME.message);
        % Return empty values
         MAPE_onesimulation = [];
        MAPE_manysimulations = [];
        MAPE_mean = [];
        MAPE_median = [];
        MAPE_mode = [];
    end
end
