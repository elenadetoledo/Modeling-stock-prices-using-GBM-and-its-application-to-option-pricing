function [S_onesimulation, S_manysimulations, S_mean, S_median, S_mode] = fmodel2(start_date_test, ndaystraining, ndayspredicting, M, callibration_option, images, Z, z)
    %MODEL 2: initial step S(i-1). (t-T0) constant, =1.muo & sigmao are constant
    %For the callibration option, choose 1 for using fitdist, choose 2 to use
    %the methods explained in Farida (sample mean, sample standard deviation)
    %Show images flag. images==0 no images, images==1 show images
    final_image=1; % Change to 0 if you dont want to see any images
    if Z == 0
        fprintf('generando randomness 1\n');
        Z = randn(M,1);
    end

    if z == 0
         fprintf('generando randomness 2\n');
        z = randn;
    end

    fprintf('\nMAKING PREDICTIONS WITH MODEL 2. USING S(i-1) AS BASIS AND CONSTANT MUO AND SIGMA\n');
    %PARAMETERS OF THE MODEL
    colors = [
        0, 174, 239;    % Vivid Cyan
        220, 20, 60;    % Crimson Red
        46, 204, 113;   % Emerald Green
        241, 196, 15;   % Sunflower Yellow
        155, 89, 182;   % Amethyst Purple
        230, 126, 34;   % Orange Hue
        191, 255, 0;    % Electric Lime
        52, 152, 219    % Deep Sky Blue
    ] / 255;  % Normalize the RGB values

    warning('off', 'stats:jbtest:PTooBig');
    warning('off', 'stats:jbtest:PTooSmall');
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    
    filename_5years = '..\data\data_AAPL_5years_27april2019_27april2024.csv';
    dataTable_total = readtable(filename_5years);
    closeData_total = dataTable_total.Close;
    datesData_total = dataTable_total.Date;
    
    % fprintf('You entered the following startdate: %s\n', start_date_test);
    test_date = datetime(start_date_test, 'InputFormat', 'dd-MMM-yyyy', 'Locale', 'en_US');
    idx_startTest = find(datesData_total == test_date);
    % fprintf('Position in the data: %d\n', idx_startTest);
   
    traindates = datesData_total(idx_startTest - ndaystraining : idx_startTest - 1);
    St = closeData_total(idx_startTest - ndaystraining : idx_startTest - 1);
    ratio=St(2:end)./St(1:end-1); %because we can fit distribution with Δt = 1
    logreturns = log(ratio);
    
    actual_data = closeData_total(idx_startTest : idx_startTest + ndayspredicting - 1); 
    preddates = datesData_total(idx_startTest : idx_startTest + ndayspredicting - 1);
    % fprintf('Train data: %s - %s (%d data points) , Test data: %s - %s (%d data points)\n', datestr(traindates(1), 'dd-mmm-yyyy'), datestr(traindates(end), 'dd-mmm-yyyy'), length(traindates), datestr(preddates(1), 'dd-mmm-yyyy'), datestr(preddates(end), 'dd-mmm-yyyy'), length(preddates));
    % fprintf('TRAINING: Hypothesis testing\n');
    % fprintf('%s - %s: Sample size Stocks: %d, Max: %f, Min: %f, Mean:%f, Std dev:%f\n', datestr(traindates(1), 'dd-mmm-yyyy'), datestr(traindates(end), 'dd-mmm-yyyy'),length(St), max(St), min(St), mean(St), std(St));
    % 
    dates_logreturns = traindates(2:end);
    if images == 1
        figure();
        plot(traindates, St);
        xlim([traindates(1), traindates(end)]);
        title(sprintf('AAPL Price %s - %s', datestr(traindates(1), 'dd-mmm-yyyy'), datestr(traindates(end), 'dd-mmm-yyyy')));
        ylabel('Closing Price ($)');
        xlabel('Date');
    end
    
    % INDEPENDENCE
    %1 - Visually assessing independence of log returns
    if images == 1
        figure();
        plot(dates_logreturns, logreturns);
        xlim([dates_logreturns(1), dates_logreturns(end)]);
        title(sprintf('Log Returns %s - %s', datestr(traindates(1), 'dd-mmm-yyyy'), datestr(traindates(end), 'dd-mmm-yyyy')));
        ylabel('Log return');
        xlabel('Date');
    
        figure;
        %  A common rule of thumb is to use a number of lags up to 10*log10(N)
        %  --> consejo de Ana: cmbiar a n/4
        autocorr(logreturns, ceil(10*log10(length(logreturns))));  % Compute ACF for 20 lags
        title(sprintf('ACF AAPL Log Returns %s - %s', datestr(traindates(1), 'dd-mmm-yyyy'), datestr(traindates(end), 'dd-mmm-yyyy')));
        xlim([1, ceil(10*log10(length(logreturns)))]);
        ylabel('ACF');
    end
    
    %Analytical tests
    [h_jb,p_jb] = jbtest(logreturns);
    [h_lbq, ~, ~, ~] = lbqtest(logreturns);
    fprintf('%s - %s: Sample size logreturns: %d, Max: %f, Min: %f, Mean:%f, Std dev:%f, Skewness:%f, Jarque-Bera Test:%d, p-value:%f, Ljung-Box Q-test:%d\n', datestr(traindates(1), 'dd-mmm-yyyy'), datestr(traindates(end), 'dd-mmm-yyyy'), length(logreturns), max(logreturns), min(logreturns), mean(logreturns), std(logreturns), skewness(logreturns), h_jb, p_jb, h_lbq);
    
    if h_jb == 0
        fprintf('** LOGRETURNS ARE NORMAL\n');
    end
    
    if h_lbq == 0
        fprintf('** LOGRETURNS ARE UNCORRELATED\n');
    end
    
    %Plotting the histogram
    if images == 1
        ceil(sqrt(length(logreturns)));
        figure;
        histogram(logreturns, ceil(sqrt(length(logreturns))));
        title(sprintf('Histogram of Log Returns %s - %s', datestr(traindates(1), 'dd-mmm-yyyy'), datestr(traindates(end), 'dd-mmm-yyyy')));
        xlabel('Value');
        ylabel('Frequency');
    end
    
    %QQplot
    if images == 1
        figure;
        qqplot(logreturns);
        title(sprintf('Q-Q Plot of Data vs. Standard Normal %s - %s', datestr(traindates(1), 'dd-mmm-yyyy'), datestr(traindates(end), 'dd-mmm-yyyy')));
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
        fprintf('TRAINING: Model calibration using fitdist (MLE)\n');
    else %callibration_option == 2
        stock_return = (St(2:end) - St(1:end-1))./St(1:end-1);
        muo = 1/length(stock_return) * sum(stock_return);
        sigmao = sqrt(1/(length(stock_return)-1)*sum((stock_return - mean(stock_return)).^2));
         fprintf('TRAINING: Model calibration using sample mean and sample std dev of arithmetic returns\n');
    end
    fprintf('Calibrated parameters: mu=%f sigma=%f\n', muo, sigmao);
    
    % PARAMETERS OF THE SIMULATION
    T = length(actual_data)-1;             % (*) Time at which we want to assess the price (in minutes)
    N = length(actual_data)-1;           %  (*) Number of time steps
    fprintf('TESTING: Model creation\n');
    fprintf('Parameters of testing: T=%d, N=%d, M=%d, S0=%f$\n\n', T, N, M, actual_data(1));
    % Time discretization
    dt = T / N;
    t = 0:dt:T;
    
    
    fprintf('Using one simulation\n');
    S = zeros(1, N+1);
    S(:,1) = actual_data(1);
    for i = 1:N
         S(:,i+1) = S(:,i) .* exp((muo - 0.5 *sigmao^2) * dt + sigmao .* sqrt(dt) .* z);
    end
    fprintf('The MAPE for one simulation is %f %%\n\n', mape(S',actual_data));
    S_onesimulation = S;
    
    
    fprintf('Using various simulations and applying mean, median and mode\n')
    S = zeros(M, N+1);
    S(:,1) = actual_data(1);
    for i = 1:N
         S(:,i+1) = S(:,i) .* exp((muo - 0.5 *sigmao^2) * dt + sigmao .* sqrt(dt) .* Z);
    end
    fprintf('The MAPE for the mean is %f %% \n', mape(mean(S)', actual_data));
    fprintf('The MAPE for the median is %f %% \n',mape(median(S)', actual_data));
    % fprintf('The MAPE for the mode is %f %% \n', mape(mode(S)', actual_data));
    % fprintf('The MAPE for the pct 55 is %f %% \n', mape(prctile(S,55)', actual_data));
    % fprintf('The MAPE for the pct 45 is %f %% \n', mape(prctile(S,45)', actual_data));
    % fprintf('The MAPE for the pct 40 is %f %% \n', mape(prctile(S,40)', actual_data));
    % fprintf('The MAPE for the pct 35 is %f %% \n', mape(prctile(S,35)', actual_data));
    % fprintf('The MAPE for the pct 30 is %f %% \n\n', mape(prctile(S,30)', actual_data));
    S_manysimulations = S;
    
    fprintf('Using closed form solution of mean, median and mode\n');
    S = zeros(1, N+1);
    S(:,1) = actual_data(1);
    for i = 1:N
         S(:, i+1) = S(:,i) .*exp((muo)*dt);
    end
    fprintf('The MAPE for the closed form of the mean (expectation) is %f %% \n', mape(S', actual_data));
    S_mean = S;
    S = zeros(1, N+1);
    S(:,1) = actual_data(1);
    for i=1:N 
        S(:, i+1) = S(:,i) .*exp((muo - 1/2*sigmao^2)*dt);
    end
    fprintf('The MAPE for the closed form of the median is %f %% \n', mape(S', actual_data));
    S_median = S;
    S = zeros(1, N+1);
    S(:,1) = actual_data(1);
    for i=1:N 
       S(:, i+1) = S(:,i) .*exp((muo - 3/2*sigmao^2)*dt);
    end
    fprintf('The MAPE for the closed form of the mode is %f %% \n', mape(S', actual_data));
    S_mode = S;
    if images == 1 || final_image==1
        figure;
        p1=plot(preddates, actual_data', 'LineWidth', 1,'Color',[0, 0, 0.5],'Marker', 's', 'MarkerFaceColor', [0, 0, 0.5] );
        hold on;
        p2=plot(preddates, S_onesimulation', 'LineWidth', 1, 'Color', colors(1,:),'Marker', 'x' , 'MarkerFaceColor', colors(1,:));
        p3=plot(preddates, mean(S_manysimulations)', 'LineWidth', 1, 'Color', colors(2,:),'Marker', 'x' , 'MarkerFaceColor', colors(2,:));
        p4=plot(preddates, median(S_manysimulations)', 'LineWidth', 1, 'Color', colors(3,:),'Marker', 'x' , 'MarkerFaceColor', colors(3,:));
        p5=plot(preddates, mode(S_manysimulations)', 'LineWidth', 1, 'Color', colors(4,:),'Marker', 'x' , 'MarkerFaceColor', colors(4,:));
        p6=plot(preddates, S_mean', 'LineWidth', 1, 'Color', colors(5,:),'Marker', 'x' , 'MarkerFaceColor', colors(5,:));
        p7=plot(preddates, S_median', 'LineWidth', 1, 'Color', colors(6,:),'Marker', 'x' , 'MarkerFaceColor', colors(6,:));
        p8=plot(preddates, S_mode', 'LineWidth', 1, 'Color', colors(7,:),'Marker', 'x' , 'MarkerFaceColor', colors(7,:));
        entry1 = 'Actual data';
        entry2 = ['One simulation, MAPE = ', num2str(mape(S_onesimulation', actual_data)),'%'];
        entry3 = ['Mean of many simulations, MAPE = ', num2str(mape(mean(S_manysimulations)', actual_data)),'%'];
        entry4 = ['Median of many simulations, MAPE = ', num2str(mape(median(S_manysimulations)', actual_data)),'%'];
        entry5 = ['Mode of many simulations, MAPE = ', num2str(mape(mode(S_manysimulations)', actual_data)),'%'];
        entry6 = ['Closed form mean, MAPE = ', num2str(mape(S_mean', actual_data)),'%'];
        entry7 = ['Closed form median, MAPE = ', num2str(mape(S_median', actual_data)),'%'];
        entry8 = ['Closed form mode, MAPE = ', num2str(mape(S_mode', actual_data)),'%'];
        legend([p1, p2, p3, p4, p5, p6, p7, p8], entry1, entry2, entry3, entry4, entry5, entry6, entry7, entry8);
        xlabel('Time (minutes)');
        ylabel('Closing price ($)');
        xlim([preddates(1), preddates(end)]);
        title(sprintf('MODEL 2: Model accuracy vs Real AAPL stock data. Training data: %s - %s Test data: %s - %s \n ', datestr(traindates(1), 'dd-mmm-yyyy'), datestr(traindates(end), 'dd-mmm-yyyy'), datestr(preddates(1), 'dd-mmm-yyyy'), datestr(preddates(end), 'dd-mmm-yyyy')));
    
        %without many simulations mode because it has always a really big mape
        figure;
        p1=plot(preddates, actual_data', 'LineWidth', 1,'Color',[0, 0, 0.5],'Marker', 's', 'MarkerFaceColor', [0, 0, 0.5] );
        hold on;
        p2=plot(preddates, S_onesimulation', 'LineWidth', 1, 'Color', colors(1,:),'Marker', 'x' , 'MarkerFaceColor', colors(1,:));
        p3=plot(preddates, mean(S_manysimulations)', 'LineWidth', 1, 'Color', colors(2,:),'Marker', 'x' , 'MarkerFaceColor', colors(2,:));
        p4=plot(preddates, median(S_manysimulations)', 'LineWidth', 1, 'Color', colors(3,:),'Marker', 'x' , 'MarkerFaceColor', colors(3,:));
        p6=plot(preddates, S_mean', 'LineWidth', 1, 'Color', colors(5,:),'Marker', 'x' , 'MarkerFaceColor', colors(5,:));
        p7=plot(preddates, S_median', 'LineWidth', 1, 'Color', colors(6,:),'Marker', 'x' , 'MarkerFaceColor', colors(6,:));
        p8=plot(preddates, S_mode', 'LineWidth', 1, 'Color', colors(7,:),'Marker', 'x' , 'MarkerFaceColor', colors(7,:));
        entry1 = 'Actual data';
        entry2 = ['One simulation, MAPE = ', num2str(mape(S_onesimulation', actual_data)),'%'];
        entry3 = ['Mean of many simulations, MAPE = ', num2str(mape(mean(S_manysimulations)', actual_data)),'%'];
        entry4 = ['Median of many simulations, MAPE = ', num2str(mape(median(S_manysimulations)', actual_data)),'%'];
        entry6 = ['Closed form mean, MAPE = ', num2str(mape(S_mean', actual_data)),'%'];
        entry7 = ['Closed form median, MAPE = ', num2str(mape(S_median', actual_data)),'%'];
        entry8 = ['Closed form mode, MAPE = ', num2str(mape(S_mode', actual_data)),'%'];
        legend([p1, p2, p3, p4, p6, p7, p8], entry1, entry2, entry3, entry4, entry6, entry7, entry8);
        xlabel('Time (minutes)');
        ylabel('Closing price ($)');
        xlim([preddates(1), preddates(end)]);
        title(sprintf('MODEL 2: Model accuracy vs Real AAPL stock data. Training data: %s - %s Test data: %s - %s \n ', datestr(traindates(1), 'dd-mmm-yyyy'), datestr(traindates(end), 'dd-mmm-yyyy'), datestr(preddates(1), 'dd-mmm-yyyy'), datestr(preddates(end), 'dd-mmm-yyyy')));
    end

end