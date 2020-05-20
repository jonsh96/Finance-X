classdef portfolio
    % About
    % Portfolio object, consists of multiple bond objects
    % 
    % Properties accessible: 
    % tickers, issue dates, maturity dates, coupon type, durations, 
    % ask prices, bid prices, prices, last prices, last yields, yields, 
    % interests (coupon rate), curve dates, zero rates, forward rates,
    % discount rates, swap rates and the current curve.
    
    properties
        ticker = [];
        issue = [];
        maturity = [];
        coupon = [];
        frequency = [];
        duration = [];
        ask = [];
        bid = [];
        price = [];
        lastPrice = [];
        lastYield = [];
        yield = [];
        interest = [];
        curveDates
        zeroRates
        forwardRates
        discountRates
        swapRates
        currentCurve = [];
    end
    
    methods
        function obj = portfolio(bond)
            % INITIALISING THE PORTFOLIO BY ADDING A SINGLE BOND
            obj.ticker = string(bond.ticker);
            obj.duration = double(bond.duration);
            obj.issue = string(bond.issue);
            obj.maturity = string(bond.maturity);
            obj.ask = bond.ask;
            obj.bid = bond.bid;
            obj.price = bond.price;
            obj.lastPrice = bond.lastPrice;
            obj.lastYield = bond.lastYield;
            obj.yield = bond.yield;
            obj.coupon = string(bond.coupon);
            
            % Most bonds are conventional coupons or bullet bonds but in
            % some cases the coupon payment is semi-annual
            if(string(bond.coupon) == "Semiannual")
                obj.frequency = 2;
            else
                obj.frequency = 1;
            end
            obj.interest = bond.interest;
            
        end
       
        function obj = addToPortfolio(obj, bond)
            % ADDING A BOND TO AN EXISTING PORTFOLIO
            obj.ticker = horzcat(obj.ticker, string(bond.ticker));
            obj.duration = horzcat(obj.duration, bond.duration);
            obj.issue = horzcat(obj.issue, string(bond.issue));
            obj.maturity = horzcat(obj.maturity, bond.maturity);
            obj.ask = horzcat(obj.ask, bond.ask);
            obj.bid = horzcat(obj.bid, bond.bid);
            obj.price = horzcat(obj.price, bond.price);
            obj.lastPrice = horzcat(obj.lastPrice, bond.lastPrice);
            obj.lastYield = horzcat(obj.lastYield, bond.lastYield);
            obj.yield = horzcat(obj.yield, bond.yield);
            obj.coupon = horzcat(obj.coupon, string(bond.coupon));
            
            % Most bonds are conventional coupons or bullet bonds but in
            % some cases the coupon payment is semi-annual
            if(string(bond.coupon) == "Semiannual")
                obj.frequency = horzcat(obj.frequency, 2);
            else
                obj.frequency = horzcat(obj.frequency, 1);
            end
            obj.interest = horzcat(obj.interest, bond.interest);
            obj = obj.calculateCurves;
        end
        
        function yieldCurve(obj)  
            % PLOTTING THE YIELD CURVE AS IT APPEARS ON WWW.BONDS.IS
            scatter(datenum(obj.maturity,'dd/mm/yyyy'),obj.yield*100,'k','filled');
            grid on
            datetick('x','dd/mm/yyyy')
            ytickformat('%.2f%%')
            xlim([today max(datenum(obj.maturity,'dd/mm/yyyy'))]);
        end
        
        function obj = zeroCurve(obj)
           % PLOTTING THE ZERO RATE CURVE FROM THE PORTFOLIO
           scatter(datenum(obj.maturity,'dd/mm/yyyy'),obj.zeroRates*100,'k','filled');
           grid on
           datetick('x','dd/mm/yyyy')
           ytickformat('%.2f%%')
           xlim([today max(datenum(obj.maturity,'dd/mm/yyyy'))]);
        end        
        
        function obj = forwardCurve(obj)
            % PLOTTING THE FORWARD RATE CURVE FROM THE PORTFOLIO
            scatter(datenum(obj.maturity,'dd/mm/yyyy'),obj.forwardRates*100,'k','filled');
            grid on
            ytickformat('%.2f%%')
            datetick('x','dd/mm/yyyy')
            xlim([today max(datenum(obj.maturity,'dd/mm/yyyy'))]);
        end        
        
        function obj = discountCurve(obj)
            % PLOTTING THE DISCOUNT RATE CURVE FROM THE PORTFOLIO
            scatter(datenum(obj.maturity,'dd/mm/yyyy'),obj.discountRates*100,'k','filled');
            grid on
            ytickformat('%.2f%%')
            datetick('x','dd/mm/yyyy')
            xlim([today max(datenum(obj.maturity,'dd/mm/yyyy'))]);
        end
        
        function obj = swapCurve(obj)
            % CALCULATING AND PLOTTING THE SWAP RATE CURVE FROM THE PORTFOLIO        
            scatter(datenum(obj.maturity,'dd/mm/yyyy'),obj.swapRates*100,'k','filled');
            grid on
            ytickformat('%.2f%%')
            datetick('x','dd/mm/yyyy')
            xlim([today max(datenum(obj.maturity,'dd/mm/yyyy'))]);
        end
        
        
        function obj = calculateCurves(obj)
           % CALCULATING THE ZERO, FORWARD, DISCOUNT AND SWAP RATE CURVES FROM THE PORTFOLIO 
            Bonds = [datenum(obj.maturity) obj.interest' 100*ones(length(obj.ticker),1) obj.frequency' 8*ones(length(obj.ticker),1)];
            Prices = obj.price;
            Settle = today();
            [zeroRates, curveDates] = zbtprice(Bonds, Prices, Settle);
            [forwardRates, curveDates] = zero2fwd(zeroRates, curveDates, Settle);
            [discRates, curveDates] = zero2disc(zeroRates, curveDates, Settle);
            for i = 1:length(obj.discountRates)
                % Finding the bond's coupon frequency 
                %   Alpha is the time period between coupons per bond
                if obj.frequency(i) == 2
                    alpha(i) = 0.5; 
                else
                    alpha(i) = 1;
                end
                
                %Forward rate vector for i time periods
                frwRate(i) = obj.forwardRates(i)*100; 
                dcRate(i) = obj.discountRates(i)*100;
                %Calculation of the Swap rate for i time periods
                obj.swapRates(i) = (sum(alpha.*dcRate.*frwRate)/sum(alpha.*dcRate))/100; 
            end
            obj.curveDates = curveDates';
            obj.zeroRates = zeroRates';
            obj.forwardRates = forwardRates';
            obj.discountRates = discRates';
        end
        
        function obj = fitMethod(obj, curve, method, polyDegree, smoothingFactor)
            % USING A FITTING METHOD ON A CURVE
            % INPUT 1 = CURVE (STRING)
            % INPUT 2 = METHOD (STRING)
            % INPUT 3 = POLYNOMIAL DEGREE (INTEGER)
            % INPUT 4 = SMOOTHING FACTOR (DOUBLE)
            % 
            %   CURVES AVAILABLE: 
            %       "Yield", "Zero rates", "Forward rates", "Discount rates", "Swap rates"
            %   METHODS AVAILABLE:
            %       "Bootstrapping", "Nelson-Siegel", "Polynomial", "Spline", "Cubic spline", "Constrained cubic spline"
            %   POLYNOMIAL DEGREE: 1 <= X <= NUMBER OF BONDS IN PORTFOLIO
            %   SMOOTHING FACTOR: 0 <= X <= 1
            
            % Choosing the curve
            obj = obj.calculateCurves;
            dates = datenum(obj.maturity,'dd/mm/yyyy');
            if curve == "Yield"
                rates = obj.yield;
            elseif curve == "Zero rates"
                rates = obj.zeroRates;
            elseif curve == "Forward rates"
                rates = obj.forwardRates;
            elseif curve == "Discount rates"
                rates = obj.discountRates;
            elseif curve == "Swap rates"
                rates = obj.swapRates;
            end
            
            % Choosing the fitting method and plotting
            if method == "Bootstrapping"
                hold on
                plot(dates, rates)
            elseif method == "Nelson-Siegel"
                % SPECIAL CASE 
                % REQUIRES FINANCIAL INSTRUMENTS TOOLBOX FOR IRFunctionCurve 
                Settle = repmat(today,[length(obj.maturity) 1]);
                Maturity = datenum(obj.maturity);
                CleanPrice = obj.price';
                CouponRate = obj.interest';
                Instruments = [Settle Maturity CleanPrice CouponRate];
                PlottingPoints = linspace(today,max(dates),max(dates)-today);
                % Type is either 'Forward' or 'Zero'
                if curve == "Zero rates" || curve == "Yield"
                    NSModel = IRFunctionCurve.fitNelsonSiegel('Zero',today,Instruments);
                    NSModel.Parameters;
                elseif curve == "Forward rates"
                    NSModel = IRFunctionCurve.fitNelsonSiegel('Forward',today,Instruments);
                    NSModel.Parameters;
                end
                plot(PlottingPoints, getParYields(NSModel, PlottingPoints)*100)
                obj.currentCurve = getParYields(NSModel, PlottingPoints)*100;
                datetick('x','dd/mm/yyyy')
                ytickformat('%.2f%%')
                xlim([today max(dates)])
            elseif method == "Nelson-Siegel-Svensson"
                % SPECIAL CASE 
                % REQUIRES FINANCIAL INSTRUMENTS TOOLBOX FOR IRFunctionCurve 
                Settle = repmat(today,[length(obj.maturity) 1]);
                Maturity = datenum(obj.maturity);
                CleanPrice = obj.price';
                CouponRate = obj.interest';
                Instruments = [Settle Maturity CleanPrice CouponRate];
                PlottingPoints = linspace(today,max(dates),max(dates)-today);
                
                % Type is either 'Forward' or 'Zero'
                if curve == "Zero rates" || curve == "Yield"
                    SvenssonModel  = IRFunctionCurve.fitSvensson('Zero',today,Instruments);
                    SvenssonModel.Parameters;
                elseif curve == "Forward rates"
                    SvenssonModel  = IRFunctionCurve.fitSvensson('Forward',today,Instruments);
                    SvenssonModel.Parameters;
                end
                plot(PlottingPoints, getParYields(SvenssonModel , PlottingPoints)*100)
                obj.currentCurve = getParYields(SvenssonModel , PlottingPoints)*100;
                datetick('x','dd/mm/yyyy')
                ytickformat('%.2f%%')
                xlim([today max(dates)])
            elseif method == "Polynomial"
                obj = obj.polynomialFit(dates, rates, polyDegree);
            elseif method == "Lagrange interpolation"
                obj = obj.lagrangeFit(dates,rates);
            elseif method == "Spline"
                obj = obj.splineFit(dates, rates);
            elseif method == "Cubic spline"
                obj = obj.cubicSplineFit(dates, rates);
            elseif method == "Constrained cubic spline"
                obj = obj.constrainedCubicSplineFit(dates, rates, smoothingFactor);
            end
        end    
        
        function obj = polynomialFit(obj,dates,rates,n)
            % APPLYING A N-TH DEGREE POLYNOMIAL FIT TO A SET OF DATA
            % INPUT 1: DATES (X VALUES)
            % INPUT 2: RATES (Y VALUES)
            % INOUT 3: N (POLYNOMIAL DEGREE)

            hold on
            degree = n;
            % Turning off warning about possible bad fit 
            ws = warning('off','all');  
            p = polyfit(dates', rates, degree);
            % Turning the warnings back on.
            warning(ws)  
            px = linspace(today,max(dates),max(dates)-today);
            py = polyval(p, px);
            plot(px, py*100)
            obj.currentCurve = py*100;
            grid on
            ytickformat('%.2f%%')
            datetick('x','dd/mm/yyyy')
            xlim([today max(dates)])
        end
        
        function obj = lagrangeFit(obj, dates, rates)            
            % APPLYING LAGRANGE INTERPOLATION FIT TO A SET OF DATA
            % INPUT 1: DATES (X VALUES)
            % INPUT 2: RATES (Y VALUES)
            
            xx = linspace(today,max(dates),max(dates)-today);
            P = lagrangepoly(dates,rates);
            plot(xx,polyval(P,xx)*100);
            obj.currentCurve = polyval(P,xx)*100;
            grid on
            ytickformat('%.2f%%');
            datetick('x','dd/mm/yyyy');
            xlim([today max(dates)]);
        end
        
        function obj = splineFit(obj, dates, rates)
            % APPLYING SPLINE FIT TO A SET OF DATA
            % INPUT 1: DATES (X VALUES)
            % INPUT 2: RATES (Y VALUES)
            
            hold on
            xsp = linspace(today,max(dates),max(dates)-today);
            sp = spline(dates,rates);
            plot(xsp,ppval(sp,xsp)*100)
            obj.currentCurve = ppval(sp,xsp)*100;
            ytickformat('%.2f%%')
            datetick('x','dd/mm/yyyy')
            xlim([today max(dates)])
        end
        
        function obj = cubicSplineFit(obj, dates, rates)
            % APPLYING CUBIC SPLINE FIT TO A SET OF DATA
            % INPUT 1: DATES (X VALUES)
            % INPUT 2: RATES (Y VALUES)

            hold on
            cs = csaps(dates,rates);
            xsp = linspace(today,max(dates),max(dates)-today);
            plot(xsp,ppval(cs,xsp)*100)
            obj.currentCurve = ppval(cs,xsp)*100;
            grid on
            ytickformat('%.2f%%')
            datetick('x','dd/mm/yyyy')
            xlim([today max(dates)])
        end
        
        function obj = constrainedCubicSplineFit(obj, dates, rates, smoothingFactor)
            % APPLYING CONSTRAINED CUBIC SPLINE FIT TO A SET OF DATA
            % INPUT 1: DATES (X VALUES)
            % INPUT 2: RATES (Y VALUES)
            % INPUT 3: SMOOTHING FACTORS [0,1] (DOUBLE)

            hold on
            x = dates;
            y = rates; 
            x_max = max(x);
            x_min = today;
            x_scaled = 10* (x - x_min) / (x_max - x_min);
            xsp = linspace(0,10,max(dates)-min(dates));
            cs = csaps(x_scaled,y,smoothingFactor);
            plot(xsp * (x_max - x_min) / 10 + x_min, ppval(cs,xsp)*100,'-');
            obj.currentCurve = ppval(cs,xsp)*100;
            grid on
            ytickformat('%.2f%%')   
            datetick('x','dd/mm/yyyy')
            xlim([today max(dates)])
        end
    end
end

% HELPER FUNCTIONS
function [P,R,S] = lagrangepoly(X,Y,XX)
    %   LAGRANGEPOLY  Lagrange interpolation polynomial fitting a set of points
    %   [P,R,S] = LAGRANGEPOLY(X,Y)  where X and Y are row vectors
    %   defining a set of N points uses Lagrange's method to find 
    %   the N-1th order polynomial in X that passes through these 
    %   points.  P returns the N coefficients defining the polynomial, 
    %   in the same order as used by POLY and POLYVAL (highest order first).
    %   Then, polyval(P,X) = Y.  R returns the x-coordinates of the N-1
    %   extrema of the resulting polynomial (roots of its derivative),
    %   and S returns the y-values  at those extrema.
    %
    %   YY = LAGRANGEPOLY(X,Y,XX) returns the values of the polynomial
    %   sampled at the points specified in XX -- the same as
    %   YY = POLYVAL(LAGRANGEPOLY(X,Y)). 
    %
    %   Example:
    %   To find the 4th-degree polynomial that oscillates between 
    %   1 and 0 across 5 points around zero, then plot the interpolation
    %   on a denser grid inbetween:
    %     X = -2:2;  Y = [1 0 1 0 1];
    %     P = lagrangepoly(X,Y);
    %     xx = -2.5:.01:2.5;
    %     plot(xx,polyval(P,xx),X,Y,'or');
    %     grid;
    %   Or simply:
    %     plot(xx,lagrangepoly(X,Y,xx));
    %
    %   Note: if you are just looking for a smooth curve passing through 
    %   a set of points, you can get a better fit with SPLINE, which 
    %   fits piecewise polynomials rather than a single polynomial.
    %
    %   See also: POLY, POLYVAL, SPLINE

    % 2006-11-20 Dan Ellis dpwe@ee.columbia.edu
    % $Header: $

    %  For more info on Lagrange interpolation, see Mathworld: 
    %  http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html

    % Make sure that X and Y are row vectors
    if size(X,1) > 1;  X = X'; end
    if size(Y,1) > 1;  Y = Y'; end
    if size(X,1) > 1 || size(Y,1) > 1 || size(X,2) ~= size(Y,2)
      error('both inputs must be equal-length vectors')
    end

    N = length(X);

    pvals = zeros(N,N);

    % Calculate the polynomial weights for each order
    for i = 1:N
      % the polynomial whose roots are all the values of X except this one
      pp = poly(X( (1:N) ~= i));
      % scale so its value is exactly 1 at this X point (and zero
      % at others, of course)
      pvals(i,:) = pp ./ polyval(pp, X(i));
    end

    % Each row gives the polynomial that is 1 at the corresponding X 
    % point and zero everywhere else, so weighting each row by the 
    % desired row and summing (in this case the polycoeffs) gives 
    % the final polynomial
    P = Y*pvals;

    if nargin==3
      % output is YY corresponding to input XX
      YY = polyval(P,XX);
      % assign to output
      P = YY;
    end

    if nargout > 1
      % Extra return arguments are values where dy/dx is zero
      % Solve for x s.t. dy/dx is zero i.e. roots of derivative polynomial
      % derivative of polynomial P scales each power by its power, downshifts
      R = roots( ((N-1):-1:1) .* P(1:(N-1)) );
      if nargout > 2
        % calculate the actual values at the points of zero derivative
        S = polyval(P,R);
      end
    end
end