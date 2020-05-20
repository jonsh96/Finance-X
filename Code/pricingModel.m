classdef pricingModel
    % About
    % Pricing Model, used to simulate zero coupon bond and price option of
    % zero coupon bonds.
    % 
    % Properties accessible: 
    % interestRateModel, maturity, optionMaturity, simulatedBonds,
    % simulatedBondPrices, calculatedBondPrices, strikePrice,
    % simulatedCalls, simulatedPuts, calculatedCall, calculatedPut, caps
    % and floor
    
    properties
        interestRateModel
        maturity
        optionMaturity
        simulatedBonds
        simulatedBondPrices
        %calculatedBondPrices held thetta se ekkert notad
        strikePrice
        simulatedCalls
        simulatedPuts
        calculatedCall
        calculatedPut
        caps
        floors
    end
    
    methods
        
        function obj = pricingModel(interestRateModel, strikePrice, optionMaturity)
            % INITIALISING THE PRICING MODEL FROM INPUTS
            obj.interestRateModel = interestRateModel;
            obj.strikePrice = strikePrice;
            obj.optionMaturity = optionMaturity;
            obj.maturity = interestRateModel.maturity; 
            
            % CALCULATING THE OPTION PRICES
            [obj.simulatedBonds, obj.simulatedBondPrices] = obj.zeroCouponBondSimulation;
            [obj.simulatedCalls, obj.simulatedPuts] = obj.optionSimulation;
            [obj.calculatedCall, obj.calculatedPut] = obj.optionPricer(0);
            %[obj.caps, obj.floors] = obj.capsAndFloor;
        end
        
        function [B, bT] = zeroCouponBondSimulation(obj)
            % ZERO COUPON SIMULATION
            % Returns bond price path and bond price at the option maturity
            interestRateModel = obj.interestRateModel;
            r0 = interestRateModel.initialRate;
            dt = interestRateModel.stepSize;
            sigma = interestRateModel.volatility;
            T = interestRateModel.maturity;
            N = interestRateModel.nrOfSimulations;
            O = obj.optionMaturity;
            
            L = T/dt;
            B = zeros(N,L); 
            bT = zeros(1,N);
            R = interestRateModel.data;
            
            switch interestRateModel.model
                case "Simple"
                    for i = 1:N 
                        for j = 1:L
                            B(i,j) = exp(-R(i,j)*((L-(j-1))*dt)+(1/6)*sigma^2*((L-(j-1))*dt)^3);
                        end
                        bT(1,i) = B(i,O/dt);
                    end 
                case "Brownian"
                    alpha = interestRateModel.longTermMeanLevel;
                    for i = 1:N 
                        for j = 1:L
                            B(i,j) = exp(-R(i,j)*((L-(j-1))*dt)-(1/2)*alpha*((L-(j-1))*dt)^2+(1/6)*sigma^2*((L-(j-1))*dt)^3);
                        end
                        bT(1,i) = B(i,O/dt);
                    end 
                case "Vasicek"
                    kappa = interestRateModel.speedOfReversion;
                    theta = interestRateModel.longTermMeanLevel;
                    
                    for i = 1:N 
                        for j = 1:L
                            BtT = (1/kappa)*(1-exp(-kappa*(L-(j-1))*dt));
                            AtT = exp((theta-(sigma^2)/(2*kappa^2))*(BtT-(L-(j-1))*dt)-((sigma^2)/(4*kappa))*BtT);
                            B(i,j) = AtT*exp(-R(i,j)*BtT);
                        end
                        bT(1,i) = B(i,O/dt);
                    end      
            end
        end
        
        function histModel(obj)
           % PLOTTING THE HISTOGRAM OF THE DISTRIBUTION OF BOND PRICES AT OPTION EXPIRY
           
           histfit(obj.simulatedBondPrices, 100,'lognormal')
           grid on
        end
        
        
        function [callPrice, putPrice] = optionSimulation(obj)
            % Pricing call and put prices from the simulated bond prices
            K = obj.strikePrice;
            T = obj.optionMaturity;
            N = T/obj.interestRateModel.stepSize; 
            callPayoff = mean(max(obj.simulatedBonds(:,N)-K,0));
            putPayoff = mean(max(K-obj.simulatedBonds(:,N),0));
            
            r0 = obj.interestRateModel.initialRate;
            callPrice = exp(-r0*T)*callPayoff;
            putPrice = exp(-r0*T)*putPayoff; 
        end
        
        function obj = optionPathSimulation(obj)
            % Shows how the price of the call and put options develop over
            % time as time reaches the option's maturity
            call = [];
            put = [];
            dt = obj.interestRateModel.stepSize;
            for i = 0:dt:obj.optionMaturity
                [call(end+1), put(end+1)] = obj.optionPricer(i);
            end
            plot(call)
            hold on
            plot(put)
        end
        
        
        function plotBonds(obj)
            % Plotting the simulated zero coupon bonds
            N = obj.interestRateModel.nrOfSimulations;
            dates = linspace(today(),today()+365*obj.maturity, 250*obj.maturity);
            for i = 1:N
                plot(dates,obj.simulatedBonds(i,:),'HandleVisibility','off')%,'Color',[0.5 0.5 0.5])
                hold on
            end
            K = obj.strikePrice;
            plot([min(dates) max(dates)],[K K],'k-','LineWidth',1.5)
            plot(today+365*obj.optionMaturity,K,'k+','LineWidth',5)
            datetick('x','dd/mm/yyyy')
            legend('Strike price','Exercise date')
            grid on
            xlim([min(dates) max(dates)])
        end

        function [C, P] = optionPricer(obj,t)    
            % Pricing the options at time t. t = 0 at the beginning
            model = obj.interestRateModel;
            dt = obj.interestRateModel.stepSize;
            sigma = obj.interestRateModel.volatility;
            
            S = obj.maturity-t; % S > T
            T = obj.optionMaturity-t;
            
            B = obj.simulatedBonds(1,1); %starting price
            R = obj.interestRateModel.data;
            K = obj.strikePrice;
            
            rt = obj.interestRateModel.initialRate;
            rT = mean(R(:,(T+t)/dt));
            rS = mean(R(:,(S+t)/dt));

            % Call
            switch model.model
                case "Simple"                 
                    [C, P] = blkprice(B,K,rT,T,sigma);

                case "Brownian"
                    alpha = model.longTermMeanLevel;

                    PtT = exp( -rT * T - (1/2) * alpha * T^2 + (1/6) * sigma^2 * T^3);
                    d1 = log(B/K + (sigma^2)*(T/2))/(sigma*sqrt(T));
                    d2 = log(B/K - (sigma^2)*(T/2))/(sigma*sqrt(T));
                    % Call
                    N1 = normcdf(d1);
                    N2 = normcdf(d2);
                    C = PtT*(B*N1-K*N2);
                    % Put
                    N1 = normcdf(-d1);
                    N2 = normcdf(-d2);
                    P = PtT*(K*N2-B*N1);
                    
                case "Vasicek"
                    Q = 1; % Principal of the bond
                    kappa = model.speedOfReversion;
                    theta = model.longTermMeanLevel;

                    BtT = (1/kappa)*(1-exp(-kappa*T));
                    BtS = (1/kappa)*(1-exp(-kappa*S));
                    AtT = exp( (((kappa^2)*theta-(sigma^2)/2)/(kappa^2))*(BtT-T) - ((sigma^2)/(4*kappa))*(BtT^2));
                    AtS = exp( (((kappa^2)*theta-(sigma^2)/2)/(kappa^2))*(BtS-S) - ((sigma^2)/(4*kappa))*(BtS^2));
                    rt = R(1,1);

                    PtT = AtT *exp(-rt*BtT);
                    PtS = AtS *exp(-rt*BtS);

                    sigmaP = (sigma/kappa)*(1-exp( -kappa*(S-T)) ) * sqrt( (1-exp(-2*kappa*T)) /(2*kappa));
                    d = (1/sigmaP)*log( (Q*PtS)/(K*PtT) ) + sigmaP/2;

                    N = @(x) normcdf(x);
                    
                    C = Q * PtS * N(d) - K * PtT * N(d-sigmaP);
                    P = K * PtT * N(-d+sigmaP) - Q * PtS * N(-d);            
            end
        end
         
        function [caps, floors] = capsAndFloor(obj,Lcc,Lcf, Q) 
            % Calculating the caps and floor
            dt = obj.interestRateModel.stepSize;            
            T = obj.optionMaturity;                         %Time2mat

            sigma = obj.interestRateModel.volatility;
            Lk = obj.interestRateModel.initialRate;                 
            N = obj.interestRateModel.nrOfSimulations;

            %Simulation of Rates and calculation of forward and discount
            %rate
            
            for n = 1:N
                for i = 2:1:T/dt-1
                    dL = sigma* Lk(i-1) *randn*sqrt(dt);
                    Lk(i) = dL + Lk(i-1);
                end
                LL(n,:) = Lk;
            end
            
            % Average of simulated rates
            Lk = [obj.interestRateModel.initialRate mean(LL)];
            
            alpha = dt;                                     %Time period            
            settle = today;                                 %Settle time
            curveDates = today+1:365/(1/dt):today+T*365;
            
            % Calculations of forward and discount rate
            [Fk,date] = zero2fwd(Lk',curveDates',settle');
            [D,cdate] = zero2disc(Lk',curveDates',settle');
            
            for i = dt:dt:T
                n = round(i/dt);
                
                d1 = real((log(Fk(n)/Lcc)) + 0.5*i*sigma^2)/(sigma*sqrt(i));
                d2 = d1-sigma*sqrt(i);
                d11 = real((log(Fk(n)/Lcf)) + 0.5*i*sigma^2)/(sigma*sqrt(i));
                d22 = d1-sigma*sqrt(i);
                
                caplet(n) = alpha*Q*D(n)*(Fk(n)*normcdf(d1)-Lcc*normcdf(d2));
                floorlet(n) = alpha*Q*D(n)*(Lcf*normcdf(-d22)-Fk(n)*normcdf(-d11));
            end
            
            caps = sum(caplet)
            floors = sum(floorlet)            

        end
    end
end

