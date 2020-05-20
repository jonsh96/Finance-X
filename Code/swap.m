classdef swap
    % About
    % Swap object, used to calculate multiple scenarios for swap
    % 
    % Properties accessible: 
    % payPortfolio, receivePortfolio, pay, receive, payments, principal, 
    % couponRate, compounding, startDate, settleDate, endDate, payCurve, 
    % receiveCurve,basisPoints, price, swapRate, AI (Accrued interest), 
    % payCF, payDates, receiveCF, receiveDates
    
    properties
        payPortfolio
        receivePortfolio
        fixedPay
        fixedReceive
        pay 
        receive
        payments % per year 1x2 array 
        principal
        couponRate % TODO
        compounding % TODO 
        startDate
        settleDate
        endDate
        payCurve
        receiveCurve
        basisPoints % 1x2 array
        swapValue
        swapRate
        AI % Accrued interest
        payCF
        payDates
        receiveCF
        receiveDates
    end
    
    methods
        function obj = swap(receivePortfolio, payPortfolio, fixedReceive, fixedPay, receive, pay, payments, principal, startDate, settleDate, endDate, payCurve, receiveCurve, basisPoints)
            % INITIALISING THE SWAP OPJECT
            obj.payPortfolio = payPortfolio;
            obj.receivePortfolio = receivePortfolio;
            obj.fixedPay = fixedPay;
            obj.fixedReceive = fixedReceive;
            obj.pay = pay;
            obj.receive = receive;
            obj.payments = payments;
            obj.principal = principal;
            obj.startDate = startDate;
            obj.settleDate = settleDate;
            obj.endDate = endDate;
            obj.payCurve = payPortfolio.currentCurve;
            obj.receiveCurve = receivePortfolio.currentCurve;
            obj.basisPoints = basisPoints;
            obj = obj.valueSwap;
        end
        
        function obj = valueSwap(obj)
            % Calculating the value of the swap and plotting the
            % (discounted) cash flow for each leg
            startDate = obj.startDate;   
            settleDate = obj.settleDate; 
            endDate = obj.endDate;
            principal = obj.principal;
            
            receiveDates = datenum(settleDate,'dd/mm/yyyy'):round(365/obj.payments(1)):datenum(endDate,'dd/mm/yyyy');
            payDates = datenum(settleDate,'dd/mm/yyyy'):round(365/obj.payments(2)):datenum(endDate,'dd/mm/yyyy');

            receiveCash = zeros(1, length(receiveDates));
            payCash = zeros(1, length(payDates));
            
            % Discount rates
            payDisc = obj.payPortfolio.discountRates;
            payMaturity = datenum(obj.payPortfolio.maturity,'dd/mm/yyyy')';
            
            receiveDisc = obj.receivePortfolio.discountRates;
            receiveMaturity = datenum(obj.receivePortfolio.maturity,'dd/mm/yyyy')';
            
            pReceive = polyfit(receiveMaturity,receiveDisc,1);
            pPay = polyfit(payMaturity,payDisc,1);
            
            discountRec = polyval(pReceive, receiveDates);
            discountPay = polyval(pPay, payDates);

            if ~isnan(obj.fixedReceive)
                % Fixed receive
                receiveCash(:) = obj.fixedReceive*principal/(100*obj.payments(1));
                receiveCash(end) = receiveCash(end)+principal;
                receiveCash = receiveCash.*discountRec;
                receiveBond = sum(receiveCash);
            else
                % Floating receive 
                receiveRates = obj.receiveCurve(receiveDates-datenum(settleDate,'dd/mm/yyyy')+1)/obj.payments(1);
                receiveCash(:) = ((receiveRates+obj.basisPoints(1))/100).*principal;
                receiveCash(end) = receiveCash(end) + principal;
                receiveCash = receiveCash.*discountRec;
                receiveBond = ((receiveRates(1)/100)*principal+principal)*discountRec(1);
            end
            
            % Fixed pay
            if ~isnan(obj.fixedPay) 
                payCash(:) = obj.fixedPay*principal/(100*obj.payments(2));
                payCash(end) = payCash(end)+principal;
                payCash = payCash.*discountPay;
                payBond = sum(payCash);
            else
                % Floating pay
                payRates = obj.payCurve(payDates-datenum(settleDate,'dd/mm/yyyy')+1)/obj.payments(2);
                payCash(:) = (((payRates+obj.basisPoints(2))/100).* principal);
                payCash(end) = payCash(end) + principal;
                payCash = payCash.*discountPay;
                payBond = ((payRates(1)/100)*principal+principal)*discountPay(1);
            end

            obj.swapValue = sum(receiveBond)-sum(payBond);

            % Joining the dates of the paying and receiving cash flows            
            totalDates = sort(unique([payDates receiveDates]));
            totalCash = zeros(1,length(totalDates));

            % Calculating the cumulative total cash received/paid 
            if totalDates(1) == payDates(1)
                totalCash(1) = totalCash(1)-payCash(1);
            end
            
            if totalDates(1) == receiveDates(1)
                totalCash(1) = totalCash(1) + receiveCash(1);
            end

            for i = 2:length(totalDates)
                totalCash(i) = totalCash(i-1);
 
                for j = 1:length(payDates)
                   if payDates(j) == totalDates(i)
                       totalCash(i) = totalCash(i) - payCash(j);
                   end
                end

                for k = 1:length(receiveDates)
                   if receiveDates(k) == totalDates(i)
                       totalCash(i) = totalCash(i) + receiveCash(k);
                   end
                end
            end
            bar(receiveDates,receiveCash)
            hold on
            bar(payDates,-payCash)
            grid on
            hold on
            plot(totalDates, totalCash,'k-','LineWidth',2)
            datetick('x','dd/mm/yyyy')
        end
        
%         function obj = valueSwap(obj)
%             % VALUING THE SWAP CONTRACT
%             alpha = obj.payments.^(-1)
%             obj.payments;
%             Settle = datenum(obj.settleDate);
%             endDatesPay = {};
%             endDatesReceive = {};
%             
%             datesPay = datenum(obj.payPortfolio.maturity,'dd/mm/yyyy');
%             payDates = datenum(obj.startDate)+1:round(365/obj.payments(1)):max(datesPay);
%             size(payDates)
%             y = obj.payCurve/100;
%             size(y)
%             ratePay = y(payDates-datenum(obj.startDate));
%             size(ratePay)
%             
%             for i = 1:length(ratePay)
%                 endDatesPay{i} = datestr(payDates(i),'dd-mmm-yyyy');
%             end
%             
%             datesReceive = datenum(obj.payPortfolio.maturity,'dd/mm/yyyy');
%             receiveDates = datenum(obj.startDate)+1:round(365/obj.payments(2)):max(datesReceive);
%             yy = obj.receiveCurve/100;
%             rateReceive = yy(receiveDates-datenum(obj.startDate)+1)
%             
%             for i = 1:length(rateReceive)
%                 endDatesReceive{i} = datestr(receiveDates(i),'dd-mmm-yyyy');
%             end
%             
%             RateSpecPay = intenvset('Rates', ratePay','StartDates',Settle, 'EndDates', endDatesPay')
%             RateSpecReceive = intenvset('Rates', rateReceive','StartDates',Settle, 'EndDates',endDatesReceive')
%             
%             Maturity = datenum(obj.endDate)
%             Spread = obj.basisPoints
% 
%             LegType = [obj.pay obj.receive];
%             
%             LegRate = [0.05 Spread(1)]; 
%             
%             if sum(LegType) == 1
%                 if(LegType(1) == 0)
%                     [obj.price, obj.swapRate, obj.AI, obj.receiveCF, obj.receiveDates, obj.payCF, obj.payDates] = ...
%                         swapbyzero(RateSpecPay, LegRate, Settle, Maturity,'LegType',LegType,...
%                         'AdjustCashFlowsBasis',true,...
%                         'BusinessDayConvention','modifiedfollow')         
%                 elseif LegType(2) == 0
%                     [obj.price, obj.swapRate, obj.AI, obj.receiveCF, obj.receiveDates, obj.payCF, obj.payDates] = ...
%                         swapbyzero(RateSpecPay, LegRate, Settle, Maturity,'LegType',LegType,...
%                         'AdjustCashFlowsBasis',true,...
%                         'BusinessDayConvention','modifiedfollow')    
%                 end
%             elseif (sum(LegType) == 0)
% %                 Rates = [ratePay' rateReceive'];
% %                 Settle = [Settle Settle];
% %                     
% %                 RateSpecPay = intenvset('Rates', ratePay','StartDates',Settle, 'EndDates', endDatesPay')
% %                 RateSpecReceive = intenvset('Rates', rateReceive','StartDates',Settle, 'EndDates',endDatesReceive')
%                 RateSpec = [RateSpecPay RateSpecReceive];
%                 [obj.price, ~ ,obj.AI,obj.receiveCF, obj.receiveDates, obj.payCF, obj.payDates] = swapbyzero(RateSpec, LegRate, Settle, Maturity, 'Principal', obj.principal);
%             end
%         end
        
%         function plotCashFlow(obj)
%            % PLOTTING THE CASH FLOW OF THE SWAP AGREEMENT
%            grid on
%            bar(obj.receiveCF)
%            hold on
%            bar(-obj.payCF)
%            hold off
%            datetick('x','dd/mm/yyyy')
% %            xlim([min(min(obj.payDates),min(obj.receiveDates)) max(max(obj.payDates),max(obj.receiveDates))])
%         end
    end
end

