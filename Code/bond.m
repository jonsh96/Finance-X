classdef bond
    % About 
    % Bond object, only used to sort and clean up the data for future use.
    % Inputs consist of structs read from www.bonds.is
    % 
    % Properties accessible: 
    % ticker, issue, maturity, coupon, duration, ask, bid, price, last price, 
    % last yield, yield, interest (coupon rate)
    properties
        ticker
        issue
        maturity
        coupon
        duration
        ask
        bid
        price
        lastPrice
        lastYield
        yield
        interest
    end
    
    methods
        function obj = bond(overview, attributes)
            % INITIALISING THE BOND
            % Assigning values from www.bonds.is tables. Changing strings
            % to double values.
            obj.ticker = overview.shortName;
            obj.issue = attributes.attributes(4).value;
            obj.maturity = attributes.attributes(5).value;
            obj.duration = str2double(overview.duration([1:strfind(overview.duration,' ')-1]));
            obj.ask = str2double(overview.askPrice);
            obj.bid = str2double(overview.bidPrice);
            obj.price = str2double(overview.price);
            obj.lastPrice = str2double(overview.lastValidPrice);
            obj.lastYield = str2double(overview.lastValidYield([1:strfind(overview.lastValidYield,'%')-1]))/100;
            obj.yield = str2double(overview.yield([1:strfind(overview.yield,'%')-1]))/100;
           
            % Adjusting for some inconsistencies in the tables
            if attributes.attributes(7).name == "Coupon"
                obj.coupon = attributes.attributes(7).value;
            else
                obj.coupon = attributes.attributes(6).value;
            end
            
            for i = 6:10
                if attributes.attributes(i).name == "Interest"
                    obj.interest = str2double(attributes.attributes(i).value([1:strfind(attributes.attributes(i).value,'%')-1]))/100;
                end
            end
        end
    end
end

