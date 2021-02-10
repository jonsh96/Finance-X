% About
% This code reads the data from www.bonds.is and sorts all the necessary 
% information about treasury bonds into either an indexed or a non-indexed 
% portfolio object (see portfolio.m) which is used later on. 
%
% This code might run into some errors which are related to internet
% connection, e.g. "Could not access server. Host not found: www.bonds.is." or
% "The connection to URL 
% 'http://www.bonds.is/api/market/LoadIndexedDetail?orderbookId=131839&lang=en' 
% timed out after 30 seconds. Set options.Timeout to a higher value."
% Running the code again usually fixes these  problems.
 
% WEB SCRAPING DATA FROM BONDS.IS (NOW LANAMAL.IS)
IndexedUrl = 'http://www.lanamal.is/api/market/LoadIndexed?lang=en&nonIndexed=false&nOrderbookId=-1';
NonIndexedUrl = 'http://www.lanamal.is/api/market/LoadIndexed?lang=en&nonIndexed=true&nOrderbookId=-1';
IndexedBonds = webread(IndexedUrl);

NonIndexedBonds = webread(NonIndexedUrl);
options = weboptions('Timeout', 30);

% CREATING A PORTFOLIO OF INDEXED BONDS
for i = 1:length(IndexedBonds)
    url = "http://www.lanamal.is/api/market/LoadIndexedDetail?orderbookId=" + IndexedBonds(i).orderbookId + "&lang=en";
    % Creating a bond object (see bond.m) to sort relevant information from
    % the market overview and attribute tables on
    % www.bonds.is/market-overview
    tempBond = bond(IndexedBonds(i), webread(url,options));
    if i == 1
        IndexedPortfolio = portfolio(tempBond);
    elseif tempBond.ticker([1:3])== 'RIK'
        % The indexed portfolio should only consist of treasury bonds
        IndexedPortfolio = IndexedPortfolio.addToPortfolio(tempBond);
    end
    clear url;
end

% CREATING A PORTFOLIO OF NONINDEXED BONDS
for i = 1:length(NonIndexedBonds)
    url = "http://www.lanamal.is/api/market/LoadIndexedDetail?orderbookId=" + NonIndexedBonds(i).orderbookId + "&lang=en";
    % Creating a bond object (see bond.m) to sort relevant information from
    % the market overview and attribute tables on
    % www.bonds.is/market-overview
    tempBond = bond(NonIndexedBonds(i), webread(url,options));
    if i == 1
        NonIndexedPortfolio = portfolio(tempBond);
    else
        NonIndexedPortfolio = NonIndexedPortfolio.addToPortfolio(tempBond);
    end
    clear url;
end


