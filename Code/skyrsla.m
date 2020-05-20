%% ICELANDIC BOND MARKET 
%% Lesa inn gögn af bonds.is
createPortfolio;    
% Núna eru til 2 portfolio object, NonIndexedPortfolio og IndexedPortfolio

%% Þetta eru öll curve og fitting methods í boði
curves = ['Yield', 'Zero rates', 'Forward rates', 'Discount rates', 'Swap rates'];
fitMethods = ['Bootstrapping','Nelson-Siegel','Nelson-Siegel-Svensson','Polynomial',... 
    'Lagrange interpolation','Spline','Cubic spline','Constrained cubic spline'];

PolynomialDegree = 2;       % Fyrir Polynomial
SmoothingFactor = 0.75;     % Fyrir Constrained cubic spline [0,1]

%% Til að fitta curve með ákveðinni aðferð + plotta
NonIndexedPortfolio = NonIndexedPortfolio.fitMethod('Forward rates', 'Lagrange interpolation',PolynomialDegree, SmoothingFactor);
IndexedPortfolio = IndexedPortfolio.fitMethod(curves(2),fitMethods(6),PolynomialDegree, SmoothingFactor);
% Núna eru þessi fittuðu curve geymd sem .currentCurve í báðum portfolioum

%%  Til að plotta
plot(datesNI, NonIndexedPortfolio.currentCurve)
% plot(datesI, IndexedPortfolio.currentCurve) % Sambærilegt fyrir indexed

%% Líka hægt að plotta gagnapunktana með:
% NonIndexedPortfolio.zeroCurve;
% NonIndexedPortfolio.yieldCurve;
 NonIndexedPortfolio.forwardCurve;
% NonIndexedPortfolio.discountCurve;
% NonIndexedPortfolio.swapCurve;

%% Til að plotta gagnapunkta OG fitted curve --> Sjá bondGUI
bondGUI

%% INTEREST RATE SIMULATIONS
% Model í boði: 'Simple', 'Brownian', 'Vasicek'
% Sjá help interestRate fyrir útskýringar
model = ['Simple', 'Brownian', 'Vasicek'];
initialRate = 0.05;
stepSize = 1/250;   
volatility = 0.03;
speedOfReversion = 0.25;    % Notað fyrir alpha í brownian og kappa í vasicek
longTermMeanLevel = 0.035;  % Notað fyrir theta í vasicek
maturity = 5;
nrOfSimulations = 1000;

IR1 = interestRate(model(1), initialRate, stepSize, volatility, speedOfReversion, longTermMeanLevel, maturity, nrOfSimulations);
IR2 = interestRate(model(2), initialRate, stepSize, volatility, speedOfReversion, longTermMeanLevel, maturity, nrOfSimulations);
IR3 = interestRate(model(3), initialRate, stepSize, volatility, speedOfReversion, longTermMeanLevel, maturity, nrOfSimulations);

%% Simulatea vextina
IR1.plotSimulations
IR2.plotSimulations
IR3.plotSimulations

%% Sýna histogram fyrir interest rate modelið
IR1.histModel
IR2.histModel
IR3.histModel

%% Sýna hvernig MLE batnar með fjölda ítrana
% Hækka maturity í IR3 til að bæta niðurstöðuna
IR3.estimationImprovement

%% Sýna hvernig OLS batnar með fjölda ítrana
% Hækka maturity í IR3 til að bæta niðurstöðuna
IR3.estimationImprovementOLS

%% Hægt að skoða þetta interactive með interestRateGUI
interestRateGUI

%% OPTION PRICER
StrikePrice = 1.0;
OptionMaturity = 4;     % Verður að vera minna en IR.maturity

% Búa til pricing model með öllum interest rate modelunum
PM1 = pricingModel(IR1, StrikePrice, OptionMaturity)
PM2 = pricingModel(IR2, StrikePrice, OptionMaturity)
PM3 = pricingModel(IR3, StrikePrice, OptionMaturity)

%% Plot zero coupon bonds með Simple
PM1.plotBonds

%% Plot zero coupon bonds með Brownian
PM2.plotBonds

%% Plot zero coupon bonds með Vasicek
PM3.plotBonds


