regionalpred = cell(12,1); 
idtranlag = cell(12,1);

NEpredictors = importdata('Predictors2022Q2NE.csv');
idtranlag{1} = NEpredictors.data(1:2,:)';
NEdata = NEpredictors.data(3:end,:);
NEdata(NEdata == -99999999) = NaN;
regionalpred{1} = NEdata;

NWpredictors = importdata('Predictors2022Q2NW.csv');
idtranlag{2} = NWpredictors.data(1:2,:)';
NWdata = NWpredictors.data(3:end,:);
NWdata(NWdata == -99999999) = NaN;
regionalpred{2} = NWdata;

Yorkpredictors = importdata('Predictors2022Q2York.csv');
idtranlag{3} = Yorkpredictors.data(1:2,:)';
Yorkdata = Yorkpredictors.data(3:end,:);
Yorkdata(Yorkdata == -99999999) = NaN;
regionalpred{3} = Yorkdata;

EMpredictors = importdata('Predictors2022Q2EM.csv');
idtranlag{4} = EMpredictors.data(1:2,:)';
EMdata = EMpredictors.data(3:end,:);
EMdata(EMdata == -99999999) = NaN;
regionalpred{4} = EMdata;

WMpredictors = importdata('Predictors2022Q2WM.csv');
idtranlag{5} = WMpredictors.data(1:2,:)';
WMdata = WMpredictors.data(3:end,:);
WMdata(WMdata == -99999999) = NaN;
regionalpred{5} = WMdata;

Eastpredictors = importdata('Predictors2022Q2East.csv');
idtranlag{6} = Eastpredictors.data(1:2,:)';
Eastdata = Eastpredictors.data(3:end,:);
Eastdata(Eastdata == -99999999) = NaN;
regionalpred{6} = Eastdata;

Londonpredictors = importdata('Predictors2022Q2London.csv');
idtranlag{7} = Londonpredictors.data(1:2,:)';
Londondata = Londonpredictors.data(3:end,:);
Londondata(Londondata == -99999999) = NaN;
regionalpred{7} = Londondata;

SEpredictors = importdata('Predictors2022Q2SE.csv');
idtranlag{8} = SEpredictors.data(1:2,:)';
SEdata = SEpredictors.data(3:end,:);
SEdata(SEdata == -99999999) = NaN;
regionalpred{8} = SEdata;

SWpredictors = importdata('Predictors2022Q2SW.csv');
idtranlag{9} = SWpredictors.data(1:2,:)';
SWdata = SWpredictors.data(3:end,:);
SWdata(SWdata == -99999999) = NaN;
regionalpred{9} = SWdata;

Walespredictors = importdata('Predictors2022Q2Wales.csv');
idtranlag{10} = Walespredictors.data(1:2,:)';
Walesdata = Walespredictors.data(3:end,:);
Walesdata(Walesdata == -99999999) = NaN;
regionalpred{10} = Walesdata;

Scotpredictors = importdata('Predictors2022Q2Scot.csv');
idtranlag{11} = Scotpredictors.data(1:2,:)';
Scotdata = Scotpredictors.data(3:end,:);
Scotdata(Scotdata == -99999999) = NaN;
regionalpred{11} = Scotdata;

NIpredictors = importdata('Predictors2022Q2NI.csv');
idtranlag{12} = NIpredictors.data(1:2,:)';
NIdata = NIpredictors.data(3:end,:);
NIdata(NIdata == -99999999) = NaN;
regionalpred{12} = NIdata;