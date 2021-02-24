// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFTemperature.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFTemperature.hh"

ClassImp(SFTemperature);

//------------------------------------------------------------------
SFTemperature::SFTemperature(int seriesNo) : fSeriesNo(seriesNo),
                                             fData(nullptr)
{

    try
    {
        fData = new SFData(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFTemperature constructor!";
    }

    TString tempFile = fData->GetTempFile();

    if (tempFile == "-")
    {
        std::cerr << "This series doesn't have temeperature measurement!" << std::endl;
        throw "##### Exception in SFTemperature!";
    }

    OpenDataBase();
}
//------------------------------------------------------------------
SFTemperature::~SFTemperature()
{

    int status = sqlite3_close(fDB);
    if (status != 0) std::cerr << "In SFData destructor. Data base corrupted!" << std::endl;
}
//------------------------------------------------------------------
bool SFTemperature::OpenDataBase(void)
{

    const char* path   = getenv("SFDATA");
    TString     DBname = std::string(path) + "/DB/ScintFib_2.db";
    int         status = sqlite3_open(DBname, &fDB);

    if (status != 0)
    {
        std::cerr << "##### Error in SFTemperature::OpenDataBase()!" << std::endl;
        std::cerr << "Could not access data base!" << std::endl;
        return false;
    }

    return true;
}
//------------------------------------------------------------------
bool SFTemperature::LoadFromDB(TString sensorID, TString name)
{

    int           status = 0;
    sqlite3_stmt* statement;

    if (!fTime.empty()) fTime.clear();
    if (!fTemperatures.empty()) fTemperatures.clear();

    TString query;

    if (name == "none")
        query = Form("SELECT TIME, TEMPERATURE FROM TEMPERATURES WHERE SERIES_ID "
                     "= %i AND SENSOR_ID = '%s'", fSeriesNo, sensorID.Data());
    else
        query = Form("SELECT TIME, TEMPERATURE FROM TEMPERATURES WHERE SERIES_ID = %i AND "
                     "SENSOR_ID = '%s' AND MEASUREMENT_NAME = '%s'",
                     fSeriesNo, sensorID.Data(), name.Data());

    status = sqlite3_prepare_v2(fDB, query, -1, &statement, nullptr);
    SFTools::CheckDBStatus(status, fDB);

    while ((status = sqlite3_step(statement)) == SQLITE_ROW)
    {
        fTime.push_back(sqlite3_column_int(statement, 0));
        fTemperatures.push_back(sqlite3_column_double(statement, 1));
    }

    SFTools::CheckDBStatus(status, fDB);
    sqlite3_finalize(statement);

    return true;
}
//------------------------------------------------------------------
SFResults* SFTemperature::CalcAverageTempMeasure(TString sensor, TString name)
{

    LoadFromDB(sensor, name);
    double mean   = SFTools::GetMean(fTemperatures);
    double stdErr = SFTools::GetStandardErr(fTemperatures);
    std::cout << mean << "\t" << stdErr << std::endl;

    SFResults* result = new SFResults(Form("TemperatureResults_S%i_%s_%s",
                                      fSeriesNo, sensor.Data(), name.Data()));

    result->AddResult(SFResultTypeNum::kTemp, mean, stdErr);
    
    return result;
}
//------------------------------------------------------------------
bool SFTemperature::CalcAverageTempSeries(TString sensor)
{

    LoadFromDB(sensor);
    double mean  = SFTools::GetMean(fTemperatures);
    double stdDev = SFTools::GetStandardErr(fTemperatures);

    SFResults* temp;
    temp->AddResult(SFResultTypeNum::kTemp, mean, stdDev);

    fAvTemp.insert(std::pair<TString, SFResults*>(sensor, temp));

    return true;
}
//------------------------------------------------------------------
bool SFTemperature::BuildTempPlotAverage(TString sensor)
{

    int                  npoints    = fData->GetNpoints();
    TString              collimator = fData->GetCollimator();
    TString              testBench  = fData->GetTestBench();
    TString              desc       = fData->GetDescription();
    std::vector<TString> names      = fData->GetNames();
    std::vector<double>  pos        = fData->GetPositions();
    SFResults*           result     = nullptr;

    TGraphErrors* gr = new TGraphErrors(npoints);
    gr->SetName(sensor);
    gr->SetMarkerStyle(4);
    
    if (desc.Contains("Regular series"))
        gr->GetXaxis()->SetTitle("source position");
    else
        gr->GetXaxis()->SetTitle("number of measurement");
    
    gr->GetYaxis()->SetTitle("average temperature during measurement [deg C]");

    for (int i = 0; i < npoints; i++)
    {
        result = CalcAverageTempMeasure(sensor, names[i]);
        if (desc.Contains("Regular series"))
        {
            gr->SetPoint(i, pos[i], result->GetValue(SFResultTypeNum::kTemp));
            gr->SetPointError(i, SFTools::GetPosError(collimator, testBench),
                              result->GetUncertainty(SFResultTypeNum::kTemp));
        }
        else
        {
            gr->SetPoint(i, i, result->GetValue(SFResultTypeNum::kTemp));
            gr->SetPointError(i, 0, result->GetUncertainty(SFResultTypeNum::kTemp));
        }
    }

    fTempPlotAv.insert(std::pair<TString, TGraphErrors*>(sensor, gr));

    return true;
}
//------------------------------------------------------------------
bool SFTemperature::BuildTempPlot(TString sensor)
{

    LoadFromDB(sensor);
    int    size = fTime.size();
    double time;
    std::cout << "size: " << size << std::endl;

    TGraphErrors* gr = new TGraphErrors(size);
    gr->SetName(sensor);
    gr->SetMarkerStyle(4);
    gr->GetXaxis()->SetTitle("minutes since beginning of measurement");
    gr->GetYaxis()->SetTitle("temperature [deg C]");

    for (int i = 0; i < size; i++)
    {
        time = (fTime[i] - fTime[0]) / 60.;
        gr->SetPoint(i, time, fTemperatures[i]);
    }

    fTempPlot.insert(std::pair<TString, TGraphErrors*>(sensor, gr));

    return true;
}
//------------------------------------------------------------------
SFResults* SFTemperature::GetAverageTempSeries(TString sensor)
{

    std::map<TString, SFResults*>::iterator itr;
    itr = fAvTemp.find(sensor);

    if (itr == fAvTemp.end())
    {
        std::cerr << "##### Error in SFTemperature::GetAverageTempSeries()!" << std::endl;
        std::cerr << "Incorrect sensor ID! Please check!" << std::endl;
        std::abort();
    }

    SFResults* temp = itr->second;

    if (temp->GetValue(SFResultTypeNum::kTemp) == -1 ||
        temp->GetUncertainty(SFResultTypeNum::kTemp) == -1)
    {
        std::cerr << "##### Error in SFTemperature::GetAverageTempSeries()!" << std::endl;
        std::cerr << "Unknown sensor requested!" << std::endl;
    }

    return temp;
}
//------------------------------------------------------------------
TGraphErrors* SFTemperature::GetTempPlot(TString sensor)
{

    std::map<TString, TGraphErrors*>::iterator itr;
    itr = fTempPlot.find(sensor);

    if (itr == fTempPlot.end())
    {
        std::cerr << "##### Error in SFTemperature::GetTempPlot()!" << std::endl;
        std::cerr << "Incorrect sensor ID! Please check!" << std::endl;
        std::abort();
    }

    TGraphErrors* temp = itr->second;

    if (temp == nullptr)
    {
        std::cerr << "##### Error in SFTemperature::GetTempPlot()!" << std::endl;
        std::cerr << "Empty graph! Please check!" << std::endl;
        std::abort();
    }

    return temp;
}
//------------------------------------------------------------------
TGraphErrors* SFTemperature::GetTempPlotAverage(TString sensor)
{

    std::map<TString, TGraphErrors*>::iterator itr;
    itr = fTempPlotAv.find(sensor);

    if (itr == fTempPlotAv.end())
    {
        std::cerr << "##### Error in SFTemperature::GetTempPlotAverage()!" << std::endl;
        std::cerr << "Incorrect sensor ID! Please check!" << std::endl;
        std::abort();
    }

    TGraphErrors* temp = itr->second;

    if (temp == nullptr)
    {
        std::cerr << "##### Error in SFTemperature::GetTempPlotAverage()!" << std::endl;
        std::cerr << "Empty graph! Please check!" << std::endl;
        std::abort();
    }

    return temp;
}
//------------------------------------------------------------------
void SFTemperature::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFTemperature class object" << std::endl;
    std::cout << "Experimental series number: " << fSeriesNo << std::endl;
    std::cout << "-------------------------------------------\n" << std::endl;
    return;
}
//------------------------------------------------------------------
