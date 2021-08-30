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
/// Standard constructor.
/// \param seriesNo - number of the experimental series
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
/// Standard destructor.
SFTemperature::~SFTemperature()
{

    int status = sqlite3_close(fDB);
    if (status != 0) std::cerr << "In SFData destructor. Data base corrupted!" << std::endl;
}
//------------------------------------------------------------------
/// Opens data base containing temperature readings.
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
/// Loads temperature readings saved in the data base for requested
/// sensor and measurement. Loaded values are saved in fTime and 
/// fTemperature vectors.
/// \param sensor - sensor ID
/// \param name - name of the measurement; if value "all" is given
/// all temperature values for this series will be loaded
bool SFTemperature::LoadFromDB(TString sensor, TString name)
{

    int           status = 0;
    sqlite3_stmt* statement;

    if (!fTime.empty()) fTime.clear();
    if (!fTemperatures.empty()) fTemperatures.clear();

    TString query;

    if (name == "all")
        query = Form("SELECT TIME, TEMPERATURE FROM TEMPERATURES WHERE SERIES_ID "
                     "= %i AND SENSOR_ID = '%s'", fSeriesNo, sensor.Data());
    else
        query = Form("SELECT TIME, TEMPERATURE FROM TEMPERATURES WHERE SERIES_ID = %i AND "
                     "SENSOR_ID = '%s' AND MEASUREMENT_NAME = '%s'",
                     fSeriesNo, sensor.Data(), name.Data());

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
/// Calculates average temperature during chosen measurement for chosen
/// temperature sensor. 
/// \param sensor - sensorID
/// \param name - measurement name
SFResults* SFTemperature::CalcAverageTempMeasure(TString sensor, TString name)
{

    LoadFromDB(sensor, name);
    double mean   = SFTools::GetMean(fTemperatures);
    double stdErr = SFTools::GetStandardErr(fTemperatures);
    std::cout << mean << "\t" << stdErr << std::endl;

    SFResults* result = new SFResults(Form("TemperatureResultsM_S%i_%s_%s",
                                      fSeriesNo, sensor.Data(), name.Data()));

    result->AddResult(SFResultTypeNum::kTemp, mean, stdErr);
    
    return result;
}
//------------------------------------------------------------------
/// Calculates average temperature during entire measurement series for 
/// chosen temperature sensor.
/// \param sensor - sensor ID
bool SFTemperature::CalcAverageTempSeries(TString sensor)
{

    LoadFromDB(sensor);
    double mean  = SFTools::GetMean(fTemperatures);
    double stdDev = SFTools::GetStandardErr(fTemperatures);

    SFResults* temp = new SFResults(Form("TemperatureResultsS_S%i_%s",
                                    fSeriesNo, sensor.Data()));
    
    temp->AddResult(SFResultTypeNum::kTemp, mean, stdDev);

    fAvTemp.insert(std::pair<TString, SFResults*>(sensor, temp));

    return true;
}
//------------------------------------------------------------------
/// Prepares plot of average temperatures during measurements for a 
/// requested sensor as a function of source position. If monitoring 
/// series is analyzed instead of source position number of measurement 
/// is used.  
/// \param sensor - sensorID
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
/// Builds plot of temperature for a requested sensor. Temperature is 
/// plotted as a function of time since measurement started.
/// \param sensor - sensor ID
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
/// Returns average temperature during the measurement series for 
/// requested sensor.
/// \param sensor - sensor ID
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
/// Returns temperature plot for requested sensor.
/// \param sensor - sensor ID
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
/// Returns average temperature plot for requested sensor.
/// \param sensor - sensor ID
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
/// Prints details of the SFTemperature class object.
void SFTemperature::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFTemperature class object" << std::endl;
    std::cout << "Experimental series number: " << fSeriesNo << std::endl;
    std::cout << "-------------------------------------------\n" << std::endl;
    return;
}
//------------------------------------------------------------------
