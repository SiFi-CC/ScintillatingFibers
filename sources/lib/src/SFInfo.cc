// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              SFInfo.cc                *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************

#include "SFInfo.hh"

//------------------------------------------------------------------
// constants
static const char*  gPath  = getenv("SFDATA"); // path to the experimental data and data base
//------------------------------------------------------------------

//------------------------------------------------------------------
SFInfo::SFInfo(int seriesNo) : fSeriesNo(seriesNo),
                               fNpoints(-1),
                               fFiberLength(-1),
                               fOvervoltage(-1),
                               fFiber("dummy"),
                               fSource("dummy"),
                               fCollimator("dummy"),
                               fDesc("dummy"),
                               fTestBench("dummy"),
                               fSiPM("dummy"),
                               fCoupling("dummy"),
                               fLogFile("dummy"),
                               fTempFile("dummy"),
                               fDAQ("dummy"),
                               fDB(nullptr)
{
    TString       query;
    sqlite3_stmt* statement;
    int           status;
    
    //----- opening data base
    TString db_name = std::string(gPath) + "/DB/ScintFib_2.db";
    std::cout << "Opening data base: " << db_name << std::endl;
    status = sqlite3_open(db_name, &fDB);
    
    if (status != 0)
    {
        std::cerr << "##### Error in SFInfo constructor" << std::endl;
        std::cerr << "Could not access data base!" << std::endl;
        exit(EXIT_FAILURE);
    }
    //----- 
    
    //-----Checking if series number is valid
    int maxSeries;
    query  = "SELECT COUNT(*) FROM SERIES";
    status = sqlite3_prepare_v2(fDB, query, -1, &statement, nullptr);

    //SFTools::CheckDBStatus(status, fDB);

    while ((status = sqlite3_step(statement)) == SQLITE_ROW)
    {
        maxSeries = sqlite3_column_int(statement, 0);
    }

    //SFTools::CheckDBStatus(status, fDB);

    sqlite3_finalize(statement);

    if (fSeriesNo < 1 || fSeriesNo > maxSeries)
    {
        std::cerr << "##### Error in SFInfo constructor! Series number out of range!" << std::endl;
        exit(EXIT_FAILURE);
    }
    //-----
    
    //----- Setting series attributes
    ///- fiber type
    ///- liber length [mm]
    ///- radioactive source type
    ///- test bench type
    ///- collimator type
    ///- SiPM type
    ///- overvoltage [V]
    ///- coupling type
    ///- number of measurements in the series
    ///- name of the measurement log file
    ///- name of the temperature log file
    ///- description of the series
    ///- analysis group number
    
    query  = Form("SELECT FIBER, FIBER_LENGTH, SOURCE, TEST_BENCH, COLLIMATOR, SIPM, OVERVOLTAGE, "
                 "COUPLING, NO_MEASUREMENTS, LOG_FILE, TEMP_FILE, DESCRIPTION, DAQ FROM "
                 "SERIES WHERE SERIES_ID = %i", fSeriesNo);
    status = sqlite3_prepare_v2(fDB, query, -1, &statement, nullptr);

    //SFTools::CheckDBStatus(status, fDB);

    while ((status = sqlite3_step(statement)) == SQLITE_ROW)
    {
        const unsigned char* fiber       = sqlite3_column_text(statement, 0);
        const unsigned char* source      = sqlite3_column_text(statement, 2);
        const unsigned char* test_bench  = sqlite3_column_text(statement, 3);
        const unsigned char* collimator  = sqlite3_column_text(statement, 4);
        const unsigned char* sipm        = sqlite3_column_text(statement, 5);
        const unsigned char* coupling    = sqlite3_column_text(statement, 7);
        const unsigned char* logfile     = sqlite3_column_text(statement, 9);
        const unsigned char* tempfile    = sqlite3_column_text(statement, 10);
        const unsigned char* description = sqlite3_column_text(statement, 11);
        const unsigned char* daq         = sqlite3_column_text(statement, 12);
        fNpoints                         = sqlite3_column_int(statement, 8);
        fFiberLength                     = sqlite3_column_double(statement, 1);
        fOvervoltage                     = sqlite3_column_double(statement, 6);
        fFiber                           = std::string(reinterpret_cast<const char*>(fiber));
        fSource                          = std::string(reinterpret_cast<const char*>(source));
        fDesc                            = std::string(reinterpret_cast<const char*>(description));
        fCollimator                      = std::string(reinterpret_cast<const char*>(collimator));
        fTestBench                       = std::string(reinterpret_cast<const char*>(test_bench));
        fSiPM                            = std::string(reinterpret_cast<const char*>(sipm));
        fCoupling                        = std::string(reinterpret_cast<const char*>(coupling));
        fLogFile                         = std::string(reinterpret_cast<const char*>(logfile));
        fTempFile                        = std::string(reinterpret_cast<const char*>(tempfile));
        fDAQ                             = std::string(reinterpret_cast<const char*>(daq));
    }

    //SFTools::CheckDBStatus(status, fDB);

    sqlite3_finalize(statement);
    //-----

    //----- Setting measurements attributes
    ///- list of measurements names
    ///- list of measurements duration times
    ///- list of source positions
    ///- list of measurements starting times
    ///- list of measurements stopping times
    ///- list of measurements IDs
    
    query  = Form("SELECT MEASUREMENT_NAME, DURATION_TIME, SOURCE_POSITION, START_TIME, STOP_TIME, "
                 "MEASUREMENT_ID FROM MEASUREMENT WHERE SERIES_ID = %i",
                 fSeriesNo);
    status = sqlite3_prepare_v2(fDB, query, -1, &statement, nullptr);

    //SFTools::CheckDBStatus(status, fDB);

    while ((status = sqlite3_step(statement)) == SQLITE_ROW)
    {
        const unsigned char* name = sqlite3_column_text(statement, 0);
        fNames.push_back(std::string(reinterpret_cast<const char*>(name)));
        fTimes.push_back(sqlite3_column_int(statement, 1));
        fPositions.push_back(sqlite3_column_double(statement, 2));
        fStart.push_back(sqlite3_column_int(statement, 3));
        fStop.push_back(sqlite3_column_int(statement, 4));
        fMeasureID.push_back(sqlite3_column_int(statement, 5));
    }

    //SFTools::CheckDBStatus(status, fDB);

    sqlite3_finalize(statement);
}
//------------------------------------------------------------------
SFInfo::SFInfo() : fSeriesNo(-1),
                   fNpoints(-1),
                   fFiberLength(-1),
                   fOvervoltage(-1),
                   fFiber("dummy"),
                   fSource("dummy"),
                   fCollimator("dummy"),
                   fDesc("dummy"),
                   fTestBench("dummy"),
                   fSiPM("dummy"),
                   fCoupling("dummy"),
                   fLogFile("dummy"),
                   fTempFile("dummy"),
                   fDAQ("dummy"),
                   fDB(nullptr)
{
    std::cout << "Warning in SFInfo::SFInfo()! You are using empty constructor!" << std::endl;
    std::cout << "Set measurement properties via available setters!" << std::endl;
}
//------------------------------------------------------------------
SFInfo::~SFInfo()
{
    int status = sqlite3_close(fDB);
    if (status != 0) 
        std::cerr << "In SFInfo destructor. Data base corrupted!" << std::endl;
}
//------------------------------------------------------------------
/// Prints details of currently analyzed experimental series.
void SFInfo::Print(void)
{
    std::cout << "\n\n------------------------------------------------" << std::endl;
    std::cout << "This is Print() for SFInfo class object" << std::endl;
    std::cout << "Number of the experimental series: " << fSeriesNo << std::endl;
    std::cout << fDesc << std::endl;
    std::cout << "Number of measurements in this series: " << fNpoints << std::endl;
    std::cout << "Fiber: " << fFiber << std::endl;
    std::cout << "Fiber length: " << fFiberLength << " mm" << std::endl;
    std::cout << "Collimator: " << fCollimator << std::endl;
    std::cout << "Test bench: " << fTestBench << std::endl;
    std::cout << "Coupling: " << fCoupling << std::endl;
    std::cout << "Radioactive source: " << fSource << std::endl;
    std::cout << "SiPM: " << fSiPM << std::endl;
    std::cout << "Overvoltage: " << fOvervoltage << " V" << std::endl;
    std::cout << "DAQ: " << fDAQ << std::endl; 
    std::cout << "Logfile: " << fLogFile << std::endl;
    std::cout << "Temperature logfile: " << fTempFile << std::endl;
    std::cout << "List of measurements in this series:" << std::endl;
    for (int i = 0; i < fNpoints; i++)
    {
        std::cout << std::setw(30);
        std::cout << fNames[i] << "\t\t" << Form("%.1f mm", fPositions[i]) << "\t\t"
                  << Form("%i s", fTimes[i]) << "\t\t" << fStart[i] << "\t\t" << fStop[i] << "\t\t"
                  << fMeasureID[i] << std::endl;
    }
    std::cout << "\n" << std::endl;
}
//------------------------------------------------------------------
