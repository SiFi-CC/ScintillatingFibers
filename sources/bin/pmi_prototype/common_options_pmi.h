#ifndef COMMON_OPTIONS_PMI_H
#define COMMON_OPTIONS_PMI_H

#include <TString.h>
#include <CmdLineConfig.hh>
#include <iostream>
#include <sys/stat.h>

int parse_common_options(int argc, char** argv, TString &data_path,
                         TString &outdir, TString &dbase, 
                         int &m, int &l, int &f)
{
    new CmdLineOption("Module", "-m", "Module number (int), default: 0", 0);
    new CmdLineOption("Layer", "-l", "Layer number (int), default: 0", 0);
    new CmdLineOption("Fiber", "-f", "Fiber number (int), default: 0", 0);
    new CmdLineOption("Output directory", "-out", "Output directory (string), default: ./", "./");
    new CmdLineOption("Database", "-db", "Database name (string), default: PMIRes.db", "PMIRes.db");
    new CmdLineArg("Log path", "Log path", CmdLineArg::kString);
    
    CmdLineConfig::instance()->ReadCmdLine(argc, argv);
    
    m = CmdLineOption::GetIntValue("Module");
    l = CmdLineOption::GetIntValue("Layer");
    f = CmdLineOption::GetIntValue("Fiber");
    outdir    = CmdLineOption::GetStringValue("Output directory");
    dbase     = CmdLineOption::GetStringValue("Database");
    data_path = CmdLineArg::GetStringValue("Log path");
    
    if (!gSystem->ChangeDirectory(outdir))
    {
        std::cout << "Creating new directory... " << std::endl;
        std::cout << outdir << std::endl;
        int stat = mkdir(outdir, 0777);
        if (stat == -1)
        {
            std::cerr << "##### Error in " << __func__ << " Unable to create new direcotry!" << std::endl;
            std::cerr << outdir << std::endl;
            return 1;
        }
    }
    
    return 0;
}

#endif /* COMMON_OPTIONS_PMI_H */
