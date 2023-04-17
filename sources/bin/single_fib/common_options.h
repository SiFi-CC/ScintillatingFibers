#ifndef COMMON_OPTIONS_H
#define COMMON_OPTIONS_H

#include <TString.h>
#include <CmdLineConfig.hh>
#include <iostream>
#include <sys/stat.h>

int parse_common_options(int argc, char** argv, TString& outdir,
                         TString& dbase, Int_t& seriesno)
{
    TString path = std::string("./");

    new CmdLineOption("Output directory", "-out", "Output directory (string), default: ./", path);
    new CmdLineOption("Database", "-db", "Database name (string), default: ScintFibRes.db", "ScintFibRes.db");
    new CmdLineArg("SeriesNo", "series number", CmdLineArg::kInt);

    CmdLineConfig::instance()->ReadCmdLine(argc, argv);

    outdir   = CmdLineOption::GetStringValue("Output directory");
    dbase    = CmdLineOption::GetStringValue("Database");
    seriesno = CmdLineArg::GetIntValue("SeriesNo");

    if (!gSystem->ChangeDirectory(outdir))
    {
        std::cout << "Creating new directory... " << std::endl;
        std::cout << outdir << std::endl;
        int stat = mkdir(outdir, 0777);
        if (stat == -1)
        {
            std::cerr << "##### Error in " << __func__ << "! Unable to create new direcotry!" << std::endl;
            return 1;
        }
    }

    return 0;
}

#endif /* COMMON_OPTIONS_H */
