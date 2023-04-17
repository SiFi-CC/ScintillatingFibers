#ifndef COMMON_OPTIONS_KRK_H
#define COMMON_OPTIONS_KRK_H

#include <TString.h>
#include <CmdLineConfig.hh>
#include <iostream>
#include <sys/stat.h>

int parse_common_options(int argc, char** argv, TString &log_path,
                         TString &data_path, TString &outdir, TString &dbase, 
                         int &m, int &l, int &f)
{
    new CmdLineOption("Module", "-m", "Module number (int), default: 0", 0);
    new CmdLineOption("Layer", "-l", "Layer number (int), default: 0", 0);
    new CmdLineOption("Fiber", "-f", "Fiber number (int), default: 0", 0);
    new CmdLineOption("Output directory", "-out", "Output directory (string), default: ./", "./");
    new CmdLineOption("Database", "-db", "Database name (string), default: KRKRes.db", "KRKRes.db");
    new CmdLineArg("Log path", "Log path", CmdLineArg::kString);

    CmdLineConfig::instance()->ReadCmdLine(argc, argv);
    
    m = CmdLineOption::GetIntValue("Module");
    l = CmdLineOption::GetIntValue("Layer");
    f = CmdLineOption::GetIntValue("Fiber");
    outdir   = CmdLineOption::GetStringValue("Output directory");
    dbase    = CmdLineOption::GetStringValue("Database");
    log_path = CmdLineArg::GetStringValue("Log path");
    
    //------ parsing logname
    std::string logname_str = std::string(log_path);
    int         nletters    = logname_str.length();
    char        letters[nletters];
    strcpy(letters, logname_str.c_str());
  
    int istop = -1;
  
    for (int i = nletters; i > 0; i--)
    {
        if (letters[i] == '/')
        {
            istop = i;
            break;
        }
    }

    int istart = 0;
   
    if (istop == -1)
    {
        std::cerr << "##### Error in " << __func__ << std::endl;
        std::cerr << "Cannot interpret log name!" << std::endl;
        return 1;
    }

    std::string data_path_str = std::string(&letters[istart], &letters[istop + 1]) + "../";
    std::string file_name_str = std::string(&letters[istop+1], &letters[nletters-4]);
    std::cout << "data path: " << data_path_str << std::endl;
    std::cout << "file name: " << file_name_str << std::endl;
    
    data_path = data_path_str;
    
    //----- creating output directory if it doesn't exist
    if (!gSystem->ChangeDirectory(outdir))
    {
        std::cout << "Creating new directory... " << std::endl;
        std::cout << outdir << std::endl;
        int stat = mkdir(outdir, 0777);
        if (stat == -1)
        {
            std::cerr << "##### Error in " << __func__ << "! Unable to create new direcotry!" << std::endl;
            std::cerr << outdir << std::endl;
            return 1;
        }
    }
    
    return 0;
}

#endif /* COMMON_OPTIONS_KRK_H */
