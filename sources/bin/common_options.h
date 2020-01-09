#ifndef COMMON_OPTIONS_H
#define COMMON_OPTIONS_H

#include <iostream>
#include <sys/stat.h> 
#include <CmdLineConfig.hh>

int parse_common_options(int argc, char ** argv, TString & outdir, TString & dbase, Int_t & seriesno)
{
  TString path = std::string(getenv("SFPATH"));

  CmdLineOption cmd_outdir("Output directory", "-out", "Output directory (string), default: $SFPATH/results", path+"results");
  
  CmdLineOption cmd_dbase("Database", "-db", "Data base name (string), default: ScintFibRes.db", "ScintFibRes.db");

  CmdLineArg serno("SeriesNo", "series number", CmdLineArg::kInt);
  
  CmdLineConfig::instance()->ReadCmdLine(argc, argv);
  
  outdir = CmdLineOption::GetStringValue("Output directory");
  dbase = CmdLineOption::GetStringValue("Database");
  seriesno = serno.GetIntValue();

  if(!gSystem->ChangeDirectory(outdir)){
    std::cout << "Creating new directory... " << std::endl;
    std::cout << outdir << std::endl;
    int stat = mkdir(outdir, 0777);
    if(stat==-1){
      std::cerr << "##### Error in parse_common_options()! Unable to create new direcotry!" << std::endl;
      return 1;
    }
  }
  
  return 0;
}

#endif /* COMMON_OPTIONS_H */
