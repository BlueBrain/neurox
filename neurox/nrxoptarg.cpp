#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "neurox/nrxoptarg.h"
#include "neurox/datatypes.h"
#include "coreneuron/nrnoc/membdef.h"

#include "tclap/CmdLine.h"

static parseCommandLine(GlobalInfo * globalInfo, char ** argv, int argc)
{
    //read command line arguments (via tclap)
    try {
        //message printed lat in help text, text delimiter, version
        TCLAP::CmdLine cmd("Welcome to NeuroX TODO whwat to write here", ' ', "0.1");

        //add all parameters
        TCLAP::ValueArg<double> tstart("s","tstart","Execution start time (msecs). The default value is 0",false, 0 ,"double");
        cmd.add(tstart);
        TCLAP::ValueArg<double> tstop("e","tstop","Execution stop time (msecs). The default value is 100",false, 100 ,"double");
        cmd.add(tstop);
        TCLAP::ValueArg<double> dt("-t","dt","Execution time step (msecs). The default value is 0.025",false, DEF_dt ,"double");
        cmd.add(dt);
        TCLAP::ValueArg<double> dt_io("-i","dt_io","I/O time step (msecs). The default value is 0.1",false, 0.1 ,"double");
        cmd.add(dt_io);
        TCLAP::ValueArg<double> celsius("-l","celsius","System temperatura (celsius degrees). The default value is 34",false, 34.0 ,"double");
        cmd.add(celsius);
        TCLAP::ValueArg<std::string> patternStimFile("-p","pattern","Apply patternstim with the spike file. No default value",false, "" ,"string");
        cmd.add(patternStimFile);
        TCLAP::ValueArg<int> prcellgid("-g","prcellgid","Output prcellstate information for given gid. The default value is -1",false, -1 ,"int");
        cmd.add(prcellgid);
        TCLAP::ValueArg<std::string> inputPath("d","inputpath","Path to input files directory",true,"./input","string");
        cmd.add(inputPath);
        TCLAP::ValueArg<std::string> outputPath("o","outputpath","Path to output directory. The default value is ./output",false,"./output","string");
        cmd.add(outputPath);
        TCLAP::ValueArg<double> forwardskip("-k","forwardskip","Set forwardskip time step (msecs). The default value is 0",false, 0 ,"double");
        cmd.add(forwardskip);

        //parse command line arguments
        cmd.parse( argc, argv );

        //copy parameters to GlobalInfo
        memcpy(globalInfo->inputPath, inputPath.getValue().c_str(), 2048);
        memcpy(globalInfo->outputPath, outputPath.getValue().c_str(), 2048);
        memcpy(globalInfo->patternStimFile, patternStimFile.getValue().c_str(), 2048);
        globalInfo->tstart 	= tstart.getValue();
        globalInfo->tstop	= tstop.getValue();
        globalInfo->dt		= dt.getValue();
        globalInfo->dt_io	= dt_io.getValue();
        globalInfo->celsius	= celsius.getValue();
        globalInfo->prcellgid = prcellgid.getValue();
        globalInfo->forwardSkip = forwardskip.getValue();
        globalInfo->voltage = DEF_vrest;
        globalInfo->maxdelay = 10;
        globalInfo->mindelay = globalInfo->dt;
        globalInfo->secondorder = DEF_secondorder;
        globalInfo->rev_dt = 1/globalInfo->dt;
        globalInfo->celsius = DEF_celsius;
    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
        { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}
