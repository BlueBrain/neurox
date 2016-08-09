#include "neurox/Neurox.h"
#include "tclap/CmdLine.h"
#include "stdio.h"
#include "string.h"

#include "coreneuron/nrnoc/membdef.h"

using namespace NeuroX;
using namespace NeuroX::Input;

InputParams::InputParams ():
  //from nrnoptarg.cpp::cn_parameters():
  tstart(0), tstop(100), dt(0.025), dt_io(0.1), celsius(34),
  voltage(-65), forwardSkip(0),  prcellgid(-1)
{
    //TODO: missing some inits
    memset(patternStim,'0',2048);
    memset(inputPath,'0',2048);
    memset(outputPath,'0',2048);
}

InputParams::~InputParams(){}

InputParams::InputParams (int argc, char** argv):
    InputParams()
{
    parseCommandLine(argc, argv);
}

void InputParams::parseCommandLine(int argc, char ** argv)
{
    //read command line arguments (via tclap)
    try {
        //message printed (help text, text delimiter, version)
        TCLAP::CmdLine cmd("NeuroX Simulator", ' ', "0.1");

        //add all parameters
        TCLAP::ValueArg<double> tstart("s","tstart","Execution start time (msecs). The default value is 0",false, 0 ,"double");
        cmd.add(tstart);
        TCLAP::ValueArg<double> tstop("e","tstop","Execution stop time (msecs). The default value is 100",false, 100 ,"double");
        cmd.add(tstop);
        TCLAP::ValueArg<double> dt("t","dt","Execution time step (msecs). The default value is 0.025",false, DEF_dt ,"double");
        cmd.add(dt);
        TCLAP::ValueArg<double> dt_io("i","dt_io","I/O time step (msecs). The default value is 0.1",false, 0.1 ,"double");
        cmd.add(dt_io);
        TCLAP::ValueArg<double> celsius("l","celsius","System temperatura (celsius degrees). The default value is 34",false, 34.0 ,"double");
        cmd.add(celsius);
        TCLAP::ValueArg<std::string> patternStim("p","pattern","Apply patternstim with the spike file. No default value",false, "" ,"string");
        cmd.add(patternStim);
        TCLAP::ValueArg<int> prcellgid("g","prcellgid","Output prcellstate information for given gid. The default value is -1",false, -1 ,"int");
        cmd.add(prcellgid);
        TCLAP::ValueArg<std::string> inputPath("d","inputpath","Path to input files directory",true,"./input","string");
        cmd.add(inputPath);
        TCLAP::ValueArg<std::string> outputPath("o","outputpath","Path to output directory. The default value is ./output",false,"./output","string");
        cmd.add(outputPath);
        TCLAP::ValueArg<double> forwardskip("k","forwardskip","Set forwardskip time step (msecs). The default value is 0",false, 0 ,"double");
        cmd.add(forwardskip);

        //parse command line arguments
        cmd.parse( argc, argv );

        //copy parameters
        memcpy(this->inputPath, inputPath.getValue().c_str(), 2048);
        memcpy(this->outputPath, outputPath.getValue().c_str(), 2048);
        memcpy(this->patternStim, patternStim.getValue().c_str(), 2048);
        this->tstart = tstart.getValue();
        this->tstop = tstop.getValue();
        this->dt = dt.getValue();
        this->dt_io = dt_io.getValue();
        this->celsius = celsius.getValue();
        this->prcellgid = prcellgid.getValue();
        this->forwardSkip = forwardskip.getValue();
        this->voltage = DEF_vrest;
        this->secondorder = DEF_secondorder;
        this->rev_dt = 1/dt.getValue();
        this->celsius = DEF_celsius;
    }
    catch (TCLAP::ArgException & e)
    {
        printf("TCLAP error: %s for argument %s\n", e.error().c_str(), e.argId().c_str());
    }
}
