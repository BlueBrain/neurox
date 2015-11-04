﻿#include "neurox/neurox.h"
#include "tclap/CmdLine.h"
#include "stdio.h"
#include "string.h"

using namespace neurox::tools;

CmdLineParser::CmdLineParser ():
  //from nrnoptarg.cpp::cn_parameters():
  tstart(0), tstop(100), dt(0.025), dt_io(0.1), celsius(34),
  voltage(-65), forwardSkip(0),  prcellgid(-1)
{
    //TODO: missing some inits
    memset(patternStim,'\0',512);
    memset(inputPath,'\0',512);
    memset(outputPath,'\0',512);
}

CmdLineParser::~CmdLineParser(){}

CmdLineParser::CmdLineParser (int argc, char** argv):
    CmdLineParser()
{
    Parse(argc, argv);
}

void CmdLineParser::Parse(int argc, char ** argv)
{
    try {
#ifndef VERSION_STRING
        TCLAP::CmdLine cmd("neurox simulator.", ' ');
#else
        TCLAP::CmdLine cmd("neurox simulator.", ' ', VERSION_STRING);
#endif

        //neurox only command line arguments (NOTE: SwitchArg does not require cmd.add())
        TCLAP::SwitchArg multiMex("M", "multimex", "activates graph-based parallelism of mechanisms.", cmd, false);
        TCLAP::SwitchArg allReduceAtLocality("L", "reduce-by-locality", "perform HPX all-reduce operation at locality level instead of neuron level (better for small cluster).", cmd, false);
        TCLAP::SwitchArg outputStatistics("1", "output-statistics", "outputs files with memory consumption and mechanism distribution.", cmd, false);
        TCLAP::SwitchArg outputMechanismsDot("2", "output-mechs", "outputs mechanisms.dot with mechanisms dependencies.", cmd, false);
        TCLAP::SwitchArg outputNetconsDot("3", "output-netcons", "outputs netcons.dot with netcons information across neurons.", cmd, false);
        TCLAP::SwitchArg outputCompartmentsDot("4", "output-compartments", "outputs compartments_*.dot files displaying neurons morpholgies.", cmd, false);

        TCLAP::ValueArg<int> branchingDepth("B","multisplix","depth tree-based parallelism of morphologies (0 = no branch parallelism, default)", false, 0, "int");
        cmd.add(branchingDepth);
        TCLAP::ValueArg<int> algorithm("A","algorithm",
                                       "BackwardEulerCoreneuron [0],\
                                        BackwardEulerWithAllReduceBarrier [1] (default),\
                                        BackwardEulerWithSlidingTimeWindow [2],\
                                        BackwardEulerWithTimeDependencyLCO [3],\
                                        BackwardEulerCoreneuronDebug [4],\
                                        All methods sequentially (NOTE: neurons data does not reset) [9]",
                                       false, 1, "int");
        cmd.add(algorithm);

        //coreneuron command line parameters
        TCLAP::SwitchArg coreneuronMpiExecution("m", "mpi", "activates coreneuron MPI based execution.", cmd, false);
        TCLAP::ValueArg<floble_t> tstart("s","tstart","Execution start time (msecs). The default value is 0",false, 0 ,"floble_t");
        cmd.add(tstart);
        TCLAP::ValueArg<floble_t> tstop("e","tstop","Execution stop time (msecs). The default value is 100",false, 100 ,"floble_t");
        cmd.add(tstop);
        TCLAP::ValueArg<floble_t> dt("t","dt","Execution time step (msecs). The default value is 0.025",false, DEF_dt ,"floble_t");
        cmd.add(dt);
        TCLAP::ValueArg<floble_t> dt_io("i","dt_io","I/O time step (msecs). The default value is 0.1",false, 0.1 ,"floble_t");
        cmd.add(dt_io);
        TCLAP::ValueArg<floble_t> celsius("l","celsius","System temperatura (celsius degrees). The default value is 34",false, 34.0 ,"floble_t");
        cmd.add(celsius);
        TCLAP::ValueArg<floble_t> forwardskip("k","forwardskip","Set forwardskip time step (msecs). The default value is 0",false, 0 ,"floble_t");
        cmd.add(forwardskip);
        TCLAP::ValueArg<neuron_id_t> prcellgid("g","prcellgid","Output prcellstate information for given gid. The default value is -1",false, -1 ,"neuron_id_t");
        cmd.add(prcellgid);
        TCLAP::ValueArg<std::string> patternStim("p","pattern","Apply patternstim with the spike file. No default value",false, "" ,"string");
        cmd.add(patternStim);
        TCLAP::ValueArg<std::string> inputPath("d","inputpath","Path to input files directory",true,"./input","string");
        cmd.add(inputPath);
        TCLAP::ValueArg<std::string> outputPath("o","outputpath","Path to output directory. The default value is ./output",false,"./output","string");
        cmd.add(outputPath);

        //parse command line arguments
        cmd.parse( argc, argv );

        //copy parameters
        sprintf(this->inputPath, "%s", inputPath.getValue().c_str());
        sprintf(this->outputPath, "%s", outputPath.getValue().c_str());
        sprintf(this->patternStim, "%s", patternStim.getValue().c_str());
        this->tstart = tstart.getValue();
        this->tstop = tstop.getValue();
        this->dt = dt.getValue();
        this->dt_io = dt_io.getValue();
        this->celsius = celsius.getValue();
        this->prcellgid = prcellgid.getValue();
        this->forwardSkip = forwardskip.getValue();
        this->voltage = DEF_vrest;
        this->secondorder = (char) DEF_secondorder;
        this->rev_dt = 1/dt.getValue();
        this->celsius = DEF_celsius;

        this->outputStatistics = outputStatistics.getValue();
        this->outputMechanismsDot = outputMechanismsDot.getValue();
        this->outputNetconsDot = outputNetconsDot.getValue();
        this->outputCompartmentsDot = outputCompartmentsDot.getValue();
        this->multiMex = multiMex.getValue();
        this->allReduceAtLocality = allReduceAtLocality.getValue();
        this->parallelDataLoading = coreneuronMpiExecution.getValue();
        this->branchingDepth = branchingDepth.getValue();
        this->algorithm = (algorithms::AlgorithmType) algorithm.getValue();
        neurox::algorithm = algorithms::Algorithm::New(this->algorithm);

        if (this->branchingDepth<0)
            throw TCLAP::ArgException("branching depth should be zero or a positive value", "branchingdepth");
        if (this->dt<=0)
            throw TCLAP::ArgException("time-step size (ms) should be a positive value", "dt");
        if (this->tstop<=0)
            throw TCLAP::ArgException("execution time (ms) should be a positive value", "tstop");
        floble_t remainder_tstop_tcomm = fmod(this->tstop, algorithms::CoreneuronAlgorithm::CommunicationBarrier::commStepSize*this->dt);
        if (! (remainder_tstop_tcomm<0.00001 || remainder_tstop_tcomm> algorithms::CoreneuronAlgorithm::CommunicationBarrier::commStepSize*this->dt-0.00001))
            throw TCLAP::ArgException("execution time " +to_string(this->tstop) + "ms  should be a multiple of the communication delay " +
                                      to_string( algorithms::CoreneuronAlgorithm::CommunicationBarrier::commStepSize * this->dt) + " ms", "tstop");
    }
    catch (TCLAP::ArgException & e)
    {
        printf("TCLAP error: %s (%s).\n", e.error().c_str(), e.argId().c_str());
        exit(1);
    }
}
