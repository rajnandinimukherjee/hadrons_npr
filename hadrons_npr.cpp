#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // parse command line //////////////////////////////////////////////////////
    std::string parameterFileName;

    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file> [Grid options]";
        std::cerr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    parameterFileName = argv[1];

    // initialise Grid /////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);

    // initialise application //////////////////////////////////////////////////
    Application            application;
    Application::GlobalPar globalPar;

    // reading parameters
    {
        XmlReader reader(parameterFileName);

        read(reader, "global", globalPar);

        // read other application-specific parameters here
    }

    // global initialisation
    application.setPar(globalPar);

    // create modules //////////////////////////////////////////////////////////

    auto geometry = GridDefaultLatt();

    if (!(geometry[Xp] == geometry[Yp] && geometry[Xp] == geometry[Zp]))
    {
        LOG(Error) << "The current implementation requires that all three spatial extents are identical." << std::endl;
        exit(EXIT_FAILURE);
    }

    int Nl = geometry[Xp];
    int Nt = geometry[Tp];
    double T_over_L = 1.0 * Nt / Nl;

    LOG(Message) << "T/L=" << T_over_L << std::endl;

    double every = 1;


    double mass = 0.1;
    double csw = 1.1;

    using MixedPrecisionSolver = MSolver::MixedPrecisionRBPrecCG;

    using FermionAction = MAction::WilsonExpClover;
    using FermionActionF = MAction::WilsonExpCloverF;

    application.createModule<MGauge::Unit>("gauge");

    MUtilities::GaugeSinglePrecisionCast::Par gaugeFPar;
    gaugeFPar.field = "gauge";
    application.createModule<MUtilities::GaugeSinglePrecisionCast>("gauge_F", gaugeFPar);

    FermionAction::Par actionDPar;
    actionDPar.gauge = "gauge";
    actionDPar.mass = mass;
    actionDPar.boundary = "1.0 1.0 1.0 1.0";
    actionDPar.twist = "0 0 0 0";
    // Wilson only
    actionDPar.cF = 1.0;
    actionDPar.csw_r = csw;
    actionDPar.csw_t = csw;
    application.createModule<FermionAction>("action", actionDPar);

    FermionActionF::Par actionFPar;
    actionFPar.gauge = "gauge_F";
    actionFPar.mass = mass;
    actionFPar.boundary = "1.0 1.0 1.0 1.0";
    actionFPar.twist = "0 0 0 0";
    // Wilson only
    actionFPar.cF = 1.0;
    actionFPar.csw_r = csw;
    actionFPar.csw_t = csw;
    application.createModule<FermionActionF>("action_F", actionFPar);

    MixedPrecisionSolver::Par solverPar;
    solverPar.outerAction = "action";
    solverPar.innerAction = "action_F";
    solverPar.residual = 1e-9;
    solverPar.maxInnerIteration = 300000;
    solverPar.maxOuterIteration = 100;
    solverPar.innerGuesser = "";
    solverPar.outerGuesser = "";
    application.createModule<MixedPrecisionSolver>("cg", solverPar);

    actionDPar.twist = actionFPar.twist = "0 0.34 0.34 0";
    application.createModule<FermionAction>("action_twisted", actionDPar);
    application.createModule<FermionActionF>("action_twisted_F", actionFPar);

    solverPar.outerAction = "action_twisted";
    solverPar.innerAction = "action_twisted_F";
    application.createModule<MixedPrecisionSolver>("cg_twisted", solverPar);
    // execution ///////////////////////////////////////////////////////////////
    try
    {
        application.run();
    }
    catch (const std::exception& e)
    {
        Exceptions::abort(e);
    }

    // epilogue ////////////////////////////////////////////////////////////////
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
