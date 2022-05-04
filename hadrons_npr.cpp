#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

/* Replaces underscores in a string by spaces */
std::string underscoreToSpace(std::string text)
{
    std::replace(text.begin(), text.end(), '_', ' ');
    return text;
}

/* Replaces . by p and - by n in a string and removes all but one trailing zeros
 * in order to create folder names based of doubles */
std::string cleanString(std::string text)
{
    std::replace(text.begin(), text.end(), '.', 'p');
    std::replace(text.begin(), text.end(), '-', 'n');
    text.erase(text.find_last_not_of("0") + 2, std::string::npos);
    return text;
}

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
    double ToverL = 1.0 * Nt / Nl;
    LOG(Message) << "T/L=" << ToverL << std::endl;

    // Input options, to be provided via xml file
    std::string outputFolder = "npr_twisted";

    bool QED = false;
    bool fourquark = false;

    if ((QED) && (fourquark))
    {
        LOG(Error) << "QED fourquark operators are not implemented." << std::endl;
        exit(EXIT_FAILURE);
    }

    double every = 2;

    double mass = 0.1;
    double csw = 1.1;

    // End of input options

    outputFolder += "/m" + cleanString(std::to_string(mass)) + "/";
    LOG(Message) << "outputFolder: " << outputFolder << std::endl;

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
    if (QED); application.createModule<FermionAction>("action", actionDPar);

    FermionActionF::Par actionFPar;
    actionFPar.gauge = "gauge_F";
    actionFPar.mass = mass;
    actionFPar.boundary = "1.0 1.0 1.0 1.0";
    actionFPar.twist = "0 0 0 0";
    // Wilson only
    actionFPar.cF = 1.0;
    actionFPar.csw_r = csw;
    actionFPar.csw_t = csw;
    if (QED); application.createModule<FermionActionF>("action_F", actionFPar);

    MixedPrecisionSolver::Par solverPar;
    solverPar.residual = 1e-9;
    solverPar.maxInnerIteration = 300000;
    solverPar.maxOuterIteration = 100;
    solverPar.innerGuesser = "";
    solverPar.outerGuesser = "";

    if (QED)
    {
        solverPar.outerAction = "action";
        solverPar.innerAction = "action_F";
        application.createModule<MixedPrecisionSolver>("cg", solverPar);
    }

    // Create base momentum source
    MSource::Momentum::Par momentumPar;
    momentumPar.mom = "0 0 0 0";
    application.createModule<MSource::Momentum>("zero_momentum_source", momentumPar);

    // Prepate propagator parameters
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.source = "zero_momentum_source";

    // Prepare ExternalLeg parameters
    MNPR::ExternalLeg::Par externalLegPar;
    externalLegPar.pIn = momentumPar.mom;

    // Prepare Bilinear parameters
    MNPR::Bilinear::Par BilinearPar;
    BilinearPar.pIn = momentumPar.mom;
    BilinearPar.pOut = momentumPar.mom;

    // Twist-2 type vertices.
    for (int i_mom=std::max(1, int(1.0 / every)); i_mom<=int(Nl / every); i_mom++)
    {
        double twist = sqrt(every * i_mom / 2.0);
        std::string p2 = std::to_string(2 * twist * twist);

        std::string momentumFolder = outputFolder + "p2_" + cleanString(p2) + "/";

        LOG(Debug) << i_mom << "\t" << twist << std::endl;
        LOG(Debug) << "p2 = " << p2 << std::endl;
        LOG(Debug) << momentumFolder << std::endl;

        std::stringstream sstream1;
        sstream1 << std::setprecision(17) << twist << "_" << twist << "_0.0_0.0";
        std::string name1 = sstream1.str();

        std::stringstream sstream2;
        sstream2 << std::setprecision(17) << twist << "_0.0_" << twist << "_0.0";
        std::string name2 = sstream2.str();

        for(const std::string& name: std::vector<std::string> {name1, name2})
        {
            LOG(Debug) << "Name: " << name << std::endl;
            LOG(Debug) << "Name space: " << underscoreToSpace(name) << std::endl;

            // Construct action and solver names
            std::string twisted_action_name = "action_twisted_" + name;
            std::string twisted_action_nameF = twisted_action_name + "F";
            std::string twisted_solver_name = "cg_twisted_" + name;

            // Create action modules for given twist
            actionDPar.twist = actionFPar.twist = underscoreToSpace(name);
            application.createModule<FermionAction>(twisted_action_name, actionDPar);
            application.createModule<FermionActionF>(twisted_action_nameF, actionFPar);

            // Create corresponding solver module
            solverPar.outerAction = twisted_action_name;
            solverPar.innerAction = twisted_action_nameF;
            application.createModule<MixedPrecisionSolver>(twisted_solver_name, solverPar);

            // Create propagator module
            std::string propagatorName = "Q_" + name;
            quarkPar.solver = twisted_solver_name;
            application.createModule<MFermion::GaugeProp>(propagatorName, quarkPar);

            // Compute and save ExternalLeg to disk
            std::string externalLegName = "ExternalLeg_" + name;
            externalLegPar.qIn = propagatorName;
            externalLegPar.output = momentumFolder + externalLegName;
            application.createModule<MNPR::ExternalLeg>(externalLegName, externalLegPar);

            // Compute and save Bilinear to disk
            std::string bilinearName = "MOM_Bilinear_" + name + "_" + name;
            BilinearPar.qIn = propagatorName;
            BilinearPar.qOut = propagatorName;
            BilinearPar.output = momentumFolder + bilinearName;
            application.createModule<MNPR::Bilinear>(bilinearName, BilinearPar);
        }
    }

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
