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

    double delta_p2 = 2.0;

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

    // Prepare FourQuarkFullyConnected parameters
    MNPR::FourQuarkFullyConnected::Par FourQuarkFullyConnectedPar;
    if (fourquark)
    {
        FourQuarkFullyConnectedPar.pIn = momentumPar.mom;
        FourQuarkFullyConnectedPar.pOut = momentumPar.mom;
        FourQuarkFullyConnectedPar.gamma_basis = "diagonal_va_sp_tt";
    }

    // Twist-2 type vertices.
    for (int i_mom=std::max(1, int(1.0 / delta_p2)); i_mom<=int(Nl / delta_p2); i_mom++)
    {

        double p2 = delta_p2 * i_mom;
        std::string momentumFolder = outputFolder + "p2_" + cleanString(std::to_string(p2)) + "/";

        LOG(Debug) << "p2 = " << p2 << std::endl;
        LOG(Debug) << momentumFolder << std::endl;

        // Calculate twist that fulfills 2*twist^2=p^2
        double twist = sqrt(delta_p2 * i_mom / 2.0);
        assert(std::abs(p2 - 2 * twist * twist) < 1e-6);

        // Calculate twist that fulfills 4*twist^2=p^2
        double sym_twist = sqrt(delta_p2 * i_mom / 4.0);
        assert(std::abs(p2 - 4 * sym_twist * sym_twist) < 1e-6);

        std::stringstream sstream;

        // 1 1 0 0 type momenta
        sstream << std::setprecision(17) << twist << "_" << twist << "_0.0_0.0";
        std::string name_1100 = sstream.str();
        sstream.str(std::string());
        sstream.clear();

        // 1 0 1 0 type momenta
        sstream << std::setprecision(17) << twist << "_0.0_" << twist << "_0.0";
        std::string name_1010 = sstream.str();
        sstream.str(std::string());
        sstream.clear();

        // 1 1 1 1 type momenta
        sstream << std::setprecision(17) << sym_twist << "_" << sym_twist << "_" << sym_twist << "_" << ToverL * sym_twist;
        std::string name_1111 = sstream.str();
        sstream.str(std::string());
        sstream.clear();

        // 1 1 1 -1 type momenta
        sstream << std::setprecision(17) << sym_twist << "_" << sym_twist << "_" << sym_twist << "_" << -ToverL * sym_twist;
        std::string name_111n1 = sstream.str();
        sstream.str(std::string());
        sstream.clear();

        // 0 0 0 2 type momenta
        sstream << std::setprecision(17) << "0.0_0.0_0.0_" << 2 * ToverL * sym_twist;
        std::string name_0002 = sstream.str();
        sstream.str(std::string());
        sstream.clear();

        // add additional names for twist-1 and 4

        for(const std::string& name: std::vector<std::string> {name_1100, name_1010, name_1111, name_111n1, name_0002})
        {

            LOG(Debug) << name << std::endl;

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
            std::string propagatorName = "Q_0_" + name;
            quarkPar.solver = twisted_solver_name;
            application.createModule<MFermion::GaugeProp>(propagatorName, quarkPar);

            // Compute and save ExternalLeg to disk
            std::string externalLegName = "ExternalLeg_0_" + name;
            externalLegPar.qIn = propagatorName;
            externalLegPar.output = momentumFolder + externalLegName;
            application.createModule<MNPR::ExternalLeg>(externalLegName, externalLegPar);

            // Compute and save Bilinear to disk
            std::string bilinearName = "MOM_Bilinear_00_" + name + "_" + name;
            BilinearPar.qIn = propagatorName;
            BilinearPar.qOut = propagatorName;
            BilinearPar.output = momentumFolder + bilinearName;
            application.createModule<MNPR::Bilinear>(bilinearName, BilinearPar);

            if (fourquark)
            {
                // Compute and save FourQuarkFullyConnected to disk
                std::string FourQuarkFullyConnectedName = "MOM_FourQuark_00_" + name + "_" + name;
                FourQuarkFullyConnectedPar.qIn = propagatorName;
                FourQuarkFullyConnectedPar.qOut = propagatorName;
                FourQuarkFullyConnectedPar.output = momentumFolder + FourQuarkFullyConnectedName;
                application.createModule<MNPR::FourQuarkFullyConnected>(FourQuarkFullyConnectedName, FourQuarkFullyConnectedPar);
            }
        }

        // Construct relevant combinations for SMOM vertices.
        std::vector<std::pair<std::string, std::string>> smom_name_pairs;
        smom_name_pairs.emplace_back(name_1100, name_1010);
        smom_name_pairs.emplace_back(name_1111, name_111n1);
        smom_name_pairs.emplace_back(name_1111, name_0002);

        for(auto const i_pair: smom_name_pairs)
        {
            std::string name_in = i_pair.first;
            std::string name_out = i_pair.second;

            LOG(Debug) << "SMOM names : " << name_in << ", " << name_out << std::endl;

            // Compute and save SMOM Bilinear to disk
            std::string bilinearName = "SMOM_Bilinear_00_" + name_in + "_" + name_out;
            BilinearPar.qIn = "Q_0_" + name_in;
            BilinearPar.qOut = "Q_0_" + name_out;
            BilinearPar.output = momentumFolder + bilinearName;
            application.createModule<MNPR::Bilinear>(bilinearName, BilinearPar);

            if (fourquark)
            {
                // Compute and save FourQuarkFullyConnected to disk
                std::string FourQuarkFullyConnectedName = "SMOM_FourQuark_00_" + name_in + "_" + name_out;
                FourQuarkFullyConnectedPar.qIn = "Q_0_" + name_in;
                FourQuarkFullyConnectedPar.qOut = "Q_0_" + name_out;
                FourQuarkFullyConnectedPar.output = momentumFolder + FourQuarkFullyConnectedName;
                application.createModule<MNPR::FourQuarkFullyConnected>(FourQuarkFullyConnectedName, FourQuarkFullyConnectedPar);
            }
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
