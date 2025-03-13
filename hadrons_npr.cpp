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
    text.erase(text.find_last_not_of("0") + 1, std::string::npos);
    return text;
}

struct ExternalLegEntry: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, qIn,
                       SqlNotNull<std::string>, momentum,
                       SqlNotNull<std::string>, photon_insertions);
};

struct VertexEntry: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, qIn,
                       SqlNotNull<std::string>, qOut,
                       SqlNotNull<std::string>, photon_insertions);
};

struct FourFermionEntry: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, qIn,
                       SqlNotNull<std::string>, qOut,
                       SqlNotNull<std::string>, lIn,
                       SqlNotNull<std::string>, lOut,
                       SqlNotNull<std::string>, photon_insertions);
};

namespace NprInputs
{
    class NprOptions : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(NprOptions,
                                        double,     min_ap2,
                                        double,     max_ap2,
                                        double,     delta_ap2,
                                        bool,       QED,
                                        bool,       fourquark,
                                        std::string, outputFolder);
    };

    class GaugeField : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeField,
                                        std::string, gaugeFieldType,
                                        std::string, gaugeFieldPath,
                                        bool,        smear,
                                        int,         steps,
                                        double,      rho);
    };

    class Action : Serializable
    {
    public:
        #ifdef MOBIUS
            GRID_SERIALIZABLE_CLASS_MEMBERS(Action,
                                            double, mass,
                                            double, M5,
                                            int, Ls);
        #else // WilsonExpClover
            GRID_SERIALIZABLE_CLASS_MEMBERS(Action,
                                            double, mass,
                                            double, csw);
        #endif
    };

    class Solver : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Solver,
                                        double, residual,
                                        int, maxInnerIteration,
                                        int, maxOuterIteration);
    };
}

struct NprPar
{
    NprInputs::NprOptions nprOptions;
    NprInputs::GaugeField gaugeField;
    NprInputs::Action action;
    NprInputs::Solver solver;
};

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

    Grid_init(&argc, &argv);

    Application            application;
    Application::GlobalPar globalPar;
    NprPar                 par;

    ExternalLegEntry       elEntry;
    VertexEntry            vertexEntry;
    FourFermionEntry       fourFermionEntry;

    // read parameters from xml file
    {
        XmlReader reader(parameterFileName);

        read(reader, "global", globalPar);
        read(reader, "nprOptions", par.nprOptions);
        read(reader, "gaugeField", par.gaugeField);
        read(reader, "action", par.action);
        read(reader, "solver", par.solver);
    }

    application.setPar(globalPar);

    // Check lattice geometry
    auto geometry = GridDefaultLatt();

    if (!(geometry[Xp] == geometry[Yp] && geometry[Xp] == geometry[Zp]))
    {
        LOG(Error) << "The current implementation requires that all three spatial extents are identical." << std::endl;
        exit(EXIT_FAILURE);
    }

    int Nl = geometry[Xp];
    int Nt = geometry[Tp];
    double ToverL = 1.0 * Nt / Nl;
    LOG(Debug) << "T/L=" << ToverL << std::endl;

    double min_ap2 = par.nprOptions.min_ap2;
    double max_ap2 = par.nprOptions.max_ap2;
    double delta_ap2 = par.nprOptions.delta_ap2;
    bool QED = par.nprOptions.QED;
    bool fourquark = par.nprOptions.fourquark;
    std::string outputFolder = par.nprOptions.outputFolder;

    std::ostringstream massStream;
    massStream.precision(17);
    massStream << std::fixed << par.action.mass;

    outputFolder += "/m" + cleanString(massStream.str()) + "/";
    LOG(Debug) << "outputFolder: " << outputFolder << std::endl;

    #ifdef BICGSTAB
        using MixedPrecisionSolver = MSolver::MixedPrecisionRBPrecBiCGSTAB;
    #else
        using MixedPrecisionSolver = MSolver::MixedPrecisionRBPrecCG;
    #endif

    #ifdef MOBIUS
        using FermionAction = MAction::MobiusDWF;
        using FermionActionF = MAction::MobiusDWFF;
    #else // WilsonExpClover
        using FermionAction = MAction::WilsonExpClover;
        using FermionActionF = MAction::WilsonExpCloverF;
    #endif

    // Create gauge field
    if (par.gaugeField.gaugeFieldType == "Unit")
    {
        application.createModule<MGauge::Unit>("gauge");
    }
    else if (par.gaugeField.gaugeFieldType == "Random")
    {
        application.createModule<MGauge::Random>("gauge");
    }
    else if (par.gaugeField.gaugeFieldType == "Nersc")
    {
        MIO::LoadNersc::Par gaugePar;
        gaugePar.file = par.gaugeField.gaugeFieldPath;
        application.createModule<MIO::LoadNersc>("gauge", gaugePar);
    }
    else if (par.gaugeField.gaugeFieldType == "openQcd")
    {
        MIO::LoadOpenQcd::Par gaugePar;
        gaugePar.file = par.gaugeField.gaugeFieldPath;
        application.createModule<MIO::LoadOpenQcd>("gauge", gaugePar);
    }
    else
    {
        LOG(Error) << "Unknown gaugeFieldType." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string gaugeFieldName = "gauge";
    if (par.gaugeField.smear)
    {
        MGauge::StoutSmearing::Par smgaugePar;
        smgaugePar.gauge = "gauge";
        smgaugePar.steps = par.gaugeField.steps;
        smgaugePar.rho = par.gaugeField.rho;
        application.createModule<MGauge::StoutSmearing>("smgauge", smgaugePar);
        gaugeFieldName = "smgauge";

    }

    MUtilities::GaugeSinglePrecisionCast::Par gaugeFPar;
    gaugeFPar.field = gaugeFieldName;
    application.createModule<MUtilities::GaugeSinglePrecisionCast>("gauge_F", gaugeFPar);

    if ((QED) && (fourquark))
    {
        // Free field required for lepton propagtors.
        application.createModule<MGauge::Unit>("free_field");
        gaugeFPar.field = "free_field";
        application.createModule<MUtilities::GaugeSinglePrecisionCast>("free_field_F", gaugeFPar);
    }

    // Set base parameters for fermion action
    FermionAction::Par actionDPar;
    actionDPar.gauge = gaugeFieldName;
    actionDPar.mass = par.action.mass;
    actionDPar.boundary = "1.0 1.0 1.0 1.0";
    actionDPar.twist = "0 0 0 0";
    #ifdef MOBIUS
        actionDPar.Ls = par.action.Ls;
        actionDPar.M5 = par.action.M5;
        actionDPar.b = 1.5;
        actionDPar.c = 0.5;
    #else // WilsonExpClover
        actionDPar.cF = 1.0;
        actionDPar.csw_r = par.action.csw;
        actionDPar.csw_t = par.action.csw;
    #endif

    FermionActionF::Par actionFPar;
    actionFPar.gauge = "gauge_F";
    actionFPar.mass = actionDPar.mass;
    actionFPar.boundary = actionDPar.boundary;
    actionFPar.twist = actionDPar.twist;
    #ifdef MOBIUS
        actionFPar.Ls = actionDPar.Ls;
        actionFPar.M5 = actionDPar.M5;
        actionFPar.b = actionDPar.b;
        actionFPar.c = actionDPar.c;
    #else // WilsonExpClover
        actionFPar.cF = actionDPar.cF;
        actionFPar.csw_r = actionDPar.csw_r;
        actionFPar.csw_t = actionDPar.csw_t;
    #endif

    if (QED)
    {
        application.createModule<FermionAction>("action", actionDPar);
        application.createModule<FermionActionF>("action_F", actionFPar);
    }
    FermionAction::Par freeActionDPar;
    FermionActionF::Par freeActionFPar;
    if ((QED) && (fourquark))
    {
        // Set base parameters for free fermion action
        // TODO: Is it correct to set the mass to zero and keep Ls and M5?
        freeActionDPar.gauge = "free_field";
        freeActionDPar.mass = 0.0;
        freeActionDPar.boundary = "1.0 1.0 1.0 1.0";
        freeActionDPar.twist = "0 0 0 0";
        #ifdef MOBIUS
            freeActionDPar.Ls = par.action.Ls;
            freeActionDPar.M5 = par.action.M5;
            freeActionDPar.b = 1.5;
            freeActionDPar.c = 0.5;
        #else // WilsonExpClover
            freeActionDPar.cF = 1.0;
            freeActionDPar.csw_r = 1.0;
            freeActionDPar.csw_t = 1.0;
        #endif

        freeActionFPar.gauge = "free_field_F";
        freeActionFPar.mass = freeActionDPar.mass;
        freeActionFPar.boundary = freeActionDPar.boundary;
        freeActionFPar.twist = freeActionDPar.twist;
        #ifdef MOBIUS
            freeActionFPar.Ls = 8;
            freeActionFPar.M5 = 1.0;
            freeActionFPar.b = freeActionDPar.b;
            freeActionFPar.c = freeActionDPar.c;
        #else // WilsonExpClover
            freeActionFPar.cF = freeActionDPar.cF;
            freeActionFPar.csw_r = freeActionDPar.csw_r;
            freeActionFPar.csw_t = freeActionDPar.csw_t;
        #endif
    }

    // Set base parameters for solver
    MixedPrecisionSolver::Par solverPar;
    solverPar.residual = par.solver.residual;
    solverPar.maxInnerIteration = par.solver.maxInnerIteration;
    solverPar.maxOuterIteration = par.solver.maxOuterIteration;
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
        FourQuarkFullyConnectedPar.gamma_basis = "va_av";
    }
    MNPR::FourFermionFullyConnected::Par FourFermionFullyConnectedPar;
    if ((QED) && (fourquark))
    {
        FourFermionFullyConnectedPar.pIn = momentumPar.mom;
        FourFermionFullyConnectedPar.pOut = momentumPar.mom;
    }

    // Prepare stochastic electromagnetic field and sequential insertion for the QED case
    MGauge::StochasticQedL::Par stochasticQedLPar;
    MSource::SeqAslash::Par seqAslashPar;
    MSource::SeqGamma::Par seqGammaPar;
    if (QED)
    {
        stochasticQedLPar.gauge = PhotonR::Gauge::feynman;
        stochasticQedLPar.improvement = "1.666666667e-1"; // QED_r
        application.createModule<MGauge::StochasticQedL>("StochEm", stochasticQedLPar);

        seqAslashPar.tA = 0;
        seqAslashPar.tB = Nt;
        seqAslashPar.emField = "StochEm";
        seqAslashPar.mom = "0 0 0 0";

        seqGammaPar.tA = 0;
        seqGammaPar.tB = Nt;
        seqGammaPar.gamma = Gamma::Algebra::Identity;
        seqGammaPar.mom = "0 0 0 0";
    }

    // Create modules that compute propagators and vertices on the three twist geometry trajectories
    for (int i_mom=0; i_mom<=int((max_ap2-min_ap2)/delta_ap2); i_mom++)
    {

        double p2 = (min_ap2+i_mom*delta_ap2)*Nl*Nl/4/M_PI/M_PI; // Convert to units of L^2/(2pi)^2
        std::string momentumFolder = outputFolder + "p2_" + cleanString(std::to_string(p2)) + "/";

        LOG(Debug) << "p2 = " << p2 << std::endl;
        LOG(Debug) << momentumFolder << std::endl;

        // Calculate twist that fulfills 2*twist^2=p^2
        double twist = sqrt(p2 / 2.0);
        assert(std::abs(p2 - 2 * twist * twist) < 1e-6);

        // Calculate twist that fulfills 4*twist^2=p^2
        double sym_twist = sqrt(p2 / 4.0);
        assert(std::abs(p2 - 4 * sym_twist * sym_twist) < 1e-6);

        std::ostringstream sstream;

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
            quarkPar.source = "zero_momentum_source";
            quarkPar.solver = twisted_solver_name; // Use solver with twisted action
            application.createModule<MFermion::GaugeProp>(propagatorName, quarkPar);

            // Compute and save ExternalLeg to disk
            std::string externalLegName = "ExternalLeg_0_" + name;
            externalLegPar.qIn = propagatorName;
            externalLegPar.output = momentumFolder + externalLegName;
            application.createModule<MNPR::ExternalLeg>(externalLegName, externalLegPar);

            elEntry.qIn = externalLegPar.qIn;
            elEntry.momentum = actionDPar.twist;
            elEntry.photon_insertions = "0";
            application.setResultMetadata(externalLegName, "externalLeg", elEntry);

            std::string propagatorName_1 = "Q_1_" + name;
            std::string propagatorName_2 = "Q_2_" + name;
            std::string propagatorName_S = "Q_S_" + name;
            if (QED)
            {
                // Prepare and calculate propagator with one photon insertion.
                std::string seqAslashName_1 = "SeqAslash_1_" + name;
                seqAslashPar.q = propagatorName;
                application.createModule<MSource::SeqAslash>(seqAslashName_1, seqAslashPar);

                quarkPar.source = seqAslashName_1;
                quarkPar.solver = "cg"; // Use untwisted solver
                application.createModule<MFermion::GaugeProp>(propagatorName_1, quarkPar);

                // Prepare and calculate propagator with two photon insertions.
                std::string seqAslashName_2 = "SeqAslash_2_" + name;
                seqAslashPar.q = propagatorName_1;
                application.createModule<MSource::SeqAslash>(seqAslashName_2, seqAslashPar);

                quarkPar.source = seqAslashName_2;
                quarkPar.solver = "cg"; // Use untwisted solver
                application.createModule<MFermion::GaugeProp>(propagatorName_2, quarkPar);

                // Compute and save ExternalLeg with two photon insertions to disk
                std::string externalLegName_2 = "ExternalLeg_2_" + name;
                externalLegPar.qIn = propagatorName_2;
                externalLegPar.output = momentumFolder + externalLegName_2;
                application.createModule<MNPR::ExternalLeg>(externalLegName_2, externalLegPar);

                elEntry.qIn = externalLegPar.qIn;
                elEntry.momentum = actionDPar.twist;
                elEntry.photon_insertions = "2";
                application.setResultMetadata(externalLegName_2, "externalLeg", elEntry);

                // Prepare and calculate propagator with one scalar insertion.
                std::string seqGammaName = "SeqGamma_" + name;
                seqGammaPar.q = propagatorName;
                application.createModule<MSource::SeqGamma>(seqGammaName, seqGammaPar);

                quarkPar.source = seqGammaName;
                quarkPar.solver = "cg"; // Use untwisted solver
                application.createModule<MFermion::GaugeProp>(propagatorName_S, quarkPar);

                // Compute and save ExternalLeg with one scalar insertion to disk
                std::string externalLegName_S = "ExternalLeg_S_" + name;
                externalLegPar.qIn = propagatorName_S;
                externalLegPar.output = momentumFolder + externalLegName_S;
                application.createModule<MNPR::ExternalLeg>(externalLegName_S, externalLegPar);

                elEntry.qIn = externalLegPar.qIn;
                elEntry.momentum = actionDPar.twist;
                elEntry.photon_insertions = "S";
                application.setResultMetadata(externalLegName_S, "externalLeg", elEntry);
            }

            // Compute and save Bilinear to disk
            std::string bilinearName = "MOM_Bilinear_00_" + name + "_" + name;
            BilinearPar.qIn = propagatorName;
            BilinearPar.qOut = propagatorName;
            BilinearPar.output = momentumFolder + bilinearName;
            application.createModule<MNPR::Bilinear>(bilinearName, BilinearPar);

            vertexEntry.qIn = BilinearPar.qIn;
            vertexEntry.qOut = BilinearPar.qOut;
            vertexEntry.photon_insertions = "00";
            application.setResultMetadata(bilinearName, "bilinear", vertexEntry);

            if (QED)
            {
                // Compute and save Bilinear_11 to disk
                std::string bilinearName_11 = "MOM_Bilinear_11_" + name + "_" + name;
                BilinearPar.qIn = propagatorName_1;
                BilinearPar.qOut = propagatorName_1;
                BilinearPar.output = momentumFolder + bilinearName_11;
                application.createModule<MNPR::Bilinear>(bilinearName_11, BilinearPar);

                vertexEntry.qIn = BilinearPar.qIn;
                vertexEntry.qOut = BilinearPar.qOut;
                vertexEntry.photon_insertions = "11";
                application.setResultMetadata(bilinearName_11, "bilinear", vertexEntry);

                // Compute and save Bilinear_20 to disk
                std::string bilinearName_20 = "MOM_Bilinear_20_" + name + "_" + name;
                BilinearPar.qIn = propagatorName_2;
                BilinearPar.qOut = propagatorName;
                BilinearPar.output = momentumFolder + bilinearName_20;
                application.createModule<MNPR::Bilinear>(bilinearName_20, BilinearPar);

                vertexEntry.qIn = BilinearPar.qIn;
                vertexEntry.qOut = BilinearPar.qOut;
                vertexEntry.photon_insertions = "20";
                application.setResultMetadata(bilinearName_20, "bilinear", vertexEntry);

                // Compute and save Bilinear_02 to disk
                std::string bilinearName_02 = "MOM_Bilinear_02_" + name + "_" + name;
                BilinearPar.qIn = propagatorName;
                BilinearPar.qOut = propagatorName_2;
                BilinearPar.output = momentumFolder + bilinearName_02;
                application.createModule<MNPR::Bilinear>(bilinearName_02, BilinearPar);

                vertexEntry.qIn = BilinearPar.qIn;
                vertexEntry.qOut = BilinearPar.qOut;
                vertexEntry.photon_insertions = "02";
                application.setResultMetadata(bilinearName_02, "bilinear", vertexEntry);

                // Compute and save Bilinear_S0 to disk
                std::string bilinearName_S0 = "MOM_Bilinear_S0_" + name + "_" + name;
                BilinearPar.qIn = propagatorName_S;
                BilinearPar.qOut = propagatorName;
                BilinearPar.output = momentumFolder + bilinearName_S0;
                application.createModule<MNPR::Bilinear>(bilinearName_S0, BilinearPar);

                vertexEntry.qIn = BilinearPar.qIn;
                vertexEntry.qOut = BilinearPar.qOut;
                vertexEntry.photon_insertions = "S0";
                application.setResultMetadata(bilinearName_S0, "bilinear", vertexEntry);

                // Compute and save Bilinear_0S to disk
                std::string bilinearName_0S = "MOM_Bilinear_0S_" + name + "_" + name;
                BilinearPar.qIn = propagatorName;
                BilinearPar.qOut = propagatorName_S;
                BilinearPar.output = momentumFolder + bilinearName_0S;
                application.createModule<MNPR::Bilinear>(bilinearName_0S, BilinearPar);

                vertexEntry.qIn = BilinearPar.qIn;
                vertexEntry.qOut = BilinearPar.qOut;
                vertexEntry.photon_insertions = "0S";
                application.setResultMetadata(bilinearName_0S, "bilinear", vertexEntry);
            }

            if (fourquark)
            {
                if (QED)
                {
                    // Construct action and solver names
                    std::string twisted_action_name = "free_action_twisted_" + name;
                    std::string twisted_action_nameF = twisted_action_name + "F";
                    std::string twisted_solver_name = "free_cg_twisted_" + name;

                    // Create action modules for given twist
                    freeActionDPar.twist = freeActionFPar.twist = underscoreToSpace(name);
                    application.createModule<FermionAction>(twisted_action_name, freeActionDPar);
                    application.createModule<FermionActionF>(twisted_action_nameF, freeActionFPar);

                    // Create corresponding solver module
                    solverPar.outerAction = twisted_action_name;
                    solverPar.innerAction = twisted_action_nameF;
                    application.createModule<MixedPrecisionSolver>(twisted_solver_name, solverPar);

                    // Create propagator module
                    std::string leptonPropagatorName = "L_0_" + name;
                    quarkPar.source = "zero_momentum_source";
                    quarkPar.solver = twisted_solver_name; // Use solver with twisted action
                    application.createModule<MFermion::GaugeProp>(leptonPropagatorName, quarkPar);

                    // Compute and save ExternalLeg to disk
                    std::string externalLegName = "LeptonExternalLeg_0_" + name;
                    externalLegPar.qIn = leptonPropagatorName;
                    externalLegPar.output = momentumFolder + externalLegName;
                    application.createModule<MNPR::ExternalLeg>(externalLegName, externalLegPar);

                    elEntry.qIn = externalLegPar.qIn;
                    elEntry.momentum = actionDPar.twist;
                    elEntry.photon_insertions = "0";
                    application.setResultMetadata(externalLegName, "externalLeg", elEntry);

                    std::string leptonPropagatorName_1 = "L_1_" + name;
                    std::string leptonPropagatorName_2 = "L_2_" + name;
                    std::string leptonPropagatorName_S = "L_S_" + name;

                    // Prepare and calculate propagator with one photon insertion.
                    std::string seqAslashName_1 = "LeptonSeqAslash_1_" + name;
                    seqAslashPar.q = leptonPropagatorName;
                    application.createModule<MSource::SeqAslash>(seqAslashName_1, seqAslashPar);

                    quarkPar.source = seqAslashName_1;
                    quarkPar.solver = "cg"; // Use untwisted solver
                    application.createModule<MFermion::GaugeProp>(leptonPropagatorName_1, quarkPar);

                    // Prepare and calculate propagator with two photon insertions.
                    std::string seqAslashName_2 = "LeptonSeqAslash_2_" + name;
                    seqAslashPar.q = leptonPropagatorName_1;
                    application.createModule<MSource::SeqAslash>(seqAslashName_2, seqAslashPar);

                    quarkPar.source = seqAslashName_2;
                    quarkPar.solver = "cg"; // Use untwisted solver
                    application.createModule<MFermion::GaugeProp>(leptonPropagatorName_2, quarkPar);

                    // Compute and save ExternalLeg with two photon insertions to disk
                    std::string externalLegName_2 = "LeptonExternalLeg_2_" + name;
                    externalLegPar.qIn = leptonPropagatorName_2;
                    externalLegPar.output = momentumFolder + externalLegName_2;
                    application.createModule<MNPR::ExternalLeg>(externalLegName_2, externalLegPar);

                    elEntry.qIn = externalLegPar.qIn;
                    elEntry.momentum = actionDPar.twist;
                    elEntry.photon_insertions = "2";
                    application.setResultMetadata(externalLegName_2, "externalLeg", elEntry);

                    // Prepare and calculate propagator with one scalar insertion.
                    std::string seqGammaName = "LeptonSeqGamma_" + name;
                    seqGammaPar.q = leptonPropagatorName;
                    application.createModule<MSource::SeqGamma>(seqGammaName, seqGammaPar);

                    quarkPar.source = seqGammaName;
                    quarkPar.solver = "cg"; // Use untwisted solver
                    application.createModule<MFermion::GaugeProp>(leptonPropagatorName_S, quarkPar);

                    // Compute and save ExternalLeg with one scalar insertion to disk
                    std::string externalLegName_S = "LeptonExternalLeg_S_" + name;
                    externalLegPar.qIn = leptonPropagatorName_S;
                    externalLegPar.output = momentumFolder + externalLegName_S;
                    application.createModule<MNPR::ExternalLeg>(externalLegName_S, externalLegPar);

                    elEntry.qIn = externalLegPar.qIn;
                    elEntry.momentum = actionDPar.twist;
                    elEntry.photon_insertions = "S";
                    application.setResultMetadata(externalLegName_S, "externalLeg", elEntry);

                    std::string FourFermionFullyConnectedName;
                    // Compute and save FourFermionFullyConnected to disk
                    FourFermionFullyConnectedName = "MOM_FourFermion_0000_" + name + "_" + name;
                    FourFermionFullyConnectedPar.qIn = propagatorName;
                    FourFermionFullyConnectedPar.qOut = propagatorName;
                    FourFermionFullyConnectedPar.lIn = leptonPropagatorName;
                    FourFermionFullyConnectedPar.lOut = leptonPropagatorName;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "0000";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_1100 to disk
                    FourFermionFullyConnectedName = "MOM_FourFermion_1100_" + name + "_" + name;
                    FourFermionFullyConnectedPar.qIn = propagatorName_1;
                    FourFermionFullyConnectedPar.qOut = propagatorName_1;
                    FourFermionFullyConnectedPar.lIn = leptonPropagatorName;
                    FourFermionFullyConnectedPar.lOut = leptonPropagatorName;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "1100";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_2000 to disk
                    FourFermionFullyConnectedName = "MOM_FourFermion_2000_" + name + "_" + name;
                    FourFermionFullyConnectedPar.qIn = propagatorName_2;
                    FourFermionFullyConnectedPar.qOut = propagatorName;
                    FourFermionFullyConnectedPar.lIn = leptonPropagatorName;
                    FourFermionFullyConnectedPar.lOut = leptonPropagatorName;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "2000";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_0200 to disk
                    FourFermionFullyConnectedName = "MOM_FourFermion_0200_" + name + "_" + name;
                    FourFermionFullyConnectedPar.qIn = propagatorName;
                    FourFermionFullyConnectedPar.qOut = propagatorName_2;
                    FourFermionFullyConnectedPar.lIn = leptonPropagatorName;
                    FourFermionFullyConnectedPar.lOut = leptonPropagatorName;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "0200";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_S000 to disk
                    FourFermionFullyConnectedName = "MOM_FourFermion_S000_" + name + "_" + name;
                    FourFermionFullyConnectedPar.qIn = propagatorName_S;
                    FourFermionFullyConnectedPar.qOut = propagatorName;
                    FourFermionFullyConnectedPar.lIn = leptonPropagatorName;
                    FourFermionFullyConnectedPar.lOut = leptonPropagatorName;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "S000";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_0S00 to disk
                    FourFermionFullyConnectedName = "MOM_FourFermion_0S00_" + name + "_" + name;
                    FourFermionFullyConnectedPar.qIn = propagatorName;
                    FourFermionFullyConnectedPar.qOut = propagatorName_S;
                    FourFermionFullyConnectedPar.lIn = leptonPropagatorName;
                    FourFermionFullyConnectedPar.lOut = leptonPropagatorName;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "0S00";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_1010 to disk
                    FourFermionFullyConnectedName = "MOM_FourFermion_1010_" + name + "_" + name;
                    FourFermionFullyConnectedPar.qIn = propagatorName_1;
                    FourFermionFullyConnectedPar.qOut = propagatorName;
                    FourFermionFullyConnectedPar.lIn = leptonPropagatorName_1;
                    FourFermionFullyConnectedPar.lOut = leptonPropagatorName;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "1010";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_0110 to disk
                    FourFermionFullyConnectedName = "MOM_FourFermion_0110_" + name + "_" + name;
                    FourFermionFullyConnectedPar.qIn = propagatorName;
                    FourFermionFullyConnectedPar.qOut = propagatorName_1;
                    FourFermionFullyConnectedPar.lIn = leptonPropagatorName_1;
                    FourFermionFullyConnectedPar.lOut = leptonPropagatorName;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "0110";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_0020 to disk
                    FourFermionFullyConnectedName = "MOM_FourFermion_0020_" + name + "_" + name;
                    FourFermionFullyConnectedPar.qIn = propagatorName;
                    FourFermionFullyConnectedPar.qOut = propagatorName;
                    FourFermionFullyConnectedPar.lIn = leptonPropagatorName_2;
                    FourFermionFullyConnectedPar.lOut = leptonPropagatorName;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "0020";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_00S0 to disk
                    FourFermionFullyConnectedName = "MOM_FourFermion_00S0_" + name + "_" + name;
                    FourFermionFullyConnectedPar.qIn = propagatorName;
                    FourFermionFullyConnectedPar.qOut = propagatorName;
                    FourFermionFullyConnectedPar.lIn = leptonPropagatorName_S;
                    FourFermionFullyConnectedPar.lOut = leptonPropagatorName;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "00S0";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);
                }
                else
                {
                    // Compute and save FourQuarkFullyConnected to disk
                    std::string FourQuarkFullyConnectedName = "MOM_FourQuark_00_" + name + "_" + name;
                    FourQuarkFullyConnectedPar.qIn = propagatorName;
                    FourQuarkFullyConnectedPar.qOut = propagatorName;
                    FourQuarkFullyConnectedPar.output = momentumFolder + FourQuarkFullyConnectedName;
                    application.createModule<MNPR::FourQuarkFullyConnected>(FourQuarkFullyConnectedName, FourQuarkFullyConnectedPar);

                    vertexEntry.qIn = FourQuarkFullyConnectedPar.qIn;
                    vertexEntry.qOut = FourQuarkFullyConnectedPar.qOut;
                    vertexEntry.photon_insertions = "00";
                    application.setResultMetadata(FourQuarkFullyConnectedName, "fourquark", vertexEntry);
                }
            }
        }

        // Construct relevant combinations for SMOM vertices.
        std::vector<std::pair<std::string, std::string>> smom_name_pairs;
        smom_name_pairs.emplace_back(name_1100, name_1010);
        smom_name_pairs.emplace_back(name_1111, name_111n1);
        smom_name_pairs.emplace_back(name_1111, name_0002);

        // Compute SMOM vertex functions
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

            vertexEntry.qIn = BilinearPar.qIn;
            vertexEntry.qOut = BilinearPar.qOut;
            vertexEntry.photon_insertions = "00";
            application.setResultMetadata(bilinearName, "bilinear", vertexEntry);

            if (QED)
            {
                // Compute and save Bilinear_11 to disk
                std::string bilinearName_11 = "SMOM_Bilinear_11_" + name_in + "_" + name_out;
                BilinearPar.qIn = "Q_1_" + name_in;
                BilinearPar.qOut = "Q_1_" + name_out;
                BilinearPar.output = momentumFolder + bilinearName_11;
                application.createModule<MNPR::Bilinear>(bilinearName_11, BilinearPar);

                vertexEntry.qIn = BilinearPar.qIn;
                vertexEntry.qOut = BilinearPar.qOut;
                vertexEntry.photon_insertions = "11";
                application.setResultMetadata(bilinearName_11, "bilinear", vertexEntry);

                // Compute and save Bilinear_20 to disk
                std::string bilinearName_20 = "SMOM_Bilinear_20_" + name_in + "_" + name_out;
                BilinearPar.qIn = "Q_2_" + name_in;
                BilinearPar.qOut = "Q_0_" + name_out;
                BilinearPar.output = momentumFolder + bilinearName_20;
                application.createModule<MNPR::Bilinear>(bilinearName_20, BilinearPar);

                vertexEntry.qIn = BilinearPar.qIn;
                vertexEntry.qOut = BilinearPar.qOut;
                vertexEntry.photon_insertions = "20";
                application.setResultMetadata(bilinearName_20, "bilinear", vertexEntry);

                // Compute and save Bilinear_02 to disk
                std::string bilinearName_02 = "SMOM_Bilinear_02_" + name_in + "_" + name_out;
                BilinearPar.qIn = "Q_0_" + name_in;
                BilinearPar.qOut = "Q_2_" + name_out;
                BilinearPar.output = momentumFolder + bilinearName_02;
                application.createModule<MNPR::Bilinear>(bilinearName_02, BilinearPar);

                vertexEntry.qIn = BilinearPar.qIn;
                vertexEntry.qOut = BilinearPar.qOut;
                vertexEntry.photon_insertions = "02";
                application.setResultMetadata(bilinearName_02, "bilinear", vertexEntry);

                // Compute and save Bilinear_S0 to disk
                std::string bilinearName_S0 = "SMOM_Bilinear_S0_" + name_in + "_" + name_out;
                BilinearPar.qIn = "Q_S_" + name_in;
                BilinearPar.qOut = "Q_0_" + name_out;
                BilinearPar.output = momentumFolder + bilinearName_S0;
                application.createModule<MNPR::Bilinear>(bilinearName_S0, BilinearPar);

                vertexEntry.qIn = BilinearPar.qIn;
                vertexEntry.qOut = BilinearPar.qOut;
                vertexEntry.photon_insertions = "S0";
                application.setResultMetadata(bilinearName_S0, "bilinear", vertexEntry);

                // Compute and save Bilinear_0S to disk
                std::string bilinearName_0S = "SMOM_Bilinear_0S_" + name_in + "_" + name_out;
                BilinearPar.qIn = "Q_0_" + name_in;
                BilinearPar.qOut = "Q_S_" + name_out;
                BilinearPar.output = momentumFolder + bilinearName_0S;
                application.createModule<MNPR::Bilinear>(bilinearName_0S, BilinearPar);

                vertexEntry.qIn = BilinearPar.qIn;
                vertexEntry.qOut = BilinearPar.qOut;
                vertexEntry.photon_insertions = "0S";
                application.setResultMetadata(bilinearName_0S, "bilinear", vertexEntry);
            }

            if (fourquark)
            {
                if (QED)
                {
                    std::string FourFermionFullyConnectedName;
                    // Compute and save FourFermionFullyConnected to disk
                    FourFermionFullyConnectedName = "SMOM_FourFermion_0000_" + name_in + "_" + name_out;
                    FourFermionFullyConnectedPar.qIn = "Q_0_" + name_in;
                    FourFermionFullyConnectedPar.qOut = "Q_0_" + name_out;
                    FourFermionFullyConnectedPar.lIn = "L_0_" + name_in;
                    FourFermionFullyConnectedPar.lOut = "L_0_" + name_out;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "0000";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_1100 to disk
                    FourFermionFullyConnectedName = "SMOM_FourFermion_1100_" + name_in + "_" + name_out;
                    FourFermionFullyConnectedPar.qIn = "Q_1_" + name_in;
                    FourFermionFullyConnectedPar.qOut = "Q_1_" + name_out;
                    FourFermionFullyConnectedPar.lIn = "L_0_" + name_in;
                    FourFermionFullyConnectedPar.lOut = "L_0_" + name_out;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "1100";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_2000 to disk
                    FourFermionFullyConnectedName = "SMOM_FourFermion_2000_" + name_in + "_" + name_out;
                    FourFermionFullyConnectedPar.qIn = "Q_2_" + name_in;
                    FourFermionFullyConnectedPar.qOut = "Q_0_" + name_out;
                    FourFermionFullyConnectedPar.lIn = "L_0_" + name_in;
                    FourFermionFullyConnectedPar.lOut = "L_0_" + name_out;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "2000";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_0200 to disk
                    FourFermionFullyConnectedName = "SMOM_FourFermion_0200_" + name_in + "_" + name_out;
                    FourFermionFullyConnectedPar.qIn = "Q_0_" + name_in;
                    FourFermionFullyConnectedPar.qOut = "Q_2_" + name_out;
                    FourFermionFullyConnectedPar.lIn = "L_0_" + name_in;
                    FourFermionFullyConnectedPar.lOut = "L_0_" + name_out;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "0200";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_S000 to disk
                    FourFermionFullyConnectedName = "SMOM_FourFermion_S000_" + name_in + "_" + name_out;
                    FourFermionFullyConnectedPar.qIn = "Q_S_" + name_in;
                    FourFermionFullyConnectedPar.qOut = "Q_0_" + name_out;
                    FourFermionFullyConnectedPar.lIn = "L_0_" + name_in;
                    FourFermionFullyConnectedPar.lOut = "L_0_" + name_out;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "S000";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_0S00 to disk
                    FourFermionFullyConnectedName = "SMOM_FourFermion_0S00_" + name_in + "_" + name_out;
                    FourFermionFullyConnectedPar.qIn = "Q_0_" + name_in;
                    FourFermionFullyConnectedPar.qOut = "Q_S_" + name_out;
                    FourFermionFullyConnectedPar.lIn = "L_0_" + name_in;
                    FourFermionFullyConnectedPar.lOut = "L_0_" + name_out;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "0S00";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_1010 to disk
                    FourFermionFullyConnectedName = "SMOM_FourFermion_1010_" + name_in + "_" + name_out;
                    FourFermionFullyConnectedPar.qIn = "Q_1_" + name_in;
                    FourFermionFullyConnectedPar.qOut = "Q_0_" + name_out;
                    FourFermionFullyConnectedPar.lIn = "L_1_" + name_in;
                    FourFermionFullyConnectedPar.lOut = "L_0_" + name_out;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "1010";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_0110 to disk
                    FourFermionFullyConnectedName = "SMOM_FourFermion_0110_" + name_in + "_" + name_out;
                    FourFermionFullyConnectedPar.qIn = "Q_0_" + name_in;
                    FourFermionFullyConnectedPar.qOut = "Q_1_" + name_out;
                    FourFermionFullyConnectedPar.lIn = "L_1_" + name_in;
                    FourFermionFullyConnectedPar.lOut = "L_0_" + name_out;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "0110";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_0020 to disk
                    FourFermionFullyConnectedName = "SMOM_FourFermion_0020_" + name_in + "_" + name_out;
                    FourFermionFullyConnectedPar.qIn = "Q_0_" + name_in;
                    FourFermionFullyConnectedPar.qOut = "Q_0_" + name_out;
                    FourFermionFullyConnectedPar.lIn = "L_2_" + name_in;
                    FourFermionFullyConnectedPar.lOut = "L_0_" + name_out;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "0020";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);

                    // Compute and save FourFermionFullyConnected_00S0 to disk
                    FourFermionFullyConnectedName = "SMOM_FourFermion_00S0_" + name_in + "_" + name_out;
                    FourFermionFullyConnectedPar.qIn = "Q_0_" + name_in;
                    FourFermionFullyConnectedPar.qOut = "Q_0_" + name_out;
                    FourFermionFullyConnectedPar.lIn = "L_S_" + name_in;
                    FourFermionFullyConnectedPar.lOut = "L_0_" + name_out;
                    FourFermionFullyConnectedPar.output = momentumFolder + FourFermionFullyConnectedName;
                    application.createModule<MNPR::FourFermionFullyConnected>(FourFermionFullyConnectedName, FourFermionFullyConnectedPar);

                    fourFermionEntry.qIn = FourFermionFullyConnectedPar.qIn;
                    fourFermionEntry.qOut = FourFermionFullyConnectedPar.qOut;
                    fourFermionEntry.lIn = FourFermionFullyConnectedPar.lIn;
                    fourFermionEntry.lOut = FourFermionFullyConnectedPar.lOut;
                    fourFermionEntry.photon_insertions = "00S0";
                    application.setResultMetadata(FourFermionFullyConnectedName, "fourfermion", fourFermionEntry);
                }

                else
                {
                    // Compute and save FourQuarkFullyConnected to disk
                    std::string FourQuarkFullyConnectedName = "SMOM_FourQuark_00_" + name_in + "_" + name_out;
                    FourQuarkFullyConnectedPar.qIn = "Q_0_" + name_in;
                    FourQuarkFullyConnectedPar.qOut = "Q_0_" + name_out;
                    FourQuarkFullyConnectedPar.output = momentumFolder + FourQuarkFullyConnectedName;
                    application.createModule<MNPR::FourQuarkFullyConnected>(FourQuarkFullyConnectedName, FourQuarkFullyConnectedPar);

                    vertexEntry.qIn = FourQuarkFullyConnectedPar.qIn;
                    vertexEntry.qOut = FourQuarkFullyConnectedPar.qOut;
                    vertexEntry.photon_insertions = "00";
                    application.setResultMetadata(FourQuarkFullyConnectedName, "fourquark", vertexEntry);
                }
            }
        }
    }

    // Save application to xml file
    Hadrons::mkdir("./hadronsxml");
    application.saveParameterFile("./hadronsxml/hadrons_npr_" + globalPar.runId + ".xml", 16);

    // Execution of the application
    try
    {
        application.run();
    }
    catch (const std::exception& e)
    {
        Exceptions::abort(e);
    }

    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
