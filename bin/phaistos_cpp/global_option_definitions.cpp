namespace module_distance_restraints {

//! Module: energy term initialization
template <typename SETTINGS_MODIFIER>
struct EnergyOptions {

     //! Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   CHAIN_TYPE *chain,
                   DBN_TYPE *dbn) {                    
     }

     //! Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   ChainFB *chain,
                   DBN_TYPE *dbn) {

          // Import namespace for make_vector
          using namespace boost::fusion;

          // Distance restraints term
          for (int counter = occurrences[prefix+"-distance-restraints"]; counter > 0; counter--) {

               // Create settings object
               typedef TermDistanceRestraints<ChainFB> EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add an options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "Distance restraints term (" + prefix + ")",
                         prefix+"-distance-restraints", settings,
                         make_vector(
                              make_vector(std::string("active-restraints"),
                                          std::string("Number of active restraints (default is 1)"),
                                          &settings->active_restraints),
                              make_vector(std::string("bias"),
                                          std::string("Force the energy difference of a change in restraints to be zero by adding a bias (default True) "),
                                          &settings->bias),
                              make_vector(std::string("upl-filename"),
                                          std::string("CYANA .UPL formatted list of distance restraints."),
                                          &settings->upl_filename)
                              )), super_group, counter==1);
          }
     }
};

}
