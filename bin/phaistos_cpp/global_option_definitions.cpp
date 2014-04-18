namespace module_noe {

//! Module: energy term initialization
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

          // noe term
          for (int counter = occurrences[prefix+"-noe"]; counter > 0; counter--) {


               // Create settings object
               typedef typename TermNoe<ChainFB>::Settings Settings;
               boost::shared_ptr<Settings> settings(new Settings());

               // Add an option for harmonic force constant
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "NOE term (" + prefix + ")",
                         prefix+"-noe", settings,
                         make_vector(
                              make_vector(std::string("force-constant"),
                                          std::string("The harmonic force constant (default is 5.0)"),
                                          &settings->force_constant),
                              make_vector(std::string("active-restraints"),
                                          std::string("Number of active restraints (default is 1)"),
                                          &settings->active_restraints),
                              make_vector(std::string("seamless"),
                                          std::string("Seamless restraint switching -- experimental! (default is False)"),
                                          &settings->seamless),
                              make_vector(std::string("contact-map-filename"),
                                          std::string("List of NOE contacts."),
                                          &settings->contact_map_filename)
                              )),
                    super_group, counter==1);

          }
    }
};

}
