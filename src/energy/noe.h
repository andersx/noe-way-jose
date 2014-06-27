// noe.h -- Simple NOE term
// Copyright (C) 2014 by Anders S. Christensen, Lars Bratholm
//
// This file is part of PHAISTOS
//
// PHAISTOS free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// PHAISTOS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with PHAISTOS.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef TERMNOE
#define TERMNOE

#include <iostream>
#include <vector>
#include <cmath>
#include <string.h>
#include <limits>

#include <boost/random/uniform_smallint.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/assign/list_of.hpp>

#include "protein/iterators/atom_iterator.h"
#include "energy/energy_term.h"
#include "protein/definitions.h"


namespace phaistos {

template <typename CHAIN_TYPE>
//! Class for the NOE energy term.
class TermNoe: public EnergyTermCommon<TermNoe<CHAIN_TYPE>, CHAIN_TYPE> {

private:

     //! For convenience, define local EnergyTermCommon.
     typedef phaistos::EnergyTermCommon<TermNoe<CHAIN_TYPE>,CHAIN_TYPE> EnergyTermCommon;

public:


     //! Class to hold information about a contact between 
     //! two residues, such as involved atoms and distance.
     class AmbiguousContact {
     public:

          //! Index in chain of residue1.
          int residue_index1;

          //! Index in chain of residue2.
          int residue_index2;

          //! Atom types of contacts in the first residue.
          std::vector<definitions::AtomEnum> atom1_types;

          //! Atom types of contacts in the second residue.
          std::vector<definitions::AtomEnum> atom2_types;

          //! Distance between atom1 and atom2.
          double distance;

          //! Default constructor.
          AmbiguousContact() {}

          //! Copy constructor.
          AmbiguousContact(int residue_index1, 
                           std::vector<definitions::AtomEnum> atom1_types,
                           int residue_index2,
                           std::vector<definitions::AtomEnum> atom2_types,
                           double distance)
              : residue_index1(residue_index1),
                residue_index2(residue_index2),
                atom1_types(atom1_types),
                atom2_types(atom2_types),
                distance(distance) {}

          //! Overload << operator for Contact (compact output).
          friend std::ostream &operator<<(std::ostream &o, 
                                          const AmbiguousContact &c) {
               o << "("
                 << c.residue_index1 + 1 << ","
                 << c.atom1_types << ","
                 << c.residue_index2 + 1 << ","
                 << c.atom2_types << ","
                 << c.distance << ")";
               return o;
          }
     };

     //! Pointer to the global random number engine.
     RandomNumberEngine *random_number_engine;

     //! Returns a random integer from {min, min+1, ... max-1, max}.
     //! \param min Lowest number in the range.
     //! \param max Highest number in the range.
     //! \return A random number.
     int rand_int(const int min, const int max, RandomNumberEngine *rne) {
          boost::uniform_smallint<> distribution(min, max);
          boost::variate_generator<RandomNumberEngine&, boost::uniform_smallint<> > generator(*rne, distribution);
          return generator();
     }

     //! Function to shuffle a std::vector in-place.
     //! \param v Vector of type AmbiguousContact to be shuffled.
     //! \param rne Pointer to the random number engine to use.
     void shuffle_vector(std::vector<AmbiguousContact> &v, RandomNumberEngine *rne) {
          boost::uniform_int<> uni_dist;
          boost::variate_generator<RandomNumberEngine&, boost::uniform_int<> > random_number(*rne, uni_dist);
          std::random_shuffle(v.begin(), v.end(), random_number);
     }


     //! List of Contact objects.
     std::vector<AmbiguousContact> contact_map;

     //! Backup list of Contact objects.
     std::vector<AmbiguousContact> contact_map_old;

     //! Variable to remember if last move involved a
     //! change in active NOE restraints.
     bool did_swap;

     //! Local setting class
     const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
     public:

         //! Number of active restraints throughout the simulation.
         unsigned int active_restraints;

         //! Whether to switch constraints with bias (=true) or without bias (=false).
         bool seamless;

         //! File containing contacts, two supported file formats:
         std::string upl_filename;

         //! Constructor and a reasonable default values.
         Settings(unsigned int active_restraints=1,
                  bool seamless=false,
                  std::string upl_filename="")
              : active_restraints(active_restraints),
                seamless(seamless),
                upl_filename(upl_filename) {}

         //! Output operator.
         friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
              o << "active-restraints:" << settings.active_restraints<< "\n";
              o << "seamless:" << settings.seamless << "\n";
              o << "upl-filename: " << settings.upl_filename << "\n";
              o << static_cast<const typename EnergyTerm<CHAIN_TYPE>::Settings>(settings);
              return o;
         }

     } settings;


     //! Converts a string name of an ambiguous atom name to a list of all possible atoms using the residue id.
     //! \param atom_name The atom name string.
     //! \param res_id The residue number.
     //! \return A std::vector of atoms enums.
     std::vector<definitions::AtomEnum> string_to_enums(const std::string atom_name, const int res_id) {

          using namespace definitions;

          std::vector<AtomEnum> return_atoms;

          std::map <std::string, std::vector<AtomEnum> > enum_map;

          // List of what an ambiguously declared restraint corresponds to
          // Add types of ambiguous restraints here
          enum_map["QB"]  = boost::assign::list_of(CB);
          enum_map["QG"]  = boost::assign::list_of(CG);
          enum_map["QG1"] = boost::assign::list_of(CG1);
          enum_map["QG2"] = boost::assign::list_of(CG2);
          enum_map["QD"]  = boost::assign::list_of(CD)(CD1)(CD2);
          enum_map["QD1"] = boost::assign::list_of(CD1);
          enum_map["QD2"] = boost::assign::list_of(CD2)(ND2);
          enum_map["QQG"] = boost::assign::list_of(CG1)(CG2);
          enum_map["QQD"] = boost::assign::list_of(CD1)(CD2);
          enum_map["QE"]  = boost::assign::list_of(CE)(NE2)(CE1)(CE2);
          enum_map["QE2"]  = boost::assign::list_of(NE2);

          // Check if atom name exists in enum_map.
          if (enum_map.count(atom_name) > 0) {

               std::vector<AtomEnum> carbon_atoms = enum_map[atom_name];

               for (std::vector<AtomEnum>::iterator carbon_atom = carbon_atoms.begin();
                    carbon_atom != carbon_atoms.end(); ++carbon_atom) {

                    if (! (*(this->chain))[res_id].has_atom((*carbon_atom))) continue;

                    std::vector<std::pair<AtomEnum, int> > neighbours =
                         (*(this->chain))[res_id][(*carbon_atom)]->covalent_neighbours;

                    for (std::vector<std::pair<AtomEnum, int> >::iterator neighbour_pair =
                         neighbours.begin(); neighbour_pair != neighbours.end(); ++neighbour_pair) {

                         AtomEnum neighbour = neighbour_pair->first;
                         double neighbour_mass = (*(this->chain))[res_id][neighbour]->mass;

                         if (! (*(this->chain))[res_id].has_atom(neighbour)) continue;

                         if (neighbour_mass == atom_h_weight) return_atoms.push_back(neighbour);
                    }
               }
          } else {

               // If it's not in the enum_map,
               // it is not an ambiguous restraints.
               // Just return the proper AtomEnum.
               return boost::assign::list_of(boost::lexical_cast<AtomEnum>(atom_name));
          }
          return return_atoms;
     }


     //! Read contact map from file, reads everything into this->contact_map
     //! \param chain Chain object for contact map.
     //! \return A list (std::vector) containing all contacts from the contact map.
     std::vector<AmbiguousContact> read_contact_map(CHAIN_TYPE *chain, std::istream &input_stream) {

          using namespace definitions;

          std::vector<AmbiguousContact> contacts_from_file;
          while (input_stream.good()) {

               std::string line;
               std::getline(input_stream, line);

               boost::trim(line);

               if (line.size() == 0 || line[0] == '#') {
                    continue;
               }

               std::vector<std::string> split_line;
               boost::split(split_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

               std::cout << "# NOE-restraints: " << split_line << " interpreted as ";

               int residue1_index = boost::lexical_cast<int>(split_line[0]) - 1;
               int residue2_index = boost::lexical_cast<int>(split_line[3]) - 1;


               // Check residue types
               int residue1_type_int_read = StrToAa()[boost::lexical_cast<std::string>(split_line[1])];
               int residue2_type_int_read = StrToAa()[boost::lexical_cast<std::string>(split_line[4])];

               int residue1_type_int_real = (int)(*(this->chain))[residue1_index].residue_type;
               int residue2_type_int_real = (int)(*(this->chain))[residue2_index].residue_type;

               if ((residue1_type_int_read != residue1_type_int_real) && 
                    (residue1_type_int_read != UNK))  {
                    std::cerr << "# ERROR: Wrong residue type read from file.\n";
                    exit(EXIT_FAILURE);
               }

               if ((residue2_type_int_read != residue2_type_int_real) &&
                    (residue2_type_int_read != UNK))  {
                    std::cerr << "# ERROR: Wrong residue type read from file.\n";
                    exit(EXIT_FAILURE);
               }

               // Generate list of atom types from ambiguous restraints.
               std::string atom1_type = boost::lexical_cast<std::string>(split_line[2]);
               std::string atom2_type = boost::lexical_cast<std::string>(split_line[5]);

               // Generate list of atoms from ambiguous restraints.
               std::vector<AtomEnum> atom1_atoms = string_to_enums(atom1_type, residue1_index);
               std::vector<AtomEnum> atom2_atoms = string_to_enums(atom2_type, residue2_index);

               // Read distance
               double distance = boost::lexical_cast<double>(split_line[6]);;

               // Create an ambiguous contact object
               AmbiguousContact contact(residue1_index, atom1_atoms,
                                        residue2_index, atom2_atoms, distance);

               // Check for residue index is less than length of protein
               if (contact.residue_index1 >= (this->chain)->size() ||
                    contact.residue_index2 >= (this->chain)->size()) {
                    std::cout << "# CONTACT ERROR: " << contact.residue_index1 + 1 << " "
                                                     << contact.atom1_types << " "
                                                     << contact.residue_index2 + 1 << " "
                                                     << contact.atom2_types << std::endl;
                    exit(EXIT_FAILURE);
               }

               // Save and print the contact
               contacts_from_file.push_back(contact);
               std::cout << contact << std::endl;

           }

          return contacts_from_file;
     }

     //! The Rosetta flat bottom potential for NOE restraints. For references see :
     //! https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/de/d50/constraint_file.html
     //! and http://www.pnas.org/content/suppl/2012/06/25/1203013109.DCSupplemental/Appendix.pdf
     //! \param r_actual The distance shortest between two restraints in the structure.
     //! \param r_equilibrium The equilibrium distance for the restrain.
     //! \return The energy corresponding to the difference in distance.
     double rosetta_bounded_potential(const double r_actual, const double r_equilibrium) {

          // Variables are named consistently with references give above.
          // Given like this for readability. Could easily be vastly optimized.
          double d0 = r_equilibrium;
          double lb = 1.5;
          double ub = d0 + 0.15;
          double sd = 0.3;
          double rswitch = ub + 0.5 * sd;

          double x = r_actual;

          double energy = 0.0;

          if (x < lb) {
               energy = std::pow((x - lb)/sd, 2.0);

          } else if (x < ub) {
               energy = 0.0;

          } else if (x < rswitch) {
               energy = std::pow((x - ub)/sd, 2.0);

          } else {
               energy = 1.0 / sd * (x - rswitch) + std::pow(rswitch, 2.0);

          }

          return energy;

     }

     //! Constructor
     TermNoe(CHAIN_TYPE *chain, const Settings &settings=Settings(),
                 RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "noe", settings, random_number_engine),
            random_number_engine(random_number_engine),
            settings(settings) {

          // Open the UPL formatted file
          std::ifstream input_stream(this->settings.upl_filename.c_str());

          if (!input_stream.is_open()) {
               std::cerr << "# Error: Cannot open UPL file " << this->settings.upl_filename << " .\n";
               exit(EXIT_FAILURE);
          }

       
          // Read data from UPL file into contact map.
          //this->contact_map = this->read_contact_map(this->chain, input_stream, this->settings);
          this->contact_map = this->read_contact_map(this->chain, input_stream);

          // Randomize the order of the contact map (since only the first N, random contact are active).
          this->shuffle_vector(this->contact_map, this->random_number_engine);

          // Create a backup for use when a change is rejected.
          this->contact_map_old = this->contact_map;

          // Check if the user set more active restraints than the length of the contact map.
          if (!(this->settings.active_restraints < this->contact_map.size())) {
               std::cerr << "# ERROR: " << this->settings.active_restraints 
                    << " active restraints out of " << this->contact_map.size() << " specified.\n";
               std::cerr << "# ERROR:  Must be less than: " << this->contact_map.size() << std::endl;
               exit(EXIT_FAILURE);
          }

          // Initially do swap of restraints is performed and this variable is thus set to false.
          this->did_swap = false;

     }


     //! Calculate the shortest distance between two sets of atoms involved in an NOE restraint.
     //! \param contact An object corresponding to an ambiguous contact.
     //! \return The shortest distance.
     double calc_distance_nearest(const AmbiguousContact &contact) {

          double r_actual = std::numeric_limits<double>::infinity();

          for (unsigned int j = 0; j < contact.atom1_types.size(); j++) {
               for (unsigned int k = 0; k < contact.atom2_types.size(); k++) {

               Vector_3D r1 = (*(this->chain))(contact.residue_index1,
                                               contact.atom1_types[j])->position;
               Vector_3D r2 = (*(this->chain))(contact.residue_index2,
                                               contact.atom2_types[k])->position;
               double r_pair = (r1 - r2).norm();

               if (r_pair < r_actual) r_actual = r_pair;

               }
          }

          return r_actual;
     }


     //! Returns the energy of a contact map.
     //! \param map A contact map.
     //! \return The energy associated with the contact map.
     double get_energy_contacts(const std::vector<AmbiguousContact> &map) {

          double energy = 0.0;

          for (unsigned int i=0; i < this->settings.active_restraints; i++){

               double r_actual = this->calc_distance_nearest(map[i]);

               energy += this->rosetta_bounded_potential(r_actual, map[i].distance);
          }

          return energy;

     }


     //! Copy constructor
     TermNoe(const TermNoe &other, RandomNumberEngine *random_number_engine,
             int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            random_number_engine(random_number_engine),
            contact_map(other.contact_map),
            contact_map_old(other.contact_map_old),
            did_swap(other.did_swap),
            settings(other.settings) {

          // Give each thread a defined random_number_engine (so runs can be replicated).
          this->random_number_engine = random_number_engine;
 
     }


     //! Get the log bias from a move in energy. If seamless is enabled, the log bias is exactly
     //! the difference in energy between the contact map before and after the move.
     //! \param move_info
     //! \return log_bias The log bias associated with the change in restraints.
     double get_log_bias(MoveInfo *move_info = NULL) {

          double log_bias = 0.0;

          if ((this->did_swap) && (this->settings.seamless)) {
               log_bias = this->get_energy_contacts(this->contact_map)
                        - this->get_energy_contacts(this->contact_map_old);
          }

          return log_bias;
     }


     //! Evaluate the repulsive energy
     //! \param move_info
     double evaluate(MoveInfo *move_info=NULL) {

          // Notify that a swap has not been performed (yet).
          this->did_swap = false;

          if (move_info) {
               if (move_info->modified_angles.empty() == true) {

                    // Index for a contact in the active part, which runs from 0 to (active_restraints - 1).
                    unsigned int disable_contact = (unsigned int)rand_int(0, 
                                                                          this->settings.active_restraints - 1,
                                                                          this->random_number_engine);

                    // Index for a contact in the inactive part, which runs from active_restraints to (contact_map.size() - 1).
                    unsigned int enable_contact = (unsigned int)rand_int(this->settings.active_restraints,
                                                                         this->contact_map.size() - 1,
                                                                         this->random_number_engine);

                    // Swap an inactive restraint for an active restraint.
                    std::swap(this->contact_map[enable_contact],
                              this->contact_map[disable_contact]);

                    // Notify accept and reject functions, that a swap was performed.
                    this->did_swap = true;
               }

          }

          // Calculate the energy of the contact map and return the energy.
          return this->get_energy_contacts(this->contact_map);

     }

     //! Accept change in restraints and backup contact map.
     void accept() {

          if (this->did_swap) {
               this->contact_map_old = this->contact_map;
          }

     }

     //! Reject change in restraints and roll back to the backup contact map.
     void reject() {

          if (this->did_swap) {
               this->contact_map = this->contact_map_old;
          }

     }

}; // End class TermNoe

//! Observable specialization for TermNoe
template <typename CHAIN_TYPE>
class Observable<TermNoe<CHAIN_TYPE> >: public TermNoe<CHAIN_TYPE>, public ObservableBase {

public:

     //! Local settings class.
     const class Settings: public TermNoe<CHAIN_TYPE>::Settings, public ObservableBase::Settings {
     public:

          //! Constructor. Defines default values for settings object.
          Settings(){}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const typename TermNoe<CHAIN_TYPE>::Settings>(settings);
               o << static_cast<const ObservableBase::Settings>(settings);
               return o;
          }
     } settings; //!< Local settings object·

     //! Constructor.
     //! \param energy_term VisibleVolume energy term object
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const TermNoe<CHAIN_TYPE> &energy_term,
                const ObservableBase::Settings &settings=ObservableBase::Settings(),
                Energy<CHAIN_TYPE> *reference_energy_function=NULL)
          : TermNoe<CHAIN_TYPE>(energy_term),
            settings(dynamic_cast<const Settings&>(settings)) {
     }

     //! Copy Constructor.
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     Observable(const Observable &other, int thread_index,
                typename TermNoe<CHAIN_TYPE>::ChainType *chain)
          : TermNoe<CHAIN_TYPE>(other, thread_index, chain),
            settings(other.settings) {
     }


     //! Clone: Corresponds to a virtual copy constructor
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermNoe<CHAIN_TYPE> *clone(int thread_index=0,
                                typename TermNoe<CHAIN_TYPE>::ChainType *chain=NULL) {
          return new Observable<TermNoe<CHAIN_TYPE> >(*this, thread_index, chain);
     }

     // https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/de/d50/constraint_file.html
     // http://www.pnas.org/content/suppl/2012/06/25/1203013109.DCSupplemental/Appendix.pdf

     //! Make observation.
     virtual std::string observe(MoveInfo *move_info=NULL,
                                 PHAISTOS_LONG_LONG current_iteration=0,
                                 bool register_only=false) {

          //! Energy to be returned
          double energy = 0.0;

          //! Number of restraints in agreement with the current structure
          unsigned int in_agreement = 0;

          //! Number of restraints within 1 angstrom from agreement with the current structure
          unsigned int violations_1 = 0;

          //! Number of restraints within 5 angstrom from agreement with the current structure
          unsigned int violations_5 = 0;

          //! Number of restraints worse than 1 angstrom from agreement with the current structure
          unsigned int violations_larger = 0;

          for (unsigned int i=0; i<this->contact_map.size(); i++){

               double r_actual = this->calc_distance_nearest(this->contact_map[i]);
               double r_equilibrium = this->contact_map[i].distance;

               if (r_actual <= r_equilibrium) {
                   in_agreement += 1;
               } else if (r_actual <= r_equilibrium + 1.0) {
                   violations_1 += 1;
               } else if (r_actual <= r_equilibrium + 5.0) {
                   violations_5 += 1;
               } else {
                   violations_larger += 1;
               }

               energy += this->rosetta_bounded_potential(r_actual, r_equilibrium);

          }

          //! Output stream
          std::stringstream s;
          s << std::fixed << std::setprecision(5) << energy;

          std::string output = "[ ";
          output += s.str();
          output += " , ";
          output += boost::lexical_cast<std::string>(in_agreement);
          output += " , ";
          output += boost::lexical_cast<std::string>(violations_1);
          output += " , ";
          output += boost::lexical_cast<std::string>(violations_5);
          output += " , ";
          output += boost::lexical_cast<std::string>(violations_larger);
          output += " ]";

          return output;
     }

}; // End class Observable


} // end namespace phaistos
#endif
