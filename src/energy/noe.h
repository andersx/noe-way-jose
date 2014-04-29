// noe.h -- Simple NOE term
// Copyright (C) 2011 by Anders Christensen.
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.

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
//! Class for the noe energy term
class TermNoe: public EnergyTermCommon<TermNoe<CHAIN_TYPE>, CHAIN_TYPE> {

private:

    //! For convenience, define local EnergyTermCommon
    typedef phaistos::EnergyTermCommon<TermNoe<CHAIN_TYPE>,CHAIN_TYPE> EnergyTermCommon;

public:

    //! Returns a random integer from {min, min+1, ... max-1, max}.
    int rand_int(int min, int max, RandomNumberEngine *rne) {
         boost::uniform_smallint<> distribution(min, max);
         boost::variate_generator<RandomNumberEngine&, boost::uniform_smallint<> > generator(*rne, distribution);
         return generator();
    }


     RandomNumberEngine *random_number_engine;

    //!
    class AmbiguousContact {
    public:

         //! Index in chain of residue1
         int residue_index1;

         //! Atom type of atom1
         std::vector<definitions::AtomEnum> atom1_types;
         int residue_index2;
         std::vector<definitions::AtomEnum> atom2_types;

         //! Distance between atom1 and atom2
         double distance;

         //! Default constructor
         AmbiguousContact() {}

         //! Copy constructor
         AmbiguousContact(int residue_index1, std::vector<definitions::AtomEnum> atom1_types,
                 int residue_index2, std::vector<definitions::AtomEnum> atom2_types,
                 double distance)
              : residue_index1(residue_index1),
                atom1_types(atom1_types),
                residue_index2(residue_index2),
                atom2_types(atom2_types),
                distance(distance) {}

          //! Overload << operator for Contact (compact output)
          friend std::ostream &operator<<(std::ostream &o, const AmbiguousContact &c) {
               o << "("
                 << c.residue_index1 + 1 << ","
                 << c.atom1_types << ","
                 << c.residue_index2 + 1 << ","
                 << c.atom2_types << ","
                 << c.distance << ")";
               return o;
          }
     };


    void shuffle_vector(std::vector<AmbiguousContact> &v, RandomNumberEngine *rne) {

        boost::uniform_int<> uni_dist;
        boost::variate_generator<RandomNumberEngine&, boost::uniform_int<> > randomNumber(*rne, uni_dist);
        std::random_shuffle(v.begin(), v.end(),randomNumber);
    }

    //! list of Contact objects
    std::vector<AmbiguousContact> contact_map;
    std::vector<AmbiguousContact> contact_map_old;

    bool did_swap;

    //! Local setting class
    const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
    public:

        //! The force constant between CA atoms
        double weight_constant;
        unsigned int active_restraints;
        bool seamless;

        //! File containing contacts, two supported file formats:
        std::string upl_filename;

        //! Constructor and a reasonable default values
        Settings(double weight_constant=1.0,
                 unsigned int active_restraints=1,
                 bool seamless=false,
                 std::string upl_filename="")
             : weight_constant(weight_constant),
               active_restraints(active_restraints),
               seamless(seamless),
               upl_filename(upl_filename) {}

        //! Output operator
        friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
             o << "weight-constant:" << settings.weight_constant << "\n";
             o << "active-restraints:" << settings.active_restraints<< "\n";
             o << "seamless:" << settings.seamless << "\n";
             o << "upl-filename: " << settings.upl_filename << "\n";
             o << static_cast<const typename EnergyTerm<CHAIN_TYPE>::Settings>(settings);
             return o;
        }

    } settings;


    std::vector<definitions::AtomEnum> string_to_enums(const std::string atom_name, const int res_id) {

        std::vector<definitions::AtomEnum> return_atoms;

        std::map <std::string, std::vector<definitions::AtomEnum> > enum_map;

        // Add types of ambiguous restraints here
        enum_map["QB"]  = boost::assign::list_of(definitions::CB);
        enum_map["QG"]  = boost::assign::list_of(definitions::CG);
        enum_map["QG1"] = boost::assign::list_of(definitions::CG1);
        enum_map["QG2"] = boost::assign::list_of(definitions::CG2);
        enum_map["QD"]  = boost::assign::list_of(definitions::CD)(definitions::CD1)(definitions::CD2);
        enum_map["QD1"] = boost::assign::list_of(definitions::CD1);
        enum_map["QD2"] = boost::assign::list_of(definitions::CD2)(definitions::ND2);
        enum_map["QQG"] = boost::assign::list_of(definitions::CG1)(definitions::CG2);
        enum_map["QQD"] = boost::assign::list_of(definitions::CD1)(definitions::CD2);
        enum_map["QE"]  = boost::assign::list_of(definitions::CE)(definitions::NE2)(definitions::CE1)(definitions::CE2);

        // Check if atom name exists in enum_map.
        if (enum_map.count(atom_name) > 0) {

            std::vector<definitions::AtomEnum> carbon_atoms = enum_map[atom_name];

            for (std::vector<definitions::AtomEnum>::iterator carbon_atom = carbon_atoms.begin();
                 carbon_atom != carbon_atoms.end(); ++carbon_atom) {

                if (! (*(this->chain))[res_id].has_atom((*carbon_atom))) continue;

                std::vector<std::pair<definitions::AtomEnum, int> > neighbours =
                        (*(this->chain))[res_id][(*carbon_atom)]->covalent_neighbours;

                for (std::vector<std::pair<definitions::AtomEnum, int> >::iterator neighbour_pair =
                     neighbours.begin(); neighbour_pair != neighbours.end(); ++neighbour_pair) {

                    definitions::AtomEnum neighbour = neighbour_pair->first;
                    double neighbour_mass = (*(this->chain))[res_id][neighbour]->mass;

                    if (! (*(this->chain))[res_id].has_atom(neighbour)) continue;

                    if (neighbour_mass == definitions::atom_h_weight)
                        return_atoms.push_back(neighbour);
                }

            }

        } else {

            // If it's not in the enum_map,
            // It is not an ambiguous restraints.
            // Just return the proper AtomEnum.
            return boost::assign::list_of(boost::lexical_cast<definitions::AtomEnum>(atom_name));

        }

        return return_atoms;

    }





    //! Read contact map from file
    //! \param chain Chain object for contact map·
    void read_contact_map(CHAIN_TYPE *chain, std::istream &input_stream, const Settings &settings=Settings()) {

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
            int residue1_type_int_read = definitions::StrToAa()[boost::lexical_cast<std::string>(split_line[1])];
            int residue2_type_int_read = definitions::StrToAa()[boost::lexical_cast<std::string>(split_line[4])];

            int residue1_type_int_real = (int)(*(this->chain))[residue1_index].residue_type;
            int residue2_type_int_real = (int)(*(this->chain))[residue2_index].residue_type;

            if (residue1_type_int_read != residue1_type_int_real) {
                std::cerr << "# ERROR: Wrong residue type read from file.\n";
                exit(EXIT_FAILURE);
            }

            if (residue2_type_int_read != residue2_type_int_real) {
                std::cerr << "# ERROR: Wrong residue type read from file.\n";
                exit(EXIT_FAILURE);
            }

            // Generate list of atoms from ambiguous restraints.
            std::string atom1_type = boost::lexical_cast<std::string>(split_line[2]);
            std::string atom2_type = boost::lexical_cast<std::string>(split_line[5]);

            std::vector<definitions::AtomEnum> atom1_atoms = string_to_enums(atom1_type, residue1_index);
            std::vector<definitions::AtomEnum> atom2_atoms = string_to_enums(atom2_type, residue2_index);

            // Read distance
            double distance = boost::lexical_cast<double>(split_line[6]);;

            AmbiguousContact contact(residue1_index,
                                     atom1_atoms,
                                     residue2_index,
                                     atom2_atoms,
                                     distance);

            // Check for residue index is less than length of protein
            if (contact.residue_index1 >= (this->chain)->size() ||
                contact.residue_index2 >= (this->chain)->size()) {
                std::cout << "# CONTACT ERROR: " << contact.residue_index1 + 1 << " "
                                                 << contact.atom1_types << " "
                                                 << contact.residue_index2 + 1 << " "
                                                 << contact.atom2_types << std::endl;
                exit(EXIT_FAILURE);
            }

            // // Check if atoms really exist within the chain.
            // for (unsigned int i = 0; i < contact.atom1_types.size(); i++) {

            //    if (!(*(this->chain))[contact.residue_index1].has_atom(contact.atom1_types[i])) {
            //        std::cout << "# CONTACT ERROR: " << contact.residue_index1 + 1<< " "
            //                                        << contact.atom1_types[i] << std::endl;
            //        exit(EXIT_FAILURE);
            //    }
            // }

            // for (unsigned int i = 0; i < contact.atom2_types.size(); i++) {

            //     if (!(*(this->chain))[contact.residue_index2].has_atom(contact.atom2_types[i])) {
            //         std::cout << "# CONTACT ERROR: " << contact.residue_index2 + 1 << " "
            //                                    << contact.atom2_types[i] << std::endl;
            //         exit(EXIT_FAILURE);
            //     }
            // }

            // Save and print the contact
            contact_map.push_back(contact);
            std::cout << contact << std::endl;

          }
     }

    // https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/de/d50/constraint_file.html
    // http://www.pnas.org/content/suppl/2012/06/25/1203013109.DCSupplemental/Appendix.pdf

    inline double rosetta_bounded_potential(double r_actual, double r_equilibrium, double weight) {

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

        return energy * weight;

    }

     //! Constructor
     TermNoe(CHAIN_TYPE *chain, const Settings &settings=Settings(),
                 RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "noe", settings, random_number_engine),
            random_number_engine(random_number_engine),
            settings(settings) {

        std::ifstream input_stream(settings.upl_filename.c_str());

        if (!input_stream.is_open()) {
            std::cerr << "# Error: Cannot open UPL file " << settings.upl_filename << " .\n";
            exit(EXIT_FAILURE);
        }

        read_contact_map(chain, input_stream, settings);
        shuffle_vector(contact_map, this->random_number_engine);

        contact_map_old = contact_map;

        if (!(settings.active_restraints < contact_map.size())) {
            std::cerr << "# ERROR: " << settings.active_restraints << " active restraints out of " << contact_map.size() << " specified.\n";
            std::cerr << "# ERROR:  Must be less than: " << contact_map.size() << std::endl;
            exit(EXIT_FAILURE);
        }



        did_swap = false;

    }



    double calc_distance_guntert(AmbiguousContact contact) {


        double r_actual = 0.0;

        for (unsigned int j = 0; j < contact.atom1_types.size(); j++) {
            for (unsigned int k = 0; k < contact.atom2_types.size(); k++) {

            Vector_3D r1 = (*(this->chain))(contact.residue_index1,
                          contact.atom1_types[j])->position;
            Vector_3D r2 = (*(this->chain))(contact.residue_index2,
                          contact.atom2_types[k])->position;
            r_actual +=  std::pow((r1 - r2).norm_squared(), -3.0);
            }
        }

        return std::pow(r_actual, (-1.0/6.0));
    }

    double get_energy_contacts(std::vector<AmbiguousContact> map) {

        double energy = 0.0;

        for (unsigned int i=0; i < settings.active_restraints; i++){

            double r_actual = calc_distance_guntert(map[i]);
            double r_equilibrium = map[i].distance;

            if (settings.debug > 0) std::cout << "# NOE Debug: " << map[i]
                                      << " r = " << r_actual
                                      << "   E = " << rosetta_bounded_potential(r_actual, r_equilibrium, settings.weight_constant) << std::endl;

            energy += rosetta_bounded_potential(r_actual, r_equilibrium, settings.weight_constant);
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
            settings(other.settings) { }


    double get_log_bias(MoveInfo *moveInfo = NULL) {

        double bias = 0.0;

        if ((did_swap) && (settings.seamless)) {
            bias = get_energy_contacts(this->contact_map)
                    - get_energy_contacts(this->contact_map_old);
        }

        return bias;
    }


     //! Evaluate the repulsive energy
     double evaluate(MoveInfo *moveInfo=NULL) {

        did_swap = false;

        if (moveInfo) {
            if (moveInfo->modified_angles.empty() == true) {

                unsigned int disable_contact = (unsigned int)rand_int(0, settings.active_restraints - 1,
                                                                      this->random_number_engine);
                unsigned int enable_contact = (unsigned int)rand_int(settings.active_restraints, this->contact_map.size() - 1,
                                                                     this->random_number_engine);

                AmbiguousContact temp_swap;
                temp_swap = this->contact_map[disable_contact];
                this->contact_map[disable_contact] = this->contact_map[enable_contact];
                this->contact_map[enable_contact] = temp_swap;

                did_swap = true;
            }

        }

        double energy = get_energy_contacts(this->contact_map);

        return energy;
     }

    void accept() {
        if (this->did_swap) this->contact_map_old = this->contact_map;
    }

    void reject() {
        if (this->did_swap) this->contact_map = this->contact_map_old;
    }


};

//! Observable specialization for TermMumu
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


     //! Make observation.
     virtual std::string observe(MoveInfo *move_info=NULL,
                                 PHAISTOS_LONG_LONG current_iteration=0,
                                 bool register_only=false) {

            // this->evaluate(move_info);

            double energy = 0.0;

            unsigned int violations_1 = 0;
            unsigned int violations_larger = 0;

            for (unsigned int i=0; i<this->contact_map.size(); i++){

                    double r_actual = calc_distance_guntert(this->contact_map[i]);
                    double r_equilibrium = this->contact_map[i].distance;

                    if (r_actual < r_equilibrium) {
                        // Do nothing
                    } else if (r_actual < r_equilibrium + 1.0) {
                        violations_1 += 1;
                    } else {
                        violations_larger += 1;
                    }


                    energy += rosetta_bounded_potential(r_actual, r_equilibrium, settings.weight_constant);

            }

          std::stringstream s;
          s << std::fixed << std::setprecision(2) << energy;

          std::string output = "[";
          output += s.str();
          output += ",";
          output += boost::lexical_cast<std::string>(violations_1);
          output += ",";
          output += boost::lexical_cast<std::string>(violations_larger);
          output += "]";



          return output;
     }

};




} // end namespace phaistos
#endif
