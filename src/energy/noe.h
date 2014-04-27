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

// #include <boost/type_traits/is_base_of.hpp>
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

    //! Contact.
    //! A contact between two atoms is specified by the residue numbers and
    //! AtomEnums, a contact distance, a width of the potential well and a
    //! weight.
    //!
    //! The energy of a contact is calculated as either
    //! \f[ -e^{-\left(\frac{r-r_0}{width}\right)^2}*weight \f]
    //! or
    //! \f[ \left(\frac{r-r_0}{width}\right)^2*weight \f]
    //! depending on the argument dist_squared. The former is the default.
    //!
    class Contact {
    public:

         //! Index in chain of residue1
         int residue_index1;

         //! Atom type of atom1
         definitions::AtomEnum atom_type1;

         //! Index in chain of residue2
         int residue_index2;

         //! Atom type of atom1
         definitions::AtomEnum atom_type2;

         //! Distance between atom1 and atom2
         double distance;

         //! Default constructor
         Contact(){}

         //! Constructor
         //! \param residue_index1 index of residue 1
         //! \param atom_type1 definitions::AtomEnum of atom in residue 1
         //! \param residue_index2 index of residue 2
         //! \param atom_type2 definitions::AtomEnum of atom in residue 2
         //! \param distance ideal distance between a1 and a2
         //! \param weight weight for this contact
         Contact(int residue_index1, definitions::AtomEnum atom_type1,
                 int residue_index2, definitions::AtomEnum atom_type2,
                 double distance)
              : residue_index1(residue_index1),
                atom_type1(atom_type1),
                residue_index2(residue_index2),
                atom_type2(atom_type2),
                distance(distance) {
         }

         //! Output as string - verbose version
         std::string to_string_verbose() const {
              std::string output = "";
              output += (boost::lexical_cast<std::string>(residue_index1) + " " +
                         boost::lexical_cast<std::string>(atom_type1) + " " +
                         boost::lexical_cast<std::string>(residue_index2) + " " +
                         boost::lexical_cast<std::string>(atom_type2) + " " +
                         boost::lexical_cast<std::string>(distance));
              return output;
         }

          //! Overload << operator for Contact (compact output)
          friend std::ostream &operator<<(std::ostream &o, const Contact &c) {
               o << "("
                 << c.residue_index1 << ","
                 << c.atom_type1 << ","
                 << c.residue_index2 << ","
                 << c.atom_type2 << ","
                 << c.distance << ")";
               return o;
          }
         // //! \param weight weight for this contact
         //  //! Overload >> operator for Contact pointer (compact input)
         //  friend std::istream &operator>>(std::istream &input, Contact &c) {
         //       input.ignore(1);
         //       input >> c.residue_index1;
         //       input.ignore(1);
         //       input >> c.atom_type1;
         //       input.ignore(1);
         //       input >> c.residue_index2;
         //       input.ignore(1);
         //       input >> c.atom_type2;
         //       input.ignore(1);
         //       input >> c.distance;
         //       input.ignore(1);
         //       return input;
         //  }

     };

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
         AmbiguousContact(){
                   // if(contact.residue_index1 - 1 >= (this->chain)->size() ||
                   //    contact.residue_index2 - 1 >= (this->chain)->size()) {
                   //     std::cout << "# CONTACT ERROR:" << contact.residue_index1 << " "
                   //                                     << contact.atom_type1 << " "
                   //                                     << contact.residue_index2 << " "
                   //                                     << contact.atom_type2 << std::endl;
                   //     exit(EXIT_FAILURE);
                   // }

                   // if(!(*(this->chain))[contact.residue_index1 - 1].has_atom(contact.atom_type1)) {
                   //     std::cout << "# CONTACT ERROR:" << contact.residue_index1 << " "
                   //                                     << contact.atom_type1 << " "
                   //                                     << contact.residue_index2 << " "
                   //                                     << contact.atom_type2 << std::endl;
                   //     exit(EXIT_FAILURE);
                   // }

                   // if(!(*(this->chain))[contact.residue_index1 - 1].has_atom(contact.atom_type1)) {
                   //     std::cout << "# CONTACT ERROR:" << contact.residue_index1 << " "
                   //                                     << contact.atom_type1 << " "
                   //                                     << contact.residue_index2 << " "
                   //                                     << contact.atom_type2 << std::endl;
                   //     exit(EXIT_FAILURE);
                   // }
         }

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

     //! vector of Contact objects
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

          //! File containing contacts, two supported file formats:·
          //! 1: One line: [(residue_index1,atom_type1,residue_index2,atom_type2,distance,weight),(...),...]
          //!    Example: [(0,CA,14,CA,8.35098,1),(0,CA,15,CA,5.72554,1),...]
          //! 2: One contact per line:
          //!     residueno1 atomname1 residueno2 atomname2 idealdistance [width [weight]].··········
          //!    Example:
          //!    5 CA 1 CA 13.029
          //!    6 CA 2 CA 13.168 2.0 1.5
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

     //! Read contact map from file
     //! \param chain Chain object for contact map·
     //! \param input_stream Stream from which data is read.
     //! \param settings settings object

     // 157 GLU  H     168 THR  H       4.5
     void read_contact_map(CHAIN_TYPE *chain, std::istream &input_stream, const Settings &settings=Settings()) {

          while (input_stream.good()) {

               std::string line;
               std::getline(input_stream, line);

               boost::trim(line);

               if (line.size()==0 || line[0] == '#') {
                    continue;
               }

               std::vector<std::string> split_line;
               boost::split(split_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

               // Parse single-line entry (compact format)
               if (split_line.size() == 1) {
                    //std::istringstream buffer(line);
                    //buffer >> contact_map;
                    //std::cout << "1: " << contact_map << "\n";
                    return;
                    //// contact_map = boost::lexical_cast<Contact>(line);
               } else {
                    int residue1_index = boost::lexical_cast<int>(split_line[0]) - 1;
                    int residue2_index = boost::lexical_cast<int>(split_line[3]) - 1;

                    std::string atom1_type = boost::lexical_cast<std::string>(split_line[2]);
                    std::string atom2_type = boost::lexical_cast<std::string>(split_line[5]);

                    std::cout << split_line << std::endl;

                    std::vector<definitions::AtomEnum> atom1_atoms;
                    std::vector<definitions::AtomEnum> atom2_atoms;

                    if (atom1_type == "QD") { atom1_atoms = boost::assign::list_of(definitions::HD1)(definitions::HD2); }

                    else if (atom1_type == "QE") { atom1_atoms = boost::assign::list_of(definitions::HE1)(definitions::HE2); }

                    else if (atom1_type == "QB") { 

                        if (boost::lexical_cast<std::string>(split_line[1]) == "ALA") {
                            atom1_atoms = boost::assign::list_of(definitions::HB1)(definitions::HB2)(definitions::HB3); 
                        } else {
                            atom1_atoms = boost::assign::list_of(definitions::HB2)(definitions::HB3); 
                        }

                    } else if (atom1_type == "QD1") { atom1_atoms = boost::assign::list_of(definitions::HD11)(definitions::HD12)(definitions::HD13); }

                    else if (atom1_type == "QD2") { atom1_atoms = boost::assign::list_of(definitions::HD21)(definitions::HD22)(definitions::HD23); }

                    else if (atom1_type == "QG1") { atom1_atoms = boost::assign::list_of(definitions::HG11)(definitions::HG12)(definitions::HG13); }

                    else if (atom1_type == "QG2") { atom1_atoms = boost::assign::list_of(definitions::HG21)(definitions::HG22)(definitions::HG23); }

                    else if (atom1_type == "QQD") { atom1_atoms = boost::assign::list_of(definitions::HD11)(definitions::HD12)(definitions::HD13)
                                                                                        (definitions::HD21)(definitions::HD22)(definitions::HD23); }

                    else if (atom1_type == "QQG") { atom1_atoms = boost::assign::list_of(definitions::HG11)(definitions::HG12)(definitions::HG13)
                                                                                        (definitions::HG21)(definitions::HG22)(definitions::HG23); }

                    else { atom1_atoms = boost::assign::list_of(boost::lexical_cast<definitions::AtomEnum>(split_line[2])); }


                    if (atom2_type == "QD") { atom2_atoms = boost::assign::list_of(definitions::HD1)(definitions::HD2); }

                    else if (atom2_type == "QE") { atom2_atoms = boost::assign::list_of(definitions::HE1)(definitions::HE2); }

                    else if (atom2_type == "QB") { 

                        if (boost::lexical_cast<std::string>(split_line[4]) == "ALA") {
                            atom2_atoms = boost::assign::list_of(definitions::HB1)(definitions::HB2)(definitions::HB3); 
                        } else {
                            atom2_atoms = boost::assign::list_of(definitions::HB2)(definitions::HB3); 
                        }

                    } else if (atom2_type == "QD1") { atom2_atoms = boost::assign::list_of(definitions::HD11)(definitions::HD12)(definitions::HD13); }

                    else if (atom2_type == "QD2") { atom2_atoms = boost::assign::list_of(definitions::HD21)(definitions::HD22)(definitions::HD23); }

                    else if (atom2_type == "QG1") { atom2_atoms = boost::assign::list_of(definitions::HG11)(definitions::HG12)(definitions::HG13); }

                    else if (atom2_type == "QG2") { atom2_atoms = boost::assign::list_of(definitions::HG21)(definitions::HG22)(definitions::HG23); }

                    else if (atom2_type == "QQD") { atom2_atoms = boost::assign::list_of(definitions::HD11)(definitions::HD12)(definitions::HD13)
                                                                                        (definitions::HD21)(definitions::HD22)(definitions::HD23); }

                    else if (atom2_type == "QQG") { atom2_atoms = boost::assign::list_of(definitions::HG11)(definitions::HG12)(definitions::HG13)
                                                                                        (definitions::HG21)(definitions::HG22)(definitions::HG23); }

                    else { atom2_atoms = boost::assign::list_of(boost::lexical_cast<definitions::AtomEnum>(split_line[5])); }

                    // std::cout << atom1_atoms << std::endl;
                    // std::cout << atom2_atoms << std::endl;




                    // std::vector<Atom> atoms1;
                    // std::vector<Atom> atoms2;


                    // std::cout << atom1_type << "  " << atom2_type << std::endl;

                    // int residue_index1 = boost::lexical_cast<int>(split_line[0]);
                    // definitions::AtomEnum atom_type1 = boost::lexical_cast<definitions::AtomEnum>(split_line[1]);
                    // int residue_index2 = boost::lexical_cast<int>(split_line[2]);
                    // definitions::AtomEnum atom_type2 = boost::lexical_cast<definitions::AtomEnum>(split_line[3]);
                    double distance = boost::lexical_cast<double>(split_line[6]);;

                    AmbiguousContact contact(residue1_index,
                                             atom1_atoms,
                                             residue2_index,
                                             atom2_atoms,
                                             distance);



                   if(contact.residue_index1 >= (this->chain)->size() ||
                      contact.residue_index2 >= (this->chain)->size()) {
                       std::cout << "# CONTACT ERROR: " << contact.residue_index1 + 1 << " "
                                                       << contact.atom1_types << " "
                                                       << contact.residue_index2 + 1 << " "
                                                       << contact.atom2_types << std::endl;
                       exit(EXIT_FAILURE);
                   }

                for(unsigned int i = 0; i < contact.atom1_types.size(); i++) {

                   if(!(*(this->chain))[contact.residue_index1].has_atom(contact.atom1_types[i])) {
                       std::cout << "# CONTACT ERROR: " << contact.residue_index1 + 1<< " "
                                                       << contact.atom1_types[i] << std::endl;
                       exit(EXIT_FAILURE);
                   }
                }

                for(unsigned int i = 0; i < contact.atom2_types.size(); i++) {

                   if(!(*(this->chain))[contact.residue_index2].has_atom(contact.atom2_types[i])) {
                       std::cout << "# CONTACT ERROR: " << contact.residue_index2 + 1 << " "
                                                       << contact.atom2_types[i] << std::endl;
                       exit(EXIT_FAILURE);
                   }
                }


                    contact_map.push_back(contact);
                    std::cout << contact << std::endl;

               }
          }
          shuffle_vector(contact_map, this->random_number_engine);
          print_active_contacts(contact_map);
     }

    // http://nmr.cit.nih.gov/xplor-nih/doc/current/xplor/node384.html
    // inline double flat_bottom_potential(double r_actual, double r_equilibrium, double epsilon) {

    //     double R = r_actual;
    //     double d = r_equilibrium;

    //     double exponent = 2.0;
    //     double dplus    = 4.0;
    //     double dminus   = 2.0;

    //     double e = 0.0;

    //     if ((d + dplus ) < R) {
    //         e = R - (d + dplus);
    //     } else if ((d - dminus) < R) {
    //         e = 0.0;
    //     } else {
    //         e = ((d - dminus) - R);
    //     }

    //     e = std::pow(e, exponent) * epsilon;

    //     return e;
    // }

    // https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/de/d50/constraint_file.html
    // http://www.pnas.org/content/suppl/2012/06/25/1203013109.DCSupplemental/Appendix.pdf

    inline double rosetta_bounded_potential(double r_actual, double r_equilibrium, double epsilon) {

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

        return energy * epsilon;

    }

    // inline double lennard_jones_potential(double r_actual, double r_equilibrium, double epsilon) {
    //     double r_ratio = r_equilibrium / r_actual;
    //     return epsilon * (std::pow(r_ratio, 12.0)
    //                         - 2 * std::pow(r_ratio, 6.0));
    // }


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

        contact_map_old = contact_map;

        did_swap = false;

    }



    void print_active_contacts(std::vector<AmbiguousContact> map) {

       //  std::cout << "# THREAD " << this->thread_index
       //        << " CONTACTS  ";

       //  for (unsigned int i=0; i < settings.active_restraints; i++){
       //      Vector_3D r1 = (*(this->chain))(map[i].residue_index1 - 1,
       //                        map[i].atom_type1)->position;
       //      Vector_3D r2 = (*(this->chain))(map[i].residue_index2 - 1,
       //                        map[i].atom_type2)->position;

       //      double r_actual = (r1 - r2).norm();

       //      std::cout << map[i] << "," << r_actual << ",";
       //  }

       //  std::cout << std::endl;
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


    double calc_distance_nearest(AmbiguousContact contact) {

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

    double get_energy_contacts(std::vector<AmbiguousContact> map) {

         double energy = 0.0;

         // for(unsigned int i=0; i<this->contact_map.size(); i++){
         for(unsigned int i=0; i < settings.active_restraints; i++){

            // double r_actual = calc_distance_guntert(map[i]);

            double r_actual = calc_distance_nearest(map[i]);

            double r_equilibrium = map[i].distance;
            double epsilon = settings.weight_constant;

            // energy += lennard_jones_potential(r_actual, r_equilibrium, epsilon);
            // energy += flat_bottom_potential(r_actual, r_equilibrium, epsilon);
            energy += rosetta_bounded_potential(r_actual, r_equilibrium, epsilon);
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
                bias = get_energy_contacts(this->contact_map) - get_energy_contacts(this->contact_map_old);
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
        // print_active_contacts(this->contact_map);
    }

    void reject() {
        if (this->did_swap) this->contact_map = this->contact_map_old;
        // print_active_contacts(this->contact_map);
    }


};

} // end namespace phaistos
#endif
