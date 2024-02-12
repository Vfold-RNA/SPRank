#include <fstream>
#include <string>
#include <cassert>
#include <map>

#include "mol.h"

namespace sprank {

const std::map<std::string,std::string> mol2_type_to_elem_map = {
    // rec
    {"H" , "H"}, {"HC", "H"}, {"H1", "H"}, {"H2", "H"}, {"H3", "H"}, {"HA", "H"}, {"H4", "H"}, {"H5", "H"}, {"HO", "H"}, {"HS", "H"},
    {"HP", "H"}, {"HZ", "H"},
    {"C" , "C"}, {"CA", "C"}, {"CB", "C"}, {"CC", "C"}, {"CD", "C"}, {"CI", "C"}, {"CK", "C"}, {"CP", "C"}, {"CM", "C"}, {"CS", "C"},
    {"CN", "C"}, {"CQ", "C"}, {"CR", "C"}, {"CT", "C"}, {"CV", "C"}, {"CW", "C"}, {"C*", "C"}, {"CX", "C"}, {"CY", "C"}, {"CZ", "C"},
    {"C5", "C"}, {"C4", "C"}, {"C0", "C"},
    {"N" , "N"}, {"NA", "N"}, {"NB", "N"}, {"NC", "N"}, {"N2", "N"}, {"N3", "N"}, {"NT", "N"}, {"N*", "N"}, {"NY", "N"},
    {"O" , "O"}, {"O2", "O"}, {"OH", "O"}, {"OS", "O"}, {"OP", "O"},
    {"P" , "P"},
    {"S" , "S"}, {"SH", "S"},
    {"F" , "F" },
    {"Cl", "CL"},
    {"Br", "BR"},
    {"I" , "I" },
    {"MG", "MG"},
    {"CU", "CU"},
    {"FE", "FE"},
    // cpd
    {"h2", "H"}, {"h3", "H"}, {"h4", "H"}, {"h5", "H"}, {"ha", "H"}, {"hc", "H"}, {"hn", "H"}, {"ho", "H"}, {"hp", "H"}, {"hs", "H"},
    {"hw", "H"}, {"hx", "H"},
    {"c" , "C"}, {"cs", "C"}, {"c1", "C"}, {"c2", "C"}, {"c3", "C"}, {"ca", "C"}, {"cp", "C"}, {"cq", "C"}, {"cc", "C"}, {"cd", "C"},
    {"ce", "C"}, {"cf", "C"}, {"cg", "C"}, {"ch", "C"}, {"cx", "C"}, {"cy", "C"}, {"cu", "C"}, {"cv", "C"}, {"cz", "C"}, {"h1", "H"},
    {"n" , "N"}, {"n1", "N"}, {"n2", "N"}, {"n3", "N"}, {"n4", "N"}, {"na", "N"}, {"nb", "N"}, {"nc", "N"}, {"nd", "N"}, {"ne", "N"},
    {"nf", "N"}, {"nh", "N"}, {"no", "N"}, {"ns", "N"}, {"nt", "N"}, {"nx", "N"}, {"ny", "N"}, {"nz", "N"}, {"n+", "N"}, {"nu", "N"},
    {"nv", "N"}, {"n7", "N"}, {"n8", "N"}, {"n9", "N"}, {"ni", "N"}, {"nj", "N"}, {"nk", "N"}, {"nl", "N"}, {"nm", "N"}, {"nn", "N"},
    {"np", "N"}, {"nq", "N"}, {"n5", "N"}, {"n6", "N"},
    {"o" , "O"}, {"oh", "O"}, {"os", "O"}, {"op", "O"}, {"oq", "O"}, {"ow", "O"},
    {"p2", "P"}, {"p3", "P"}, {"p4", "P"}, {"p5", "P"}, {"pb", "P"}, {"pc", "P"}, {"pd", "P"}, {"pe", "P"}, {"pf", "P"}, {"px", "P"},
    {"py", "P"},
    {"s" , "S"}, {"s2", "S"}, {"s4", "S"}, {"s6", "S"}, {"sh", "S"}, {"ss", "S"}, {"sp", "S"}, {"sq", "S"}, {"sx", "S"}, {"sy", "S"},
    {"f" , "F"},
    {"cl", "CL"},
    {"br", "BR"},
    {"i" , "I"},
    {"Cu", "Cu"}
};

void parse_single_pose_from_stream_each_time(std::ifstream& pose_stream, std::string& name, std::vector<Atom>& atoms) {
    std::string sline;
    Size_Type number_atoms, number_bonds;
    std::vector<Atom> atoms_temp;
    std::string name_temp;
    while(std::getline(pose_stream, sline)) {
        if(sline.find("@<TRIPOS>MOLECULE") != std::string::npos) {
            std::getline(pose_stream, sline);
            name_temp = trim(sline);
            std::getline(pose_stream, sline);
            const std::vector<std::string> token = string2vector(sline);
            assert(token.size() >= 2);
            number_atoms = std::stoi(token[0]);
            number_bonds = std::stoi(token[1]);
        }
        if(sline.find("@<TRIPOS>ATOM") != std::string::npos) {
            Size_Type heavy_index = 0;
            for(Size_Type index = 0; index < number_atoms; ++index) {
                std::getline(pose_stream, sline);
                const std::vector<std::string> token = string2vector(sline);
                //if(token.size() != 9) {
                //    std::cout << "Error: failed to parse the following mol2 line (9 columns expected):" << std::endl;
                //    std::cout << sline << std::endl;
                //    exit(2);
                //}
                Atom mol2_atom;
                // mol2_atom.index = index;
                mol2_atom.serial = std::stoi(token[0]);
                mol2_atom.name = token[1];
                mol2_atom.xyz = Vec3d(std::stod(token[2]), std::stod(token[3]), std::stod(token[4]));
                //mol2_atom.mol2_type = token[5];
                //mol2_atom.elem_type = mol2_type_to_elem_map.at(mol2_atom.mol2_type);
                //mol2_atom.res_seq = std::stoi(token[6]);
                //mol2_atom.res_name = token[7];
                //mol2_atom.charge = std::stod(token[8]);
                //if(mol2_atom.elem_type != "H") {
                mol2_atom.index = heavy_index;
                atoms_temp.push_back(mol2_atom);
                heavy_index++;
                //}
            }// for atom
        }

        if(atoms_temp.size() != 0) {
            name = name_temp;
            atoms = atoms_temp;
            // assert(atoms.size() == number_atoms);
            break; // break the while for this read
        }
    }
}

// read heavy atoms from mol2
void parse_single_mol2_file(const Parameter& param, const std::string file_path, std::string& name, std::vector<Atom>& atoms, const std::string mol_type) {
    // read all atom part
    std::ifstream inmol2(file_path, std::ios::in);
    std::string sline;
    Size_Type number_atoms, number_bonds;
    std::vector<Atom> atoms_temp;
    std::map<Size_Type,Size_Type> index_to_heavy_index;
    while(std::getline(inmol2, sline)) {
        if(sline.find("@<TRIPOS>MOLECULE") != std::string::npos) {
            std::getline(inmol2, sline);
            name = trim(sline);
            std::getline(inmol2, sline);
            const std::vector<std::string> token = string2vector(sline);
            assert(token.size() >= 2);
            number_atoms = std::stoi(token[0]);
            number_bonds = std::stoi(token[1]);
        }
        if(sline.find("@<TRIPOS>ATOM") != std::string::npos) {
            Size_Type heavy_index = 0;
            for(Size_Type index = 0; index < number_atoms; ++index) {
                std::getline(inmol2, sline);
                const std::vector<std::string> token = string2vector(sline);
                if(mol_type == "cpd") {
                    if(token.size() != 9) {
                        std::cout << "Error: failed to parse the following mol2 line (9 columns expected):" << std::endl;
                        std::cout << sline << std::endl;
                        exit(2);
                    }
                } else if(mol_type == "rec") {
                    if(token.size() != 10) {
                        std::cout << "Error: failed to parse the following mol2 line (10 columns expected):" << std::endl;
                        std::cout << sline << std::endl;
                        exit(2);
                    }
                } else {
                    std::cout << "Error: unknown mol_type: " << mol_type << std::endl;
                    exit(2);
                }
                Atom mol2_atom;
                mol2_atom.index = index;
                mol2_atom.serial = std::stoi(token[0]);
                mol2_atom.name = token[1];
                mol2_atom.xyz = Vec3d(std::stod(token[2]), std::stod(token[3]), std::stod(token[4]));
                mol2_atom.mol2_type = token[5];
                mol2_atom.elem_type = mol2_type_to_elem_map.at(mol2_atom.mol2_type);
                mol2_atom.res_seq = std::stoi(token[6]);
                mol2_atom.res_name = token[7];
                mol2_atom.charge = std::stod(token[8]);
                if(mol2_atom.elem_type != "H") {
                    index_to_heavy_index.insert({index, heavy_index});
                    heavy_index++;
                }
                atoms_temp.push_back(mol2_atom);
            }
        }
        // read bond, update degree, update num_h and merge hydrogen charge
        if(sline.find("@<TRIPOS>BOND") != std::string::npos) {
            for(Size_Type ib = 0; ib != number_bonds; ++ib) {
                std::getline(inmol2, sline);
                const std::vector<std::string> token = string2vector(sline);
                const Size_Type left_idx = std::stoi(token[1])-1;
                const Size_Type right_idx = std::stoi(token[2])-1;

                // get bonded atoms and heavy atom degree
                if(atoms_temp[left_idx].elem_type != "H" && atoms_temp[right_idx].elem_type != "H") {
                    atoms_temp[right_idx].degree++;
                    atoms_temp[left_idx].degree++;
                    atoms_temp[right_idx].bonded_atom_ids.push_back(left_idx);
                    atoms_temp[left_idx].bonded_atom_ids.push_back(right_idx);
                    atoms_temp[right_idx].bond_types.push_back(token[3]);
                    atoms_temp[left_idx].bond_types.push_back(token[3]);
                }
                // merge h
                if(atoms_temp[left_idx].elem_type == "H" && atoms_temp[right_idx].elem_type != "H") {
                    atoms_temp[right_idx].charge += atoms_temp[left_idx].charge;
                    atoms_temp[left_idx].charge = 0.0;
                    atoms_temp[left_idx].degree++;
                    atoms_temp[right_idx].num_h++;
                } else if(atoms_temp[left_idx].elem_type != "H" && atoms_temp[right_idx].elem_type == "H") {
                    atoms_temp[left_idx].charge += atoms_temp[right_idx].charge;
                    atoms_temp[right_idx].charge = 0.0;
                    atoms_temp[right_idx].degree++;
                    atoms_temp[left_idx].num_h++;
                }
            }
        }
    }// end reading
    inmol2.close();

    //// update aromatic info
    //for(const auto& aromatic_idx : aromatic_indices) {
    //    atoms_temp[aromatic_idx].is_aromatic = true;
    //}
    //// update h bond info
    //for(const auto& h_don_idx : h_don_indices) {
    //    atoms_temp[h_don_idx].is_don = true;
    //}
    //for(const auto& h_acc_idx : h_acc_indices) {
    //    atoms_temp[h_acc_idx].is_acc = true;
    //}
    
    // assign elem_enum and radius
    for(auto& a : atoms_temp) {
        if(a.elem_type == "H") { a.elem_enum = EE_H; a.radius = elem_radii[EE_H]; }
        else if(a.elem_type == "C") { a.elem_enum = EE_C; a.radius = elem_radii[EE_C]; }
        else if(a.elem_type == "N") { a.elem_enum = EE_N; a.radius = elem_radii[EE_N]; }
        else if(a.elem_type == "O") { a.elem_enum = EE_O; a.radius = elem_radii[EE_O]; }
        else if(a.elem_type == "P") { a.elem_enum = EE_P; a.radius = elem_radii[EE_P]; }
        else if(a.elem_type == "S") { a.elem_enum = EE_S; a.radius = elem_radii[EE_S]; }
        else if(a.elem_type == "F" || a.elem_type == "CL" || a.elem_type == "BR" || a.elem_type == "I") { a.elem_enum = EE_HA; a.radius = elem_radii[EE_HA]; }
        else { a.elem_enum = EE_UNK; a.radius = elem_radii[EE_UNK]; }
    }

    // assign atom type
    if(mol_type == "rec") {
        for(auto& a : atoms_temp) {
            a.atom_type = a.mol2_type;
        }
    } else if(mol_type == "cpd") {
        for(auto& a : atoms_temp) {
            const std::string &et = a.elem_type;
            const std::string &mt = a.mol2_type;
            std::string &at = a.atom_type;
            at = "unk";
            if(et == "H") {
                at = "H.*";
            } else if(et == "C") {
                if(mt == "c3" || mt == "cx" || mt == "cy") {
                    at = "C.c3";
                } else if(mt == "ca") {
                    at = "C.ca";
                } else if(mt == "cc" || mt == "cd") {
                    at = "C.cc";
                } else if(mt == "c" || mt == "cs") {
                    at = "C.c";
                } else if(mt == "c2" || mt == "cp" || mt == "cq" || mt == "ce" || mt == "cf" || mt == "cu" || mt == "cv" || mt == "cz" ) {
                    at = "C.c2";
                } else if(mt == "c1" || mt == "cg" || mt == "ch") {
                    //at = "C.c1";
                    at = "R.cn";
                }
            } else if(et == "N") {
                if(mt == "n4" || mt == "nx" || mt == "ny" || mt == "nz" || mt == "n+" || mt == "nk" || mt == "nl") {
                    at = "N.n4";
                } else if(mt == "nb") {
                    at = "N.nb";
                } else if(mt == "nh" || mt == "nv" || mt == "nu" || mt == "nm" || mt == "nn") {
                    //assert((a.degree+a.num_h) == 3 || (a.degree+a.num_h) == 4);
                    //if((a.degree+a.num_h) == 3) {
                    //    at = "N.nh3";
                    //} else if((a.degree+a.num_h) == 4) {
                    //    at = "N.nh4";
                    //}
                    at = "N.nh";
                } else if(mt == "na") {
                    at = "N.na";
                } else if(mt == "nc" || mt == "nd") {
                    at = "N.nc";
                } else if(mt == "n" || mt == "ns" || mt == "nt" || mt == "ni" || mt == "nj") {
                    at = "N.n";
                } else if(mt == "n3" || mt == "n5" || mt == "n6" || mt == "n7" || mt == "n8" || mt == "n9" || mt == "np" || mt == "nq") {
                    at = "N.n3";
                } else if(mt == "n2" || mt == "ne" || mt == "nf") {
                    at = "N.n2";
                } else if(mt == "n1") {
                    //at = "N.n1";
                    at = "R.cn";
                }
            } else if(et == "O") {
                if(mt == "os" || mt == "op" || mt == "oq") {
                    at = "O.os";
                } else if(mt == "oh") {
                    at = "O.oh";
                } else if(mt == "o") {
                    //assert(a.bond_types.size() == 1);
                    if(a.bond_types.size() == 1) {
                        if(a.bond_types[0] == "1") {
                            at = "O.o1";
                        } else if(a.bond_types[0] == "2") {
                            at = "O.o2";
                        }
                    }
                }
            } else if(et == "P") {
                at = "P.*";
            } else if(et == "S") {
                //if(mt == "ss" | mt == "sh" | mt == "sp" | mt == "sq") {
                //    at = "O.os";
                //} else if(mt == "s") {
                //    assert(a.bond_types.size() == 1);
                //    if(a.bond_types.size() == 1) {
                //        if(a.bond_types[0] == "1") {
                //            at = "O.o1";
                //        } else if(a.bond_types[0] == "2") {
                //            at = "O.o2";
                //        }
                //    }
                //} else if(mt == "s4" || mt == "sx") {
                //    at = "S.s4";
                //} else if(mt == "s2" || mt == "s6" || mt == "sy") {
                //    at = "S.*";
                //}
                at = "S.*";
            } else if(mt == "f" || mt == "cl" || mt == "br" || mt == "i") {
                at = "HA.*";
            }
        }
    }

    Size_Type updated_index = 0;
    // extract heavy atoms from atoms_temp, and update index, and bonded ids
    for(const auto& a : atoms_temp) {
        if(a.elem_type != "H") {
            atoms.push_back(a);
            auto& new_atom = atoms[updated_index];
            new_atom.index = index_to_heavy_index[a.index];
            new_atom.bonded_atom_ids.clear();
            for(const Size_Type b_id : a.bonded_atom_ids) {
                new_atom.bonded_atom_ids.push_back(index_to_heavy_index[b_id]);
            }
            new_atom.bond_types.clear();
            new_atom.bond_types = a.bond_types;
            updated_index++;
        }
    }
}//end read mol2

void Receptor::info(bool verbose) const {
    std::cout << "************************************************" << std::endl;
    std::cout << "*                   REC INFO                   *" << std::endl;
    std::cout << "************************************************" << std::endl;
    std::cout << "name: " << this->name << std::endl;
    std::cout << "file_path: " << this->file_path << std::endl;
    std::cout << "atoms num: " << this->atoms.size() << std::endl;

    if(verbose) {
        for(const auto a : this->atoms) {
            std::cout << "-----------------------------------" << std::endl;
            std::cout << "serial: " << a.serial << " | index: " << a.index << " | name: " << a.name << " | res_seq: " << a.res_seq << " | res_name: " << a.res_name << " | charge: " << a.charge << " | radius: " << a.radius << std::endl;
            std::cout << "xyz: "; print(a.xyz); std::cout << std::endl;
            std::cout << "mol2_type: " << a.mol2_type << " | elem_type: " << a.elem_type << " | atom_type: " << a.atom_type << " | num_h: " << a.num_h << " | degree: " << a.degree << std::endl;
            std::cout << "bonded_atom_ids: "; print(a.bonded_atom_ids); std::cout << std::endl;
            std::cout << "bond_types: "; print(a.bond_types); std::cout << std::endl;
            std::vector<std::string> bonded_atom_names;
            for(const auto &ba_id : a.bonded_atom_ids) {
                bonded_atom_names.push_back(this->atoms[ba_id].name);
            }
            std::cout << "bonded_atom_names: "; print(bonded_atom_names); std::cout << std::endl;
        }
    }
}

void Compound::info(bool verbose) const {
    std::cout << "************************************************" << std::endl;
    std::cout << "*                   CPD INFO                   *" << std::endl;
    std::cout << "************************************************" << std::endl;
    std::cout << "name: " << this->name << std::endl;
    std::cout << "file_path: " << this->file_path << std::endl;
    std::cout << "atoms num: " << this->atoms.size() << std::endl;

    if(verbose) {
        for(const auto a : this->atoms) {
            std::cout << "-----------------------------------" << std::endl;
            std::cout << "serial: " << a.serial << " | index: " << a.index << " | name: " << a.name << " | res_seq: " << a.res_seq << " | res_name: " << a.res_name << " | charge: " << a.charge << " | radius: " << a.radius << std::endl;
            std::cout << "xyz: "; print(a.xyz); std::cout << std::endl;
            std::cout << "mol2_type: " << a.mol2_type << " | elem_type: " << a.elem_type << " | atom_type: " << a.atom_type << " | num_h: " << a.num_h << " | degree: " << a.degree << std::endl;
            std::cout << "bonded_atom_ids: "; print(a.bonded_atom_ids); std::cout << std::endl;
            std::cout << "bond_types: "; print(a.bond_types); std::cout << std::endl;
            std::vector<std::string> bonded_atom_names;
            for(const auto &ba_id : a.bonded_atom_ids) {
                bonded_atom_names.push_back(this->atoms[ba_id].name);
            }
            std::cout << "bonded_atom_names: "; print(bonded_atom_names); std::cout << std::endl;
        }
    }
}

}
