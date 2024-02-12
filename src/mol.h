#pragma once

#include <fstream>
#include <string>
#include <cassert>
#include <map>
#include <vector>
#include <set>
#include "common.h"
#include "parameter.h"
#include "str.h"
#include "vec3d.h"
#include "print.h"

namespace sprank {

struct Atom {
    Size_Type index;
    int serial;
    std::string name;
    Vec3d xyz;
    std::string mol2_type;
    std::string elem_type;
    ELEM_ENUM elem_enum;
    int res_seq;
    std::string res_name;
    Float charge;
    Float radius;
    std::string atom_type;
    Size_Type num_h = 0;
    Size_Type degree = 0; //heavy atom degree
    std::vector<Size_Type> bonded_atom_ids; // bonded heavy atom ids
    std::vector<std::string> bond_types; // bond types

    //bool is_acc = false;
    //bool is_don = false;
    //bool is_aromatic = false;
};

void parse_single_pose_from_stream_each_time(std::ifstream& pose_stream, std::string& name, std::vector<Atom>& atoms);
// read atoms from mol2
void parse_single_mol2_file(const Parameter& param, const std::string file_path, std::string& name, std::vector<Atom>& atoms, const std::string mol_type);

struct Receptor {
    std::string name;
    std::string file_path;
    std::vector<Atom> atoms;
    Receptor(const Parameter& param, const std::string fp) : file_path(fp) {
        parse_single_mol2_file(param, this->file_path, this->name, this->atoms, "rec");
    }
    void info(bool verbose) const;
};

struct Compound {
    std::string name;
    std::string file_path;
    std::vector<Atom> atoms;
    //Float torsion_num = 0.0; // assigned in initialization
    Compound(const Parameter& param, const std::string fp) : file_path(fp) {
        parse_single_mol2_file(param, this->file_path, this->name, this->atoms, "cpd");
    }
    void info(bool verbose) const;
};

}
