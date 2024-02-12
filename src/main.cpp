#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <iomanip>
#include <cmath>
#include <unistd.h> // for getopt
#include <getopt.h> // for getopt_long()
#include <sys/stat.h> // for checking if path exists

#include "config.h"
#include "common.h"
#include "parameter.h"
#include "mol.h"
#include "sasa.h"

// inline bool file_exists(const std::string& name) {
//     std::ifstream f(name);
//     return f.good();
// }

inline bool path_exists(const std::string& path_name, const char mode = 'p') {
    struct stat sb;
    if(stat(path_name.c_str(), &sb) == 0) {
        switch(mode) {
            case 'p':
                return true;
                break;
            case 'f':
                return S_ISREG(sb.st_mode);
                break;
            case 'd':
                return S_ISDIR(sb.st_mode);
                break;
            default:
                std::cout << "Unkown path mode: " << mode << std::endl;
                exit(2);
                break;
        }
    }
    return false;
}

//start main//////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
    const std::string usage_hint = "sprank -r in_rec -c in_cpd -p in_poses [-e in_potentials] [-s] [-v] [-h]";
    int stats = 0; // 0 for false
    int verbose = 0; // 0 for false
    int help = 0; // 0 for false

    std::string in_rec = "";
    std::string in_cpd = "";
    std::string in_poses = "";
    std::string in_potentials = "";

    sprank::Float LR_a = 1.0; // linear regression coefficient a
    sprank::Float LR_b = 0.0; // linear regression coefficient b

    /////////////////////////////////////////

    int option_index = 0;
    struct option longopts[] = {
        { "receptor",     required_argument,  NULL,      'r' },
        { "compound",     required_argument,  NULL,      'c' },
        { "poses",        required_argument,  NULL,      'p' },

        { "potentials",   required_argument,  NULL,      'e' },

        { "stats",        no_argument,        & stats,    1  },
        { "verbose",      no_argument,        & verbose,  1  },
        { "help",         no_argument,        & help,     1  },

        {0, 0, 0, 0}
    };

    int c; // switch
    bool args_correct = true;
    // opterr, optopt, optind, optarg are reserved variables for getopt
    while( (c = getopt_long(argc, argv, ":r:c:p:e:svh", longopts, &option_index)) != -1 ) {
        switch(c) {
            case 'r':
                in_rec = std::string(optarg);
                break;
            case 'c':
                in_cpd = std::string(optarg);
                break;
            case 'p':
                in_poses = std::string(optarg);
                break;
            case 'e':
                in_potentials = std::string(optarg);
                break;
            case 's':
                stats = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                help = true;
                break;
            case 0: {
                /* getopt_long() set a variable, just keep going */
                const std::string arg_name = std::string(longopts[option_index].name);
                // std::cout << arg_name << std::endl;
                // std::cout << longopts[option_index].flag << std::endl;
                if(arg_name == "stats" || arg_name == "verbose" || arg_name == "help") {}
                break;
            }
            case ':':
                if (optopt == 'r' || optopt == 'c' || optopt == 'p' || optopt == 'e') {
                    std::cout << "Option -" << char(optopt) << " requires an argument." << std::endl;
                }
                args_correct = false;
                break;
            case '?':
            default:
                if (std::isprint(optopt)) {
                    std::cout << "Invalid option `-" << char(optopt) << "'." << std::endl;
                } else {
                    std::cout << "Invalid option character `\\x" << char(optopt) << "'." << std::endl;
                }
                args_correct = false;
                break;
        }
    }

    // no options are passed, or args are not correct
    if(optind == 1 || !args_correct) {
        std::cout << "Usage: " << usage_hint << std::endl;
        return 2;
    }

    std::cout << "################################################" << std::endl;
    std::cout << "##                   SPRank                   ##" << std::endl;
    std::cout << "################################################" << std::endl;

    if(help) {
        std::cout << "Required:" << std::endl;
        std::cout << "  --r arg    receptor with amber RNA.OL3 atom types (mol2)" << std::endl;
        std::cout << "  --l arg    compound with amber gaff2 atom types and bond table (mol2)" << std::endl;
        std::cout << "  --p arg    compound poses to be scored (mol2)" << std::endl;
        std::cout << "               ->This file shouldn't contain any hydrogens. Only heavy atoms' coordinates are used." << std::endl;
        std::cout << "               ->The order of heavy atoms must be same as the order in input compound." << std::endl;
        std::cout << std::endl;

        std::cout << "Optional:" << std::endl;
        std::cout << "  --e arg    user-provided pairwise potentials file" << std::endl;
        std::cout << "  --s        flag to output contact statistics" << std::endl;
        std::cout << "  --v        flag to output verbose info" << std::endl;

        std::cout << std::endl;
        std::cout << "Example folder structure:" << std::endl;
        std::cout << "  --|work_dir/" << std::endl;
        std::cout << "  ----|rec.mol2" << std::endl;
        std::cout << "  ----|cpd.mol2" << std::endl;
        std::cout << "  ----|poses.mol2" << std::endl;
        return 2;
    }

    for (int index = optind; index < argc; index++) {
        std::cout << "Not used argument: " << argv[index] << std::endl;
    }

    if(!path_exists(in_rec))   { std::cout << "in_rec: "    << in_rec      << " not exists!" << std::endl; exit(2); }
    if(!path_exists(in_cpd))   { std::cout << "in_cpd: "    << in_cpd      << " not exists!" << std::endl; exit(2); }
    if(!path_exists(in_poses)) { std::cout << "in_poses: "  << in_poses    << " not exists!" << std::endl; exit(2); }
    if(in_potentials != "" && !path_exists(in_potentials))  { std::cout << "in_potentials: "  << in_potentials  << " not exists!" << std::endl; exit(2); }

    std::cout << "in_rec:          " << in_rec          << std::endl;
    std::cout << "in_cpd:          " << in_cpd          << std::endl;
    std::cout << "in_poses:        " << in_poses        << std::endl;
    if(in_potentials != "") {
        std::cout << "in_potentials:   " << in_potentials   << std::endl;
    }
    std::cout << "stats:           " << stats           << std::endl;
    std::cout << "verbose:         " << verbose         << std::endl;

    sprank::Parameter param(in_potentials);

    if(in_potentials != "") {
        std::cout << "Load potentials from: " << in_potentials << std::endl;
    }

    sprank::Receptor rec(param, in_rec);
    std::cout << "Load receptor: "      << rec.name << " from " << in_rec << std::endl;
    std::cout << "--> heavy atom num: " << rec.atoms.size() << std::endl;

    sprank::Compound cpd(param, in_cpd);
    std::cout << "Load compound: "      << cpd.name << " from " << in_cpd << std::endl;
    std::cout << "--> heavy atom num: " << cpd.atoms.size() << std::endl;

    param.info();
    rec.info(verbose);
    cpd.info(verbose);

    std::vector<std::string> score_strs;
    std::vector<std::string> stats_strs;
    std::ostringstream score_oss;
    std::ostringstream stats_oss;

    score_oss << "************************************************" << std::endl;
    score_oss << "*                    Scoring                   *" << std::endl;
    score_oss << "************************************************" << std::endl;
    score_oss << "#Scoring ---------------------------------------" << std::endl;
    score_oss << "#pose pair sasa stack" << std::endl;
    score_strs.push_back(score_oss.str());

    stats_oss << "************************************************" << std::endl;
    stats_oss << "*                  Statistics                  *" << std::endl;
    stats_oss << "************************************************" << std::endl;
    stats_oss << "#Statistics ------------------------------------" << std::endl;
    stats_oss << "#pose rec_atom_type cpd_atom_type radial_dist..." << std::endl;
    stats_strs.push_back(stats_oss.str());

    const std::map<std::string,std::string> RBBRBA_map = {
        {"CT", "RBB"}, {"CI", "RBB"}, {"OS", "RBB"}, {"O2", "RBB"}, {"P" , "RBB"}, {"OH", "RBB"},
        {"O" , "RBA"}, {"C" , "RBA"}, {"CB", "RBA"}, {"N2", "RBA"}, {"NC", "RBA"}, {"CA", "RBA"},
        {"N*", "RBA"}, {"NA", "RBA"}, {"C4", "RBA"}, {"NB", "RBA"}, {"CS", "RBA"}, {"CP", "RBA"},
        {"CQ", "RBA"}, {"C5", "RBA"}
    };
    // const std::set<std::string> rec_aroma_set = {"C", "CA", "CB", "CP", "CQ", "CS", "C4", "C5", "NA", "NB", "NC", "N*"};
    const std::set<std::string> rec_aromaC_set = {"C", "CA", "CB", "CP", "CQ", "CS", "C4", "C5"};
    // const std::set<std::string> cpd_aroma_set = {"C.ca", "C.cc", "N.nb", "N.nc"};
    const std::set<std::string> cpd_aromaC_set = {"C.ca", "C.cc"};

    std::ifstream in_poses_f(in_poses, std::ios::in);
    while(true) {
        score_oss.str(std::string());//clear content
        stats_oss.str(std::string());//clear content
        score_oss.clear();//clear state
        stats_oss.clear();//clear state
        std::string pose_name;
        std::vector<sprank::Atom> pose_atoms;
        sprank::parse_single_pose_from_stream_each_time(in_poses_f, pose_name, pose_atoms);

        if(pose_name == "" && pose_atoms.size() == 0 && in_poses_f.eof()) { break; }

        assert(cpd.atoms.size() == pose_atoms.size());
        for(sprank::Size_Type i = 0; i < cpd.atoms.size(); ++i) {
            cpd.atoms[i].xyz = pose_atoms[i].xyz;
        }

        //cal intermolecular pairwise energy terms and sasa
        std::map<std::string,sprank::Floats> pairwise_statistics_map;
        if(stats == true) {
            for(const auto& ppm : param.pairwise_potentials_map) {
                pairwise_statistics_map[ppm.first] = sprank::Floats(param.pairwise_bin_num, 0.0);
            }
        }

        std::vector<sprank::Atom> pocket_atoms;
        sprank::Float inter_pairwise_energy = 0.0;
        // sprank::Float inter_pairwise_aroma_energy = 0.0;
        sprank::Float inter_pairwise_aromaC_energy = 0.0;
        for(const auto& a : rec.atoms) {
            bool pocket_flag = false;
            for(const auto& b : cpd.atoms) {
                const sprank::Float dis = std::sqrt(vec3d_distance_square(a.xyz, b.xyz));
                if(dis < param.pairwise_dis_upper_bound) {
                    pocket_flag = true;
                    const sprank::Size_Type i_bin = dis/param.pairwise_bin_size;
                    const sprank::Float& bin_vol = param.bin_volumes[i_bin];

                    std::string rec_at = a.atom_type;
                    if(b.atom_type == "P.*" || b.atom_type == "S.*" || b.atom_type == "R.cn" || b.atom_type == "HA.*") {
                        const auto it = RBBRBA_map.find(rec_at);
                        if(it != RBBRBA_map.end()) {
                            rec_at = it->second;
                        }
                    }
                    const std::string p_name = rec_at+"-"+b.atom_type;

                    if(param.pairwise_potentials_map.find(p_name) != param.pairwise_potentials_map.end()) {
                        inter_pairwise_energy += param.pairwise_potentials_map[p_name][i_bin];

                        // if((rec_aroma_set.find(rec_at) != rec_aroma_set.end()) and (cpd_aroma_set.find(b.atom_type) != cpd_aroma_set.end())) {
                        //     inter_pairwise_aroma_energy += param.pairwise_potentials_map[p_name][i_bin];
                        // }
                        if((rec_aromaC_set.find(rec_at) != rec_aromaC_set.end()) and (cpd_aromaC_set.find(b.atom_type) != cpd_aromaC_set.end())) {
                            inter_pairwise_aromaC_energy += param.pairwise_potentials_map[p_name][i_bin];
                        }

                        if(stats == true) {
                            pairwise_statistics_map[p_name][i_bin] += 1.0;
                        }
                    }
                    //// deal with S
                    //if(b.atom_type == "S.*") {
                    //    const std::string &mt = b.mol2_type;
                    //    if(mt == "ss" or mt == "sh" or mt == "sp" or mt == "sq") {
                    //        p_name = a.atom_type+"-O.os";
                    //    } else if(mt == "s") {}
                    //}
                    //if(param.pairwise_potentials_map.find(p_name) != param.pairwise_potentials_map.end()) {
                    //    inter_pairwise_energy += param.pairwise_potentials_map[p_name][i_bin];
                    //    if(stats == true) {
                    //        pairwise_statistics_map[p_name][i_bin] += 1.0;
                    //    }
                    //}
                }
            }
            if(pocket_flag == true) { pocket_atoms.push_back(a); }
        }
        inter_pairwise_energy = param.rt_to_kcal_per_mol*inter_pairwise_energy; // 1 RT = 0.593 kcal/mol at T=298K
        // inter_pairwise_aroma_energy = param.rt_to_kcal_per_mol*inter_pairwise_aroma_energy; // 1 RT = 0.593 kcal/mol at T=298K
        inter_pairwise_aromaC_energy = param.rt_to_kcal_per_mol*inter_pairwise_aromaC_energy; // 1 RT = 0.593 kcal/mol at T=298K

        // cal sasa
        std::vector<std::vector<sprank::Size_Type>> sasa_neighboring_rec_indices_for_rec_atoms = std::vector<std::vector<sprank::Size_Type>>(rec.atoms.size());
        std::vector<std::vector<sprank::Size_Type>> sasa_neighboring_cpd_indices_for_cpd_atoms = std::vector<std::vector<sprank::Size_Type>>(cpd.atoms.size());
        std::vector<std::vector<sprank::Size_Type>> sasa_neighboring_cpd_indices_for_rec_atoms = std::vector<std::vector<sprank::Size_Type>>(rec.atoms.size());
        std::vector<std::vector<sprank::Size_Type>> sasa_neighboring_rec_indices_for_cpd_atoms = std::vector<std::vector<sprank::Size_Type>>(cpd.atoms.size());
        // for rec
        for(const auto& a : pocket_atoms) {
            for(const auto& b : pocket_atoms) {
                if(a.index < b.index) {
                    const sprank::Float dis = std::sqrt(vec3d_distance_square(a.xyz, b.xyz));
                    if(dis < (a.radius+b.radius+2*param.r_water)) {
                        sasa_neighboring_rec_indices_for_rec_atoms[a.index].push_back(b.index);
                        sasa_neighboring_rec_indices_for_rec_atoms[b.index].push_back(a.index);
                    }
                }
            }
        }
        // for cpd
        for(const auto& a : cpd.atoms) {
            for(sprank::Size_Type b_index = 0; b_index < a.index; b_index++) {
                const sprank::Atom& b = cpd.atoms[b_index];
                const sprank::Float dis = std::sqrt(vec3d_distance_square(a.xyz, b.xyz));
                if(dis < (a.radius+b.radius+2*param.r_water)) {
                    sasa_neighboring_cpd_indices_for_cpd_atoms[a.index].push_back(b.index);
                    sasa_neighboring_cpd_indices_for_cpd_atoms[b.index].push_back(a.index);
                }
            }
        }
        // for rec+cpd
        for(const auto& a : pocket_atoms) {
            for(const auto& b : cpd.atoms) {
                const sprank::Float dis = std::sqrt(vec3d_distance_square(a.xyz, b.xyz));
                if(dis < (a.radius+b.radius+2*param.r_water)) {
                    sasa_neighboring_cpd_indices_for_rec_atoms[a.index].push_back(b.index);
                    sasa_neighboring_rec_indices_for_cpd_atoms[b.index].push_back(a.index);
                }
            }
        }

        // cpd sasa before
        const sprank::Float cpd_sasa_before = sprank::cal_mol_sasa(param, cpd, sasa_neighboring_cpd_indices_for_cpd_atoms);
        // rec delta sasa upon binding
        const sprank::Float delta_rec_sasa = sprank::cal_delta_rec_sasa(param, rec, cpd, sasa_neighboring_rec_indices_for_rec_atoms, sasa_neighboring_cpd_indices_for_rec_atoms);
        //cal sasa after
        const sprank::Float cpd_sasa_after = sprank::cal_mol_sasa(param, cpd, rec, sasa_neighboring_cpd_indices_for_cpd_atoms, sasa_neighboring_rec_indices_for_cpd_atoms);
        const sprank::Float delta_cpd_sasa = cpd_sasa_after - cpd_sasa_before;
        const sprank::Float delta_sasa_energy = param.gamma_SASA*(delta_rec_sasa + delta_cpd_sasa);

        const sprank::Float total_energy = LR_a*inter_pairwise_energy + LR_b;

        score_oss << pose_name << " ";
        score_oss << std::fixed << std::setprecision(3) << total_energy << " " << delta_sasa_energy << " " << inter_pairwise_aromaC_energy << std::endl;
        score_strs.push_back(score_oss.str());

        if(stats == true) {
            stats_oss << std::fixed << std::setprecision(3);
            for(const auto &psm : pairwise_statistics_map) {
                const std::string &p_name = psm.first;
                const sprank::Floats &p_stats = psm.second;
                stats_oss << pose_name << " " << std::setw(5) << p_name;
                for(const auto &s_v : p_stats) {
                    stats_oss << " " << s_v;
                }
                stats_oss << std::endl;
            }
            stats_strs.push_back(stats_oss.str());
        }

        if(verbose) {
            std::cout << "************************************************" << std::endl;
            std::cout << "*                 VERBOSE INFO                 *" << std::endl;
            std::cout << "************************************************" << std::endl;
            std::cout  << "name                 " << pose_name    << std::endl;
            std::cout  << "total_energy:        " << total_energy << std::endl;
        }
    }
    in_poses_f.close();

    for(const auto& out_str : score_strs) {
        std::cout << out_str;
    }
    if(stats == true) {
        for(const auto& out_str : stats_strs) {
            std::cout << out_str;
        }
    }

    return 0;
}//end main
