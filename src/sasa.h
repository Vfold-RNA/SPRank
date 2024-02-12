#pragma once

#include <iostream>
#include <vector>
//#include <array>
//#include <map>
#include <string>
//#include <cmath>
#include "common.h"
#include "parameter.h"
//#include "matrix.h"
#include "mol.h"

namespace sprank {

    template<typename T>
    const Float cal_mol_sasa(const Parameter& param, const T& target_mol, const std::vector<std::vector<Size_Type>>& intra_near_indices_list) {
        Float SumSASAtmp = 0.0;
        for(const auto& a : target_mol.atoms) {
            const auto& intra_near_indices = intra_near_indices_list[a.index];
            if(intra_near_indices.size() != 0) {
                Size_Type effective_count = 0;
                for(const auto& point : param.sasa_points[a.elem_enum]) {
                    const Vec3d surface_point = point + a.xyz;
                    bool buried = false;
                    for(Size_Type l=0; l != intra_near_indices.size(); l++) {
                        const Atom& intra_near_atom = target_mol.atoms[intra_near_indices[l]];
                        const Float dis_square = vec3d_distance_square(surface_point, intra_near_atom.xyz);
                        if(dis_square < ((intra_near_atom.radius+param.r_water)*(intra_near_atom.radius+param.r_water))) {
                            buried = true;
                            break;
                        }
                    }
                    if(buried == true) {
                        continue;
                    }
                    effective_count++;
                }
                SumSASAtmp += 4.0 * k_pi * (a.radius+param.r_water) * (a.radius+param.r_water) * static_cast<Float>(effective_count) / static_cast<Float>(param.sasa_points[a.elem_enum].size());
            } else {
                SumSASAtmp += 4.0 * k_pi * (a.radius+param.r_water) * (a.radius+param.r_water);
            }
        }
        return SumSASAtmp;
    }

    const Float cal_delta_rec_sasa(const Parameter& param, const Receptor& rec, const Compound& cpd, const std::vector<std::vector<Size_Type>>& intra_near_indices_list, const std::vector<std::vector<Size_Type>>& inter_near_indices_list) {
        Float SumSASAtmp = 0.0;
        for(const auto& a : rec.atoms) {
            const auto& inter_near_indices = inter_near_indices_list[a.index];
            const auto& intra_near_indices = intra_near_indices_list[a.index];
            if(inter_near_indices.size() == 0) {
                continue;
            }
            Size_Type effective_count = 0;
            for(const auto& point : param.sasa_points[a.elem_enum]) {
                const Vec3d surface_point = point + a.xyz;
                bool buried = false;
                for(Size_Type l=0; l != intra_near_indices.size(); l++) {
                    const Atom& intra_near_atom = rec.atoms[intra_near_indices[l]];
                    const Float dis_square = vec3d_distance_square(surface_point, intra_near_atom.xyz);
                    if(dis_square < ((intra_near_atom.radius+param.r_water)*(intra_near_atom.radius+param.r_water))) {
                        buried = true;
                        break;
                    }
                }
                if(buried == true) {
                    continue;
                }
                for(Size_Type l=0; l != inter_near_indices.size(); l++) {
                    const Atom& inter_near_atom = cpd.atoms[inter_near_indices[l]];
                    const Float dis_square = vec3d_distance_square(surface_point, inter_near_atom.xyz);
                    if(dis_square < ((inter_near_atom.radius+param.r_water)*(inter_near_atom.radius+param.r_water))) {
                        effective_count++;
                        break;
                    }
                }
            }
            SumSASAtmp += 4.0 * k_pi * (a.radius+param.r_water) * (a.radius+param.r_water) * static_cast<Float>(effective_count) / static_cast<Float>(param.sasa_points[a.elem_enum].size());
        }
        return -SumSASAtmp;
    }

    template<typename T1, typename T2>
    const Float cal_mol_sasa(const Parameter& param, const T1& target_mol, const T2& near_mol, const std::vector<std::vector<Size_Type>>& intra_near_indices_list, const std::vector<std::vector<Size_Type>>& inter_near_indices_list) {
        Float SumSASAtmp = 0.0;
        for(const auto& a : target_mol.atoms) {
            const auto& inter_near_indices = inter_near_indices_list[a.index];
            const auto& intra_near_indices = intra_near_indices_list[a.index];
            if(inter_near_indices.size() == 0 and intra_near_indices.size() == 0) {
                SumSASAtmp += 4.0 * k_pi * (a.radius+param.r_water) * (a.radius+param.r_water);
            } else {
                Size_Type effective_count = 0;
                for(const auto& point : param.sasa_points[a.elem_enum]) {
                    const Vec3d surface_point = point + a.xyz;
                    bool intra_buried = false;
                    for(Size_Type l=0; l != intra_near_indices.size(); l++) {
                        const Atom& intra_near_atom = target_mol.atoms[intra_near_indices[l]];
                        const Float dis_square = vec3d_distance_square(surface_point, intra_near_atom.xyz);
                        if(dis_square < ((intra_near_atom.radius+param.r_water)*(intra_near_atom.radius+param.r_water))) {
                            intra_buried = true;
                            break;
                        }
                    }
                    if(intra_buried == true) {
                        continue;
                    }
                    bool inter_buried = false;
                    for(Size_Type l=0; l != inter_near_indices.size(); l++) {
                        const Atom& inter_near_atom = near_mol.atoms[inter_near_indices[l]];
                        const Float dis_square = vec3d_distance_square(surface_point, inter_near_atom.xyz);
                        if(dis_square < ((inter_near_atom.radius+param.r_water)*(inter_near_atom.radius+param.r_water))) {
                            inter_buried = true;
                            break;
                        }
                    }
                    if(inter_buried) {
                        continue;
                    }
                    effective_count++;
                }
                SumSASAtmp += 4.0 * k_pi * (a.radius+param.r_water) * (a.radius+param.r_water) * static_cast<Float>(effective_count) / static_cast<Float>(param.sasa_points[a.elem_enum].size());
            }
        }
        return SumSASAtmp;
    }
}
