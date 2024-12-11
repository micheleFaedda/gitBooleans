/*****************************************************************************************
 *              MIT License                                                              *
 *                                                                                       *
 * Copyright (c) 2022 G. Cherchi, F. Pellacini, M. Attene and M. Livesu                  *
 *                                                                                       *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this  *
 * software and associated documentation files (the "Software"), to deal in the Software *
 * without restriction, including without limitation the rights to use, copy, modify,    *
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to    *
 * permit persons to whom the Software is furnished to do so, subject to the following   *
 * conditions:                                                                           *
 *                                                                                       *
 * The above copyright notice and this permission notice shall be included in all copies *
 * or substantial portions of the Software.                                              *
 *                                                                                       *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,   *
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A         *
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT    *
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION     *
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE        *
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                *
 *                                                                                       *
 * Authors:                                                                              *
 *      Gianmarco Cherchi (g.cherchi@unica.it)                                           *
 *      https://www.gianmarcocherchi.com                                                 *
 *                                                                                       *
 *      Fabio Pellacini (fabio.pellacini@uniroma1.it)                                    *
 *      https://pellacini.di.uniroma1.it                                                 *
 *                                                                                       *
 *      Marco Attene (marco.attene@ge.imati.cnr.it)                                      *
 *      https://www.cnr.it/en/people/marco.attene/                                       *
 *                                                                                       *
 *      Marco Livesu (marco.livesu@ge.imati.cnr.it)                                      *
 *      http://pers.ge.imati.cnr.it/livesu/                                              *
 *                                                                                       *
 * ***************************************************************************************/
#include "booleans.h"
#include "debug.h"
#include "io_functions.h"
#include <tbb/tbb.h>
#include <cinolib/rationals.h>
#include <intersect_point_rationals.h>
#include <ranges> // For std::ranges::contains
inline void customBooleanPipeline(std::vector<genericPoint*>& arr_verts, std::vector<uint>& arr_in_tris,
                                  std::vector<uint>& arr_out_tris, std::vector<std::bitset<NBIT>>& arr_in_labels,
                                  std::vector<DuplTriInfo>& dupl_triangles, Labels& labels,
                                  std::vector<phmap::flat_hash_set<uint>>& patches, cinolib::Octree& octree,
                                  const BoolOp &op, std::vector<double> &bool_coords, std::vector<uint> &bool_tris,
                                  std::vector< std::bitset<NBIT>> &bool_labels)
{
    FastTrimesh tm(arr_verts, arr_out_tris, false);

    computeAllPatches(tm, labels, patches, false);

    // the informations about duplicated triangles (removed in arrangements) are restored in the original structures
    addDuplicateTrisInfoInStructures(dupl_triangles, arr_in_tris, arr_in_labels, octree);

    // parse patches with octree and rays
    cinolib::vec3d max_coords(octree.root->bbox.max.x() +0.5, octree.root->bbox.max.y() +0.5, octree.root->bbox.max.z() +0.5);
    computeInsideOut(tm, patches, octree, arr_verts, arr_in_tris, arr_in_labels, max_coords, labels);

    // boolean operations
    uint num_tris_in_final_solution;


    if(op == INTERSECTION)
        num_tris_in_final_solution = boolIntersection(tm, labels);
    else if(op == UNION)
        num_tris_in_final_solution = boolUnion(tm, labels);
    else if(op == SUBTRACTION)
        num_tris_in_final_solution = boolSubtraction(tm, labels);
    else if(op == XOR)
        num_tris_in_final_solution = boolXOR(tm, labels);
    else
    {
        std::cerr << "boolean operation not implemented yet" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    computeFinalExplicitResult(tm, labels, num_tris_in_final_solution, bool_coords, bool_tris, bool_labels, true);
/*

    uint n_patch = 1;
    std::vector<double> patch_coords;
    patch_coords.resize(patches.at(n_patch).size() * 3 * 3);
    std::vector<uint> patch_tris;
    patch_tris.resize(patches.at(n_patch).size() * 3);

    std::cout<<"Patch size: "<<patches.at(n_patch).size()<<std::endl;
    auto p = patches.at(n_patch);

    for(auto t : p){
        std::cout<<"Triangle: "<<t<<std::endl;
        for(int j = 0; j < 3; ++j){
            uint v_id = tm.tri(t)[j];
            double x, y, z;
            tm.vert(v_id)->getApproxXYZCoordinates(x, y, z);
            patch_coords.push_back(x);
            patch_coords.push_back(y);
            patch_coords.push_back(z);
            patch_tris.push_back(v_id);
        }
    }


    cinolib::write_OBJ("patch1.obj", patch_coords, patch_tris, {});*/
}

extern int arr_time;
extern int bool_time;
extern std::vector<std::string> files;

inline void booleanPipeline(const std::vector<double> &in_coords, const std::vector<uint> &in_tris,
                            const std::vector<uint> &in_labels, const BoolOp &op, std::vector<double> &bool_coords,
                            std::vector<uint> &bool_tris, std::vector< std::bitset<NBIT> > &bool_labels)
{
    initFPU();

    point_arena arena;
    std::vector<genericPoint*> arr_verts; // <- it contains the original expl verts + the new_impl verts
    std::vector<uint> arr_in_tris, arr_out_tris;
    std::vector<std::bitset<NBIT>> arr_in_labels;
    std::vector<DuplTriInfo> dupl_triangles;
    Labels labels;
    std::vector<phmap::flat_hash_set<uint>> patches;
    cinolib::Octree octree; // built with arr_in_tris and arr_in_labels

    customArrangementPipeline(in_coords, in_tris, in_labels, arr_in_tris, arr_in_labels, arena, arr_verts,
                              arr_out_tris, labels, octree, dupl_triangles);

    customBooleanPipeline(arr_verts, arr_in_tris, arr_out_tris, arr_in_labels, dupl_triangles, labels,
                          patches, octree, op, bool_coords, bool_tris, bool_labels);
}



//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

/* a custom arrangement pipeline in witch we can expose the octree used to find the starting intersection list */
inline void customArrangementPipeline(const std::vector<double> &in_coords, const std::vector<uint> &in_tris, const std::vector<uint> &in_labels,
                                      std::vector<uint> &arr_in_tris, std::vector< std::bitset<NBIT>> &arr_in_labels,
                                      point_arena& arena, std::vector<genericPoint *> &vertices, std::vector<uint> &arr_out_tris, Labels &labels,
                                      cinolib::Octree &octree, std::vector<DuplTriInfo> &dupl_triangles)
{
    arr_in_labels.resize(in_labels.size());
    std::bitset<NBIT> mask;
    for(uint i = 0; i < in_labels.size(); i++)
    {
        arr_in_labels[i][in_labels[i]] = true;
        mask[in_labels[i]] = true;
    }

    labels.num = mask.count();
    initFPU();
    double multiplier = computeMultiplier(in_coords);

    mergeDuplicatedVertices(in_coords, in_tris, arena, vertices, arr_in_tris, false);

    customRemoveDegenerateAndDuplicatedTriangles(vertices, arr_in_tris, arr_in_labels, dupl_triangles, false);

    TriangleSoup ts(arena, vertices, arr_in_tris, arr_in_labels, multiplier, false);

    AuxiliaryStructure g;

    customDetectIntersections(ts, g.intersectionList(), octree);

    g.initFromTriangleSoup(ts);

    classifyIntersections(ts, arena, g);

    triangulation(ts, arena, g, arr_out_tris, labels.surface);

    //print labels.surface to check if it is correct
   /*for (int i = 0; i < labels.surface.size(); i++){
         std::cout <<"Triangle "<< i << " labels.surface: " <<labels.surface[i] << std::endl;
    }*/


    ts.appendJollyPoints();

    labels.inside.resize(arr_out_tris.size() / 3);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void customRemoveDegenerateAndDuplicatedTriangles(const std::vector<genericPoint *> &verts, std::vector<uint> &tris,
                                                  std::vector<std::bitset<NBIT>> &labels, std::vector<DuplTriInfo> &dupl_triangles,
                                                  bool parallel)
{
    if(parallel)
    {
        using vec3i = std::array<uint, 3>;
        uint num_orig_tris = static_cast<uint>(tris.size() / 3);
        vec3i* data_orig_tris = (vec3i*)tris.data();

        // compute colinear
        auto colinear = vector<bool>(num_orig_tris, false);
        tbb::parallel_for((uint)0, num_orig_tris, [data_orig_tris, &colinear, &verts](uint t_id) {
            auto& t = data_orig_tris[t_id];
            colinear[t_id] = cinolib::points_are_colinear_3d(
                verts[t[0]]->toExplicit3D().ptr(),
                verts[t[1]]->toExplicit3D().ptr(),
                verts[t[2]]->toExplicit3D().ptr());
        });

        // loop as before by use simpler way
        uint t_off = 0, l_off = 0;

        phmap::flat_hash_map < std::array<uint, 3>, std::pair<uint, uint> > tris_map; // tri_vertices -> <l_off, t_off>
        tris_map.reserve(num_orig_tris);

        for(uint t_id = 0; t_id < num_orig_tris; t_id++)
        {
            uint v0_id = tris[(3 * t_id)];
            uint v1_id = tris[(3 * t_id) +1];
            uint v2_id = tris[(3 * t_id) +2];
            std::bitset<NBIT> l = labels[t_id];

            if(!colinear[t_id]) // good triangle
            {
                std::array<uint, 3> tri = {v0_id, v1_id, v2_id};
                std::sort(tri.begin(), tri.end());

                auto ins= tris_map.insert({tri, std::make_pair(l_off, t_off)});

                if(ins.second) // first time for tri v0, v1, v2
                {
                    labels[l_off] = l;
                    l_off++;

                    tris[t_off] = v0_id, tris[t_off +1] = v1_id, tris[t_off +2] = v2_id;
                    t_off += 3;
                }
                else // triangle already present
                {
                    uint pos = ins.first->second.first;
                    labels[pos] |= l; // label for duplicates
                }

                if(!ins.second) // triangle already present -> save info about duplicates
                {
                    uint orig_tri_off = ins.first->second.second;

                    uint mesh_l = bitsetToUint(l);
                    assert(mesh_l >= 0);

                    uint curr_tri_verts[] = {v0_id, v1_id, v2_id};
                    uint orig_tri_verts[] ={tris[orig_tri_off], tris[orig_tri_off +1], tris[orig_tri_off +2]};

                    bool w = consistentWinding(curr_tri_verts, orig_tri_verts);

                    dupl_triangles.push_back({orig_tri_off / 3, // original triangle id
                                            static_cast<uint>(mesh_l), // label of the actual triangle
                                            w}); // winding with respect to the triangle stored in mesh (true -> same, false -> opposite)
                }
            }
        }

        tris.resize(t_off);
        labels.resize(l_off);
    } else {
        uint num_orig_tris = static_cast<uint>(tris.size() / 3);
        uint t_off = 0, l_off = 0;

        phmap::flat_hash_map < std::array<uint, 3>, std::pair<uint, uint> > tris_map; // tri_vertices -> <l_off, t_off>
        tris_map.reserve(num_orig_tris);

        for(uint t_id = 0; t_id < num_orig_tris; t_id++)
        {
            uint v0_id = tris[(3 * t_id)];
            uint v1_id = tris[(3 * t_id) +1];
            uint v2_id = tris[(3 * t_id) +2];
            std::bitset<NBIT> l = labels[t_id];

            if(!cinolib::points_are_colinear_3d(verts[v0_id]->toExplicit3D().ptr(),
                                                verts[v1_id]->toExplicit3D().ptr(),
                                                verts[v2_id]->toExplicit3D().ptr())) // good triangle
            {
                std::array<uint, 3> tri = {v0_id, v1_id, v2_id};
                std::sort(tri.begin(), tri.end());

                auto ins= tris_map.insert({tri, std::make_pair(l_off, t_off)});

                if(ins.second) // first time for tri v0, v1, v2
                {
                    labels[l_off] = l;
                    l_off++;

                    tris[t_off] = v0_id, tris[t_off +1] = v1_id, tris[t_off +2] = v2_id;
                    t_off += 3;
                }
                else // triangle already present
                {
                    uint pos = ins.first->second.first;
                    labels[pos] |= l; // label for duplicates
                }

                if(!ins.second) // triangle already present -> save info about duplicates
                {
                    uint orig_tri_off = ins.first->second.second;

                    uint mesh_l = bitsetToUint(l);
                    assert(mesh_l >= 0);

                    uint curr_tri_verts[] = {v0_id, v1_id, v2_id};
                    uint orig_tri_verts[] ={tris[orig_tri_off], tris[orig_tri_off +1], tris[orig_tri_off +2]};

                    bool w = consistentWinding(curr_tri_verts, orig_tri_verts);

                    dupl_triangles.push_back({orig_tri_off / 3, // original triangle id
                                            static_cast<uint>(mesh_l), // label of the actual triangle
                                            w}); // winding with respect to the triangle stored in mesh (true -> same, false -> opposite)
                }
            }
        }

        tris.resize(t_off);
        labels.resize(l_off);
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void customDetectIntersections(const TriangleSoup &ts, std::vector<std::pair<uint, uint> > &intersection_list, cinolib::Octree &o)
{
    std::vector<cinolib::vec3d> verts(ts.numVerts());



    for(uint v_id = 0; v_id < ts.numVerts(); v_id++)
        verts[v_id] = cinolib::vec3d(ts.vertX(v_id), ts.vertY(v_id), ts.vertZ(v_id));


    o.build_from_vectors(verts, ts.trisVector());

    intersection_list.reserve(ts.numTris());

    tbb::spin_mutex mutex;
    tbb::parallel_for((uint)0, (uint)o.leaves.size(), [&](uint i)
    {
        auto & leaf = o.leaves[i];
        if(leaf->item_indices.empty()) return;
        for(uint j=0; j<leaf->item_indices.size()-1; ++j)
            for(uint k=j+1; k<leaf->item_indices.size();   ++k)
            {
                uint tid0 = leaf->item_indices[j];
                uint tid1 = leaf->item_indices[k];
                auto T0 = o.items[tid0];
                auto T1 = o.items[tid1];
                if(T0->aabb.intersects_box(T1->aabb)) // early reject based on AABB intersection
                {
                    const cinolib::Triangle *t0 = reinterpret_cast<cinolib::Triangle*>(T0);
                    const cinolib::Triangle *t1 = reinterpret_cast<cinolib::Triangle*>(T1);
                    if(t0->intersects_triangle(t1->v,true)) // precise check (exact if CINOLIB_USES_EXACT_PREDICATES is defined)
                    {
                        std::lock_guard<tbb::spin_mutex> guard(mutex);
                        intersection_list.push_back(cinolib::unique_pair(tid0,tid1));
                    }
                }
            }
    });
    remove_duplicates(intersection_list);
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void addDuplicateTrisInfoInStructures(const std::vector<DuplTriInfo> &dupl_tris, std::vector<uint> &in_tris,
                                             std::vector<std::bitset<NBIT>> &in_labels, cinolib::Octree &octree)
{
    for(auto &item : dupl_tris)
    {
        uint v0_id = in_tris[3 * item.t_id];
        uint v1_id = in_tris[3 * item.t_id + 1];
        uint v2_id = in_tris[3 * item.t_id + 2];

        std::bitset<NBIT> new_label;
        new_label[item.l_id] = true;

        const cinolib::Triangle *orig_tri = dynamic_cast<const cinolib::Triangle*>(octree.items.at(item.t_id));

        uint new_t_id = in_tris.size() / 3;

        if(item.w)
        {
            in_tris.push_back(v0_id);
            in_tris.push_back(v1_id);
            in_tris.push_back(v2_id);
            octree.items.push_back(new cinolib::Triangle(new_t_id, orig_tri->v[0], orig_tri->v[1], orig_tri->v[2]));
        }
        else
        {
            in_tris.push_back(v0_id);
            in_tris.push_back(v2_id);
            in_tris.push_back(v1_id);
            octree.items.push_back(new cinolib::Triangle(new_t_id, orig_tri->v[0], orig_tri->v[1], orig_tri->v[2]));
        }

        in_labels.push_back(new_label); // we add the new_label to the new_triangle
        in_labels[item.t_id][item.l_id] = false; // we remove the dupl label from the orig triangle


    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void computeAllPatches(FastTrimesh &tm, const Labels &labels, std::vector<phmap::flat_hash_set<uint>> &patches, bool parallel)
{
    if(parallel) {
        tm.resetVerticesInfo();
        auto adjT2E = tm.adjT2EAll(parallel);

        for(uint t_id = 0; t_id < tm.numTris(); t_id++)
        {
            if(tm.triInfo(t_id) != 1)
            {
                patches.emplace_back();
                computeSinglePatch(tm, t_id, labels, patches.back(), adjT2E);
            }
        }
    } else {
        tm.resetVerticesInfo();

        for(uint t_id = 0; t_id < tm.numTris(); t_id++)
        {
            if(tm.triInfo(t_id) != 1)
            {
                patches.emplace_back();
                computeSinglePatch(tm, t_id, labels, patches.back());
            }
        }
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void computeSinglePatch(FastTrimesh &tm, uint seed_t, const Labels &labels, phmap::flat_hash_set<uint> &patch)
{
    std::bitset<NBIT> ref_l = labels.surface[seed_t];

    std::stack<uint> tris_stack;
    tris_stack.push(seed_t);

    while(!tris_stack.empty())
    {
        uint curr_t = tris_stack.top();
        tris_stack.pop();

        tm.setTriInfo(curr_t, 1); // set triangle as visited
        patch.insert(curr_t);

        for(uint e_id : tm.adjT2E(curr_t))
        {
            if(tm.edgeIsManifold(e_id))
            {
                for(uint t_id : tm.adjE2T(e_id))
                {
                    if(t_id != curr_t && tm.triInfo(t_id) != 1)
                    {
                        assert(labels.surface[t_id] == ref_l);
                        tris_stack.push(t_id);
                    }
                }
            }
            else // e_id is not manifold -> stop flooding
            {
                // we set the vertices in the patch border with 1 (useful for ray computation funcion)
                tm.setVertInfo(tm.edgeVertID(e_id, 0), 1);
                tm.setVertInfo(tm.edgeVertID(e_id, 1), 1);
            }
        }
    }
}

inline void computeSinglePatch(FastTrimesh &tm, uint seed_t, const Labels &labels, phmap::flat_hash_set<uint> &patch, const std::vector<std::array<uint, 3>>& adjT2E)
{
    std::bitset<NBIT> ref_l = labels.surface[seed_t];

    std::stack<uint> tris_stack;
    tris_stack.push(seed_t);

    while(!tris_stack.empty())
    {
        uint curr_t = tris_stack.top();
        tris_stack.pop();

        tm.setTriInfo(curr_t, 1); // set triangle as visited
        patch.insert(curr_t);

        for(uint e_id : adjT2E[curr_t])
        {
            if(tm.edgeIsManifold(e_id))
            {
                for(uint t_id : tm.adjE2T(e_id))
                {
                    if(t_id != curr_t && tm.triInfo(t_id) != 1)
                    {
                        assert(labels.surface[t_id] == ref_l);
                        tris_stack.push(t_id);
                    }
                }
            }
            else // e_id is not manifold -> stop flooding
            {
                // we set the vertices in the patch border with 1 (useful for ray computation funcion)
                tm.setVertInfo(tm.edgeVertID(e_id, 0), 1);
                tm.setVertInfo(tm.edgeVertID(e_id, 1), 1);
            }
        }
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void findRayEndpoints(const FastTrimesh &tm, const phmap::flat_hash_set<uint> &patch, const cinolib::vec3d &max_coords, Ray &ray)
{
    // check for an explicit point (all operations with explicits are faster)
    int v_id = -1;
    for(uint t_id : patch)
    {
        const uint tv[3] = {tm.triVertID(t_id, 0), tm.triVertID(t_id, 1), tm.triVertID(t_id, 2)};

        if (tm.vert(tv[0])->isExplicit3D() && tm.vertInfo(tv[0]) == 0)      v_id = static_cast<int>(tv[0]);
        else if (tm.vert(tv[1])->isExplicit3D() && tm.vertInfo(tv[1]) == 0) v_id = static_cast<int>(tv[1]);
        else if (tm.vert(tv[2])->isExplicit3D() && tm.vertInfo(tv[2]) == 0) v_id = static_cast<int>(tv[2]);

        if (v_id != -1)
        {
            const explicitPoint3D &v = tm.vert(v_id)->toExplicit3D();
            ray.v0  = explicitPoint3D(v.X(), v.Y(), v.Z());
            ray.v1 = explicitPoint3D(max_coords.x(), v.Y(), v.Z());
            return;
        }
    }

    int tri_counter = 0;
    // parse triangles with all implicit points
    for(uint t_id : patch)
    {
        double x0, x1, x2, y0, y1, y2, z0, z1, z2;
        tm.triVert(t_id, 0)->getApproxXYZCoordinates(x0, y0, z0);
        tm.triVert(t_id, 1)->getApproxXYZCoordinates(x1, y1, z1);
        tm.triVert(t_id, 2)->getApproxXYZCoordinates(x2, y2, z2);

        explicitPoint3D tv0(x0, y0, z0), tv1(x1, y1, z1), tv2(x2, y2, z2);
        if(!genericPoint::misaligned(tv0, tv1, tv2)) continue;

        int dir = genericPoint::maxComponentInTriangleNormal(x0, y0, z0, x1, y1, z1, x2, y2, z2);
        if(dir == 0) // dir = X
        {
            ray.v0 = explicitPoint3D(((x0 + x1 + x2) / 3.0) - 0.1, (y0 + y1 + y2) / 3.0, (z0 + z1 + z2) / 3.0);
            ray.v1 = explicitPoint3D(max_coords.x(), ray.v0.Y(), ray.v0.Z());
            ray.dir = 'X';
        }
        else if(dir == 1) // dir = Y
        {
            ray.v0 = explicitPoint3D((x0 + x1 + x2) / 3.0, ((y0 + y1 + y2) / 3.0) -0.1, (z0 + z1 + z2) / 3.0);
            ray.v1 = explicitPoint3D(ray.v0.X(), max_coords.y(), ray.v0.Z());
            ray.dir = 'Y';
        }
        else // dir = Z
        {
            ray.v0 = explicitPoint3D((x0 + x1 + x2) / 3.0, (y0 + y1 + y2) / 3.0, ((z0 + z1 + z2) / 3.0)  -0.1);
            ray.v1 = explicitPoint3D(ray.v0.X(), ray.v0.Y(), max_coords.z());
            ray.dir = 'Z';
        }

        int orf = genericPoint::orient3D(*tm.triVert(t_id, 0), *tm.triVert(t_id, 1), *tm.triVert(t_id, 2), ray.v0);
        int ors = genericPoint::orient3D(*tm.triVert(t_id, 0), *tm.triVert(t_id, 1), *tm.triVert(t_id, 2), ray.v1);

        if((orf < 0 && ors > 0) || (orf > 0 && ors < 0)) // the ray passes through the triangle
        {
            if(checkIntersectionInsideTriangle3DImplPoints(ray, tm.triVert(t_id, 0), tm.triVert(t_id, 1), tm.triVert(t_id, 2))) // the ray passes inside the triangle
            {
                ray.tv[0] = static_cast<int>(tm.triVertID(t_id, 0));
                ray.tv[1] = static_cast<int>(tm.triVertID(t_id, 1));
                ray.tv[2] = static_cast<int>(tm.triVertID(t_id, 2));
                return;
            }
        }
    }

   // std::cout << "WARNING: the arrangement contains a fully implicit patch that requires exact rationals for evaluation. This version of the code does not support rationals, therefore the output result may contain open boundaries." << std::endl;

    //std::exit(EXIT_FAILURE);
}


inline void findRayEndpointsCustom(const FastTrimesh &tm, const phmap::flat_hash_set<uint> &patch, const cinolib::vec3d &max_coords, Ray &ray, RationalRay &rational_ray)
{
    // check for an explicit point (all operations with explicits are faster)
    int v_id = -1;
    for(uint t_id : patch)
    {
        const uint tv[3] = {tm.triVertID(t_id, 0), tm.triVertID(t_id, 1), tm.triVertID(t_id, 2)};

        if (tm.vert(tv[0])->isExplicit3D() && tm.vertInfo(tv[0]) == 0)      v_id = static_cast<int>(tv[0]);
        else if (tm.vert(tv[1])->isExplicit3D() && tm.vertInfo(tv[1]) == 0) v_id = static_cast<int>(tv[1]);
        else if (tm.vert(tv[2])->isExplicit3D() && tm.vertInfo(tv[2]) == 0) v_id = static_cast<int>(tv[2]);

        if (v_id != -1)
        {
            const explicitPoint3D &v = tm.vert(v_id)->toExplicit3D();
            ray.v0  = explicitPoint3D(v.X(), v.Y(), v.Z());
            ray.v1 = explicitPoint3D(max_coords.x(), v.Y(), v.Z());
            return;
        }
    }
    //commenta
    /*
    int tri_counter = 0;
    bool flag_impl = false;
    // parse triangles with all implicit points
    uint t_id_aux;

    for(uint t_id : patch)
    {
        double x0, x1, x2, y0, y1, y2, z0, z1, z2;
        tm.triVert(t_id, 0)->getApproxXYZCoordinates(x0, y0, z0);
        tm.triVert(t_id, 1)->getApproxXYZCoordinates(x1, y1, z1);
        tm.triVert(t_id, 2)->getApproxXYZCoordinates(x2, y2, z2);

        explicitPoint3D tv0(x0, y0, z0), tv1(x1, y1, z1), tv2(x2, y2, z2);
        if(!genericPoint::misaligned(tv0, tv1, tv2)) continue;

        int dir = genericPoint::maxComponentInTriangleNormal(x0, y0, z0, x1, y1, z1, x2, y2, z2);
        if(dir == 0) // dir = X
        {
            ray.v0 = explicitPoint3D(((x0 + x1 + x2) / 3.0) - 0.1, (y0 + y1 + y2) / 3.0, (z0 + z1 + z2) / 3.0);
            ray.v1 = explicitPoint3D(max_coords.x(), ray.v0.Y(), ray.v0.Z());
            ray.dir = 'X';
        }
        else if(dir == 1) // dir = Y
        {
            ray.v0 = explicitPoint3D((x0 + x1 + x2) / 3.0, ((y0 + y1 + y2) / 3.0) -0.1, (z0 + z1 + z2) / 3.0);
            ray.v1 = explicitPoint3D(ray.v0.X(), max_coords.y(), ray.v0.Z());
            ray.dir = 'Y';
        }
        else // dir = Z
        {
            ray.v0 = explicitPoint3D((x0 + x1 + x2) / 3.0, (y0 + y1 + y2) / 3.0, ((z0 + z1 + z2) / 3.0)  -0.1);
            ray.v1 = explicitPoint3D(ray.v0.X(), ray.v0.Y(), max_coords.z());
            ray.dir = 'Z';
        }

        int orf = genericPoint::orient3D(*tm.triVert(t_id, 0), *tm.triVert(t_id, 1), *tm.triVert(t_id, 2), ray.v0);
        int ors = genericPoint::orient3D(*tm.triVert(t_id, 0), *tm.triVert(t_id, 1), *tm.triVert(t_id, 2), ray.v1);

        if((orf < 0 && ors > 0) || (orf > 0 && ors < 0)) // the ray passes through the triangle
        {
            if(checkIntersectionInsideTriangle3DImplPoints(ray, tm.triVert(t_id, 0), tm.triVert(t_id, 1), tm.triVert(t_id, 2))) // the ray passes inside the triangle
            {
                ray.tv[0] = static_cast<int>(tm.triVertID(t_id, 0));
                ray.tv[1] = static_cast<int>(tm.triVertID(t_id, 1));
                ray.tv[2] = static_cast<int>(tm.triVertID(t_id, 2));
                return;
            }
        }
    }

    //take the patch with triangle with t_id_aux and iterate through the triangles
        */
    for(uint t_id : patch){
        std::cout << "Triangle: " << t_id << std::endl;
        bigrational x0_rat, x1_rat, x2_rat, y0_rat, y1_rat, y2_rat, z0_rat, z1_rat, z2_rat;
        tm.triVert(t_id, 0)->getExactXYZCoordinates(x0_rat, y0_rat, z0_rat);
        tm.triVert(t_id, 1)->getExactXYZCoordinates(x1_rat, y1_rat, z1_rat);
        tm.triVert(t_id, 2)->getExactXYZCoordinates(x2_rat, y2_rat, z2_rat);

        std::vector<bigrational> tv0_rat = {x0_rat, y0_rat, z0_rat};
        std::vector<bigrational> tv1_rat = {x1_rat, y1_rat, z1_rat};
        std::vector<bigrational> tv2_rat = {x2_rat, y2_rat, z2_rat};

        //if(!genericPoint::misaligned(tv0_rat, tv1_rat, tv2_rat)) continue;
        int dir = maxComponentInTriangleNormalRationals(x0_rat, y0_rat, z0_rat, x1_rat, y1_rat, z1_rat, x2_rat, y2_rat, z2_rat);

        if(dir == 0) // dir = X
        {
            rational_ray.v0 = {(x0_rat + x1_rat + x2_rat) / bigrational(3.0),
                               (y0_rat + y1_rat + y2_rat) / bigrational(3.0),
                               (z0_rat + z1_rat + z2_rat) / bigrational(3.0)};


            rational_ray.v1 = {bigrational(max_coords.x()), rational_ray.v0[1], rational_ray.v0[2]};
            rational_ray.dir = 'X';
        }
        else if(dir == 1) // dir = Y
        {

            rational_ray.v0 = {(x0_rat + x1_rat + x2_rat) / bigrational(3.0),
                               (y0_rat + y1_rat + y2_rat) / bigrational(3.0),
                               (z0_rat + z1_rat + z2_rat) / bigrational(3.0)};

            rational_ray.v1  = {rational_ray.v0[0], bigrational(max_coords.y()), rational_ray.v0[2]};
            rational_ray.dir = 'Y';
        }
        else if(dir == 2) // dir = Z
        {
            rational_ray.v0  = {(x0_rat + x1_rat + x2_rat) / bigrational(3.0),
                               (y0_rat + y1_rat + y2_rat) / bigrational(3.0),
                               (z0_rat + z1_rat + z2_rat) / bigrational(3.0)};

            rational_ray.v1 = {rational_ray.v0[0], rational_ray.v0[1], bigrational(max_coords.z())};
            rational_ray.dir = 'Z';
        }else{
            std::cout<<"Error in direction"<<std::endl;
            std::exit(EXIT_FAILURE);
        }

        //TODO: check se centroide sta sullo stesso piano dei vertici del triangolo
        if(cinolib::orient3d(&tv0_rat[0], &tv1_rat[0], &tv2_rat[0], &rational_ray.v0[0]).sgn() != 0){
            std::cout << "The centroid is not on the triangle or the orient3d doesn't work properly" << std::endl;
            std::cout << "The value of orient3d is: " << cinolib::orient3d(&tv0_rat[0], &tv1_rat[0], &tv2_rat[0], &rational_ray.v0[0]) << std::endl;
            std::exit(EXIT_FAILURE);
        } //the centroid is exaclty on the plane of the triangle


        //TODO: controlla che le tre orient siano di segno concorde per controllare che il centroide sta all'interno
        bigrational e0_rat = cinolib::orient3d(&tv0_rat[0], &tv1_rat[0], &rational_ray.v1[0], &rational_ray.v0[0]);
        bigrational e1_rat = cinolib::orient3d(&tv1_rat[0], &tv2_rat[0], &rational_ray.v1[0], &rational_ray.v0[0]);
        bigrational e2_rat = cinolib::orient3d(&tv2_rat[0], &tv0_rat[0], &rational_ray.v1[0], &rational_ray.v0[0]);

        //print the coords of the triangle
        //std::cout << "Triangle coords: " << std::endl;
        //std::cout << "v0: " << tv0_rat[0] << " " << tv0_rat[1] << " " << tv0_rat[2] << "In double: " << tv0_rat[0].get_d() << " " << tv0_rat[1].get_d() << " " << tv0_rat[2].get_d() << std::endl;
        //std::cout << "v1: " << tv1_rat[0] << " " << tv1_rat[1] << " " << tv1_rat[2] << "In double: " << tv1_rat[0].get_d() << " " << tv1_rat[1].get_d() << " " << tv1_rat[2].get_d() << std::endl;
        //std::cout << "v2: " << tv2_rat[0] << " " << tv2_rat[1] << " " << tv2_rat[2] << "In double: " << tv2_rat[0].get_d() << " " << tv2_rat[1].get_d() << " " << tv2_rat[2].get_d() << std::endl;

        //print the coords of the ray
        //std::cout << "Ray coords: " << std::endl;
        //std::cout << "v0: " << rational_ray.v0[0] << " " << rational_ray.v0[1] << " " << rational_ray.v0[2] << "In double: "<< rational_ray.v0[0].get_d() <<" " << rational_ray.v0[1].get_d() << " " << rational_ray.v0[2].get_d() << std::endl;
        //std::cout << "v1: " << rational_ray.v1[0] << " " << rational_ray.v1[1] << " " << rational_ray.v1[2] << "In double: "<< rational_ray.v1[0].get_d() <<" " << rational_ray.v1[1].get_d() << " " << rational_ray.v1[2].get_d() << std::endl;

        //print the dir
        //std::cout << "Direction: " << rational_ray.dir << std::endl;


        if((e0_rat > bigrational(0,0,0) && e1_rat > bigrational(0,0,0) && e2_rat > bigrational(0,0,0)) ||
           (e0_rat < bigrational(0,0,0) && e1_rat < bigrational(0,0,0) && e2_rat < bigrational(0,0,0))){
            rational_ray.tv[0] = static_cast<int>(tm.triVertID(t_id, 0));
            rational_ray.tv[1] = static_cast<int>(tm.triVertID(t_id, 1));
            rational_ray.tv[2] = static_cast<int>(tm.triVertID(t_id, 2));
            return;
        } else{

            std::cout << "The centroid is not inside the triangle or the orient3D doesn't work properly" << std::endl;
            std::cout << "The value of orient3d for the first edge is: " << e0_rat << std::endl;
            std::cout << "The value of orient3d for the second edge is: " << e1_rat << std::endl;
            std::cout << "The value of orient3d for the third edge is: " << e2_rat << std::endl;
            std::cout << "Direction : " << rational_ray.dir << std::endl;

            //print the coords of the triangle
            std::cout << "Triangle coords: " << std::endl;
            std::cout << "v0: " << tv0_rat[0] << " " << tv0_rat[1] << " " << tv0_rat[2] << " In double: " << tv0_rat[0].get_d() << " " << tv0_rat[1].get_d() << " " << tv0_rat[2].get_d() << std::endl;
            std::cout << "v1: " << tv1_rat[0] << " " << tv1_rat[1] << " " << tv1_rat[2] << " In double: " << tv1_rat[0].get_d() << " " << tv1_rat[1].get_d() << " " << tv1_rat[2].get_d() << std::endl;
            std::cout << "v2: " << tv2_rat[0] << " " << tv2_rat[1] << " " << tv2_rat[2] << " In double: " << tv2_rat[0].get_d() << " " << tv2_rat[1].get_d() << " " << tv2_rat[2].get_d() << std::endl;

            //print the coords of the ray
            std::cout << "Ray coords: " << std::endl;
            std::cout << "v0: " << rational_ray.v0[0] << " " << rational_ray.v0[1] << " " << rational_ray.v0[2] << " In double: "<< rational_ray.v0[0].get_d() <<" " << rational_ray.v0[1].get_d() << " " << rational_ray.v0[2].get_d() << std::endl;
            std::cout << "v1: " << rational_ray.v1[0] << " " << rational_ray.v1[1] << " " << rational_ray.v1[2] << " In double: "<< rational_ray.v1[0].get_d() <<" " << rational_ray.v1[1].get_d() << " " << rational_ray.v1[2].get_d() << std::endl;


            std::exit(EXIT_FAILURE);
        }

        //if((orf_rat < 0 && ors_rat > 0) || (orf_rat > 0 && ors_rat < 0)) // the ray passes through the triangle
        //{
        //if(checkIntersectionInsideTriangle3DImplPoints(ray, tv0_rat, tv1_rat, tv2_rat)) // the ray passes inside the triangle
          /*  {
                ray.tv[0] = static_cast<int>(tm.triVertID(t_id, 0));
                ray.tv[1] = static_cast<int>(tm.triVertID(t_id, 1));
                ray.tv[2] = static_cast<int>(tm.triVertID(t_id, 2));
                return;
            }*/
        //}

    }
    /*
    std::cout <<"t_id :" << t_id_aux << std::endl;
    std::cout << "WARNING: the arrangement contains a fully implicit patch that requires exact rationals for evaluation. This version of the code does not support rationals, therefore the output result may contain open boundaries." << std::endl;

    std::exit(EXIT_FAILURE);*/
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline bool intersects_box(const cinolib::Octree& tree, const cinolib::AABB & b, phmap::flat_hash_set<uint> & ids)
{
    auto root = tree.root;
    auto& items = tree.items;

    std::stack<cinolib::OctreeNode*> lifo;
    if(root && root->bbox.intersects_box(b))
    {
        lifo.push(root);
    }

    while(!lifo.empty())
    {
        cinolib::OctreeNode *node = lifo.top();
        lifo.pop();
        assert(node->bbox.intersects_box(b));

        if(node->is_inner())
        {            
            for(int i = 0; i < 8; ++i)
            {
                if(node->children[i]->bbox.intersects_box(b))
                {
                    lifo.push(node->children[i]);
                }
            }
        }
        else
        {
            for(uint i : node->item_indices)
            {
                if(items.at(i)->aabb.intersects_box(b))
                {
                    ids.insert(items.at(i)->id);
                }
            }
        }
    }

    return !ids.empty();
}


inline void findIntersectionsAlongRayRationals(const FastTrimesh &tm, const std::vector<phmap::flat_hash_set<uint>> &patches, const cinolib::Octree& tree ,const std::vector<genericPoint *> &in_verts,
                                               const std::vector<std::bitset<NBIT>> &in_labels, Labels &labels, const RationalRay &rational_ray, uint curr_p_id,
                                               phmap::flat_hash_set<uint> &tmp_inters,
                                               std::vector<IntersectionPointRationals> &inter_rat)
{
    tbb::spin_mutex mutex;
    //tbb::parallel_for((uint)0, (uint)patches.size(), [&](uint p_id)
    std::vector<uint> inters_tris_rat;

    //for (uint p_id_other = 0; p_id_other < patches.size(); p_id_other++) {
       //  if (p_id_other == curr_p_id) continue;

    //I take the coordinates of the end points of the ray
    std::vector<bigrational> ray_v0 = {rational_ray.v0[0], rational_ray.v0[1], rational_ray.v0[2]};
    std::vector<bigrational> ray_v1 = {rational_ray.v1[0], rational_ray.v1[1], rational_ray.v1[2]};

    //for (uint t_id: patches.at(p_id_other)){
    for(uint t_id = 0; t_id < tm.numTris(); t_id++){

            const uint tv[3] = {tm.triVertID(t_id, 0), tm.triVertID(t_id, 1), tm.triVertID(t_id, 2)};

            //I take the coordinates of the vertices of the triangle
            bigrational x0, x1, x2, y0, y1, y2, z0, z1, z2;
            in_verts[tv[0]]->getExactXYZCoordinates(x0, y0, z0);
            in_verts[tv[1]]->getExactXYZCoordinates(x1, y1, z1);
            in_verts[tv[2]]->getExactXYZCoordinates(x2, y2, z2);

            //I put them into a vector that represents the vertices of the triangle
            std::vector<bigrational> tv0 = {x0, y0, z0};
            std::vector<bigrational> tv1 = {x1, y1, z1};
            std::vector<bigrational> tv2 = {x2, y2, z2};

            bool intersection = segment_triangle_intersect_3d(&ray_v0[0], &ray_v1[0], &tv0[0], &tv1[0], &tv2[0]);

            //I check if the ray intersects the triangle
            if (segment_triangle_intersect_3d(&ray_v0[0], &ray_v1[0], &tv0[0], &tv1[0], &tv2[0])) {
                tmp_inters.insert(t_id);

                //fai un vettore di pair di <coordinate, id triangolo>
                std::vector<bigrational> p_int(3);

                plane_line_intersection(&tv0[0], &tv1[0], &tv2[0], &ray_v0[0], &ray_v1[0], &p_int[0]);

                //search the patch id
                uint p_id_other = 0;
                for(uint p_id = 0; p_id < patches.size(); p_id++){
                    if (patches.at(p_id).contains(t_id)){
                        p_id_other = p_id;
                        break;
                    }
                }
                if(t_id ==11){
                    std::cout << "Triangle 11: " << std::endl;
                    std::cout << "Intersection point coordinates: " << p_int[0].get_d() << " " << p_int[1].get_d() << " " << p_int[2].get_d() << std::endl;
                    std::cout << "Intersection point coordinates: " << p_int[0] << " " << p_int[1] << " " << p_int[2] << std::endl;
                }
                IntersectionPointRationals p_int_rat(p_int[0], p_int[1], p_int[2], t_id, p_id_other);
                inter_rat.push_back(p_int_rat);
            }
        }

        //print tmp_inters
        std::cout << "tmp_inters: " << std::endl;
        for (uint t_id: tmp_inters){
            std::cout << t_id << std::endl;
        }

  }
//}

inline void computeInsideOut(const FastTrimesh &tm, const std::vector<phmap::flat_hash_set<uint>> &patches, const cinolib::Octree &octree,
                             const std::vector<genericPoint *> &in_verts, const std::vector<uint> &in_tris,
                             const std::vector<std::bitset<NBIT>> &in_labels, const cinolib::vec3d &max_coords, Labels &labels)
{

    // print the vertex 15 of the mesh
    bigrational x, y, z;
    bigrational x1, y1, z1;
    in_verts[15]->getExactXYZCoordinates(x, y, z);
    in_verts[14]->getExactXYZCoordinates(x1, y1, z1);
    std::cout << "Vertex 15 of the mesh: " << tm.vert(15)->toExplicit3D().X() << " " << tm.vert(15)->toExplicit3D().Y() << " " << tm.vert(15)->toExplicit3D().Z() << std::endl;
    std::cout << "Vertex 15 of the mesh: " << x<< " " << y << " " << z << std::endl;
    std::cout << "Vertex 14 of the mesh: " << tm.vert(14)->toExplicit3D().X() << " " << tm.vert(14)->toExplicit3D().Y() << " " << tm.vert(14)->toExplicit3D().Z() << std::endl;
    std::cout << "Vertex 14 of the mesh: " << x1<< " " << y1 << " " << z1 << std::endl;

    tbb::spin_mutex mutex;
    //tbb::parallel_for((uint)0, (uint)patches.size(), [&](uint p_id)

    std::bitset <NBIT> labels_patch [patches.size()];

   /* for (uint i = 0; i < patches.size(); i++) {
        for(uint t_id : patches[i]) {
                labels_patch[i] = labels.surface.at(t_id);
                break;
        }
    }

    for (uint pp_id = 0; pp_id < patches.size(); pp_id++) {
        std::cout << "patch n° " << pp_id << " of mesh  " << labels_patch[pp_id] << " contains triangles: ";
        for (uint t_id : patches[pp_id])
            std::cout << t_id << " ";
        std::cout << std::endl;
    }*/

    for(uint p_id = 0; p_id < patches.size(); p_id++) //For each patch
    {
        const phmap::flat_hash_set<uint> &patch_tris = patches[p_id];
        const std::bitset<NBIT> &patch_surface_label = labels.surface[*patch_tris.begin()]; // label of the first triangle of the patch

        Ray ray;
        RationalRay rational_ray;
        //findRayEndpoints(tm, patch_tris, max_coords, ray);
        findRayEndpointsCustom(tm, patch_tris, max_coords, ray, rational_ray);


        phmap::flat_hash_set<uint> tmp_inters;

        //std::cout << "rational_ray.tv[0]: " << rational_ray.tv[0]  <<  std::endl;
        if(rational_ray.tv[0] != -1) {//is defined

            std::cout << "Direction:  " << rational_ray.dir <<  std::endl;
            std::vector<IntersectionPointRationals> inter_rat;

            findIntersectionsAlongRayRationals(tm, patches, octree, in_verts, in_labels, labels, rational_ray, p_id, tmp_inters, inter_rat);
            std::vector<uint> inters_tris_rat;

            //pruneIntersectionPart
            /************PRUNE INTERSECTIONS ***********/

            pruneIntersectionsAndSortAlongRayRationals(rational_ray, tm, in_verts, in_tris, in_labels, tmp_inters, patch_surface_label, inter_rat, inters_tris_rat,labels);
           /* phmap::flat_hash_set<uint> visited_tri;
            visited_tri.reserve(inter_rat.size() / 6);
            std::pair<phmap::flat_hash_set<uint>::iterator, bool> ins;


            for (uint t_id_int: tmp_inters) {
                //uint t_id_app = t_id_int;
                ins = visited_tri.insert(t_id_int);

                if (!ins.second) continue; // triangle already analyzed or in the one ring of a vert or in the adj of an edge
*/
                /*int t_id_father = -1;


                if (t_id_int >= in_labels.size()){
                   uint patch_id= IntersectionPointRationals::findPatchIdByTriId(inter_rat, t_id_int);
                   phmap::flat_hash_set<uint> patch_app = patches.at(patch_id);
                    for(uint t_id : patch_app){
                        if (t_id < in_labels.size()){
                            t_id_father = t_id;
                            break;
                        }
                    }
                    t_id_int = t_id_father;
                }*/


               /* const std::bitset<NBIT> tested_tri_label = labels.surface.at(t_id_int);
                uint uint_tri_label = bitsetToUint(tested_tri_label);
                //t_id_int = t_id_app;

                if (patch_surface_label[uint_tri_label]) continue; // <-- triangle of the same label of the tested patch

                bigrational tv0_x, tv0_y, tv0_z,
                            tv1_x, tv1_y, tv1_z,
                            tv2_x, tv2_y, tv2_z;

                const uint tv[3] = {tm.triVertID(t_id_int, 0), tm.triVertID(t_id_int, 1), tm.triVertID(t_id_int, 2)};

                in_verts[tv[0]]->getExactXYZCoordinates(tv0_x, tv0_y, tv0_z);
                in_verts[tv[1]]->getExactXYZCoordinates(tv1_x, tv1_y, tv1_z);
                in_verts[tv[2]]->getExactXYZCoordinates(tv2_x, tv2_y, tv2_z);


                //print size in tris
                //in_verts[in_tris[3 * t_id_int]]->getExactXYZCoordinates(tv0_x, tv0_y, tv0_z);
                //in_verts[in_tris[3 * t_id_int + 1]]->getExactXYZCoordinates(tv1_x, tv1_y, tv1_z);
                //in_verts[in_tris[3 * t_id_int + 2]]->getExactXYZCoordinates(tv2_x, tv2_y, tv2_z);

                const std::vector<bigrational> tv0_exact = {tv0_x, tv0_y, tv0_z};
                const std::vector<bigrational> tv1_exact = {tv1_x, tv1_y, tv1_z};
                const std::vector<bigrational> tv2_exact = {tv2_x, tv2_y, tv2_z};

                //std::cout <<"t_id_int before fast2DIntersectionOnRayRationals: " << t_id_int << std::endl;

                IntersInfo ii = fast2DCheckIntersectionOnRayRationals(rational_ray, tv0_exact, tv1_exact, tv2_exact);

                std::cout<< "ii: " << ii << std::endl;

                if (ii == DISCARD){
                    std::cout<<"DISCARD" << std::endl;
                     continue;}
                if (ii == NO_INT){
                    std::cout<<"NO_INT" << std::endl;
                    continue;}

                if (ii == INT_IN_TRI) {
                    std::cout << "INT_IN_TRI" << " t_id-> "<<t_id_int<<std::endl;

                    inters_tris_rat.push_back(t_id_int);

                } else if (ii == INT_IN_V0 || ii == INT_IN_V1 || ii == INT_IN_V2) {
                    std::cout << "INT_IN_V" << std::endl;
                    uint v_id;
                    if (ii == INT_IN_V0) v_id = in_tris[3 * t_id_int];
                    else if (ii == INT_IN_V1) v_id = in_tris[3 * t_id_int + 1];
                    else v_id = in_tris[3 * t_id_int + 2];

                    std::vector<uint> vert_one_ring;
                    findVertRingTris(v_id, tested_tri_label, tmp_inters, in_tris, in_labels, vert_one_ring);

                    for (uint t: vert_one_ring)
                        visited_tri.insert(t); // mark all the one ring as visited

                    int winner_tri = -1;
                    winner_tri = perturbRayAndFindIntersTri(rational_ray, in_verts, in_tris,
                                                            vert_one_ring); // the first inters triangle after ray perturbation

                    if (winner_tri != -1)
                        inters_tris_rat.push_back(winner_tri);

                } else if (ii == INT_IN_EDGE01 || ii == INT_IN_EDGE12 || ii == INT_IN_EDGE20) {
                    std::cout << "INT_IN_EDGE" << std::endl;

                   int ev0_id, ev1_id;
                   if (ii == INT_IN_EDGE01) {
                       ev0_id = in_tris[3 * t_id_int];
                       ev1_id = in_tris[3 * t_id_int + 1];
                   } else if (ii == INT_IN_EDGE12) {
                       ev0_id = in_tris[3 * t_id_int + 1];
                       ev1_id = in_tris[3 * t_id_int + 2];
                   } else {
                       ev0_id = in_tris[3 * t_id_int + 2];
                       ev1_id = in_tris[3 * t_id_int];
                   }

                   std::cout << "t_id_int: " << t_id_int<<std::endl;
                   std::cout << "ev0_id: " << ev0_id << std::endl;
                   std::cout << "ev1_id: " << ev1_id << std::endl;

                   std::vector<uint> edge_tris;

                    findEdgeTrisCustom(ev0_id, ev1_id, tested_tri_label, tmp_inters, in_tris, in_labels, edge_tris);

                    for (uint t: edge_tris)
                        visited_tri.insert(t); // mark all the one ring as visited

                    int winner_tri = -1;
                    winner_tri = perturbRayAndFindIntersTri(rational_ray, in_verts, in_tris, edge_tris);

                    if (winner_tri != -1)
                        inters_tris_rat.push_back(winner_tri);

                }

                for (auto &t_id_tmp : inters_tris_rat) {
                    bool found = false;
                    for (auto &p_int: inter_rat) {
                        if (t_id_tmp == p_int.getTriId()) {
                            found = true;
                            break;
                        }
                    }
                    assert(found && "The triangle is not in the intersection points");
                }

                std::cout << "inter_rat size: " << inter_rat.size() << std::endl;
                if (rational_ray.dir == 'X')
                    //sort the inter_rat along the x coordinate
                    sort(inter_rat.begin(), inter_rat.end(),
                         [](const IntersectionPointRationals &a, const IntersectionPointRationals &b) {
                             return a.lessThanX(b);
                         });
                else if (rational_ray.dir == 'Y')
                    sort(inter_rat.begin(), inter_rat.end(),
                         [](const IntersectionPointRationals &a, const IntersectionPointRationals &b) {
                             return a.lessThanY(b);
                         });
                else
                    sort(inter_rat.begin(), inter_rat.end(),
                         [](const IntersectionPointRationals &a, const IntersectionPointRationals &b) {
                             return a.lessThanZ(b);
                         });
*/
                /*********ANALYZE INTERSECITONS **/
                //analyse the sorted intersection points
                std::bitset<NBIT> patch_inner_label;
                analyzeSortedIntersectionsRationals(rational_ray, tm, in_verts, inter_rat, patch_inner_label, labels);
               /* std::bitset<NBIT> visited_labels;


                //copy the tri_id into inter_rat_tris_ordered from inter_rat
                std::vector<uint> inter_rat_tris_ordered;

                for (auto &p: inter_rat) {
                    inter_rat_tris_ordered.push_back(p.getTriId());
                }
                bool print_ray = true;
                std::vector<uint> list_t_id;

                for (int i = 0; i < inter_rat_tris_ordered.size(); i++) {

                    uint t_id = inter_rat_tris_ordered[i];
                    uint t_label = bitsetToUint(labels.surface.at(t_id));

                    //if (visited_labels[t_label]) continue;

                    for (auto &p: inter_rat) {
                        if (p.getTriId() == t_id) {
                            //get the three points of the triangle
                            bigrational x0_aux, x1_aux, x2_aux, y0_aux, y1_aux, y2_aux, z0_aux, z1_aux, z2_aux;

                            uint tv_aux[3] = {tm.triVertID(t_id, 0), tm.triVertID(t_id, 1), tm.triVertID(t_id, 2)};

                            in_verts[tv_aux[0]]->getExactXYZCoordinates(x0_aux, y0_aux, z0_aux);
                            in_verts[tv_aux[1]]->getExactXYZCoordinates(x1_aux, y1_aux, z1_aux);
                            in_verts[tv_aux[2]]->getExactXYZCoordinates(x2_aux, y2_aux, z2_aux);

                            //in_verts[in_tris[3 * t_id]]->getExactXYZCoordinates(x0_aux, y0_aux, z0_aux);
                            //in_verts[in_tris[3 * t_id + 1]]->getExactXYZCoordinates(x1_aux, y1_aux, z1_aux);
                            //in_verts[in_tris[3 * t_id + 2]]->getExactXYZCoordinates(x2_aux, y2_aux, z2_aux);

                            std::vector<bigrational> tv0_aux = {x0_aux, y0_aux, z0_aux};
                            std::vector<bigrational> tv1_aux = {x1_aux, y1_aux, z1_aux};
                            std::vector<bigrational> tv2_aux = {x2_aux, y2_aux, z2_aux};


                            //printInfoTriangleRationals(rational_ray, tv0_aux, tv1_aux, tv2_aux, tv_aux, t_id, print_ray);


                            if (checkTriangleOrientationRationals(rational_ray, tv0_aux, tv1_aux, tv2_aux) == 1) {
                                patch_inner_label[t_label] = true;
                            }
                            visited_labels[t_label] = true;
                        }
                    }
                }*/

                propagateInnerLabelsOnPatch(patch_tris, patch_inner_label, labels);

                /*for (uint t_id_app: patch_tris){
                    labels.inside[t_id_app] = patch_inner_label;
                    std::cout<<"label.inside["<< t_id_app <<"]: "<< patch_inner_label << std::endl;
                }*/

            //}
        }else {

            cinolib::AABB rayAABB(cinolib::vec3d(ray.v0.X(), ray.v0.Y(), ray.v0.Z()),
                                  cinolib::vec3d(ray.v1.X(), ray.v1.Y(), ray.v1.Z()));

            intersects_box(octree, rayAABB, tmp_inters);
            std::vector<uint> sorted_inters;
            pruneIntersectionsAndSortAlongRay(ray, in_verts, in_tris, in_labels, tmp_inters, patch_surface_label,
                                              sorted_inters);

            std::bitset<NBIT> patch_inner_label;
            analyzeSortedIntersections(ray, in_verts, in_tris, in_labels, sorted_inters, patch_inner_label);
            propagateInnerLabelsOnPatch(patch_tris, patch_inner_label, labels);
        }
    }
    //);
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void pruneIntersectionsAndSortAlongRay(const Ray &ray, const std::vector<genericPoint*> &in_verts,
                                              const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                                              const phmap::flat_hash_set<uint> &tmp_inters, const std::bitset<NBIT> &patch_surface_label,
                                              std::vector<uint> &inters_tris)
{
    phmap::flat_hash_set<uint> visited_tri;
    visited_tri.reserve(tmp_inters.size()/6);
    std::pair<phmap::flat_hash_set<uint>::iterator, bool> ins;

    for(uint t_id : tmp_inters)
    {
        ins = visited_tri.insert(t_id);
        if(!ins.second) continue; // triangle already analyzed or in the one ring of a vert or in the adj of an edge

        const std::bitset<NBIT> tested_tri_label = in_labels[t_id];
        uint uint_tri_label = bitsetToUint(tested_tri_label);
        if(patch_surface_label[uint_tri_label]) continue; // <-- triangle of the same label of the tested patch

        //genericPoints in floating point

        const explicitPoint3D &tv0 = in_verts[in_tris[3 * t_id]]->toExplicit3D();
        const explicitPoint3D &tv1 = in_verts[in_tris[3 * t_id +1]]->toExplicit3D();
        const explicitPoint3D &tv2 = in_verts[in_tris[3 * t_id +2]]->toExplicit3D();

        IntersInfo ii = fast2DCheckIntersectionOnRay(ray, tv0, tv1, tv2);

        if(ii == DISCARD || ii == NO_INT) continue;

        if(ii == INT_IN_TRI)
        {
            inters_tris.push_back(t_id);
        }
        else if(ii == INT_IN_V0 || ii == INT_IN_V1 || ii == INT_IN_V2)
        {
            uint v_id;
            if(ii == INT_IN_V0) v_id = in_tris[3 * t_id];
            else if(ii == INT_IN_V1) v_id = in_tris[3 * t_id +1];
            else v_id = in_tris[3 * t_id +2];

            std::vector<uint> vert_one_ring;
            findVertRingTris(v_id, tested_tri_label, tmp_inters, in_tris, in_labels, vert_one_ring);

            for(uint t : vert_one_ring)
                visited_tri.insert(t); // mark all the one ring as visited

            int winner_tri = -1;
            winner_tri = perturbRayAndFindIntersTri(ray, in_verts, in_tris, vert_one_ring); // the first inters triangle after ray perturbation

            if(winner_tri != -1)
                inters_tris.push_back(winner_tri);
        }
        else if(ii == INT_IN_EDGE01 || ii == INT_IN_EDGE12 || ii == INT_IN_EDGE20)
        {
            uint ev0_id, ev1_id;
            if(ii == INT_IN_EDGE01)
            {
                ev0_id = in_tris[3 * t_id];
                ev1_id = in_tris[3 * t_id +1];
            }
            else if(ii == INT_IN_EDGE12)
            {
                ev0_id = in_tris[3 * t_id +1];
                ev1_id = in_tris[3 * t_id +2];
            }
            else
            {
                ev0_id = in_tris[3 * t_id +2];
                ev1_id = in_tris[3 * t_id];
            }

            std::vector<uint> edge_tris;
            findEdgeTris(ev0_id, ev1_id, tested_tri_label, tmp_inters, in_tris, in_labels, edge_tris);

            for(uint t : edge_tris)
                visited_tri.insert(t); // mark all the one ring as visited

            int winner_tri = -1;
            winner_tri = perturbRayAndFindIntersTri(ray, in_verts, in_tris, edge_tris);

            if(winner_tri != -1)
                inters_tris.push_back(winner_tri);
        }
    }

    if(ray.dir == 'X')
        sortIntersectedTrisAlongX(ray, in_verts, in_tris, inters_tris);
    else if(ray.dir == 'Y')
        sortIntersectedTrisAlongY(ray, in_verts, in_tris, inters_tris);
    else
        sortIntersectedTrisAlongZ(ray, in_verts, in_tris, inters_tris);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void pruneIntersectionsAndSortAlongRayRationals(const RationalRay &ray, const FastTrimesh &tm, const std::vector<genericPoint*> &in_verts,
                                                       const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                                                       const phmap::flat_hash_set<uint> &tmp_inters, const std::bitset<NBIT> &patch_surface_label,
                                                       std::vector<IntersectionPointRationals> &inter_rat, std::vector<uint> &inters_tris_rat, Labels &labels){

    phmap::flat_hash_set<uint> visited_tri;
    visited_tri.reserve(inter_rat.size() / 6);
    std::pair<phmap::flat_hash_set<uint>::iterator, bool> ins;


    for (uint t_id_int: tmp_inters) {
        //uint t_id_app = t_id_int;
        ins = visited_tri.insert(t_id_int);

        if (!ins.second) continue; // triangle already analyzed or in the one ring of a vert or in the adj of an edge

        /*int t_id_father = -1;


        if (t_id_int >= in_labels.size()){
           uint patch_id= IntersectionPointRationals::findPatchIdByTriId(inter_rat, t_id_int);
           phmap::flat_hash_set<uint> patch_app = patches.at(patch_id);
            for(uint t_id : patch_app){
                if (t_id < in_labels.size()){
                    t_id_father = t_id;
                    break;
                }
            }
            t_id_int = t_id_father;
        }*/


        const std::bitset<NBIT> tested_tri_label = labels.surface.at(t_id_int);
        uint uint_tri_label = bitsetToUint(tested_tri_label);
        //t_id_int = t_id_app;

        if (patch_surface_label[uint_tri_label]) continue; // <-- triangle of the same label of the tested patch

        bigrational tv0_x, tv0_y, tv0_z,
                tv1_x, tv1_y, tv1_z,
                tv2_x, tv2_y, tv2_z;

        const uint tv[3] = {tm.triVertID(t_id_int, 0), tm.triVertID(t_id_int, 1), tm.triVertID(t_id_int, 2)};

        in_verts[tv[0]]->getExactXYZCoordinates(tv0_x, tv0_y, tv0_z);
        in_verts[tv[1]]->getExactXYZCoordinates(tv1_x, tv1_y, tv1_z);
        in_verts[tv[2]]->getExactXYZCoordinates(tv2_x, tv2_y, tv2_z);


        //print size in tris
        //in_verts[in_tris[3 * t_id_int]]->getExactXYZCoordinates(tv0_x, tv0_y, tv0_z);
        //in_verts[in_tris[3 * t_id_int + 1]]->getExactXYZCoordinates(tv1_x, tv1_y, tv1_z);
        //in_verts[in_tris[3 * t_id_int + 2]]->getExactXYZCoordinates(tv2_x, tv2_y, tv2_z);

        const std::vector<bigrational> tv0_exact = {tv0_x, tv0_y, tv0_z};
        const std::vector<bigrational> tv1_exact = {tv1_x, tv1_y, tv1_z};
        const std::vector<bigrational> tv2_exact = {tv2_x, tv2_y, tv2_z};

        std::cout <<"t_id_int before fast2DIntersectionOnRayRationals: " << t_id_int << std::endl;

        IntersInfo ii = fast2DCheckIntersectionOnRayRationals(ray, tv0_exact, tv1_exact, tv2_exact);

        std::cout<< "ii: " << ii << std::endl;

        if (ii == DISCARD){
            std::cout<<"DISCARD" << std::endl;
            continue;}
        if (ii == NO_INT){
            std::cout<<"NO_INT" << std::endl;
            continue;}

        if (ii == INT_IN_TRI) {
            std::cout << "INT_IN_TRI" << " t_id-> "<<t_id_int<<std::endl;

            inters_tris_rat.push_back(t_id_int);

        } else if (ii == INT_IN_V0 || ii == INT_IN_V1 || ii == INT_IN_V2) {
            std::cout << "INT_IN_V" << std::endl;
            uint v_id;
            if (ii == INT_IN_V0) v_id = in_tris[3 * t_id_int];
            else if (ii == INT_IN_V1) v_id = in_tris[3 * t_id_int + 1];
            else v_id = in_tris[3 * t_id_int + 2];

            std::vector<uint> vert_one_ring;
            findVertRingTris(v_id, tested_tri_label, tmp_inters, in_tris, in_labels, vert_one_ring);

            for (uint t: vert_one_ring)
                visited_tri.insert(t); // mark all the one ring as visited

            int winner_tri = -1;
            winner_tri = perturbRayAndFindIntersTri(ray, in_verts, in_tris,
                                                    vert_one_ring); // the first inters triangle after ray perturbation

            if (winner_tri != -1)
                inters_tris_rat.push_back(winner_tri);

        } else if (ii == INT_IN_EDGE01 || ii == INT_IN_EDGE12 || ii == INT_IN_EDGE20) {
            std::cout << "INT_IN_EDGE" << std::endl;

            int ev0_id, ev1_id;
            if (ii == INT_IN_EDGE01) {
                ev0_id = in_tris[3 * t_id_int];
                ev1_id = in_tris[3 * t_id_int + 1];
            } else if (ii == INT_IN_EDGE12) {
                ev0_id = in_tris[3 * t_id_int + 1];
                ev1_id = in_tris[3 * t_id_int + 2];
            } else {
                ev0_id = in_tris[3 * t_id_int + 2];
                ev1_id = in_tris[3 * t_id_int];
            }

            std::cout << "t_id_int: " << t_id_int<<std::endl;
            std::cout << "ev0_id: " << ev0_id << std::endl;
            std::cout << "ev1_id: " << ev1_id << std::endl;

            std::vector<uint> edge_tris;

            findEdgeTrisCustom(ev0_id, ev1_id, tested_tri_label, tmp_inters, in_tris, in_labels, edge_tris);

            for (uint t: edge_tris)
                visited_tri.insert(t); // mark all the one ring as visited

            int winner_tri = -1;
            winner_tri = perturbRayAndFindIntersTri(ray, in_verts, in_tris, edge_tris);

            if (winner_tri != -1)
                inters_tris_rat.push_back(winner_tri);

        }

        //Sanity check if there is an intersection point that it's found before
        for (auto &t_id_tmp : inters_tris_rat) {
            bool found = false;
            for (auto &p_int: inter_rat) {
                if (t_id_tmp == p_int.getTriId()) {
                    found = true;
                    break;
                }
            }
            assert(found && "The triangle is not in the intersection points");
        }

        std::cout << "inter_rat size: " << inter_rat.size() << std::endl;
        if (ray.dir == 'X')
            //sort the inter_rat along the x coordinate
            sort(inter_rat.begin(), inter_rat.end(),
                 [](const IntersectionPointRationals &a, const IntersectionPointRationals &b) {
                     return a.lessThanX(b);
                 });
        else if (ray.dir == 'Y')
            sort(inter_rat.begin(), inter_rat.end(),
                 [](const IntersectionPointRationals &a, const IntersectionPointRationals &b) {
                     return a.lessThanY(b);
                 });
        else
            sort(inter_rat.begin(), inter_rat.end(),
                 [](const IntersectionPointRationals &a, const IntersectionPointRationals &b) {
                     return a.lessThanZ(b);
                 });

    }
}




//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::....

inline void analyzeSortedIntersections(const Ray &ray, const std::vector<genericPoint*> &in_verts, const std::vector<uint> &in_tris,
                                       const std::vector<std::bitset<NBIT>> &in_labels, const std::vector<uint> &sorted_inters,
                                       std::bitset<NBIT> &patch_inner_label)
{
    std::bitset<NBIT> visited_labels;

    for(uint t_id : sorted_inters)
    {
        uint t_label = bitsetToUint(in_labels[t_id]);
        if(visited_labels[t_label]) continue; // already visited patch

        const explicitPoint3D &tv0 = in_verts[in_tris[3 * t_id]]->toExplicit3D();
        const explicitPoint3D &tv1 = in_verts[in_tris[3 * t_id +1]]->toExplicit3D();
        const explicitPoint3D &tv2 = in_verts[in_tris[3 * t_id +2]]->toExplicit3D();

        if(checkTriangleOrientation(ray, tv0, tv1, tv2) == 1) // checkOrientation -> 1 if inside, 0 if outside
            patch_inner_label[t_label] = true;

        visited_labels[t_label] = true;
    }
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void analyzeSortedIntersectionsRationals(const RationalRay &rational_ray, const FastTrimesh &tm, const std::vector<genericPoint*> &in_verts,
                                               std::vector<IntersectionPointRationals> &inter_rat, std::bitset<NBIT> &patch_inner_label, Labels &labels){

    std::bitset<NBIT> visited_labels;


    //copy the tri_id into inter_rat_tris_ordered from inter_rat
    /*std::vector<uint> inter_rat_tris_ordered;

    for (auto &p: inter_rat) {
        inter_rat_tris_ordered.push_back(p.getTriId());
    }*/
    bool print_ray = true;

    for (int i = 0; i < inter_rat.size(); i++) {

        uint t_id = inter_rat[i].getTriId();
        uint t_label = bitsetToUint(labels.surface.at(t_id));

        if (visited_labels[t_label]) continue;

        for (auto &p: inter_rat) {
            if (p.getTriId() == t_id) {
                //get the three points of the triangle
                bigrational x0_aux, x1_aux, x2_aux, y0_aux, y1_aux, y2_aux, z0_aux, z1_aux, z2_aux;

                uint tv_aux[3] = {tm.triVertID(t_id, 0), tm.triVertID(t_id, 1), tm.triVertID(t_id, 2)};

                in_verts[tv_aux[0]]->getExactXYZCoordinates(x0_aux, y0_aux, z0_aux);
                in_verts[tv_aux[1]]->getExactXYZCoordinates(x1_aux, y1_aux, z1_aux);
                in_verts[tv_aux[2]]->getExactXYZCoordinates(x2_aux, y2_aux, z2_aux);

                std::vector<bigrational> tv0_aux = {x0_aux, y0_aux, z0_aux};
                std::vector<bigrational> tv1_aux = {x1_aux, y1_aux, z1_aux};
                std::vector<bigrational> tv2_aux = {x2_aux, y2_aux, z2_aux};

                //printInfoTriangleRationals(rational_ray, tv0_aux, tv1_aux, tv2_aux, tv_aux, t_id, print_ray);

                if (checkTriangleOrientationRationals(rational_ray, tv0_aux, tv1_aux, tv2_aux) == 1) {
                    patch_inner_label[t_label] = true;
                }
                visited_labels[t_label] = true;
            }
        }
    }

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline bool triContainsVert(uint t_id, uint v_id, const std::vector<uint> &in_tris)
{
    if(in_tris[3 * t_id]     == v_id) return true;
    if(in_tris[3 * t_id + 1] == v_id) return true;
    if(in_tris[3 * t_id + 2] == v_id) return true;
    return false;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void findVertRingTris(uint v_id, const std::bitset<NBIT> &ref_label, const phmap::flat_hash_set<uint> &inters_tris,
                             const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                             std::vector<uint> &one_ring)
{
    for(uint t_id : inters_tris)
    {
        if(in_labels[t_id] == ref_label && triContainsVert(t_id, v_id, in_tris))
            one_ring.push_back(t_id);
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void findEdgeTris(uint ev0_id, uint ev1_id, const std::bitset<NBIT> &ref_label, const phmap::flat_hash_set<uint> &inters_tris,
                         const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                         std::vector<uint> &edge_tris)
{
    for(auto t_id : inters_tris)
    {
        if(in_labels[t_id] == ref_label  && triContainsVert(t_id, ev0_id, in_tris) && triContainsVert(t_id, ev1_id, in_tris))
            edge_tris.push_back(t_id);
    }

    assert(edge_tris.size() == 2 && "problem in finding edge triangles"); // always true in closed and manifold meshes
}

inline void findEdgeTrisCustom(uint ev0_id, uint ev1_id, const std::bitset<NBIT> &ref_label, const phmap::flat_hash_set<uint> &inters_tris,
                         const std::vector<uint> &in_tris, const std::vector<std::bitset<NBIT>> &in_labels,
                         std::vector<uint> &edge_tris)
{
    for(auto t_id : inters_tris)
    {
        std::cout << t_id << std::endl;
        std::cout << in_labels[t_id] << std::endl;
        std::cout << ref_label << std::endl;

        if(in_labels[t_id] == ref_label  && triContainsVert(t_id, ev0_id, in_tris) && triContainsVert(t_id, ev1_id, in_tris))
            edge_tris.push_back(t_id);
    }
    std::cout << "edge_tris.size(): " << edge_tris.size() << std::endl;

    assert(edge_tris.size() == 2 && "problem in finding edge triangles"); // always true in closed and manifold meshes
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline Ray perturbXRay(const Ray &ray, uint offset) // offset is used to perturb the ray in all the possible directions
{
    Ray new_ray = ray;

    switch (offset)
    {
        case 0: // -> +y
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, ray.v1.Z());
        } break;

        case 1: // -> +y+z
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, new_z);
        } break;

        case 2: // -> +z
        {
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), ray.v1.Y(), new_z);
        } break;

        case 3: //-> -y+z
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, new_z);
        } break;

        case 4: // -> -y
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, ray.v1.Z());
        } break;

        case 5: // -> -y-z
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, new_z);
        } break;

        case 6: // -> -z
        {
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), ray.v1.Y(), new_z);
        } break;

        case 7: // -> +y-z
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, new_z);
        } break;

        default:
        {
            assert(false && "non-valid offset value");
        } break;
    }

    return new_ray;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline Ray perturbYRay(const Ray &ray, uint offset)
{
    Ray new_ray = ray;

    switch (offset)
    {
        case 0: // -> +x
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() + 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), ray.v1.Z());
        } break;

        case 1: // -> +x+z
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() + 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), new_z);
        } break;

        case 2: // -> +z
        {
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), ray.v1.Y(), new_z);
        } break;

        case 3: //-> -x+z
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() - 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), new_z);
        } break;

        case 4: // -> -x
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() - 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), ray.v1.Z());
        } break;

        case 5: // -> -x-z
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() - 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), new_z);
        } break;

        case 6: // -> -z
        {
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), ray.v1.Y(), new_z);
        } break;

        case 7: // -> +x-z
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() + 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), new_z);
        } break;

        default:
        {
            assert(false && "non-valid offset value");
        } break;
    }

    return new_ray;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline Ray perturbZRay(const Ray &ray, uint offset)
{
    Ray new_ray = ray;

    switch (offset)
    {
        case 0: // -> +x
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() + 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), ray.v1.Z());
        } break;

        case 1: // -> +x+y
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() + 1.0));
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            new_ray.v1 = explicitPoint3D(new_x, new_y, ray.v1.Z());
        } break;

        case 2: // -> +y
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, ray.v1.Z());
        } break;

        case 3: //-> -x+y
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() - 1.0));
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            new_ray.v1 = explicitPoint3D(new_x, new_y, ray.v1.Z());
        } break;

        case 4: // -> -x
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() - 1.0));
            new_ray.v1 = explicitPoint3D(new_x, ray.v1.Y(), ray.v1.Z());
        } break;

        case 5: // -> -x-y
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() - 1.0));
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            new_ray.v1 = explicitPoint3D(new_x, new_y, ray.v1.Z());
        } break;

        case 6: // -> -y
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, ray.v1.Z());
        } break;

        case 7: // -> +x-y
        {
            double new_x = std::nextafter(ray.v1.X(), (ray.v1.X() + 1.0));
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            new_ray.v1 = explicitPoint3D(new_x, new_y, ray.v1.Z());
        } break;

        default:
        {
            assert(false && "non-valid offset value");
        } break;
    }

    return new_ray;

}


/*******
inline Ray perturbXRay(const RationalRay &ray, uint offset) // offset is used to perturb the ray in all the possible directions
{
    RationalRay new_ray = ray;

    switch (offset)
    {
        case 0: // -> +y
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, ray.v1.Z());
        } break;

        case 1: // -> +y+z
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, new_z);
        } break;

        case 2: // -> +z
        {
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), ray.v1.Y(), new_z);
        } break;

        case 3: //-> -y+z
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() + 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, new_z);
        } break;

        case 4: // -> -y
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, ray.v1.Z());
        } break;

        case 5: // -> -y-z
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() - 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, new_z);
        } break;

        case 6: // -> -z
        {
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), ray.v1.Y(), new_z);
        } break;

        case 7: // -> +y-z
        {
            double new_y = std::nextafter(ray.v1.Y(), (ray.v1.Y() + 1.0));
            double new_z = std::nextafter(ray.v1.Z(), (ray.v1.Z() - 1.0));
            new_ray.v1 = explicitPoint3D(ray.v1.X(), new_y, new_z);
        } break;

        default:
        {
            assert(false && "non-valid offset value");
        } break;
    }

    return new_ray;
}
********/

inline int perturbRayAndFindIntersTri(const Ray &ray, const std::vector<genericPoint*> &in_verts, const std::vector<uint> &in_tris,
                                       const std::vector<uint> &tris_to_test)
{
    std::vector<uint> inters_tris;
    Ray p_ray;

    for(uint i = 0; i <= 7; i++)
    {
        if(ray.dir == 'X')      p_ray = perturbXRay(ray, i);
        else if(ray.dir == 'Y') p_ray = perturbYRay(ray, i);
        else if(ray.dir == 'Z') p_ray = perturbZRay(ray, i);

        for(uint t_id : tris_to_test)
        {
            const explicitPoint3D &tv0 = in_verts[in_tris[3 * t_id]]->toExplicit3D();
            const explicitPoint3D &tv1 = in_verts[in_tris[3 * t_id +1]]->toExplicit3D();
            const explicitPoint3D &tv2 = in_verts[in_tris[3 * t_id +2]]->toExplicit3D();

            if(checkIntersectionInsideTriangle3D(p_ray, tv0, tv1, tv2))
                inters_tris.push_back(t_id);

            if(!inters_tris.empty()) break;
        }
    }

    if(inters_tris.empty())
        return -1;

    if(ray.dir == 'X')      sortIntersectedTrisAlongX(p_ray, in_verts, in_tris, inters_tris);
    else if(ray.dir == 'Y') sortIntersectedTrisAlongY(p_ray, in_verts, in_tris, inters_tris);
    else                    sortIntersectedTrisAlongZ(p_ray, in_verts, in_tris, inters_tris);

    return static_cast<int>(inters_tris[0]); // return the first triangle intersected
}




inline int perturbRayAndFindIntersTri(const RationalRay &ray, const std::vector<genericPoint*> &in_verts, const std::vector<uint> &in_tris,
                                      const std::vector<uint> &tris_to_test)
{
    std::vector<uint> inters_tris;
    RationalRay p_ray;

    //perturb the ray using a centroid of one of the triangle adjacent to the edge
    bigrational x0, y0, z0, x1, y1, z1, x2, y2, z2;
    in_verts[in_tris[3 * tris_to_test[0]]]->getExactXYZCoordinates(x0, y0, z0);
    in_verts[in_tris[3 * tris_to_test[0] + 1]]->getExactXYZCoordinates(x1, y1, z1);
    in_verts[in_tris[3 * tris_to_test[0] + 2]]->getExactXYZCoordinates(x2, y2, z2);

    std::vector<bigrational> p0 = {x0, y0, z0};
    std::vector<bigrational> p1 = {x1, y1, z1};
    std::vector<bigrational> p2 = {x2, y2, z2};


    p_ray.v1 = {(x0 + x1 + x2) / bigrational(3.0),
                (y0 + y1 + y2) / bigrational(3.0) ,
                (z0 + z1 + z2) / bigrational(3.0)};

    for(uint t_id : tris_to_test)
    {

        bigrational x0_aux, x1_aux, x2_aux, y0_aux, y1_aux, y2_aux, z0_aux, z1_aux, z2_aux;
        in_verts[in_tris[3 * t_id]]->getExactXYZCoordinates(x0_aux, y0_aux, z0_aux);
        in_verts[in_tris[3 * t_id +1]]->getExactXYZCoordinates(x1_aux, y1_aux, z1_aux);
        in_verts[in_tris[3 * t_id +2]]->getExactXYZCoordinates(x2_aux, y2_aux, z2_aux);

        std::vector<bigrational> tv0 = {x0_aux, y0_aux, z0_aux};
        std::vector<bigrational> tv1 = {x1_aux, y1_aux, z1_aux};
        std::vector<bigrational> tv2 = {x2_aux, y2_aux, z2_aux};


        if(checkIntersectionInsideTriangle3DRationals(p_ray, tv0, tv1, tv2))
            inters_tris.push_back(t_id);

        if(!inters_tris.empty()) break;
    }


    //TODO: check se centroide sta sullo stesso piano dei vertici del triangolo
    if(cinolib::orient3d(&p0[0], &p1[0], &p2[0], &p_ray.v1[0]).sgn() != 0){
        std::cout << "The centroid is not on the triangle or the orient3d doesn't work properly" << std::endl;
        std::exit(EXIT_FAILURE);
    } //the centroid is exaclty on the plane of the triangle


    //TODO: controlla che le tre orient siano di segno concorde per controllare che il centroide sta all'interno
    double e0_rat = cinolib::orient3d(&p0[0], &p1[0], &ray.v1[0], &p_ray.v1[0]).get_d();
    double e1_rat = cinolib::orient3d(&p1[0], &p2[0], &ray.v1[0], &p_ray.v1[0]).get_d();
    double e2_rat = cinolib::orient3d(&p2[0], &p0[0], &ray.v1[0], &p_ray.v1[0]).get_d();

    if(!(e0_rat > 0 && e1_rat > 0 && e2_rat > 0) && !(e0_rat < 0 && e1_rat < 0 && e2_rat < 0)){
        std::cout << "The centroid is not inside the triangle or the orient3D doesn't work properly" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    inters_tris.push_back(tris_to_test[0]);

    //TODO: Ordina i triangoli in base alla distanza dal punto di intersezione
    if(inters_tris.empty())
        return -1;

    std::vector<IntersectionPointRationals> inter_rat;
    for(uint t_id : tris_to_test)
    {
        bigrational x0_aux, x1_aux, x2_aux, y0_aux, y1_aux, y2_aux, z0_aux, z1_aux, z2_aux;
        in_verts[in_tris[3 * t_id]]->getExactXYZCoordinates(x0_aux, y0_aux, z0_aux);
        in_verts[in_tris[3 * t_id +1]]->getExactXYZCoordinates(x1_aux, y1_aux, z1_aux);
        in_verts[in_tris[3 * t_id +2]]->getExactXYZCoordinates(x2_aux, y2_aux, z2_aux);

        std::vector<bigrational> tv0 = {x0_aux, y0_aux, z0_aux};
        std::vector<bigrational> tv1 = {x1_aux, y1_aux, z1_aux};
        std::vector<bigrational> tv2 = {x2_aux, y2_aux, z2_aux};

        if(segment_triangle_intersect_3d(&ray.v0[0], &ray.v1[0], &tv0[0], &tv1[0], &tv2[0])){
            std::vector<bigrational> p_int(3);
            plane_line_intersection(&tv0[0], &tv1[0], &tv2[0], &ray.v0[0], &ray.v1[0], &p_int[0]);
            IntersectionPointRationals p_int_rat(p_int[0], p_int[1], p_int[2], t_id,0);
            inter_rat.push_back(p_int_rat);
        }
    }

    if(ray.dir == 'X'){
        sort(inter_rat.begin(), inter_rat.end(), [](const IntersectionPointRationals &a, const IntersectionPointRationals &b){
            return a.lessThanX(b);
        });
    }
    else if(ray.dir == 'Y'){
        sort(inter_rat.begin(), inter_rat.end(), [](const IntersectionPointRationals &a, const IntersectionPointRationals &b){
            return a.lessThanY(b);
        });
    }
    else{
        sort(inter_rat.begin(), inter_rat.end(), [](const IntersectionPointRationals &a, const IntersectionPointRationals &b){
            return a.lessThanZ(b);
        });
    }

    return static_cast<int>(inter_rat[0].getTriId()); // return the first triangle intersected
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline IntersInfo fast2DCheckIntersectionOnRay(const Ray &ray, const explicitPoint3D &tv0, const explicitPoint3D &tv1, const explicitPoint3D &tv2)
{
    double v0[2], v1[2], v2[2], vq[2];

    switch (ray.dir)
    {
        case 'X': // only YZ coordinates
        {
            v0[0] = tv0.Y(); v0[1] = tv0.Z();
            v1[0] = tv1.Y(); v1[1] = tv1.Z();
            v2[0] = tv2.Y(); v2[1] = tv2.Z();
            vq[0] = ray.v1.Y(); vq[1] = ray.v1.Z();
        } break;

        case 'Y': //only XZ coordinates
        {
            v0[0] = tv0.X(); v0[1] = tv0.Z();
            v1[0] = tv1.X(); v1[1] = tv1.Z();
            v2[0] = tv2.X(); v2[1] = tv2.Z();
            vq[0] = ray.v1.X(); vq[1] = ray.v1.Z();
        } break;

        case 'Z': //only XY coordinates
        {
            v0[0] = tv0.X(); v0[1] = tv0.Y();
            v1[0] = tv1.X(); v1[1] = tv1.Y();
            v2[0] = tv2.X(); v2[1] = tv2.Y();
            vq[0] = ray.v1.X(); vq[1] = ray.v1.Y();
        } break;

    }

    double or01 = cinolib::orient2d(v0, v1, vq);
    double or12 = cinolib::orient2d(v1, v2, vq);
    double or20 = cinolib::orient2d(v2, v0, vq);

    if((or01 >= 0 && or12 >= 0 && or20 >= 0) || (or01 <= 0 && or12 <= 0 && or20 <= 0))
    {
        // check if the the ray passes through a vert
        if(v0[0] == vq[0] && v0[1] == vq[1]) return INT_IN_V0;
        if(v1[0] == vq[0] && v1[1] == vq[1]) return INT_IN_V1;
        if(v2[0] == vq[0] && v2[1] == vq[1]) return INT_IN_V2;

        // check if the triangle is coplanar with the ray
        if(or01 == 0 && or12 == 0) return DISCARD;
        if(or12 == 0 && or20 == 0) return DISCARD;
        if(or20 == 0 && or01 == 0) return DISCARD;

        // check if the ray passes through an edge
        if(or01 == 0) return INT_IN_EDGE01;
        if(or12 == 0) return INT_IN_EDGE12;
        if(or20 == 0) return INT_IN_EDGE20;

        return INT_IN_TRI; // so the triangle intersect insede the triangle area
    }

    return NO_INT;
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline IntersInfo fast2DCheckIntersectionOnRayRationals(const RationalRay &ray, const std::vector<bigrational> &tv0, const std::vector<bigrational> &tv1, const std::vector<bigrational> &tv2)
{
    bigrational v0_rat[2], v1_rat[2], v2_rat[2], vq_rat[2];

    bigrational ray_v1_x, ray_v1_y, ray_v1_z;

    ray_v1_x = ray.v1[0];
    ray_v1_y = ray.v1[1];
    ray_v1_z = ray.v1[2];

    switch (ray.dir)
    {
        case 'X': // only YZ coordinates
        {
            v0_rat[0] = tv0[1]; v0_rat[1] = tv0[2];
            v1_rat[0] = tv1[1]; v1_rat[1] = tv1[2];
            v2_rat[0] = tv2[1]; v2_rat[1] = tv2[2];

            vq_rat[0] = ray_v1_y; vq_rat[1] = ray_v1_z;
        } break;

        case 'Y': //only XZ coordinates
        {
            v0_rat[0] = tv0[0]; v0_rat[1] = tv0[2];
            v1_rat[0] = tv1[0]; v1_rat[1] = tv1[2];
            v2_rat[0] = tv2[0]; v2_rat[1] = tv2[2];
            vq_rat[0] = ray_v1_x; vq_rat[1] = ray_v1_z;
        } break;

        case 'Z': //only XY coordinates
        {
            v0_rat[0] = tv0[0]; v0_rat[1] = tv0[1];
            v1_rat[0] = tv1[0]; v1_rat[1] = tv1[1];
            v2_rat[0] = tv2[0]; v2_rat[1] = tv2[1];
            vq_rat[0] = ray_v1_x; vq_rat[1] = ray_v1_y;
        } break;

    }

    bigrational or01_rat = cinolib::orient2d(&v0_rat[0], &v1_rat[0], &vq_rat[0]);
    bigrational or12_rat = cinolib::orient2d(&v1_rat[0], &v2_rat[0], &vq_rat[0]);
    bigrational or20_rat = cinolib::orient2d(&v2_rat[0], &v0_rat[0], &vq_rat[0]);
    bigrational zero_rat = bigrational(0,0,0);

    std::cout << or01_rat << " " << or12_rat << " " << or20_rat << std::endl;
    //TODO: IMPORTANT!! If the orientation is not correct we need to invert the boolean test like >= to <=
    if(((or01_rat < zero_rat || or01_rat.sgn() == 0) && (or12_rat < zero_rat || or12_rat.sgn() == 0) && (or20_rat < zero_rat || or20_rat.sgn() == 0)) ||
       ((or01_rat > zero_rat || or01_rat.sgn() == 0) && (or12_rat > zero_rat || or12_rat.sgn() == 0) && (or20_rat > zero_rat || or20_rat.sgn() == 0)))
    {
        if(v0_rat[0] == vq_rat[0] && v0_rat[1] == vq_rat[1]) return INT_IN_V0;
        if(v1_rat[0] == vq_rat[0] && v1_rat[1] == vq_rat[1]) return INT_IN_V1;
        if(v2_rat[0] == vq_rat[0] && v2_rat[1] == vq_rat[1]) return INT_IN_V2;

        if(or01_rat.sgn() == 0 && or12_rat.sgn() == 0) return DISCARD;
        if(or12_rat.sgn() == 0 && or20_rat.sgn() == 0) return DISCARD;
        if(or20_rat.sgn() == 0 && or01_rat.sgn() == 0) return DISCARD;

        if(or01_rat.sgn() == 0) return INT_IN_EDGE01;
        if(or12_rat.sgn() == 0) return INT_IN_EDGE12;
        if(or20_rat.sgn() == 0) return INT_IN_EDGE20;

        return INT_IN_TRI; // so the triangle intersect insede the triangle area
    }

    return NO_INT;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline bool checkIntersectionInsideTriangle3D(const Ray &ray, const explicitPoint3D &tv0, const explicitPoint3D &tv1, const explicitPoint3D &tv2)
{
    // we check the orientation of ray.v1 with respct to the planes v0-v1-ray.v0, v1-v2-ray.v0, v2-v0-ray.v0
    double or01f = cinolib::orient3d(tv0.ptr(), tv1.ptr(), ray.v0.ptr(), ray.v1.ptr());
    double or12f = cinolib::orient3d(tv1.ptr(), tv2.ptr(), ray.v0.ptr(), ray.v1.ptr());
    double or20f = cinolib::orient3d(tv2.ptr(), tv0.ptr(), ray.v0.ptr(), ray.v1.ptr());

    if(or01f > 0 && or12f > 0 && or20f > 0) return true;
    if(or01f < 0 && or12f < 0 && or20f < 0) return true;

    return false;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline bool checkIntersectionInsideTriangle3DRationals(const RationalRay &ray, const std::vector<bigrational> &tv0, const std::vector<bigrational> &tv1, const std::vector<bigrational> &tv2)
{


    bigrational or01f_rat = cinolib::orient3d(&tv0[0], &tv1[0], &ray.v0[0], &ray.v1[0]);
    bigrational or12f_rat = cinolib::orient3d(&tv1[0], &tv2[0], &ray.v0[0], &ray.v1[0]);
    bigrational or20f_rat = cinolib::orient3d(&tv2[0], &tv0[0], &ray.v0[0], &ray.v1[0]);

    if(or01f_rat > bigrational(0,0,0) && or12f_rat > bigrational(0,0,0) && or20f_rat > bigrational(0,0,0)) return true;
    if(or01f_rat < bigrational(0,0,0) && or12f_rat < bigrational(0,0,0) && or20f_rat < bigrational(0,0,0)) return true;

    return false;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline bool checkIntersectionInsideTriangle3DImplPoints(const Ray &ray, const genericPoint *tv0, const genericPoint *tv1, const genericPoint *tv2)
{
    // we check the orientation of ray.v1 with respct to the planes v0-v1-ray.v0, v1-v2-ray.v0, v2-v0-ray.v0
    double or01f = genericPoint::orient3D(*tv0, *tv1, ray.v0, ray.v1);
    double or12f = genericPoint::orient3D(*tv1, *tv2, ray.v0, ray.v1);
    double or20f = genericPoint::orient3D(*tv2, *tv0, ray.v0, ray.v1);

    if(or01f > 0 && or12f > 0 && or20f > 0) return true;
    if(or01f < 0 && or12f < 0 && or20f < 0) return true;

    return false;

}



inline bool checkIntersectionInsideTriangle3DImplPoints(const RationalRay &ray, const std::vector<bigrational> &tv0, const std::vector<bigrational> &tv1, const std::vector<bigrational> &tv2)
{
    // we check the orientation of ray.v1 with respct to the planes v0-v1-ray.v0, v1-v2-ray.v0, v2-v0-ray.v0
    double or01f = cinolib::orient3d(&tv0[0], &tv1[0], &ray.v0[0], &ray.v1[0]).get_d();
    double or12f = cinolib::orient3d(&tv1[0], &tv2[0], &ray.v0[0], &ray.v1[0]).get_d();
    double or20f = cinolib::orient3d(&tv2[0], &tv0[0], &ray.v0[0], &ray.v1[0]).get_d();

    if(or01f > 0 && or12f > 0 && or20f > 0) return true;
    if(or01f < 0 && or12f < 0 && or20f < 0) return true;

    return false;

}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// sort all intersected triangles from ray.v0 to ray.v1 (intersections before ray.v0 are discarded)
inline void sortIntersectedTrisAlongX(const Ray &ray, const std::vector<genericPoint*> &in_verts,
                                      const std::vector<uint> &in_tris, std::vector<uint> &inters_tris)
{
    phmap::btree_set<std::pair<genericPoint*, uint>, less_than_GP_on_X > inters_set; // <- <t_id, impl_point>
    std::vector<implicitPoint3D_LPI> arena;
    arena.reserve(inters_tris.size());

    for(uint t_id : inters_tris)
    {
        uint v0_id = in_tris[3 * t_id];
        uint v1_id = in_tris[3 * t_id +1];
        uint v2_id = in_tris[3 * t_id +2];

        std::pair<genericPoint*, uint> pair;
        pair.first = &arena.emplace_back(ray.v0, ray.v1,
                                         in_verts[v0_id]->toExplicit3D(),
                                         in_verts[v1_id]->toExplicit3D(),
                                         in_verts[v2_id]->toExplicit3D());
        pair.second = t_id;

        inters_set.insert(pair);
    }

    inters_tris.clear();
    auto curr_int = inters_set.begin();

    // we discard the intersection before ray.first along X
    if(ray.tv[0] != -1) // the ray is generated
    {
        const genericPoint *tv0 = in_verts[ray.tv[0]];
        const genericPoint *tv1 = in_verts[ray.tv[1]];
        const genericPoint *tv2 = in_verts[ray.tv[2]];

        if(genericPoint::orient3D(*tv0, *tv1, *tv2, ray.v1) > 0)
        {
            while(curr_int != inters_set.end() && genericPoint::orient3D(*tv0, *tv1, *tv2, *curr_int->first) < 0)
                curr_int++;
        }
        else
        {
            while(curr_int != inters_set.end() && genericPoint::orient3D(*tv0, *tv1, *tv2, *curr_int->first) > 0)
                curr_int++;
        }
    }
    else // the ray is composed of 2 real explicit points
    {
        while(curr_int != inters_set.end() && genericPoint::lessThanOnX(*curr_int->first, ray.v0) < 0)
            curr_int++;
    }

    // we save all the intersecting triangles from ray.first to ray.second
    while(curr_int != inters_set.end())
    {
        inters_tris.push_back(curr_int->second);
        curr_int++;
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void sortIntersectedTrisAlongY(const Ray &ray, const std::vector<genericPoint*> &in_verts,
                                      const std::vector<uint> &in_tris, std::vector<uint> &inters_tris)
{
    phmap::btree_set< std::pair<genericPoint*, uint>, less_than_GP_on_Y > inters_set; // <- <t_id, impl_point>
    std::vector<implicitPoint3D_LPI> arena;
    arena.reserve(inters_tris.size());

    for(uint t_id : inters_tris)
    {
        uint v0_id = in_tris[3 * t_id];
        uint v1_id = in_tris[3 * t_id +1];
        uint v2_id = in_tris[3 * t_id +2];

        std::pair<genericPoint*, uint> pair;
        pair.first = &arena.emplace_back(ray.v0, ray.v1,in_verts[v0_id]->toExplicit3D(), in_verts[v1_id]->toExplicit3D(), in_verts[v2_id]->toExplicit3D());
        pair.second = t_id;

        inters_set.insert(pair);
    }

    inters_tris.clear();
    auto curr_int = inters_set.begin();

    // we discard the intersection before ray.first along Y
    if(ray.tv[0] != -1) // the ray is generated
    {
        const genericPoint *tv0 = in_verts[ray.tv[0]];
        const genericPoint *tv1 = in_verts[ray.tv[1]];
        const genericPoint *tv2 = in_verts[ray.tv[2]];

        if(genericPoint::orient3D(*tv0, *tv1, *tv2, ray.v1) > 0)
        {
            while(curr_int != inters_set.end() && genericPoint::orient3D(*tv0, *tv1, *tv2, *curr_int->first) < 0)
                curr_int++;
        }
        else
        {
            while(curr_int != inters_set.end() && genericPoint::orient3D(*tv0, *tv1, *tv2, *curr_int->first) > 0)
                curr_int++;
        }
    }
    else // the ray is composed of 2 real explicit points
    {
        while(curr_int != inters_set.end() && genericPoint::lessThanOnY(*curr_int->first, ray.v0) < 0)
            curr_int++;
    }

    // we save all the intersecting triangles from ray.first to ray.second
    while(curr_int != inters_set.end())
    {
        inters_tris.push_back(curr_int->second);
        curr_int++;
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void sortIntersectedTrisAlongZ(const Ray &ray, const std::vector<genericPoint*> &in_verts,
                                      const std::vector<uint> &in_tris, std::vector<uint> &inters_tris)
{
    phmap::btree_set< std::pair<genericPoint*, uint>, less_than_GP_on_Z> inters_set; // <- <t_id, impl_point>
    std::vector<implicitPoint3D_LPI> arena;
    arena.reserve(inters_tris.size());

    for(uint t_id : inters_tris)
    {
        uint v0_id = in_tris[3 * t_id];
        uint v1_id = in_tris[3 * t_id +1];
        uint v2_id = in_tris[3 * t_id +2];

        std::pair<genericPoint*, uint> pair;
        pair.first = &arena.emplace_back(ray.v0, ray.v1,in_verts[v0_id]->toExplicit3D(), in_verts[v1_id]->toExplicit3D(), in_verts[v2_id]->toExplicit3D());
        pair.second = t_id;

        inters_set.insert(pair);
    }

    inters_tris.clear();
    auto curr_int = inters_set.begin();

    // we discard the intersection before ray.first along Z
    if(ray.tv[0] != -1) // the ray is generated
    {
        const genericPoint *tv0 = in_verts[ray.tv[0]];
        const genericPoint *tv1 = in_verts[ray.tv[1]];
        const genericPoint *tv2 = in_verts[ray.tv[2]];

        if(genericPoint::orient3D(*tv0, *tv1, *tv2, ray.v1) > 0)
        {
            while(curr_int != inters_set.end() && genericPoint::orient3D(*tv0, *tv1, *tv2, *curr_int->first) < 0)
                curr_int++;
        }
        else
        {
            while(curr_int != inters_set.end() && genericPoint::orient3D(*tv0, *tv1, *tv2, *curr_int->first) > 0)
                curr_int++;
        }
    }
    else // the ray is composed of 2 real explicit points
    {
        while(curr_int != inters_set.end() && genericPoint::lessThanOnZ(*curr_int->first, ray.v0) < 0)
            curr_int++;
    }

    // we save all the intersecting triangles from ray.first to ray.second
    while(curr_int != inters_set.end())
    {
        inters_tris.push_back(curr_int->second);
        curr_int++;
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// return 1 if inside, 0 if outside
inline uint checkTriangleOrientation(const Ray &ray, const explicitPoint3D &tv0, const explicitPoint3D &tv1, const explicitPoint3D &tv2)
{
    double res = cinolib::orient3d(tv0.ptr(), tv1.ptr(), tv2.ptr(), ray.v1.ptr());

    assert(res != 0 && "Problem in PointOrientation(...)");

    /* in res we have sign(area(v0, v1, v2, ray.second))
     * if the area is >0 the ray is doing INSIDE -> OUTSIDE, so the patch is INSIDE
     * else the ray is doing OUTSIDE -> INSIDE so the patch is OUTSIDE */
    return (res < 0) ? 1 : 0;
}


// return 1 if inside, 0 if outside
inline uint checkTriangleOrientationRationals(const RationalRay &ray, const std::vector<bigrational> &tv0, const std::vector<bigrational> &tv1, const std::vector<bigrational> &tv2)
{
    bigrational res = cinolib::orient3d(&tv0[0],&tv1[0],&tv2[0], &ray.v1[0]);

    assert(res.sgn() != 0 && "Problem in PointOrientation(...)");

    /* in res we have sign(area(v0, v1, v2, ray.second))
     * if the area is >0 the ray is doing INSIDE -> OUTSIDE, so the patch is INSIDE
     * else the ray is doing OUTSIDE -> INSIDE so the patch is OUTSIDE */
    return (res < bigrational(0,0,0)) ? 1 : 0;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void propagateInnerLabelsOnPatch(const phmap::flat_hash_set<uint> &patch_tris, const std::bitset<NBIT> &patch_inner_label, Labels &labels)
{
    for(uint t_id : patch_tris)
        labels.inside[t_id] = patch_inner_label;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void computeFinalExplicitResult(const FastTrimesh &tm, const Labels &labels, uint num_tris_in_final_res,
                                       std::vector<double> &out_coords, std::vector<uint> &out_tris, 
                                       std::vector<std::bitset<NBIT>> &out_label, bool flat_array)
{
    if(flat_array)
    {
        // loop over triangles and fix vertex indices
        uint num_vertices = 0;
        std::vector<int>  vertex_index(tm.numVerts(), -1);
        out_tris.resize(3 * num_tris_in_final_res);
        out_label.resize(num_tris_in_final_res);
        uint tri_offset = 0;
        for(uint t_id = 0; t_id < tm.numTris(); t_id++)
        {
            if(tm.triInfo(t_id) == 0) continue; // triangle not included in final version
            const uint *triangle = tm.tri(t_id);
            for(uint i = 0; i < 3; i++)
            {
                uint old_vertex = triangle[i];
                if (vertex_index[old_vertex] == -1) {
                    vertex_index[old_vertex] = num_vertices ++;

                }
                out_tris[3 * tri_offset + i] = vertex_index[old_vertex];
            }
            out_label[tri_offset] = labels.surface[t_id];
            tri_offset++;
        }

        // loop over vertices
        out_coords.resize(3 * num_vertices);
        for(uint v_id = 0; v_id < (uint)tm.numVerts(); v_id++) {
            if (vertex_index[v_id] == -1) continue;
            double* v = out_coords.data() + (3 * vertex_index[v_id]);
            tm.vert(v_id)->getApproxXYZCoordinates(v[0], v[1], v[2]);
        }

        // rescale output
        double multiplier = tm.vert(tm.numVerts() - 1)->toExplicit3D().X();
        for(double &c : out_coords) c /= multiplier;
    } else
    {
        out_coords.reserve(3 * 3 * num_tris_in_final_res);
        out_tris.resize(3 * num_tris_in_final_res);
        out_label.resize(num_tris_in_final_res);
        phmap::flat_hash_map<uint, uint> v_map;
        uint tri_offset = 0;

        double multiplier = tm.vert(tm.numVerts() - 1)->toExplicit3D().X();

        for(uint t_id = 0; t_id < tm.numTris(); t_id++)
        {
            if(tm.triInfo(t_id) == 0) continue; // triangle not included in final version

            const uint *v_id = tm.tri(t_id);

            for(uint i = 0; i < 3; i++)
            {
                uint fresh_v_id = (uint)(out_coords.size() / 3);
                auto ins = v_map.insert({v_id[i], fresh_v_id});

                if(ins.second) // vert added
                {
                    double x, y, z;
                    tm.vert(v_id[i])->getApproxXYZCoordinates(x, y, z);
                    out_coords.push_back(x);
                    out_coords.push_back(y);
                    out_coords.push_back(z);
                }

                out_tris[3 * tri_offset + i] = ins.first->second;
            }
            out_label[tri_offset] = labels.surface[t_id];
            tri_offset++;
        }

        for(double &c : out_coords)
            c /= multiplier;
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline uint boolIntersection(FastTrimesh &tm, const Labels &labels)
{
    uint num_tris_in_final_solution = 0;
    tm.resetTrianglesInfo();

    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {
        if((labels.surface[t_id] ^ labels.inside[t_id]).count() == labels.num) // triangle to keep
        {
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_solution++;
        }
    }

    return num_tris_in_final_solution;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline uint boolUnion(FastTrimesh &tm, const Labels &labels)
{
    uint num_tris_in_final_solution = 0;
    tm.resetTrianglesInfo();

    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {   tm.setTriInfo(t_id, 1);
        num_tris_in_final_solution++;
        /*if(labels.inside[t_id].count() == 0) // triangle to keep
        {
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_solution++;
        }*/
    }

    return num_tris_in_final_solution;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// if more than 2 models -> model 0 - all the others
inline uint boolSubtraction(FastTrimesh &tm, const Labels &labels)
{
    uint num_tris_in_final_solution = 0;
    tm.resetTrianglesInfo();

    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {
        if(labels.surface[t_id][0] && labels.inside[t_id].count() == 0) // triangle to keep
        {
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_solution++;
        }
        else if(!labels.surface[t_id][0] && labels.inside[t_id][0] && labels.inside[t_id].count() == 1)
        {
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_solution++;
        }
    }

    //fix triangles orientation
    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {
        if(tm.triInfo(t_id) == 1 && labels.surface[t_id][0] != 1)
            tm.flipTri(t_id);
    }

    return num_tris_in_final_solution;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline uint boolXOR(FastTrimesh &tm, const Labels &labels)
{
    uint num_tris_in_final_solution = 0;
    tm.resetTrianglesInfo();

    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {
        if((labels.inside[t_id].count() == 0) || ((labels.surface[t_id] ^ labels.inside[t_id]).count() == labels.num)) // triangle to keep
        {
            tm.setTriInfo(t_id, 1);
            num_tris_in_final_solution++;
        }
    }

    // fix triangles orientation
    for(uint t_id = 0; t_id < tm.numTris(); t_id++)
    {
        if(tm.triInfo(t_id) == 1 && labels.inside[t_id].count() > 0)
            tm.flipTri(t_id);
    }

    return num_tris_in_final_solution;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline uint bitsetToUint(const bitset<NBIT> &b)
{
    assert(b.count() == 1 && "more than 1 bit set to 1");

    for(uint i = 0; i < NBIT; i++)
        if (b[i]) return i;

    return 0; //warning killer
}


inline uint bitsetToUintCustom(const bitset<NBIT> &b)
{
    //assert(b.count() == 1 && "more than 1 bit set to 1");

    for(uint i = 0; i < NBIT; i++)
        if (b[i]) return i;

    return 0; //warning killer
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool consistentWinding(const uint *t0, const uint *t1) // t0 -> vertex ids of triangle t0, t1 -> vertex ids of triangle t1
{
    int j = 0;

    while(t0[0] != t1[j] && j < 3) j++;
    assert(j < 3 && "not same triangle");

    if (t0[1] == t1[(j+1)%3] && t0[2] == t1[(j+2)%3]) return true;
    else return false;
}

/// :::::::::::::: DEBUG :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline std::string printBitset(const std::bitset<NBIT> &b, uint num_label)
{
    std::string s = b.to_string();
    s = s.substr(b.size()-num_label, b.size());
    std::cerr << s << std::endl;

    return s;
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void saveOutputWithLabels(const std::string &filename, cinolib::Trimesh<> &m, const std::vector<std::bitset<NBIT> > &labels)
{
    std::vector<double> coords(3 * m.num_verts());
    for(uint v_id = 0; v_id < m.num_verts(); v_id++)
    {
        coords[3 * v_id] = m.vert(v_id).x();
        coords[3 * v_id +1] = m.vert(v_id).y();
        coords[3 * v_id +2] = m.vert(v_id).z();
    }

    std::vector<int> int_labels(m.num_polys());
    for(uint t_id = 0; t_id < m.num_polys(); t_id++)
    {
        int l = static_cast<int>(labels[t_id].to_ulong());
        int_labels[t_id] = l;
    }

    cinolib::write_OBJ(filename.c_str(), coords, m.vector_polys(), int_labels);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void loadInputWithLabels(const string &filename, std::vector<double> &coords, std::vector<uint> &tris, std::vector< std::bitset<NBIT> > &labels)
{
    cinolib::Trimesh<> m(filename.c_str());
    m.poly_label_wrt_color();

    coords.reserve(3 * m.num_verts());
    tris.reserve(3 * m.num_polys());
    labels.reserve(m.num_polys());

    for(uint v_id = 0; v_id < m.num_polys(); v_id++)
    {
        coords.push_back(m.vert(v_id).x());
        coords.push_back(m.vert(v_id).y());
        coords.push_back(m.vert(v_id).z());
    }

    for(uint t_id = 0; t_id < m.num_polys(); t_id++)
    {
        tris.push_back(m.poly_vert_id(t_id, 0));
        tris.push_back(m.poly_vert_id(t_id, 1));
        tris.push_back(m.poly_vert_id(t_id, 2));

        std::bitset<NBIT> l(static_cast<unsigned long>(m.poly_data(t_id).label));
        labels.push_back(l);
    }
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

inline void loadInputWithLabels(const string &filename, std::vector<double> &coords, std::vector<uint> &tris, std::vector<uint> &labels)
{
    cinolib::Trimesh<> m(filename.c_str());
    m.poly_label_wrt_color();

    coords.reserve(3 * m.num_verts());
    tris.reserve(3 * m.num_polys());
    labels.reserve(m.num_polys());

    for(uint v_id = 0; v_id < m.num_polys(); v_id++)
    {
        coords.push_back(m.vert(v_id).x());
        coords.push_back(m.vert(v_id).y());
        coords.push_back(m.vert(v_id).z());
    }
    for(uint t_id = 0; t_id < m.num_polys(); t_id++)
    {
        tris.push_back(m.poly_vert_id(t_id, 0));
        tris.push_back(m.poly_vert_id(t_id, 1));
        tris.push_back(m.poly_vert_id(t_id, 2));

        labels.push_back(static_cast<uint>(m.poly_data(t_id).label));
    }
}


inline void printInfoTriangleRationals(RationalRay &rational_ray, std::vector <bigrational> &tv0_aux, std::vector <bigrational> &tv1_aux, std::vector <bigrational> &tv2_aux,
                                       uint *tv_aux, uint &t_id, bool &print_ray){

    if(print_ray){
        std::cout << "Info about ray:" << std::endl;
        std::cout << "Start point: " << "v0: " << rational_ray.v0[0].get_d() << " v1: " << rational_ray.v0[1].get_d() << " v2: " << rational_ray.v0[2].get_d() << std::endl;
        std::cout << "End point: "   << "v0: " << rational_ray.v1[0].get_d() << " v1: " << rational_ray.v1[1].get_d() << " v2: " << rational_ray.v1[2].get_d() << std::endl;
        print_ray = false;
    }

    std::cout << "Info about triangle whith t_id: " << t_id << std::endl;
    std::cout << "v_id: " <<tv_aux[0]<< " v0: " << tv0_aux[0].get_d() << " v1: " << tv0_aux[1].get_d() << " v2: " << tv0_aux[2].get_d() << std::endl;
    std::cout << "v_id: " <<tv_aux[1]<< " v0: " << tv1_aux[0].get_d() << " v1: " << tv1_aux[1].get_d() << " v2: " << tv1_aux[2].get_d() << std::endl;
    std::cout << "v_id: " <<tv_aux[2]<< " v0: " << tv2_aux[0].get_d() << " v1: " << tv2_aux[1].get_d() << " v2: " << tv2_aux[2].get_d() << std::endl;


}

