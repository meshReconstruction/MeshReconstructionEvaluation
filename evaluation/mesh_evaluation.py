import glob
import os
import sys

import pymeshlab as ml
import numpy as np


def compute_distances(ms: ml.MeshSet):
    ms.apply_filter("compute_scalar_by_distance_from_another_mesh_per_vertex", measuremesh=ms.current_mesh_id(),
                    refmesh=0)
    ms.apply_filter("compute_scalar_by_distance_from_another_mesh_per_vertex", measuremesh=0,
                    refmesh=ms.current_mesh_id())


def chamfer_distance_from_reference_mesh(ms: ml.MeshSet, mesh_id):
    ms.set_current_mesh(mesh_id)
    m1 = ms.current_mesh()
    ms.set_current_mesh(0)
    m2 = ms.current_mesh()
    dist1 = np.mean(np.abs(m1.vertex_scalar_array()))
    dist2 = np.mean(np.abs(m2.vertex_scalar_array()))
    return 0.5 * (dist1 + dist2), dist1, dist2


def hausdorff_score(ms: ml.MeshSet, face_num, mesh_id):
    res1 = ms.apply_filter("get_hausdorff_distance", sampledmesh=mesh_id, targetmesh=0,
                           samplenum=face_num * 3)
    res2 = ms.apply_filter("get_hausdorff_distance", sampledmesh=0, targetmesh=mesh_id,
                           samplenum=face_num * 3)
    return max(res1['max'], res2['max'])


def f_score_from_reference_mesh(ms: ml.MeshSet, threshold: float, mesh_id):
    ms.set_current_mesh(mesh_id)
    m1 = ms.current_mesh()
    ms.set_current_mesh(0)
    m2 = ms.current_mesh()
    dist1 = m1.vertex_scalar_array()
    dist2 = m2.vertex_scalar_array()
    precision = (dist1 < threshold).sum() / len(dist1)
    recall = (dist2 < threshold).sum() / len(dist2)
    return (2 * precision * recall) / (precision + recall), precision, recall


def sub_folder_eval(directory):
    scores = []
    all_path = glob.glob(os.path.join(directory, "*.ply"))
    print(all_path)
    g_path = glob.glob(os.path.join(directory, "*_gt.ply"))
    if len(g_path) == 0:
        print("No ground-truth mesh was detected!")
        sys.exit()

    ms = ml.MeshSet()
    ms.load_new_mesh(g_path[0])
    ms.apply_filter("meshing_remove_duplicate_vertices")
    face_num = ms.current_mesh().face_number()
    if face_num < 8: face_num = 8

    for a in all_path:
        if g_path[0] != a:
            ms.load_new_mesh(a)
            mesh_id = ms.current_mesh_id()
            ms.apply_filter("meshing_remove_duplicate_vertices")
            compute_distances(ms)
            chamfer = chamfer_distance_from_reference_mesh(ms, mesh_id)
            fscore = f_score_from_reference_mesh(ms, 0.01, mesh_id)
            hausforff = hausdorff_score(ms, face_num, mesh_id)
            scores.append([a, hausforff, chamfer[0], chamfer[1], chamfer[2], fscore[0], fscore[1], fscore[2]])
    return scores


if __name__ == '__main__':
    sub_folder_eval("PATH/TO/MESHES/TO/EVALUATE")