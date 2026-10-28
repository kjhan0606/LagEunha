import copy
from unittest.mock import patch
import numpy as np
from scipy.spatial import cKDTree
from quadratic_velocity_reconstruction import radial_neighbors,prepare
from test_quadratic_velocity_reconstruction import fixture


def test_candidates_preserve_original_indices_and_distances():
    rng=np.random.default_rng(71)
    for shift in (0.,1e6):
        centers=rng.uniform(-2,2,(384,2))+shift
        tree=cKDTree(centers)
        for i in (0,27,383):
            for scale in (.001,.1,1.,10.):
                radius=np.linalg.norm(centers-centers[i],axis=1)/(3*scale)
                ids=np.flatnonzero((radius<1)&(np.arange(len(centers))!=i))
                actual,distances=radial_neighbors(centers,centers[i],scale,i,tree)
                np.testing.assert_array_equal(actual,ids)
                np.testing.assert_array_equal(distances,radius[ids])


def test_support_boundary_uses_original_strict_predicate():
    points=np.array([[0.,0.],[1.,0.],[np.nextafter(1.,0.),0.],
                     [np.nextafter(1.,2.),0.],[0.,-1.],[.6,.6]])
    old=radial_neighbors(points,points[0],1/3,0)
    new=radial_neighbors(points,points[0],1/3,0,cKDTree(points))
    for a,b in zip(old,new):np.testing.assert_array_equal(a,b)


def test_full_reconstruction_matches_dense_search():
    mass,_,factors,cells,edges,states=fixture()
    mass=np.tile(mass,16);factors=factors*16;cells=cells*16
    dense=copy.deepcopy(states*16);indexed=copy.deepcopy(states*16)
    # Repeated cells exercise the large-N route, including coincident centers.
    # deepcopy each element separately so reconstruction metadata is not shared.
    dense=[copy.deepcopy(s) for s in dense];indexed=[copy.deepcopy(s) for s in indexed]
    with patch.dict('os.environ',{'LAGFORCE_VELOCITY_FIT':'radial_fixed_mean',
            'LAGFORCE_VECTOR_QUADRATIC_BUDGET':'0','LAGFORCE_QUADRATIC_BUDGET':'dominant'}):
        with patch('scipy.spatial.cKDTree',return_value=None):
            reference=prepare(dense,cells,edges,mass,factors)
        actual=prepare(indexed,cells,edges,mass,factors)
    assert actual==reference
    for old,new in zip(dense,indexed):
        for key,value in old['velocity_quadratic'].items():
            np.testing.assert_array_equal(value,new['velocity_quadratic'][key])
