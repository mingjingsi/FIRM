// #include <RcppArmadillo.h>
// #include <Rcpp.h>
// #include <ANN/ANN.h>
#include "nn.hpp"
// using namespace arma;

umat nn2_cpp(arma::mat data, const int k) {

  mat query = data;
  uword d = data.n_cols;
  uword nd = data.n_rows;
  uword nq = query.n_rows;

  ANNkd_tree	*the_tree;	// Search structure

  ANNpointArray data_pts 	= annAllocPts(nd,d);		// Allocate data points
  ANNidxArray nn_idx 		= new ANNidx[k];		// Allocate near neigh indices
  ANNdistArray dists 		= new ANNdist[k];		// Allocate near neighbor dists

  // now construct the points
  for(int i = 0; i < nd; i++)
  {
    for(int j = 0; j < d; j++)
    {
      data_pts[i][j]=data(i,j);
    }
  }
  the_tree = new ANNkd_tree( data_pts, nd, d);

  // return values here
  umat ridx(nq, k);

  //now iterate over query points
  ANNpoint pq = annAllocPt(d);

  for(int i = 0; i < nq; i++)	// Run all query points against tree
  {
    // read coords of current query point
    for(int j = 0; j < d; j++)
    {
      pq[j]=query(i,j);
    }

    the_tree->annkSearch(	// search
        pq,	// query point
        k,		// number of near neighbors
        nn_idx,		// nearest neighbors (returned)
        dists,		// distance (returned)
        0);	// error bound
    for (int j = 0; j < k; j++)
    {
      ridx(i,j) = nn_idx[j];	// put indices in returned array
    }
  }

  annDeallocPt(pq);
  annDeallocPts(data_pts);
  delete the_tree;
  delete [] nn_idx;
  delete [] dists;

  return ridx ;
}

